/**
 Given FASTA files for ancestral population sequences and a "query" individual
 sequence, reads through the sequences in windows and computes scores to be used
 by AD-LIBS. Output should be piped to ad-libs-hmm.py.
**/
#include <getopt.h>
#include <argp.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "kseq.h"

// Initialize FASTA parser
KSEQ_INIT(gzFile, gzread);

/**
 * Function to parse FASTA, including gzipped files.
 */
kseq_t* parse_fasta(char *filename){
    gzFile fp = gzopen(filename, "r");
    if (!fp){
        fprintf(stderr, "ERROR: unable to open file %s\n", filename);
        exit(1);
    } else{
        kseq_t* seq;
        seq = kseq_init(fp);
        return seq;
    }
}

const char capitalize(char base){
    if (base == 'a'){
        return 'A';
    } else if (base == 'c'){
        return 'C';
    } else if (base == 'g'){
        return 'G';
    } else if (base == 't'){
        return 'T';
    } else if (base == 'n'){
        return 'N';
    }
    return base;
}

// Define growable array to store IBS tract lengths.
// http://stackoverflow.com/questions/3536153/c-dynamically-growing-array
struct flex_array {
    size_t len;
    int index;
    int *arr;
};

void init_array(struct flex_array *a){
    a->arr = (int *) malloc(11 * sizeof(int));
    a->index = 0;
    a->len = 10;
}

void add_item(struct flex_array *a, int item){
    if (a->index == a->len){
        a->len *= 2;
        a->arr = (int *) realloc(a->arr, (a->len) * sizeof(int) + 1);
        if (! a->arr){
            fprintf(stderr, "ERROR reallocating memory.\n");
            exit(1);
        }
    }
    a->arr[a->index] = item;
    a->index++;
}

void free_array(struct flex_array *a){
    // First zero it out to make sure it doesn't mess things up in the
    // future.
    free(a->arr);
    a->arr = NULL;
    a->index = 0;
    a->len = 0;
}

/**
 * Function to generate random numbers between 0 and 1
 */
float get_rand(){
    float r = (float) rand() / (float) RAND_MAX;
    return r;
}

/**
 * Function to compute a score in a given window.
 */
const float calc_score_window(kseq_t* pop1[], int num_pop1, \
    kseq_t* pop2[], int num_pop2, kseq_t* hybrid, long int win_start, long int win_end, \
    float skip, int skipscore, int mask_cpg){
    
    long int baseIndex = win_start;
    int pop1_index = 0;
    int pop2_index = 0;
        
    // These arrays will keep track of whether or not we are currently in 
    // an IBS tract with each individual from each ancestral population.
    // If so, they will store the index at which that match began.
    long int match_1[num_pop1];
    for (pop1_index; pop1_index < num_pop1; pop1_index++){
        match_1[pop1_index] = -1;
    }
    pop1_index = 0;
    long int match_2[num_pop2];
    for (pop2_index; pop2_index < num_pop2; pop2_index++){
        match_2[pop2_index] = -1;
    }
    pop2_index = 0;
    
    // These arrays will store an array of each IBS tract length found with
    // each individual in both of the ancestral populations.
    struct flex_array *ibs_1 = malloc(num_pop1 * sizeof(struct flex_array));
    for (pop1_index; pop1_index < num_pop1; pop1_index++){
        struct flex_array arr;
        init_array(&arr);
        memcpy(&ibs_1[pop1_index], &arr, sizeof(struct flex_array));
    }
    pop1_index = 0;
    
    struct flex_array *ibs_2 = malloc(num_pop2 * sizeof(struct flex_array));
    for (pop2_index; pop2_index < num_pop2; pop2_index++){
        struct flex_array arr;
        init_array(&arr);
        memcpy(&ibs_2[pop2_index], &arr, sizeof(struct flex_array));
    }
    pop2_index = 0;
    
    // Store the number of "N" bases found in query sequence
    int n_count_query = 0;
    
    // Store number of "N" bases found in all pop 1 sequences
    float n_count_1 = 0;
    // Store number of "N" bases found in all pop 2 sequences
    float n_count_2 = 0;
    
    int added_count = 0;
    
    for (baseIndex; baseIndex < win_end; baseIndex++){
        
        unsigned char pop1chr;
        unsigned char pop2chr;
        unsigned char hchr = capitalize(hybrid->seq.s[baseIndex]);
        if (mask_cpg && ((hchr == 'G' && baseIndex > 0 && 
            capitalize(hybrid->seq.s[baseIndex-1]) == 'C') ||
            (hchr == 'C' && baseIndex < win_end-1 && 
            capitalize(hybrid->seq.s[baseIndex+1]) == 'G'))){
            hchr = 'N';
        }
        if (hchr == 'N'){
            n_count_query++;
        }
        
        for (pop1_index; pop1_index < num_pop1; pop1_index++){
            pop1chr = capitalize(pop1[pop1_index]->seq.s[baseIndex]);
            if (mask_cpg && ((pop1chr == 'G' && baseIndex > 0 && 
                capitalize(pop1[pop1_index]->seq.s[baseIndex-1]) == 'C') ||
                (pop1chr == 'C' && baseIndex < win_end-1 && 
                capitalize(pop1[pop1_index]->seq.s[baseIndex+1]) == 'G'))){
                pop1chr = 'N';
            }
                    
            if (pop1chr == 'N' || hchr == 'N'){
                // This is a tough decision to make. We can't always count an N as
                // a match or a mismatch, as that will mess things up.
                
                // We could always break here with a certain probability, but then
                // we would be biasing toward looking heterozygous.
                
                // Instead, we will look to see if any IBS tracts have already been found
                // with this individual (without seeing Ns). If they have, we will break
                // on an N with a probability determined by the length of previously-seen
                // IBS tracts.
                                
                if (ibs_1[pop1_index].index > 0){
                    float prev_ibs_sum = 0;
                    float prev_ibs_tot = 0;
                    int prev_ibs_index = 0;
                    for (prev_ibs_index; prev_ibs_index < ibs_1[pop1_index].index; \
                        prev_ibs_index++){
                        prev_ibs_sum += ibs_1[pop1_index].arr[prev_ibs_index];
                        prev_ibs_tot++;   
                    }
                    float prev_ibs_avg = prev_ibs_sum/prev_ibs_tot;
                    float r = get_rand();
                    if (r < 1.0/prev_ibs_avg){
                        
                        if (match_1[pop1_index] > -1){
                            // If we're currently in an IBS block, end it.
                            int ibs_len = baseIndex-1-match_1[pop1_index];
                            if (ibs_len > 1){
                                add_item(&ibs_1[pop1_index], ibs_len);
                            }
                            match_1[pop1_index] = -1;
                        }
                        else{
                            // If we're not currently in an IBS block, begin one.
                            match_1[pop1_index] = baseIndex;
                        }
                    }
                }
                // Note the fact that we saw an N here.
                if (pop1chr == 'N'){
                    n_count_1++;
                }
            }
            else if (pop1chr == hchr){
                if (match_1[pop1_index] == -1){
                    match_1[pop1_index] = baseIndex;
                }
                // If we're already tracking an IBS tract, do nothing until it
                // ends.
            }
            else if (match_1[pop1_index] > -1){
                // If the bases do not match and we're currently in an IBS
                // tract, end that IBS tract and store its length.
                int ibs_len = baseIndex-1-match_1[pop1_index];
                
                if (ibs_len > 1){
                    add_item(&ibs_1[pop1_index], ibs_len);
                }
                match_1[pop1_index] = -1;
            }
        }
        pop1_index = 0;
        
        for (pop2_index; pop2_index < num_pop2; pop2_index++){
            pop2chr = capitalize(pop2[pop2_index]->seq.s[baseIndex]);
            if (mask_cpg && ((pop1chr == 'G' && baseIndex > 0 && 
                capitalize(pop1[pop1_index]->seq.s[baseIndex-1]) == 'C') ||
                (pop1chr == 'C' && baseIndex < win_end-1 && 
                capitalize(pop1[pop1_index]->seq.s[baseIndex+1]) == 'G'))){
                pop1chr = 'N';
            }
            if (pop2chr == 'N' || hchr == 'N'){
                // (see above)
                if (ibs_2[pop2_index].index > 0){
                    float prev_ibs_sum = 0;
                    float prev_ibs_tot = 0;
                    int prev_ibs_index = 0;
                    for (prev_ibs_index; prev_ibs_index < ibs_2[pop2_index].index; \
                        prev_ibs_index++){
                        prev_ibs_sum += ibs_2[pop2_index].arr[prev_ibs_index];
                        prev_ibs_tot ++;   
                    }
                    float prev_ibs_avg = prev_ibs_sum/prev_ibs_tot;
                    float r = get_rand();
                    if (r < 1.0/prev_ibs_avg){
                        
                        if (match_2[pop2_index] > -1){
                            // If we're currently in an IBS block, end it.
                            int ibs_len = baseIndex-1-match_2[pop2_index];
                            if (ibs_len > 1){
                                add_item(&ibs_2[pop2_index], ibs_len);
                            }
                            match_2[pop2_index] = -1;
                        }
                        else{
                            // If we're not currently in an IBS block, begin one.
                            match_2[pop2_index] = baseIndex;
                        }
                    }
                }
                // Note the fact that we saw an N here.
                if (pop2chr == 'N'){
                    n_count_2++;
                }
            }
            else if (pop2chr == hchr){
                if (match_2[pop2_index] == -1){
                    match_2[pop2_index] = baseIndex;
                }
                // If we're already tracking an IBS tract, do nothing until it
                // ends.
            }
            else if (match_2[pop2_index] > -1){
                // If the bases do not match and we're currently in an IBS
                // tract, end that IBS tract and store its length.
                int ibs_len = baseIndex-1-match_2[pop2_index];
                if (ibs_len > 1){
                    add_item(&ibs_2[pop2_index], ibs_len);
                }
                match_2[pop2_index] = -1;
            }
  
        }
        pop2_index = 0;
    }
    
    
    // The code below will add any "final" IBS tracts for which a beginning, but not
    // end, was found. Doing this will probably bias things downward (since these are
    // incomplete IBS tracts).
    /**
    for (pop1_index; pop1_index < num_pop1; pop1_index++){
        //if (match_1[pop1_index] > -1){
        if (match_1[pop1_index] > -1 && ibs_1[pop1_index].index == 0){
            add_item(&ibs_1[pop1_index], win_end-match_1[pop1_index]);
        }
    }
    pop1_index = 0;
    for (pop2_index; pop2_index < num_pop2; pop2_index++){
        //if (match_2[pop2_index] > -1){
        if (match_2[pop2_index] > -1 && ibs_2[pop2_index].index == 0){
            add_item(&ibs_2[pop2_index], win_end-match_2[pop2_index]);
        }
    }
    pop2_index = 0;
    **/
    
    // Calculate score.
    
    // First, determine whether to skip. If so, we can return early.
    float winsize = win_end - win_start;
    if ((float)n_count_query/winsize >= skip ||
        (float)n_count_1/winsize/(float)num_pop1 >= skip || 
        (float)n_count_2/winsize/(float)num_pop2 >= skip){
        // Free everything.
        for (pop1_index; pop1_index < num_pop1; pop1_index++){
            free_array(&ibs_1[pop1_index]);
        }
        pop1_index = 0;
        for (pop2_index; pop2_index < num_pop2; pop2_index++){
            free_array(&ibs_2[pop2_index]);
        }
        pop2_index = 0;
        return skipscore;
    }
    
    // Calculate average IBS tract lengths for both ancestral populations.
    
    float ibs_sum_1 = 0;
    float ibs_tot_1 = 0;
    float ibs_sum_2 = 0;
    float ibs_tot_2 = 0;
    
    pop1_index = 0; 
    for (pop1_index; pop1_index < num_pop1; pop1_index++){
        int pop_total = (int) ibs_1[pop1_index].index;
        int ibs_index = 0;
        for (ibs_index; ibs_index < pop_total; ibs_index++){
            ibs_sum_1 += (float) ibs_1[pop1_index].arr[ibs_index];
            ibs_tot_1++;
        }
    }
    
    if (ibs_tot_1 == 0){
        ibs_tot_1 = 1;
    }
    
    pop2_index = 0;
    for (pop2_index; pop2_index < num_pop2; pop2_index++){
        int ibs_index = 0;
        for (ibs_index; ibs_index < ibs_2[pop2_index].index; ibs_index++){
            ibs_sum_2 += (float) ibs_2[pop2_index].arr[ibs_index];
            ibs_tot_2++;
        }
    }
    if (ibs_tot_2 == 0){
        ibs_tot_2 = 1;
    }
    
    float ibs_avg_1 = ibs_sum_1/ibs_tot_1;
    float ibs_avg_2 = ibs_sum_2/ibs_tot_2;
    // Free everything.
    for (pop1_index; pop1_index < num_pop1; pop1_index++){
        free_array(&ibs_1[pop1_index]);
    }
    pop1_index = 0;
    for (pop2_index; pop2_index < num_pop2; pop2_index++){
        free_array(&ibs_2[pop2_index]);
    }
    pop2_index = 0;
    if (ibs_avg_1 == 0 || ibs_avg_2 == 0){
        return skipscore;
    }
    else{
        return log(ibs_avg_1) - log(ibs_avg_2);
    }
}

// Set a variable to tell all functions whether or not to print out
// progress messages, etc.
static int verbose_flag;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
   fprintf(stderr, "adlibs_score [OPTIONS]\n");
   fprintf(stderr, "Given multiple FASTA files representing two populations, \
a window size, and a FASTA file for a query hybrid individual, divides the \
sequences into windows and computes scores on each. Output should then be \
piped to ad-libs-hmm.py.\n");
   fprintf(stderr, "[OPTIONS]:\n");
   fprintf(stderr, "    --pop1 -1 a FASTA file for an individual from the first \
population (specify multiple times)\n");
   fprintf(stderr, "    --pop2 -2 a FASTA file for an individual from the \
second population (specify multiple times)\n");
   fprintf(stderr, "    --window -w a window size (integer # of bases) to use\n");
   fprintf(stderr, "    --hybrid -h a FASTA file for the hybrid individual\n");
   fprintf(stderr, "    --skip, -s a value between 0 and 1 representing what \
percent of all bases across all sequences may be \"N\" before a window is considered \
invalid (default 0.25).\n");
   fprintf(stderr, "    --skipscore, -z an integer representing scores for skipped \
windows. (default 999).\n");
    fprintf(stderr, "   --mask_cpg, -c a flag to mask out (ignore/treat as Ns) CpG sites \
in all sequences. (default=False).\n");
   fprintf(stderr, "    --verbose -v <OPTIONAL> print status messages to stderr\n");
   exit(code); 
}

int main(int argc, char *argv[]) {    
    // Define arguments 
    static struct option long_options[] = {
       {"pop1", required_argument, 0, '1'},
       {"pop2", required_argument, 0, '2'},
       {"window", required_argument, 0, 'w'},
       {"hybrid", required_argument, 0, 'h'},
       {"skip", optional_argument, 0, 's'},
       {"skipscore", optional_argument, 0, 'z'},
       {"mask_cpg", optional_argument, 0, 'm'},
       {"verbose", no_argument, &verbose_flag, 1},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    int num_pop1 = 0;
    int num_pop2 = 0;
    kseq_t* pop1[10];
    kseq_t* pop2[10];
    kseq_t* hybrid;
    long int window = 0;
    float skip = 0.25;
    int skipscore = 999;
    int mask_cpg = 0;
    
    int option_index = 0;
    int ch;
    
    // Pointer to end character of strings when converting to floats(?)
    //char *end;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "1:2:w:h:s:z:m:v", long_options, &option_index )) 
        != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case '1':
                pop1[num_pop1] = parse_fasta(optarg);
                num_pop1++;
                break;
            case '2':
                pop2[num_pop2] = parse_fasta(optarg);
                num_pop2++;
                break;
            case 'w':
                window = atoi(optarg);
                break;
            case 'h':
                hybrid = parse_fasta(optarg);
                break;
            case 's':
                skip = atof(optarg);
                break;
            case 'z':
                skipscore = atoi(optarg);
                break;
            case 'm':
                mask_cpg = 1;
                break;
            case 'v':
                verbose_flag = 1;
                break;
            case '?':
                //help(0);
                break;
            default:
                help(0);
        }    
    }
    
    // Error check.
    if (window == 0){
        fprintf(stderr, "ERROR: please provide a valid window width in base pairs.\n");
        exit(1);
    }
    if (skip < 0 || skip > 1){
        fprintf(stderr, "ERROR: please provide a valid skip threshold between 0 and 1.\n");
        exit(1);
    }
    
    // Seed random number generator
    srand(time(NULL));
    
    int progress_h;
    int progress_1[num_pop1];
    int progress_2[num_pop2];
    
    int pop1_index = 0;
    int pop2_index = 0;
    
    long int win_start = 0;
    long int win_end = 0;
    
    //float score;
    
    while(progress_h = kseq_read(hybrid) >= 0){
        if (verbose_flag){
            fprintf(stderr, "Processing seq %s\n", hybrid->name.s);
        }
        
        // Advance all files by one sequence.
        for (pop1_index; pop1_index < num_pop1; pop1_index++){
            progress_1[pop1_index] = kseq_read(pop1[pop1_index]);
            if (strcmp(pop1[pop1_index]->name.s, hybrid->name.s) != 0){
                fprintf(stderr, "ERROR: sequence IDs %s and %s do not match. \
Please make sure sequences are sorted in the same order.\n", \
pop1[pop1_index]->name.s, hybrid->name.s);
                exit(1);
            }
        }
        pop1_index = 0;
        for (pop2_index; pop2_index < num_pop2; pop2_index++){
            progress_2[pop2_index] = kseq_read(pop2[pop2_index]);
            if (strcmp(pop2[pop2_index]->name.s, hybrid->name.s) != 0){
                fprintf(stderr, "ERROR: sequence IDs %s and %s do not match. \
Please make sure sequences are sorted in the same order.\n", \
pop2[pop2_index]->name.s, hybrid->name.s);
                exit(1);
            }
        }
        pop2_index = 0;

        // Process this sequence.
        
        // First, determine shortest sequence.
        long int shortest = pop1[num_pop1-1]->seq.l;
        for (pop1_index; pop1_index < num_pop1-1; pop1_index++){
            if (pop1[pop1_index]->seq.l < shortest){
                shortest = pop1[pop1_index]->seq.l;
            }
        }
        pop1_index = 0;
        for (pop2_index; pop2_index < num_pop2; pop2_index++){
            if (pop2[pop2_index]->seq.l < shortest){
                shortest = pop2[pop2_index]->seq.l;
            }
        }
        pop2_index = 0;
        
        // Skip this sequence if ANY individual is missing it.
        if (shortest < window){
            continue;
        }
        else{
            // Calculate score in every window.
            for (win_start; win_start < shortest; win_start = win_start + window){
                win_end = win_start + window;
                if (win_end >= shortest){
                    win_end = shortest-1;
                }
                float score = calc_score_window(pop1, num_pop1, pop2, num_pop2, \
                    hybrid, win_start, win_end, skip, skipscore, mask_cpg);
            
                // Print score to stdout (in BED format).
                if (score == skipscore){
                    fprintf(stdout, "%s\t%ld\t%ld\tinf\n", hybrid->name.s, win_start, 
                        win_end);
                }
                else{
                    fprintf(stdout, "%s\t%ld\t%ld\t%g\n", hybrid->name.s, win_start, 
                        win_end, score);
                }
            }
            win_start = 0;
            win_end = 0;
        }
        
    }
    
    // Clean up.
    kseq_destroy(hybrid);
    
    for (pop1_index; pop1_index < num_pop1; pop1_index++){
        kseq_destroy(pop1[pop1_index]);
    }
    for (pop2_index; pop2_index < num_pop2; pop2_index++){
        kseq_destroy(pop2[pop2_index]);
    }
    
    return 0;
}

