/*
Computes estimated % ancestry a hybrid individual derives from an ancestral
population of interest.

*/
#include <getopt.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kseq.h"

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

// Data structure to hold statistics for a single scaffold
struct pi_stats {
    long int diffs_pop1;
    long int diffs_pop2;
    long int diffs_between;
    long int seqlen_pop1;
    long int seqlen_pop2;
    long int seqlen_between;
};

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

/**
 * Function to compute numbers needed for computation of pi, on a single scaffold.
 */
const struct pi_stats pi_seq(kseq_t* pop1[], kseq_t* pop2[], int num_pop1, 
    int num_pop2, long int shortest){
    struct pi_stats stats;
    stats.diffs_pop1 = 0;
    stats.diffs_pop2 = 0;
    stats.diffs_between = 0;
    stats.seqlen_pop1 = shortest;
    stats.seqlen_pop2 = shortest;
    stats.seqlen_between = shortest;
    
    for (int baseIndex = 0; baseIndex < shortest; baseIndex++){
        
        long int this_pop1_diffs = 0;
        long int this_pop2_diffs = 0;
        long int this_between_diffs = 0;
        
        int pop1_n = 0;
        int pop2_n = 0;
        
        for (int pop1_index = 0; pop1_index < num_pop1; pop1_index++){
            char pop1chr = capitalize(pop1[pop1_index]->seq.s[baseIndex]);
            if (pop1chr == 'N'){
                pop1_n = 1;
                break;
            }
            if (pop1_index > 0 && pop1_n == 0){
                for (int pop1_index2 = 0; pop1_index2 < pop1_index; pop1_index2++){
                    char pop1bchr = capitalize(pop1[pop1_index2]->seq.s[baseIndex]);
                    // We already checked this base as pop1chr, so we don't
                    // need to check to see if it's N.
                    if (pop1chr != pop1bchr){
                        this_pop1_diffs++;
                    }
                }
            }
            for (int pop2_index = 0; pop2_index < num_pop2; pop2_index++){
                char pop2chr = capitalize(pop2[pop2_index]->seq.s[baseIndex]);
                if (pop2chr == 'N'){
                    pop2_n = 1;
                    break;
                }
                if (pop2chr != pop1chr){
                    this_between_diffs++;
                }
            }
        }
        
        for (int pop2_index = 0; pop2_index < num_pop2; pop2_index++){
            char pop2chr = capitalize(pop2[pop2_index]->seq.s[baseIndex]);
            if (pop2chr == 'N' || pop2_n == 1){
                pop2_n = 1;
                break;
            }
            if (pop2_index > 0 && pop2_n == 0){
                for (int pop2_index2 = 0; pop2_index2 < pop2_index; pop2_index2++){
                    char pop2bchr = capitalize(pop2[pop2_index2]->seq.s[baseIndex]);
                    // We already checked this base as pop2chr, so we don't
                    // need to check to see if it's N.
                    if (pop2chr != pop2bchr){
                        this_pop2_diffs++;
                    }
                }
            }
        }
        
        // Skip site if there are any Ns.
        if (pop1_n == 0 && pop2_n == 0){
            stats.diffs_pop1 += this_pop1_diffs;
            stats.diffs_pop2 += this_pop2_diffs;
            stats.diffs_between += this_between_diffs;
        } else{
            stats.seqlen_pop1--;
            stats.seqlen_pop2--;
            stats.seqlen_between--;
        }
        /*
        if (pop1_n == 0){
            stats.diffs_pop1 += this_pop1_diffs;
            if (pop2_n == 0){
                stats.diffs_between += this_between_diffs;
            }
        } else{
            stats.seqlen_pop1--;
        } 
        if (pop2_n == 0){
            stats.diffs_pop2 += this_pop2_diffs;
        } else{
            stats.seqlen_pop2--;
            if (pop1_n == 1){
                stats.seqlen_between--;
            }
        }
        */
        
    }
    
    return stats;
}

const int binomial_coeff(int n, int k){
    int numerator = 1;
    int denominator = 1;
    
    for (int term = n; term >= n-(k-1); term--){
        numerator = numerator*term;
    }
    for (int term = k; term >= 1; term--){
        denominator = denominator*term;
    }
    return numerator/denominator;
}

// Set a variable to tell all functions whether or not to print out
// progress messages, etc.
static int verbose_flag;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
   fprintf(stderr, "pi_pops [OPTIONS]\n");
   fprintf(stderr, "Given multiple FASTA files representing two populations, \
calculates nucleotide diversity per site (pi) within both populations and \
between them.\n");
   fprintf(stderr, "[OPTIONS]:\n");
   fprintf(stderr, "    -1 a FASTA file for an individual from the first \
population (specify multiple times)\n");
   fprintf(stderr, "    -2 a FASTA file for an individual from the \
second population (specify multiple times)\n");
   fprintf(stderr, "    -v <OPTIONAL> print status messages to stderr\n");
   fprintf(stderr, "    -n <OPTIONAL> the number of bases to scan. \
If this is set, once enough sequences have been read from the FASTA files to \
just exceed this number, file I/O will quit and numbers will be calculated just \
using the bases already read.\n");
   exit(code); 
}

int main(int argc, char *argv[]) {
    
    // Set default argument values
    int num_pop1 = 0;
    int num_pop2 = 0;
    kseq_t* pop1[10];
    kseq_t* pop2[10];
    long int base_limit = -1;
    
    int option_index = 0;
    int ch;
    
    // Pointer to end character of strings when converting to floats(?)
    char *end;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt(argc, argv, "1:2:n:v")) != -1){
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
            case 'v':
                verbose_flag = 1;
                break;
            case 'n':
                base_limit = strtol(optarg, &end, 10);
            case '?':
                //help(0);
                break;
            default:
                help(0);
        }    
    }

    // Validate arguments.
    if (num_pop1 < 2 || num_pop2 < 2){
        fprintf(stderr, "ERROR: You must provide at least two individuals from each \
population.\n");
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "\n");
        help(1);   
    }
    
    /**
     * Initialize variables
     */
    float diffs_pop1_sum = 0;
    float diffs_pop2_sum = 0;
    float diffs_between_sum = 0;
    float seqlen_pop1 = 0;
    float seqlen_pop2 = 0;
    float seqlen_between = 0;
    
    long long int seqlen_tot = 0;
    
    // Number of comparisons can be computed via the binomial coefficient
    float comps_pop1 = binomial_coeff(num_pop1, 2);
    float comps_pop2 = binomial_coeff(num_pop2, 2);
    float comps_between = num_pop1*num_pop2;
    
    // Iterate through every sequence in the FASTA file; for each,
    // compute relevant statistics for that sequence.
    
    // Set progress through all files to 0 initially
    int progress_1[num_pop1];
    int progress_2[num_pop2];
    
    for (int pop1_index = 0; pop1_index < num_pop1; pop1_index++){
        progress_1[pop1_index] = 0;
    }
    
    for (int pop2_index = 0; pop2_index < num_pop2; pop2_index++){
        progress_2[pop2_index] = 0;
    }
    
    while((progress_1[0] = kseq_read(pop1[0])) >= 0){
        if (verbose_flag){
            fprintf(stderr, "Processing seq %s\n", pop1[0]->name.s);
        }
        
        // Advance all files by one sequence.
        for (int pop1_index = 1; pop1_index < num_pop1; pop1_index++){
            progress_1[pop1_index] = kseq_read(pop1[pop1_index]);
        }
        for (int pop2_index = 0; pop2_index < num_pop2; pop2_index++){
            progress_2[pop2_index] = kseq_read(pop2[pop2_index]);
        }

        // Process this sequence.
        
        // First, determine shortest sequence.
        long int shortest = pop1[num_pop1-1]->seq.l;
        for (int pop1_index = 0; pop1_index < num_pop1-1; pop1_index++){
            if (pop1[pop1_index]->seq.l < shortest){
                shortest = pop1[pop1_index]->seq.l;
            }
        }
        for (int pop2_index = 0; pop2_index < num_pop2; pop2_index++){
            if (pop2[pop2_index]->seq.l < shortest){
                shortest = pop2[pop2_index]->seq.l;
            }
        }
        
        seqlen_tot += shortest;
        
        // Calculate statistics on this scaffold.
        struct pi_stats stats = pi_seq(pop1, pop2, num_pop1, num_pop2, shortest);
        
        // Add statistics to running totals.
        diffs_pop1_sum += stats.diffs_pop1;
        diffs_pop2_sum += stats.diffs_pop2;
        diffs_between_sum += stats.diffs_between;
        seqlen_pop1 += stats.seqlen_pop1;
        seqlen_pop2 += stats.seqlen_pop2;
        seqlen_between += stats.seqlen_between;
        
        // If any sequence has reached the end of the file, bail out here.
        for (int pop1_index = 0; pop1_index < num_pop1; pop1_index++){
            if (progress_1[pop1_index] < 0){
                break;
            }
        }
        for (int pop2_index = 0; pop2_index < num_pop2; pop2_index++){
            if (progress_2[pop2_index] < 0){
                break;
            }
        }
        
        // If we are set to stop scanning after a given number of bases,
        // bail out here.
        if (base_limit > -1 && seqlen_tot >= base_limit){
            if (verbose_flag){
                fprintf(stderr, "Given limit of %ld bases reached. Exiting.\n", 
                    base_limit);
            }
            break;
        }
    }
    
    
    // Calculate pi.
    float pi1 = (float)diffs_pop1_sum/(float)comps_pop1/(float)seqlen_pop1;
    float pi2 = (float)diffs_pop2_sum/(float)comps_pop2/(float)seqlen_pop2;
    float pibetween = (float)diffs_between_sum/(float)comps_between/(float)seqlen_between;
    
    // Output format: pi1 pi2 pibetween fst (tab separated)
    printf("%f\t%f\t%f\n", pi1, pi2, pibetween);
    
    // Clean up.
    for (int pop1_index = 0; pop1_index < num_pop1; pop1_index++){
        kseq_destroy(pop1[pop1_index]);
    }
    for (int pop2_index = 0; pop2_index < num_pop2; pop2_index++){
        kseq_destroy(pop2[pop2_index]);
    }

    return 0;
}

