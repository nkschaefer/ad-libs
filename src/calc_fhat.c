/*
Computes estimated % ancestry a hybrid individual derives from an ancestral
population of interest.
*/
#include <getopt.h>
#include <argp.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
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
struct fhat_stats {
    long int abxba;
    long int baxba;
    long int axbba;
    long int bxaba;
};

/**
 * Function to compute f_hat statistics (numbers of different types of sites)
 * on a given scaffold.
 */
inline const struct fhat_stats fhat_seq(char* anc2, char* hybrid, char* anc1a, char* anc1b, 
    char* out, long int shortest){
    struct fhat_stats stats;
    stats.abxba = 0;
    stats.baxba = 0;
    stats.axbba = 0;
    stats.bxaba = 0;
    
    int baseIndex = 0;
    for (baseIndex; baseIndex < shortest; baseIndex++){
        if (anc2[baseIndex] != 'N' && hybrid[baseIndex] != 'N' &&
            anc1a[baseIndex] != 'N' && anc1b[baseIndex] != 'N' &&
            out[baseIndex] != 'N'){
            
            if (anc2[baseIndex] == out[baseIndex]){
                if (hybrid[baseIndex] != anc2[baseIndex] && 
                    hybrid[baseIndex] == anc1b[baseIndex]){
                    // ABBA site.
                    stats.abxba++;
                }
                if (anc1b[baseIndex] != anc2[baseIndex] &&
                    anc1a[baseIndex] == anc1b[baseIndex]){
                    // ABBA site for ancestral individual.
                    stats.axbba++;   
                }
            } else if (anc2[baseIndex] == anc1b[baseIndex]){
                if (hybrid[baseIndex] != anc2[baseIndex] &&
                    hybrid[baseIndex] == out[baseIndex]){
                    // BABA site.
                    stats.baxba++;
                }
                if (anc1a[baseIndex] != anc2[baseIndex] &&
                    anc1a[baseIndex] == out[baseIndex]){
                    // BABA site for ancestral individual.
                    stats.bxaba++;   
                }
            }   
            
        }
    }
    
    return stats;
}

// Set a variable to tell all functions whether or not to print out
// progress messages, etc.
static int verbose_flag;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
   fprintf(stderr, "calc_fhat [OPTIONS]\n");
   fprintf(stderr, "Estimates the percent ancestry an admixed individual \
derives from a particular ancestral population, given FASTA genomes \
for the admixed individual, two individuals from the ancestral \
population of interest, an individual from the other ancestral/admixing \
population, and an outgroup. If split time between the two ancestral \
populations and the admixture time are known, they can be provided \
in any units to help correct the value given to account for missing \
gene flow between the time of admixture and now. See Green et al 2010, \
Durand et al 2011.\n");
   fprintf(stderr, "[OPTIONS]:\n");
   fprintf(stderr, "    --pop1a -a a FASTA file for an individual from the \
admixing population of interest\n");
   fprintf(stderr, "    --pop1b -b a FASTA file for a second individual from \
the admixing population of interest\n");
   fprintf(stderr, "    --pop2 -2 a FASTA file for an individual from the \
second admixing population\n");
   fprintf(stderr, "    --hybrid -h a FASTA file for the admixed individual \
of interest\n");
   fprintf(stderr, "    --outgroup -o a FASTA file for an individual from \
an outgroup population to all the others\n");
   fprintf(stderr, "    --time_split -s <OPTIONAL> the split time (ago) for \
the two admixing populations, if known, in any units\n");
   fprintf(stderr, "    --time_admixture -t <OPTIONAL> the time of admixture \
(ago) of the two admixing populations, if known, in any units\n");
   fprintf(stderr, "    --verbose -v <OPTIONAL> print status messages to stderr\n");
   fprintf(stderr, "    --num_bases -n <OPTIONAL> the number of bases to scan. \
If this is set, once enough sequences have been read from the FASTA files to \
just exceed this number, file I/O will quit and numbers will be calculated just \
using the bases already read.\n");
   exit(code); 
}

void main(int argc, char *argv[]) {    
    // Define arguments 
    static struct option long_options[] = {
        {"pop1a", required_argument, 0, 'a'},
        {"pop1b", required_argument, 0, 'b'},
        {"pop2", required_argument, 0, '2'},
        {"hybrid", required_argument, 0, 'h'},
        {"outgroup", required_argument, 0, 'o'},
        {"time_split", required_argument, 0, 's'},
        {"time_admixture", required_argument, 0, 't'},
        {"verbose", no_argument, &verbose_flag, 1},
        {"num_bases", required_argument, 0, 'n'},
        {0, 0, 0, 0}
    };

    // Set default values
    kseq_t* pop1a = NULL;
    kseq_t* pop1b = NULL;
    kseq_t* pop2 = NULL;
    kseq_t* hybrid = NULL;
    kseq_t* outgroup = NULL;
    float time_split = -1;
    float time_admixture = -1;
    long int base_limit = -1;
    
    int option_index = 0;
    int ch;
    
    // Pointer to end character of strings when converting to floats(?)
    char *end;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "a:b:2:h:o:s:t:n:v", long_options, 
        &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'a':
                pop1a = parse_fasta(optarg);
                break;
            case 'b':
                pop1b = parse_fasta(optarg);
                break;
            case '2':
                pop2 = parse_fasta(optarg);
                break;
            case 'h':
                hybrid = parse_fasta(optarg);
                break;
            case 'o':
                outgroup = parse_fasta(optarg);
                break;
            case 's':
                time_split = strtof(optarg, &end);
                break;
            case 't':
                time_admixture = strtof(optarg, &end);
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
    if (pop1a == NULL || pop1b == NULL || pop2 == NULL || hybrid == NULL || 
        outgroup == NULL){
        fprintf(stderr, 
            "ERROR: one or more required arguments were not provided. Aborting.\n");
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "\n");
        help(1);   
    }
    
    /**
     * Initialize variables
     */
    float abxba = 0;
    float baxba = 0;
    float axbba = 0;
    float bxaba = 0;
    
    long int num_bases = 0;
    
    // Iterate through every sequence in the FASTA file; for each,
    // compute relevant statistics for that sequence.
    int progress_1a = 0;
    int progress_1b = 0;
    int progress_2 = 0;
    int progress_hybrid = 0;
    int progress_out = 0;
    
    while(progress_1a = kseq_read(pop1a) >= 0){
        if (verbose_flag){
            fprintf(stderr, "Processing seq %s\n", pop1a->name.s);
        }

        // Advance all files by one sequence.
        progress_1b = kseq_read(pop1b);
        progress_2 = kseq_read(pop2);
        progress_hybrid = kseq_read(hybrid);
        progress_out = kseq_read(outgroup);
        
        // Process this sequence.
        
        // First make sure we are reading the same sequence in each individual.
        if (strcmp(pop1a->name.s, pop1b->name.s) != 0 || 
            strcmp(pop1a->name.s, pop2->name.s) != 0 
            || strcmp(pop1a->name.s, hybrid->name.s) != 0 || 
            strcmp(pop1a->name.s, outgroup->name.s) != 0){
            fprintf(stderr, "%s\n", pop1a->name.s);
            fprintf(stderr, "%s\n", pop1b->name.s);
            fprintf(stderr, "%s\n", pop2->name.s);
            fprintf(stderr, "%s\n", hybrid->name.s);
            fprintf(stderr, "%s\n", outgroup->name.s);
            fprintf(stderr, 
                "ERROR: sequence IDs do not match; files not sorted in the same order.\n");
            exit(1);
        }
        
        // Determine shortest sequence.
        long int shortest = pop1a->seq.l;
        if (pop1b->seq.l < shortest){
            shortest = pop1b->seq.l;
        }
        if (pop2->seq.l < shortest){
            shortest = pop2->seq.l;
        }
        if (hybrid->seq.l < shortest){
            shortest = hybrid->seq.l;
        }
        if (outgroup->seq.l < shortest){
            shortest = outgroup->seq.l;
        }
        
        // Calculate statistics on this scaffold.
        struct fhat_stats stats = fhat_seq(pop2->seq.s, hybrid->seq.s,
            pop1a->seq.s, pop1b->seq.s, outgroup->seq.s, shortest);
        
        // Add statistics to running totals.
        abxba += stats.abxba;
        baxba += stats.baxba;
        axbba += stats.axbba;
        bxaba += stats.bxaba;
                   
        // If any sequence has reached the end of the file, bail out here.
        if (progress_1b < 0 || progress_2 < 0 || progress_hybrid < 0 || 
            progress_out < 0){
            break;   
        }
        
        // If we are set to stop scanning after a given number of bases,
        // bail out here.
        num_bases += shortest;
        if (base_limit > -1 && num_bases >= base_limit){
            if (verbose_flag){
                fprintf(stderr, "Given limit of %ld bases reached. Exiting.\n", 
                    base_limit);
            }
            break;
        }
    }
    
    // Perform genome-wide fhat calculation.
    
    float fhat = (float)(abxba-baxba)/(float)(axbba-bxaba);
    
    if (time_split > 0 && time_admixture > 0){
        float fhat_upper = (time_split/(time_split-time_admixture)) * fhat;
        fhat = (fhat + fhat_upper)/2;
    }
    
    printf("%f\n", fhat);
    
    kseq_destroy(pop1a);
    kseq_destroy(pop1b);
    kseq_destroy(pop2);
    kseq_destroy(hybrid);
    kseq_destroy(outgroup);
}

