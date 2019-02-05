#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
import argparse
import pomegranate
from adlibs_dists import *
from adlibs_math import *
import subprocess
import random
from collections import Counter, defaultdict
import numpy
from networkx import DiGraph
import math
import gzip
from allele_freq_drift import recomb_resample_prob, recomb_resample_prob_aux,\
    phi_diffusion_t, phi_diffusion
import os
from multiprocessing import Pool
import time
"""
adlibs.py

This is the main program for running AD-LIBS. This program can be used to estimate
parameters or to run AD-LIBS. Running can be accomplished either through a 
configuration file or via command-line arguments. One run can include multiple
hybrid individuals.
"""

# Get directory of this script
ADLIBS_DIR = os.path.dirname(os.path.realpath(__file__))

overlap_thresh = 0.25
prob_lim = 0.1

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument("--ancestral1", "-1", help=\
        "The names of one or more FASTA files representing ancestral population 1",
        nargs="+",
        required=True)
    parser.add_argument("--ancestral2", "-2", help=\
        "The names of one or more FASTA files representing ancestral population 2",
        nargs="+",
        required=True)
    parser.add_argument("--hybrid", "-h", help=\
        "The name of a FASTA file representing the query sequence.",
        required=True,
        nargs="+")
    parser.add_argument("--window", "-w", help=\
        "The size (in bp) of windows used to compute BED file",
        type=int,
        required=False)
    parser.add_argument("--names", "-n", help=\
        "The names of all hybrid individuals, to be used when creating filenames "
        "(optional; default = use input filename)",
        required=False,
        nargs="+")
        
    parser.add_argument("--outgroup", "-og", help=\
        "If you want to estimate admixture proportions using f_hat, you must provide "
        "an outgroup sequence in the same format as all other sequences.",
        required=False)
    parser.add_argument("--f", "-f", help=\
        "The percent ancestry each query individual derives from ancestral population 1, "
        "if available. If you do not provide this, it will be calculated. You must provide "
        "a value for each hybrid individual; if you do not know the admixture proportion "
        "for one or more individuals, use 0 for those.",
        type=float,
        nargs="+",
        required=False)
    
    parser.add_argument("--pi", "-p", help=\
        "Three per-site nucleotide diversity values: for ancestral population 1, "
        "ancestral population 2, and between the two ancestral populations. If not "
        "provided, it will be calculated.",
        type=float,
        required=False,
        nargs=3)
    
    parser.add_argument("--skip", "-s", help=\
        "The percent of a window that must be made up of N bases for it to be skipped",
        type=float,
        default=0.25,
        required=False)

    parser.add_argument("--debug", "-d", help=\
        "Specify to print out debug messages to stderr.",
        action="store_true",
        default=False)
    parser.add_argument("--x_seqs", "-x", help=\
        "The names of all sequences that are likely to be long to the X chromosome. A "
        "different model will be used to scan these sequences.",
        required=False,
        nargs="+")
    parser.add_argument("--sex", help=\
        "For every hybrid individual provided, indicate its sex (M or F). You must "
        "provide a sex for every hybrid individual, or leave this blank, in which case "
        "it will be assumed that all individuals are female.",
        nargs="+",
        required=False)
    parser.add_argument("--resample_prob", "-rs", help=\
        "The probability of sampling the same ancestral recombination event twice "
        "in a single individual since the admixture event happened. If you do not "
        "have this value, it can either be generated by running allele_freq_drift.py "
        "or by this program. Since calculating it can be slow, watching stderr and "
        "providing the number calculated in subsequent runs can speed up performance.",
        type=float,
        required=False,
        default=None)
    parser.add_argument("--resample_prob_x", "-rsx", help=\
        "The probability of sampling the same ancestral recombination event twice "
        "in a single individual since the admixture event happened on the X chromosome. If you do not "
        "have this value, it can either be generated by running allele_freq_drift.py "
        "or by this program. Since calculating it can be slow, watching stderr and "
        "providing the number calculated in subsequent runs can speed up performance.",
        type=float,
        required=False,
        default=None)
    parser.add_argument("--gens", "-g", help=\
        "The estimated number of generations since admixture.",
        required=False,
        type=int)
    parser.add_argument("--pop_size", "-N", help=\
        "The estimated population size of the admixed population.",
        required=True,
        type=int)
    
    parser.add_argument("--print_scores", "-ps", help=\
        "Specify to print out scores in windows to stderr.",
        action="store_true",
        required=False,
        default=False)
    
    parser.add_argument("--reduce_het", "-rh", help=\
        "Build into the model transition probabilities the expected reduction in "
        "heterozygosity due to drift in the time since admixture. Leave out this "
        "option if you want to see if the model detects a reduction in heterozygosity "
        "on its own, without building in this assumption.",
        action="store_true",
        default=False,
        required=False)
    
    outparams = parser.add_mutually_exclusive_group(required=True)
    outparams.add_argument("--output_prefix", "-o", help=\
        "The prefix to use for all output files. Other filenames will be based on "
        "input filenames.",
        required=False)
    outparams.add_argument("--est_params", "-ep", help=\
        "Instead of running AD-LIBS, calculates parameters and ",
        action="store_true",
        required=False,
        default=False)
    
    parser.add_argument("--num_procs", "-np", help=\
        "The number of processes to use. Note that each process needs to read in on the "
        "order of the size of the largest genomic scaffold * number of individuals bytes "
        "into memory, so spawning too many processes can cause issues. The maximum number "
        "of processes allowed is the number of hybrid individuals to process (you can only "
        "use one process per individual).",
        type=int, 
        default=1,
        required=False)
        
    parsed = parser.parse_args()
    
    if parsed.sex is not None and len(parsed.sex) != len(parsed.hybrid):
        print("ERROR: you must either leave sex blank (assume female) or provide a value "
            "for every hybrid individual", file=sys.stderr)
        exit(1)
    elif parsed.f is not None and len(parsed.f) != len(parsed.hybrid):
        print("ERROR: you must either leave f blank (calculate admixture proportion for "
            "every individual using f-hat technique) or provide an estimate for every "
            "hybrid individual", file=sys.stderr)
        exit(1)
    if parsed.f is None:
        parsed.f = [0] * len(parsed.hybrid)
    
    # Check to see if the admixture proportion will need to be calculated
    # using f-hat. If so, we will require an outgroup sequence.
    calc_admixprop = False
    for index, value in enumerate(parsed.f):
        if value < 0 or value > 1:
            print("ERROR: invalid admixture proportion given for hybrid individual {}".format(index+1), \
                file=sys.stderr)
            exit(1)
        if value == 0:
            parsed.f[index] = None
            calc_admixprop = True
    if calc_admixprop and parsed.outgroup is None:
        print("ERROR: you have not provided an outgroup sequence. An outgroup is "
            "required if admixture proportion is to be calculated for any individuals. "
            "If you do not have an outgroup, please provide an estimate of admixture "
            "proportion for each hybrid individual.", file=sys.stderr)
        exit(1)
    if parsed.sex is not None:
        for sex in parsed.sex:
            if sex != "M" and sex != 'm' and sex != 'F' and sex != 'f':
                print("ERROR: invalid sex value {} provided. Valid values are M/m, F/f.".\
                    format(sex), file=sys.stderr)
                exit(1)
    else:
        parsed.sex = ['F'] * len(parsed.hybrid)
    if parsed.names is not None and len(parsed.names) != len(parsed.hybrid):
        print("ERROR: you have provided names for hybrid individuals that are not "
            "the same in number as the number of hybrid sequences.", file=sys.stderr)
        exit(1)
    elif parsed.names is None:
        parsed.names = []
        for index, hybrid in enumerate(parsed.hybrid):
            # Take "base" filename as the name of each individual.
            parsed.names.append(hybrid.split('/')[-1].split('.')[0])
        
    return parsed

def parse_config_file(filename):
    """
    Given the path to an appropriately-formatted configuration file, parses it and returns
    a list representing the new argument string. The argparse module will later be told
    to parse this string rather than the actual argument string.
    
    Arguments:
        filename -- the full path to an appropriately formatted configuration file
    Returns:
        a list representing arguments from the config file, which will be interpreted
            by the argparse module
    """
    dirname = os.path.dirname(filename)
    argstr = defaultdict(list)
    f = open(filename, 'r')
    has_newline = False
    for line in f:
        has_newline = line[-1] == "\n"
        line = line.rstrip()
        if len(line) == 0:
            continue
        if line[0] == "#":
            continue
        data = line.split()
        key = data[0].upper()
        if key == "A":
            argstr['-1'] = data[1:]
        elif key == "B":
            argstr['-2'] = data[1:]
        elif key == "PI":
            argstr['-p'] = data[1:]
        elif key == "W" or key == "WINDOW":
            argstr['-w'] = [data[1]]
        elif key == "RESAMPLE":
            argstr['-rs'] = [data[1]]
        elif key == "RESAMPLE_X":
            argstr['-rsx'] = [data[1]]
        elif key == "X":
            argstr['-x'] = [data[1]]
        elif key == "N" or key == "POPSIZE" or key == "POP_SIZE":
            argstr['-N'] = [data[1]]
        elif key == "G" or key == "GENERATIONS" or key == "GENS":
            argstr['-g'] = [data[1]]
        elif key == "OUTPUT":
            argstr['-o'] = [data[1]]
        elif key == "DEBUG":
            argstr['-d'] = []
        elif key == "SKIP":
            argstr['-s'] = [data[1]]
        elif key == "PRINT_SCORES":
            argstr['-ps'] = []
        elif key == "REDUCE_HET":
            argstr['-rh'] = []
        elif key == "OUTGROUP":
            argstr['-og'] = [data[1]]
        elif key == "PROCS" or key == "NUM_PROCS":
            argstr['-np'] = [data[1]]
        elif key == "HYBRID" or key == "H":
            # There can be multiple lines.
            # Each line MUST contain the name of the file, and can optionally
            # be followed by name, sex, and admixture proportion.
            argstr['-h'].append(data[1])
            if len(data) > 2:
                argstr['-n'].append(data[2])
            else:
                argstr['-n'].append(data[1].split('/')[-1].split('.')[0])
            if len(data) > 3:
                argstr['--sex'].append(data[3])
            else:
                # Default to female.
                argstr['--sex'].append('F')
            if len(data) > 4:
                argstr['-f'].append(data[4])
            else:
                # Ensure that admixture proportions have the same indices as
                # file names.
                argstr['-f'].append("0")
    f.close()
    
    # Fix file paths, if necessary.
    file_keys = ['-1', '-2', '-h', '-og']
    for file_key in file_keys:
        if file_key in argstr:
            for k, v in enumerate(argstr[file_key]):
                if (os.path.basename(v) == v):
                    # Path to file was not given. Assume files are in config
                    # file directory.
                    argstr[file_key][k] = dirname + '/' + v

    if not has_newline:
        f = open(filename, 'a')
        print("\n", sep="", end="", file=f)
        f.close()
    argstr_list = []
    for key in argstr:
        argstr_list.append(key)
        for item in argstr[key]:
            argstr_list.append(item)
    return argstr_list

def reconcile_args(argstr, otherargs):
    """
    Given a list of arguments parsed from a config file, and any other arguments 
    supplied via command line, add anything from the command line not already set via
    the config file to the list.
    
    Arguments:
        argstr -- the arguments obtained from the config file, as a list (returned by
            parse_config_file())
        otherargs -- the arguments actually provided to the program on the command line
            (as a list)
    
    Returns:
        A list of arguments containing both the items from the config file and those 
            provided via command line. We can then set sys.argv to this and have
            argparse handle it.
    """
    argstr_dict = defaultdict(list)
    item_counts = Counter()
    key = None
    for item in argstr:
        if item[0] == '-':
            key = item
            item_counts[key] = 0
        elif key is not None:
            argstr_dict[key].append(item)
            item_counts[key] += 1
    for key in item_counts:
        if item_counts[key] == 0:
            argstr_dict[key] = []
    
    key = None
    otherargs_dict = defaultdict(list)
    item_counts = Counter()
    for item in otherargs:
        if item[0] == '-':
            key = item
            item_counts[key] = 0
        elif key is not None:
            otherargs_dict[key].append(item)
            item_counts[key] += 1
    for key in item_counts:
        if item_counts[key] == 0:
            otherargs_dict[key] = []
    
    for key in otherargs_dict:
        # Try this:
        # If this was not in the config file, add it to the arguments.
        # If it was in the config file, replace it.
        argstr_dict[key] = otherargs_dict[key]
    
    # Convert back to list so sys.argv can handle it.
    newargs = []
    for key in argstr_dict:
        newargs.append(key)
        if len(argstr_dict[key]) > 0:
            for item in argstr_dict[key]:
                newargs.append(item)
    return newargs

def calc_pi(options):
    global ADLIBS_DIR
    
    if options.pi is None or len(options.pi) < 3:
        if options.debug:
            print("# Calculating nucleotide diversity values...", file=sys.stderr)
            
        cmd = ['{}/pi_pops'.format(ADLIBS_DIR)]
        for filename in options.ancestral1:
            cmd.append('-1')
            cmd.append(filename)
        for filename in options.ancestral2:
            cmd.append('-2')
            cmd.append(filename)
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        out, err = proc.communicate()

        pi1, pi2, pi_between = out.rstrip("\n\r\t\s").split("\t")
        pi1 = float(pi1)
        pi2 = float(pi2)
        pi_between = float(pi_between)
    
    else:
        pi1, pi2, pi_between = options.pi
        
    # Choose the ancestral population with lower nucleotide diversity as A, to ensure
    # that the model works.
    # Keep track of the fact that we did this so that, later, we can re-convert the state
    # names back to the ones the user specified.
    pops_flipped = False
    if pi1 > pi2:
        tmp = pi2
        pi2 = pi1
        pi1 = tmp
        tmp = options.ancestral2
        options.ancestral2 = options.ancestral1
        options.ancestral1 = tmp
        pops_flipped = True
        
    return (pi1, pi2, pi_between, pops_flipped)

def calc_f(options, pops_flipped):
    """
    Calculates admixture proportion. This is done via the c program calc_fhat and called
    via subprocess module (if necessary). If the populations have been flipped (A is B
    while the program runs and vice versa), then the given admixture proportion
    will be 1 - whatever is calculated.
    
    Arguments:
        options -- the parsed command-line (and config file) options, as returned by
            argparse
        pops_flipped -- True or False: have the ancestral populations been flipped for
            this run?
    
    Returns:
        A list, where indices are indices of hybrid individuals as provided to the program
            and values are floating point admixture proportions.
    """
    global ADLIBS_DIR
    # Now see if we need to calculate admixture proportions. If so, do that for
    # every individual and print them, either to the config file (if it exists), or
    # to stderr.
    
    # Determine which hybrid individuals need admixture proportion calculated.
    need_f = []

    for hybrid_index, admix_prop in enumerate(options.f):
        if admix_prop is None:
            need_f.append(hybrid_index)
            
    if len(need_f) > 0:
        # Make sure we have at least two individuals from at least one of the ancestral
        # populations.
        if len(options.ancestral1) < 2 and len(options.ancestral2) < 2:
            print("ERROR: you must provide at least two individuals from at least one "
                "of the ancestral populations to estimate admixture proportions, or else "
                "provide your own estimate for each individual (using -f at the command "
                "line or provide after the hybrid individuals' filenames on H/HYBRID lines "
                "in the config file.", file=sys.stderr)
            if configfile is not None:
                configfile.close()
            exit(1)
        
        f = options.f[:]
        
        for hybrid_index in need_f:
            hybrid = options.hybrid[hybrid_index]
            if options.debug:
                print("# Calculating admixture proportion of {}".format(hybrid), \
                    file=sys.stderr)
            a1cpy = options.ancestral1[:]
            if hybrid in a1cpy:
                a1cpy.remove(hybrid)
            a2cpy = options.ancestral2[:]
            if hybrid in a2cpy:
                a2cpy.remove(hybrid)
            if len(a1cpy) < 2:
                if len(a2cpy) < 2:
                    print("ERROR: not enough ancestral individuals to calculate admixture "
                        "proportion.", file=sys.stderr)
                    exit(1)
                a1file1 = a2cpy[0]
                a1file2 = a2cpy[1]
                a2file = a1cpy[0]
                f_flipped = True
            else:
                a1file1 = a1cpy[0]
                a1file2 = a1cpy[1]
                a2file = a2cpy[0]
                f_flipped = False
            cmd = ['{}/calc_fhat'.format(ADLIBS_DIR), \
                '--pop1a', a1file1, \
                '--pop1b', a1file2, \
                '--pop2', a2file, \
                '--hybrid', hybrid, \
                '--outgroup', options.outgroup]

            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            out, err = proc.communicate()
            p = float(out.rstrip("\n\r\t\s"))
            
            # If we flipped the ancestral populations for the sake of being able to calculate
            # admixture proportion, flip it back (this is the value of f we will use for
            # real stuff).
            if f_flipped:
                p = 1-p
            
            # Ensure p is positive and non-zero.
            if p <= 0:
                p = 0.001
                
            # If we flipped ancestral populations from how the user specified them, we 
            # keep the value of p we calculated, but print the opposite to stderr/config
            # file for the sake of the user/future runs.
            if pops_flipped:
                p_print = 1-p
            else:
                p_print = p
            
            sex = 'F'
            if options.sex is not None:
                sex = options.sex[hybrid_index]
            f[hybrid_index] = p      
    else:
        f = options.f[:]
        for p_index, p in enumerate(f):
            # Ensure p is positive and nonzero.
            if p <= 0:
                p = 0.001
            if pops_flipped:
                p = 1-p
            f[p_index] = p
    return f

def calc_resample_prob(options, r):
    """
    Calculates z, the probability of resampling the same ancestral recombination
    event twice in a given individual. This is hard and is done by an external program,
    which uses a diffusion approximation suggested by McKane and Waxman. 
    
    Arguments:
        options -- the command line (and config file) options, as supplied by argparse
        r -- the per site, per generation recombination probability (float)
    Returns:
        a tuple where the first item is the probability of resampling the same ancestral
            recombination event twice in an individual, per site, per generation (z), 
            and the second item is the same value on the X chromosome.
    """
    # Determine probability of resampling recombination events twice in the same
    # individual.
    if options.resample_prob is not None:
        resample_prob = options.resample_prob
    else:
        if options.debug:
            print("# Calculating resample probability", file=sys.stderr)
        resample_prob = recomb_resample_prob(r, options.gens, \
            options.pop_size, nthreads=4)
    if options.x_seqs is not None and len(options.x_seqs) > 0:
        if options.resample_prob_x is not None:
            resample_prob_x = options.resample_prob_x
        else:
            if options.debug:
                print("# Calculating resample probability on X chromosome", file=sys.stderr)
            resample_prob_x = recomb_resample_prob(r * (2/3), options.gens, \
                options.pop_size * (3/4), nthreads=4)
    else:
        resample_prob_x = -1
    return (resample_prob, resample_prob_x)

def calc_winsize(options, pi1, pi2, pi_between, f, r, resample_prob, \
    resample_prob_x):
    """
    Determines an optimal window size given other parameters, if the user opted to have
    AD-LIBS choose one.
    
    Arguments:
        options -- command line (and config file) arguments, as returned by argparse
        pi1 -- (float) nucleotide diversity in ancestral population A
        pi2 -- (float) nucleotide diversity in ancestral population B
        pi_between -- (float) nucleotide diversity when one haplotype is sampled from
            population A and another from population B
        f -- (float) estimated percent ancestry hybrid individuals derive from population A
        r -- (float) per site, per generation recombination probability
        resample_prob -- (float) probability of resampling the same ancestral recombination
            event twice in an individual after the set number of generations since admixture
            (referred to as z in the paper)
        resample_prob_x -- (float) same quantity as above (z), but calculated on X 
            chromosome instead of autosomes.
    
    Returns:
        (int) estimated optimal window size, rounded to the nearest 1000 bp
    """
    global overlap_thresh
    global prob_lim
    
    if options.reduce_het:
        # Calculate reduction in heterozygosity due to drift.
        # This number is different on the X chromosome due to its smaller
        # effective population size (3/4 that of the autosomes).
        h_t = math.exp(-options.gens/(2*options.pop_size))
        h_t_x = math.exp(-options.gens/(2*options.pop_size*0.75))
    else:
        h_t = 1
        h_t_x = 1
            
    # Determine optimal window size.
    # Determine maximum and minimum bounds on window size.
    pi_min = min(pi1, pi2, pi_between)
    pi_max = max(pi1, pi2, pi_between)

    # Minimum window size: depends on standard deviations of emission probability
    # distributions
    if options.x_seqs is not None and len(options.x_seqs) > 0:
        # NOTE: use 0.75 * pi instead of pi because on X-chromosome, pi is only
        # 3/4 the autosomal value.
        win_min = int(math.ceil((1/(0.75 * pi_min) + 1)/1000)*1000)
    else:
        win_min = int(math.ceil((1/pi_min + 1)/1000)*1000)
            
    # Determine maximum window based on transition probabilities.
    # NOTE: don't need to worry about X-chromosome here, since recombination rate
    # on X is only 2/3 the autosomal value; transition probabilities are therefore
    # lower.
    win_max = 999999
    # Test for every admixture proportion value.
    for p in f:
        probs = get_trans_probs(r, options.gens, p, resample_prob)
        prob_sum_max = max((probs['aa_ab'] + probs['aa_bb'], \
            probs['ab_aa'] + probs['ab_bb'], \
            probs['bb_aa'] + probs['bb_ab']))
        
        #win_lim_indv_nonround = int(math.floor(((h_t_x-2*prob_lim)/prob_sum_max)))

        win_lim_indv = int(math.floor(((h_t_x-2*prob_lim)/prob_sum_max)/1000)*1000)
        if win_lim_indv < win_max:
            win_max = win_lim_indv
            
            win_max_prev = win_max
            while win_max < win_min:
                # Iteratively reduce max prob threshold until we fix it.
                prob_lim *= 0.5
                win_max = int(math.floor(((h_t_x-2*prob_lim)/prob_sum_max)))
                if win_max < win_min and win_max_prev == win_max:
                    # Infinite loop.
                    print("ERROR: the parameters you have provided make choosing a window size impossible. "
                        "Please enter a different number of generations since admixture and try again.", file=sys.stderr)
                    exit(1)
                win_max_prev = win_max
    window = None
    
    if options.window is not None:
        # Ensure the given window size is within the bounds.
        window = options.window
        if window < win_min:
            print("# Given window size of {} is below minimum value of {}; using "
                "this value instead.".format(window, win_min), file=sys.stderr)
            window = win_min
        elif window > win_max:
            print("# Given window size of {} is greater than maximum value of {}; using "
                "this value instead.".format(window, win_max), file=sys.stderr)
            window = win_max
    else:
        # Choose an optimal window size based on distribution overlap.
        if options.debug:
            print("# Testing window sizes for distribution overlap...", file=sys.stderr)
            print("# size\toverlap", file=sys.stderr)
        for win_test in xrange(win_min, win_max, 1000):
            window = win_test
            # NOTE: debug mode will cause these distributions to be plotted 
            # at every step.
            overlap = get_dist_params(pi1, pi2, pi_between, win_test, x_chr=False, \
                overlap_only=True, debug=options.debug)
            if options.x_seqs is not None:
                overlap_x = get_dist_params(pi1, pi2, pi_between, win_test, x_chr=True, \
                    overlap_only=True, debug=options.debug)
                if overlap_x > overlap:
                    overlap = overlap_x
            if options.debug:
                print("# {}\t{}".format(win_test, overlap), file=sys.stderr)
            if overlap < overlap_thresh:
                break
    
    return window

def print_configfile(options, pi1, pi2, pi_between, pops_flipped, f, \
        resample_prob, resample_prob_x, window, configfile):
    """
    Writes out a new configuration file containing all the parameters used in this run,
    so another run can easily be done.
    
    Arguments:
        options -- command line (and config file) options, as returned by argparse
        pi1 -- (float) nucleotide diversity in ancestral population A
        pi2 -- (float) nucleotide diversity in ancestral population B
        pi_between -- (float) nucleotide diversity when one haplotype is sampled from
            population A and another from population B
        pops_flipped (boolean) -- whether or not ancestral population A and B have 
            been flipped during this run
        f -- (float) estimated percent ancestry hybrid individuals derive from population A
        resample_prob -- (float) probability of resampling the same ancestral recombination
            event twice in an individual after the set number of generations since admixture
            (referred to as z in the paper)
        resample_prob_x -- (float) same quantity as above (z), but calculated on X 
            chromosome instead of autosomes.
        window -- (int) the number of base pairs per window in this run
        configfile -- the file handle to write to (this should have already been opened
            and writeable)
    """
    print("# AD-LIBS config file created on {} at {}".format(\
        time.strftime("%Y/%m/%d"), time.strftime("%H:%M:%S")), file=configfile)
   
    if pops_flipped:
        # Flip back before printing
        anc1print = options.ancestral2
        anc2print = options.ancestral1
        f_print = f[:]
        for index, admixprop in enumerate(f):
            f_print[index] = 1-admixprop
        pi1print = pi2
        pi2print = pi1
    else:
        anc1print = options.ancestral1
        anc2print = options.ancestral2
        f_print = f[:]
        pi1print = pi1
        pi2print = pi2
        
    print("A\t{}".format(" ".join(anc1print)), file=configfile)
    print("B\t{}".format(" ".join(anc2print)), file=configfile)
    print("PI\t{}\t{}\t{}".format(pi1print, pi2print, pi_between), file=configfile)
    print("WINDOW\t{}".format(window), file=configfile)
    print("POPSIZE\t{}".format(options.pop_size), file=configfile)
    print("GENS\t{}".format(options.gens), file=configfile)
    print("RESAMPLE\t{}".format(resample_prob), file=configfile)
    if resample_prob_x is not None and resample_prob_x >= 0:
        print("RESAMPLE_X\t{}".format(resample_prob_x), file=configfile)
    if options.x_seqs is not None and len(options.x_seqs) > 0:
        print("X\t{}".format(" ".join(options.x_seqs)), file=configfile)
    if options.outgroup is not None:
        print("OUTGROUP\t{}".format(options.outgroup), file=configfile)
    print("SKIP\t{}".format(options.skip), file=configfile)
    print("PROCS\t{}".format(options.num_procs), file=configfile)
    if options.debug:
        print("DEBUG", file=configfile)
    if options.print_scores:
        print("PRINT_SCORES", file=configfile)
    if options.reduce_het:
        print("REDUCE_HET", file=configfile)
    if options.output_prefix is not None:
        print("OUTPUT\t{}".format(options.output_prefix), file=configfile)
    for index, hybrid in enumerate(options.hybrid):
        print("HYBRID\t{}\t{}\t{}\t{}".format(hybrid, options.names[index], \
            options.sex[index], f_print[index]), file=configfile)
    
def worker(options, index, pops_flipped, window, f, pi1, pi2, pi_between, resample_prob, \
    resample_prob_x, skip_score, outprefix):
    """
    Performs the actual work of an AD-LIBS run for a single individual. This makes
    it multi-process friendly: each hybrid individual gets its own process, so many 
    can be run at the same time.
    
    Arguments:
        options -- command-line (and config file) options, as returned by argparse
        index -- (int) the index of the hybrid individual for this AD-LIBS run. This is the
            key into the list of hybrid individuals stored in options.
        pops_flipped -- (boolean) have the ancestral populations A and B been flipped 
            for this run?
        window -- (int) the number of base pairs per window in this run
        f -- (float) estimated percent ancestry hybrid individuals derive from population A
        pi1 -- (float) nucleotide diversity in ancestral population A
        pi2 -- (float) nucleotide diversity in ancestral population B
        pi_between -- (float) nucleotide diversity when one haplotype is sampled from
            population A and another from population B
        resample_prob -- (float) probability of resampling the same ancestral recombination
            event twice in an individual after the set number of generations since admixture
            (referred to as z in the paper)
        resample_prob_x -- (float) same quantity as above (z), but calculated on X 
            chromosome instead of autosomes.
        skip_score -- (float) the numeric score that the window-scoring part of the program
            will use to represent skipped windows
        outprefix -- (string) the prefix (including directory) for output files
    """
    global ADLIBS_DIR
    
    # Determine base filename of individual
    if options.names is not None:
        fnbase = options.names[index]
    else:
        fnbase = options.hybrid[index].split('/')[-1].split('.')[0]
    
    # Ensure that the hybrid individual is not in one of the ancestral populations.
    # this makes it possible to specify one of the ancestral individuals as a hybrid
    # individual too, and to remove it from the pool of ancestral individuals when
    # processing.
    a1cpy = options.ancestral1[:]
    if options.hybrid[index] in a1cpy:
        a1cpy.remove(options.hybrid[index])
        if len(a1cpy) == 0:
            print("ERROR: not enough ancestral individuals left after removing hybrid "
                "individual from ancestral population.", file=sys.stderr)
            exit(1)
    a2cpy = options.ancestral2[:]
    if options.hybrid[index] in a2cpy:
        a2cpy.remove(options.hybrid[index])
        if len(a2cpy) == 0:
            print("ERROR: not enough ancestral individuals left after removing hybrid "
                "individual from ancestral population.", file=sys.stderr)
            exit(1)
    
    # Create a worker (full pipe for running AD-LIBS), given the index of a hybrid
    # individual. Wait for it to finish.
    cmd1 = ['{}/adlibs_score'.format(ADLIBS_DIR), '-z', str(skip_score), \
        '-w', str(window), '-s', str(options.skip)]
    for filename in a1cpy:
        cmd1.append('-1')
        cmd1.append(filename)
    for filename in a2cpy:
        cmd1.append('-2')
        cmd1.append(filename)
    cmd1.append('-h')
    cmd1.append(options.hybrid[index])
    
    cmd2 = ['{}/adlibs_hmm.py'.format(ADLIBS_DIR), \
        '-w', str(window), \
        '-p', str(pi1), str(pi2), str(pi_between), \
        '-rs', str(resample_prob), \
        '-rsx', str(resample_prob_x), \
        '-N', str(options.pop_size), \
        '-g', str(options.gens)]
    if options.sex is not None and (options.sex[index] == 'm' or options.sex[index] == "M"):
        cmd2.append('-m')
    if options.x_seqs is not None:
        cmd2.append('-x')
        for x_seq in options.x_seqs:
            cmd2.append(x_seq)
    if options.reduce_het is not None:
        cmd2.append('-rh')
    cmd2.append('-f')
    cmd2.append(str(f[index]))
    if pops_flipped:
        cmd2.append("-F")
        
    scoresfile = None
    
    outfilename = "{}{}.bed".format(outprefix, fnbase)
    outfile = open(outfilename, 'w')
    
    final_proc = None
    
    if options.print_scores:
        scoresfilename = "{}{}.scores".format(outprefix, fnbase)
        scoresfile = open(scoresfilename, 'w')
        p1 = subprocess.Popen(cmd1, stdout=scoresfile)
        subprocess.Popen.wait(p1)
        scoresfile.close()
        scoresfile = open(scoresfilename, 'r')        
        p2 = subprocess.Popen(cmd2, stdin=scoresfile, stdout=outfile)
        final_proc = p2
    else:
        p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=outfile)
        final_proc = p2
    
    subprocess.Popen.wait(final_proc)
    out, err = final_proc.communicate()
    outfile.close()
    if options.print_scores:
        scoresfile.close()
    
def main(args):
    """
    Main method
    
    Parses command line arguments and/or configuration file, sets up parameters, and
    uses the multiprocessing module to run AD-LIBS for each individual in its own 
    process. The work of this program is mainly done by external programs: 
    bin/adlibs_score takes FASTA files and outputs scores in windows, then bin/adlibs_hmm.py
    reads those scores into a hidden Markov model to do ancestry inference.
    """

    # Look for config file and attempt to parse it. If it's not there, parse
    # args normally.
    configfile_exists = False
    
    configfile_dir = None
    
    if len(args) > 1 and os.path.isfile(args[1]):
        configfile_exists = True
        argstr = parse_config_file(args[1])
        # Look for any arguments passed normally that weren't provided in the config file; 
        # Insert them into the arg string.
        if len(args) > 2:
            argstr = reconcile_args(argstr, args[2:])
        sys.argv = [args[0]] + argstr

    options = parse_args()
            
    # Calculate nucleotide diversity in ancestral populations.
    pi1, pi2, pi_between, pops_flipped = calc_pi(options)
    
    # Calculate admixture proportions for all hybrid individuals.
    f = calc_f(options, pops_flipped)

    # Set recombination rate
    r = 0.01/1000000
    
    # Calculate probability of resampling the same ancestral recombination event
    # twice in an individual.
    resample_prob, resample_prob_x = calc_resample_prob(options, r)
        
    # Calculate window size.
    window = calc_winsize(options, pi1, pi2, pi_between, f, r, resample_prob, \
        resample_prob_x)
    
        
    # Create a new config file that will contain all parameters, including
    # those that need to be calculated, for use in future runs.
    if options.est_params:
        configfile = sys.stderr
        outprefix = None
    else:
        outprefix = options.output_prefix
        if os.path.isdir(outprefix) and outprefix[-1] != '/':
            outprefix += '/'
        else:
            outprefix += '.'
        configfile = open("{}adlibs.config".format(outprefix), 'w')
    
    print_configfile(options, pi1, pi2, pi_between, pops_flipped, f, \
        resample_prob, resample_prob_x, window, configfile)
    
    if not options.est_params:
        configfile.close()
        
    if options.est_params:
        print("# All parameters calculated. You can find them either appended to the "
            "config file you provided, or printed to stdout. Exiting.", file=sys.stderr)
        exit(0)
    
    # We can now run the program.    
    
    # Define a numeric score to represent windows in which there was insufficient data
    # to calculate a score.
    skip_score = 999
        
    pool = Pool(processes=options.num_procs)
    for hybrid_index in range(0, len(options.hybrid)):
        #worker(options, hybrid_index, pops_flipped, window, f, pi1, pi2, pi_between, resample_prob, resample_prob_x, skip_score, outprefix)
        proc = pool.apply_async(worker, [options, hybrid_index, pops_flipped, window, \
            f, pi1, pi2, pi_between, resample_prob, \
            resample_prob_x, skip_score, outprefix])
    pool.close()
    pool.join()
    
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))

