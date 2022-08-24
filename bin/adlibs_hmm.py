#! /usr/bin/env python3
import sys
import argparse
import subprocess
import random
from collections import Counter, defaultdict
import numpy
from adlibs_dists import *
from adlibs_math import *
import pomegranate
"""
adlibs_hmm.py

Part of AD-LIBS, a program for inferring ancestry of genomic windows of hybrid
individuals from low-coverage shotgun sequence data.

adlibs.py is the main program that runs everything; adlibs_score takes in sequences
and computes scores on windows. Output from adlibs_score should be piped to this
program, which builds and runs a hidden Markov model on these scores to infer
ancestry. Output of this program is in BED format: chromosome/scaffold name, start
position (0-indexed), end-position + 1 (0-indexed), ancestry state
where ancestry state is AA (homozygous population A), AB (heterozygous), or
BB (homozygous population B). Everything is printed to stdout.
"""

# The cap on any individual transition probability
prob_lim = 0.01

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
        
    parser.add_argument("--window", "-w", help=\
        "The size (in bp) of windows used to compute BED file",
        type=int,
        required=True)
        
    parser.add_argument("--pi", "-p", help=\
        "Three per-site nucleotide diversity values: for ancestral population 1, "
        "ancestral population 2, and between the two ancestral populations. If not "
        "provided, it will be calculated.",
        type=float,
        nargs=3,
        required=True)
        
    parser.add_argument("--f", "-f", help=\
        "The percent ancestry the query individual derives from ancestral population 1",
        type=float,
        required=True)
    parser.add_argument("--x_seqs", "-x", help=\
        "The names of all sequences that are likely to be long to the X chromosome. A "
        "different model will be used to scan these sequences.",
        required=False,
        nargs="+")
    parser.add_argument("--male", "-m", help=\
        "Use this to specify that the individual is a male. If so, a different, haploid "
        "model will be used on any scaffolds determined to be X chromomsome (since the X "
        "is haploid in males.",
        action="store_true",
        default=False)
    parser.add_argument("--resample_prob", "-rs", help=\
        "The probability of sampling the same ancestral recombination event twice "
        "in a single individual since the admixture event happened.",
        type=float,
        required=True,
        default=None)
    parser.add_argument("--resample_prob_x", "-rsx", help=\
        "The probability of sampling the same ancestral recombination event twice "
        "in a single individual since the admixture event happened on the X chromosome.",
        type=float,
        required=True,
        default=None)
    parser.add_argument("--gens", "-g", help=\
        "The estimated number of generations since admixture.",
        required=False,
        type=int)
    parser.add_argument("--pop_size", "-N", help=\
        "The estimated population size of the admixed population.",
        required=True,
        type=int)
    parser.add_argument("--flip_pops", "-F", help=\
        "Specify to flip populations A and B -- that is, state AA will become BB and \
        state BB will become AA in output.", \
        required=False,\
        action="store_true")
    
    parser.add_argument("--reduce_het", "-rh", help=\
        "Build into the model transition probabilities the expected reduction in "
        "heterozygosity due to drift in the time since admixture. Leave out this "
        "option if you want to see if the model detects a reduction in heterozygosity "
        "on its own, without building in this assumption.",
        action="store_true",
        default=False,
        required=False)
        
    parsed = parser.parse_args()
    
    return parsed
    
def get_model(r, params, window_size, num_skipped, seq_len, p, \
    g, resample_prob, x_chr=False, haploid=False, debug=False, h_t=1, skip_score=float("-Inf")):
    """
    Builds the hidden Markov model for a given chromosome or scaffold, using the
    Pomegranate module.
    
    Arguments:
        r -- (float) the per site, per generation recombination probability
        params -- a dict where keys are names of states (AA, AB, and BB) and values
            are dicts where values are mu and sd, which are floats representing
            means and standard deviations of emission probability distributions
        window_size -- (int) the window size for this run, in bp
        num_skipped -- (int) the number of windows that were skipped due to not passing
            criteria
        seq_len -- (int) the number of windows in the current chromosome/scaffold
        p -- (float) the percent ancestry the admixed population derives from ancestral
            population A (estimated beforehand)
        g -- (int) the number of generations since admixture (estimated beforehand)
        resample_prob -- (float) probability of resampling the same ancestral recombination
            event twice in an individual after the set number of generations since admixture
            (referred to as z in the paper)
        x_chr -- (boolean) does this chromosome/scaffold belong to a hemizygous sex
            chromosome?
        haploid -- (boolean) is this individual haploid along this chromosome/scaffold?
        debug -- (boolean) should debugging messages be printed to the screen?
        h_t -- (float) if the user has specified that expected reduction in heterozygosity
            given the number of generations since admixture should be incorporated into
            the model, this is the expected fraction of the initial heterozygosity that
            remains after g generations.
        skip_score -- (float) the number emitted by adlibs_score when "skipped" windows
            are encountered
    
    Returns:
        a Pomegranate HMM object for the current chromosome/scaffold
    """
    global prob_lim
    
    model = pomegranate.HiddenMarkovModel(name='ancestry')
    
    # Compute probabilities of transitioning to a skip state or the end. Cap these
    # both at the specified probability limit.
    skip_prob = num_skipped/seq_len
    if skip_prob > prob_lim:
        skip_prob = prob_lim
    state_end = 1/seq_len
    if state_end > prob_lim:
        state_end = prob_lim
    
    if x_chr:
        r *= (2/3)
        
    # Determine probabilities of transitions
    if haploid:
        # Should 2 be 1.5? I don't think so -- we already multiplied r by (2/3)
        # so that's in here already.
        aa_bb = g*r*(1-p)
        bb_aa = g*r*p
        # Eliminate the heterozygous state.
        aa_ab = 0
        ab_aa = 0
        bb_ab = 0
        ab_bb = 0
    else:
        probs = get_trans_probs(r, g, p, resample_prob)
        aa_ab = probs['aa_ab']
        ab_aa = probs['ab_aa']
        aa_bb = probs['aa_bb']
        bb_ab = probs['bb_ab']
        ab_bb = probs['ab_bb']
        bb_aa = probs['bb_aa']
    
    aa_ab *= window_size
    ab_aa *= window_size
    aa_bb *= window_size
    bb_ab *= window_size
    ab_bb *= window_size
    bb_aa *= window_size
    
    aa_aa = 1-(aa_ab + aa_bb + state_end + skip_prob)
    ab_ab = 1-(ab_aa + ab_bb + state_end + skip_prob)
    bb_bb = 1-(bb_aa + bb_ab + state_end + skip_prob)
    
    # Account for reduction in heterozygosity due to genetic drift
    
    if haploid:
        pass
        #aa_aa += (aa_bb - aa_bb*h_t)
        #aa_bb *= h_t
        #bb_bb += (bb_aa - bb_aa*h_t)
        #bb_aa *= h_t
    else:
        aa_aa += (aa_aa/(aa_aa+aa_bb)) * (aa_ab - aa_ab*h_t)
        aa_bb += (aa_bb/(aa_aa+aa_bb)) * (aa_ab - aa_ab*h_t)
        bb_aa += (bb_aa/(bb_aa+bb_bb)) * (bb_ab - bb_ab*h_t)
        bb_bb += (bb_bb/(bb_aa+bb_bb)) * (bb_ab - bb_ab*h_t)
        aa_ab *= h_t
        bb_ab *= h_t
        ab_aa += (ab_aa/(ab_aa+ab_bb)) * (ab_ab - ab_ab*h_t)
        ab_bb += (ab_bb/(ab_aa+ab_bb)) * (ab_ab - ab_ab*h_t)
        ab_ab *= h_t
        
    
    if debug:
        print("# AA -> AA {}".format(aa_aa), file=sys.stderr)
        print("# AA -> AB {}".format(aa_ab), file=sys.stderr)
        print("# AA -> BB {}".format(aa_bb), file=sys.stderr)
        print("# AB -> AA {}".format(ab_aa), file=sys.stderr)
        print("# AB -> AB {}".format(ab_ab), file=sys.stderr)
        print("# AB -> BB {}".format(ab_bb), file=sys.stderr)
        print("# BB -> AA {}".format(bb_aa), file=sys.stderr)
        print("# BB -> AB {}".format(bb_ab), file=sys.stderr)
        print("# BB -> BB {}".format(bb_bb), file=sys.stderr)
        print("# SKIP {}".format(skip_prob), file=sys.stderr)

    aaDist = pomegranate.NormalDistribution(params['AA']['mu'], params['AA']['sd'])
    abDist = pomegranate.NormalDistribution(params['AB']['mu'], params['AB']['sd'])
    bbDist = pomegranate.NormalDistribution(params['BB']['mu'], params['BB']['sd'])
                
    aaState = pomegranate.State(aaDist, name="AA")
    abState = pomegranate.State(abDist, name="AB")
    bbState = pomegranate.State(bbDist, name="BB")
    
    model.add_state(aaState)
    if not haploid:
        model.add_state(abState)
    model.add_state(bbState)
    
    #### ADD skip states
    
    skip_dist = pomegranate.UniformDistribution(skip_score-0.01, skip_score)
    
    aa_skip_state = pomegranate.State(skip_dist, name="skip-AA")
    ab_skip_state = pomegranate.State(skip_dist, name="skip-AB")
    bb_skip_state = pomegranate.State(skip_dist, name="skip-BB")
    
    model.add_state(aa_skip_state)
    if not haploid:
        model.add_state(ab_skip_state)
    model.add_state(bb_skip_state)

    if haploid:
        model.add_transition(model.start, aaState, p * (1-skip_prob))
        model.add_transition(model.start, aa_skip_state, p*skip_prob)
        model.add_transition(model.start, bbState, (1-p) * (1-skip_prob))
        model.add_transition(model.start, bb_skip_state, (1-p)*skip_prob)
    else:
        model.add_transition(model.start, aaState, p**2 * (1-skip_prob))
        model.add_transition(model.start, aa_skip_state, p**2 * skip_prob)
        model.add_transition(model.start, abState, 2*p*(1-p) * (1-skip_prob))
        model.add_transition(model.start, ab_skip_state, 2*p*(1-p)*skip_prob)
        model.add_transition(model.start, bbState, (1-p)**2 * (1-skip_prob))
        model.add_transition(model.start, bb_skip_state, (1-p)**2 * skip_prob)
    
    model.add_transition(aaState, model.end, 1/seq_len)
    if not haploid:
        model.add_transition(abState, model.end, 1/seq_len)
    model.add_transition(bbState, model.end, 1/seq_len)
    
    model.add_transition(aaState, bbState, aa_bb)
    model.add_transition(aaState, aaState, aa_aa)
    model.add_transition(bbState, aaState, bb_aa)
    model.add_transition(bbState, bbState, bb_bb)
    
    if not haploid:
        model.add_transition(aaState, abState, aa_ab)
        model.add_transition(abState, aaState, ab_aa)
        model.add_transition(abState, bbState, ab_bb)
        model.add_transition(abState, abState, ab_ab)
        model.add_transition(bbState, abState, bb_ab)
    
    ### Add skip state transitions
    model.add_transition(aaState, aa_skip_state, skip_prob)
    if not haploid:
        model.add_transition(abState, ab_skip_state, skip_prob)
    model.add_transition(bbState, bb_skip_state, skip_prob)
    
    model.add_transition(aa_skip_state, aa_skip_state, skip_prob)
    if not haploid:
        model.add_transition(ab_skip_state, ab_skip_state, skip_prob)
    model.add_transition(bb_skip_state, bb_skip_state, skip_prob)
    
    model.add_transition(aa_skip_state, bbState, aa_bb)
    model.add_transition(bb_skip_state, aaState, bb_aa)
    
    if not haploid:
        model.add_transition(aa_skip_state, abState, aa_ab)
        model.add_transition(ab_skip_state, aaState, ab_aa)
        model.add_transition(ab_skip_state, bbState, ab_bb)
        model.add_transition(bb_skip_state, abState, bb_ab)
    
    model.add_transition(aa_skip_state, model.end, 1/seq_len)
    if not haploid:
        model.add_transition(ab_skip_state, model.end, 1/seq_len)
    model.add_transition(bb_skip_state, model.end, 1/seq_len)
    
    model.add_transition(aa_skip_state, aaState, 1-skip_prob-aa_ab-aa_bb-1/seq_len)
    if not haploid:
        model.add_transition(ab_skip_state, abState, 1-skip_prob-ab_aa-ab_bb-1/seq_len)
    model.add_transition(bb_skip_state, bbState, 1-skip_prob-bb_aa-bb_ab-1/seq_len)
    ###    
    
    model.bake()
    
    return model

def compute_scores_file(score_file):
    """
    Given a file containing scores in windows (in BED format), reads the file and returns
    a numpy.array of scores.
    
    Arguments:
        score_file -- an open file handle for a file containing scores in BED format.
            Each line should be tab separated, of the form
            seqname start_index end_index score
    Returns:
        a numpy array where columns are seqname, start, end, score -- yielded unique for
            every chromosome/scaffold
    """
    scores = []
    prev_chr = None
    for line in score_file:
        line = line.rstrip()
        chrname, start, end, score = line.split('\t')
        start = int(start)
        end = int(end)
        score = float(score)
        if chrname != prev_chr and prev_chr is not None:
            yield numpy.array(scores)
            scores = []
            
        scores.append((chrname, start, end, score))
        prev_chr = chrname
        
    score_file.close()
    yield numpy.array(scores)

def path_to_bed(path, win_scores, flip_pops=False):
    """
    Given a Pomegranate HMM path returned by the viterbi method, prints the path
    in the form of seqname start end state in BED format, where state is the name of
    the ancestry state used in the model (AA, AB, or BB).
    
    Arguments:
        path -- the path the Pomegranate module returns from the viterbi() method
        win_scores -- a numpy.array of scores in windows, as returned by compute_scores_file()
    """
    for state_index, state_tup in enumerate(path):
        state_num, state = state_tup
        if state.name == "ancestry-start" or state.name == "ancestry-end":
            continue
        elif state.name[0:5] == "skip-":
            continue
        else:
            sname = state.name
            if flip_pops:
                if sname == "AA":
                    sname = "BB"
                elif sname == "BB":
                    sname = "AA"
            chrom, start_index, end_index = win_scores[state_index-1,0:3]
            print("{}\t{}\t{}\t{}".format(chrom, start_index, end_index, sname))
                
def main(args):
    """
    Parses command-line arguments, reads scores from stdin, builds a Pomegranate
    HMM object for every unique chromosome/scaffold, runs Pomegranate's Viterbi 
    algorithm to find the most likely path through states along that chromosome/scaffold,
    then prints the result in BED format to stdout.
    """
    
    options = parse_args()
    
    p = options.f
    # Ensure that p is positive and nonzero. If p were zero, there would be no
    # point in running this program.
    if p <= 0:
        p = 0.001
    
    # Score to represent skipped windows/missing data
    skip_score = float("inf")
    
    # Background recombination rate
    r = 0.01/1000000
    
    resample_prob = options.resample_prob
    resample_prob_x = options.resample_prob_x
    
    x_seqs = []

    if options.x_seqs is not None:
        x_seqs = options.x_seqs
    
    if options.reduce_het:
        # Calculate reduction in heterozygosity due to drift.
        # This number is different on the X chromosome due to its smaller
        # effective population size (3/4 that of the autosomes).
        h_t = math.exp(-options.gens/(2*options.pop_size))
        h_t_x = math.exp(-options.gens/(2*options.pop_size*0.75))
    else:
        h_t = 1
        h_t_x = 1
            
    probs = get_trans_probs(r, options.gens, p, resample_prob)
    
    # These will store emission probability distribution means and standard
    # deviations, once calculated (they depend on window size).
    params = None
    params_x = None
    
    pi1, pi2, pi_between = options.pi
    
    # Get distribution parameters.
    params = get_dist_params(pi1, pi2, pi_between, options.window, x_chr=False, debug=False)
    params_x = None
    if options.x_seqs is not None:
        params_x = get_dist_params(pi1, pi2, pi_between, options.window, x_chr=True, debug=False)
       
    for win_scores in compute_scores_file(sys.stdin):
        
        query_data = win_scores[:,3].astype('float')
        num_skipped = sum(query_data==skip_score)
        
        # Determine whether this scaffold is on the X chromosome.
        on_x_chr = win_scores[0][0] in x_seqs
        if on_x_chr:
            params_scaf = params_x
            rs = resample_prob_x
            h_t_scaff = h_t_x
            
        else:
            params_scaf = params
            rs = resample_prob
            h_t_scaff = h_t
                
        # Determine if haploid model should be used
        is_haploid = on_x_chr and options.male
        
        model = get_model(r, params_scaf, \
            options.window, num_skipped, len(query_data), p, \
            options.gens, rs, \
            x_chr=on_x_chr, haploid=is_haploid, debug=False, h_t=h_t_scaff, \
            skip_score=skip_score)
        
        likelihood, path = model.viterbi(query_data)
    
        path_to_bed(path, win_scores, flip_pops=options.flip_pops)
    

if __name__ == "__main__" :
    sys.exit(main(sys.argv))

