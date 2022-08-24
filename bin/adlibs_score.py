#! /usr/bin/env python3
import sys
import argparse
import random
from collections import Counter, defaultdict
import numpy
import math
import gzip
from Bio import SeqIO
"""
adlibs_score.py

This program reads FASTA sequences into memory, runs through them, computes scores
in windows, and prints to stdout. This is also done by adlibs_score (C program).
This is a backup in case that one doesn't work or compile properly. Output of this
program should be piped to adlibs_hmm.py. This should not need to be used.
"""
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
        required=True)
    parser.add_argument("--window", "-w", help=\
        "The size (in bp) of windows used to compute BED file",
        type=int,
        required=True)
    
    parser.add_argument("--skip", "-s", help=\
        "The percent of a window that must be made up of N bases for it to be skipped",
        type=float,
        default=0.25)

    parser.add_argument("--skip_score", "-z", help=\
        "The score to emit when a window is skipped (default 999).",
        type=float,
        default=999)
        
    return parser.parse_args()

def parse_fasta(filename):
    """
    Parse a FASTA (return a SeqIO obejct of it) given a filename.
    Handles gzipped files properly.
    """
    if filename[-3:] == '.gz':
        f = gzip.open(filename, 'r')
    else:
        f = open(filename, 'r')
    return SeqIO.parse(f, 'fasta')
    
def compute_scores(ancestral1files, ancestral2files, query, winsize, skip_threshold=0.25,
    skip_score=999):
    """
    Divides sequences into windows, calls compute_scores_window to compute the score
    on every window, then prints results to stdout (in BED format).
    """
    
    # Transform everything from filenames into SeqIO objects.
    ancestral1 = [None] * len(ancestral1files)
    ancestral2 = [None] * len(ancestral2files)
    
    for index, anc1file in enumerate(ancestral1files):
        ancestral1[index] = parse_fasta(anc1file)
    for index, anc2file in enumerate(ancestral2files):
        ancestral2[index] = parse_fasta(anc2file)

    query = parse_fasta(query)
    for qseq in query:
        shortest = len(qseq.seq)
        
        # If this sequence is too small to use, just skip it.
        if shortest < winsize:
            continue
            
        a1seqs = []
        a2seqs = []
        try:
            for a1indv in ancestral1:
                a1seq = a1indv.next()
                if a1seq.id != qseq.id:
                    print("ERROR: sequence IDs {} and {} do not match.".\
                        format(qseq.id, a1seq.id), file=sys.stderr)
                    exit(1)
                
                if len(a1seq.seq) > 0:
                    if len(a1seq.seq) < shortest:
                        shortest = len(a1seq.seq)
                    a1seqs.append(str(a1seq.seq).upper())
            
            for a2indv in ancestral2:
                a2seq = a2indv.next()
                if a2seq.id != qseq.id:
                    print("ERROR: sequence IDs {} and {} do not match.".\
                        format(qseq.id, a2seq.id), file=sys.stderr)
                    exit(1)
                
                if len(a2seq.seq) > 0:
                    if len(a2seq.seq) < shortest:
                        shortest = len(a2seq.seq)
                    a2seqs.append(str(a2seq.seq).upper())
                        
        except StopIteration:
            print("ERROR: sequence {} not found in one or more ancestral individual files.".\
                format(qseq.id), file=sys.stderr)
            exit(1)
        
        # Skip this sequence if ANY individual is missing it.
        if shortest < winsize:
            continue
            
        seqId = qseq.id
        qseq = str(qseq.seq).upper()
        
        win = winsize
            
        # Define windows.
        win_scores = []
        
        winStarts = range(0, shortest, win)
        
        for winStart in winStarts:
            winEnd = winStart + win
            if winEnd > shortest:
                winEnd = shortest
            qWin = qseq[winStart:winEnd]
            a1Win = []
            for a1seq in a1seqs:
                a1Win.append(a1seq[winStart:winEnd])
            a2Win = []
            for a2seq in a2seqs:
                a2Win.append(a2seq[winStart:winEnd])
            win_score = compute_scores_window(a1Win, a2Win, qWin, \
                skip_threshold=skip_threshold, skip_score=skip_score)
            
            print("{}\t{}\t{}\t{}".format(seqId, winStart, winEnd, win_score))
    
def compute_scores_window(a1seqs, a2seqs, qseq, skip_threshold=0.25, skip_score=999):
    """
    Compute and return a single score on a given window.
    """
    numA1 = len(a1seqs)
    numA2 = len(a2seqs)
    
    winLen = len(qseq)
    for a1seq in a1seqs:
        if len(a1seq) < winLen:
            winLen = len(a1seq)
    for a2seq in a2seqs:
        if len(a2seq) < winLen:
            winLen = len(a2seq)
            
    ibs_a1 = []
    for i in range(0, len(a1seqs)):
        ibs_a1.append([])
    ibs_a2 = []
    for i in range(0, len(a2seqs)):
        ibs_a2.append([])
    
    match_a1 = []
    for i in range(0, len(a1seqs)):
        match_a1.append(False)
    match_a2 = []
    for i in range(0, len(a2seqs)):
        match_a2.append(False)
    
    
    # Skip a window if the query sequence or at least the threshold % of at least one
    # ancestral population is over the threshold %N.
    n_count_a1 = Counter()
    n_count_a2 = Counter()
    n_count_query = 0
    
    for siteIndex in range(0, winLen):
        # Determine how useful site is for determining ancestry.
        qChr = qseq[siteIndex]
        if qChr == "N":
            n_count_query += 1
        
        for a1Index, a1seq in enumerate(a1seqs):
            # Determine if we should break the IBS tract here.
            if qChr == "N" or a1seq[siteIndex] == "N":
            
                # This is a tough decision to make. We can't always count an N as
                # a match or a mismatch, as that will mess things up.
                
                # We could always break here with a certain probability, but then
                # we would be biasing toward looking heterozygous.
                
                # Instead, we will look to see if any IBS tracts have already been found
                # with this individual (without seeing Ns). If they have, we will break
                # on an N with a probability determined by the length of previously-seen
                # IBS tracts.
                
                if match_a1[a1Index] != False and len(ibs_a1[a1Index]) > 0:
                    prev_ibs_sum = 0
                    prev_ibs_tot = 0
                    for prev_ibs in ibs_a1[a1Index]:
                        prev_ibs_sum += prev_ibs
                        prev_ibs_tot += 1
                    prev_ibs_avg = prev_ibs_sum/prev_ibs_tot
                    if random.random() < 1/prev_ibs_avg:
                        ibs_len = siteIndex-1-match_a1[a1Index]
                        if ibs_len > 1:
                            ibs_a1[a1Index].append(ibs_len)
                        match_a1[a1Index] = False
                if a1seq[siteIndex]:
                    n_count_a1[siteIndex] += 1
                    
            elif a1seq[siteIndex] == qChr:
                if match_a1[a1Index] == False:
                    match_a1[a1Index] = siteIndex
                # If we already have a haplotype going, do nothing until it ends.
            elif match_a1[a1Index] != False:
                # end the haplotype.
                ibs_len = siteIndex-1-match_a1[a1Index]
                if ibs_len > 1:
                    ibs_a1[a1Index].append(ibs_len)
                match_a1[a1Index] = False
            
        for a2Index, a2seq in enumerate(a2seqs):
            # Determine if we should break the IBS tract here.
            if qChr == "N" or a2seq[siteIndex] == "N":
                if match_a2[a2Index] != False and len(ibs_a2[a2Index]) > 0:
                    prev_ibs_sum = 0
                    prev_ibs_tot = 0
                    for prev_ibs in ibs_a2[a2Index]:
                        prev_ibs_sum += prev_ibs
                        prev_ibs_tot += 1
                    prev_ibs_avg = prev_ibs_sum/prev_ibs_tot
                    if random.random() < 1/prev_ibs_avg:
                        ibs_len = siteIndex-1-match_a2[a2Index]
                        if ibs_len > 1:
                            ibs_a2[a2Index].append(ibs_len)
                        match_a2[a2Index] = False
                if a2seq[siteIndex]:
                    n_count_a2[siteIndex] += 1
            elif a2seq[siteIndex] == qChr:
                if match_a2[a2Index] == False:
                    match_a2[a2Index] = siteIndex
                # If we already have a haplotype going, do nothing until it ends.
            elif match_a2[a2Index] != False:
                # end the haplotype.
                ibs_len = siteIndex-1-match_a2[a2Index]
                if ibs_len > 1:
                    ibs_a2[a2Index].append(ibs_len)
                match_a2[a2Index] = False
    
    # Determine whether to skip (if so, we can break early).
    n_thresh = skip_threshold * winLen
    if n_count_query > n_thresh:
        return skip_score
    else:
        a1_skipped = 0
        for a1_ncount in n_count_a1.values():
            if a1_ncount > n_thresh:
                a1_skipped += 1
        if a1_skipped/numA1 > skip_threshold:
            return skip_score
        a2_skipped = 0
        for a2_ncount in n_count_a2.values():
            if a2_ncount > n_thresh:
                a2_skipped += 1
        if a2_skipped/numA2 > skip_threshold:
            return skip_score
    
    # Add any final IBS tracts.
    # should we do this???
    #for a1Index, match_a1_indv in enumerate(match_a1):
    #    if match_a1_indv != False:
    #        ibs_a1[a1Index].append(winLen-match_a1_indv)
    #for a2Index, match_a2_indv in enumerate(match_a2):
    #    if match_a2_indv != False:
    #        ibs_a2[a2Index].append(winLen-match_a2_indv)
    
    # Compute average IBS tract len
    ibs_sum_a1 = 0
    ibs_tot_a1 = 0
    ibs_sum_a2 = 0
    ibs_tot_a2 = 0
    
    for indv in ibs_a1:
        for ibs_len in indv:
            ibs_sum_a1 += ibs_len/winLen
            ibs_tot_a1 += 1
    if ibs_tot_a1 == 0:
        ibs_tot_a1 = 1
    ibs_avg_a1 = ibs_sum_a1/ibs_tot_a1
    
    for indv in ibs_a2:
        for ibs_len in indv:
            ibs_sum_a2 += ibs_len/winLen
            ibs_tot_a2 += 1
    if ibs_tot_a2 == 0:
        ibs_tot_a2 = 1
    ibs_avg_a2 = ibs_sum_a2/ibs_tot_a2
    
    if ibs_avg_a1 == 0 or ibs_avg_a2 == 0:
        return skip_score
        
    return math.log(ibs_avg_a1) - math.log(ibs_avg_a2)

def main(args):
    """
    Main method.
    """
    options = parse_args()
    compute_scores(options.ancestral1, options.ancestral2, options.hybrid, \
        options.window, skip_threshold=options.skip, skip_score=options.skip_score)

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
