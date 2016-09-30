#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
from scipy.special import gegenbauer, hyp2f1, binom
import scipy
import numpy
import math
from multiprocessing import Pool
import argparse
"""
allele_freq_drift.py
Created by Nathan Schaefer on 10/26/15 at 10:41

This file contains functions used to calculate the probability of resampling
the same ancestral recombination event twice in an individual, given a 
population size, recombination probability per site, and a number of generations
spanned by a time period of interest.

This was motivated by the need to set transition probabilities in a hidden
Markov model used to infer ancestry of individuals. In order to fully account
for the probability of a given locus representing a transition from 
homozygous ancestry of one type to homozygous ancestry of the other, we need
to know the probability that the individual has inherited the same ancestral
recombination event on both chromosomes.

To solve this problem, recombination events can be thought of as alleles that
arise at frequency 1/(2N) in a diploid population. Since selection does not 
act on recombination events, we are interested in how drift affects their
frequency in a population. The Wright-Fisher model of drift allows one to build
a Markov chain with a transition probability matrix in which each matrix entry
represents the probability of transitioning from i to j copies of an allele
in one generation. The gth power of the matrix then contains, in each [i,j] 
position, the probability of transitioning from i to j copies of an allele over
the course of g generations.

This problem is computationally intractable for large N or g, however, and
continuous approximations exist. We used a solution worked out by McKane and
Waxman (2007) that expands upon the work of Kimura (1955), applying a diffusion
approximation to the genetic drift problem, but including solutions where
the allele frequency is 0 or 2N (loss and fixation). The article is here:
http://www.sciencedirect.com/science/article/pii/S0022519307001981

Furthermore, to make computation tractable, we allow multiprocessing and bin
allele frequencies, performing no more than 1000 calculations. This is useful
when population sizes are very large: there is no need to hold all possible
allele frequencies and their probabilities in memory.
"""
def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gens", "-g", help=\
        "The number of generations within the time span of interest.",
        type=int,
        required=True)
    parser.add_argument("--pop_size", "-N", help=\
        "The number of diploid individuals in the population.",
        type=int,
        required=True)
    parser.add_argument("--recomb", "-r", help=\
        "The probability of recombination (per site)",
        type=float,
        required=False,
        default=1e-8)
    parser.add_argument("--threads", "-t", help=\
        "The number of threads to use for calculations.",
        type=int,
        required=False,
        default=4)
    return parser.parse_args()
    
def phi_diffusion(x, t, N, p):
    """
    Uses the diffusion approximation to the genetic drift equation to 
    calculate the probability of all possible allele frequencies in a population
    of size N.
    
    Reference for full solution including loss and fixation:
    http://www.sciencedirect.com/science/article/pii/S0022519307001981
    Reference for original solution (excluding loss & fixation) by Kimura, 1955:
    http://www.pnas.org/content/41/3/144.full.pdf
    
    Parameters:
        x: the allele frequency for which to calculate probability (eg 5/2N)
        t: the number of generations for which drift has been allowed to
            take place (eg 10)
        N: the population size
        p: the starting allele frequency (should be 1/2N)
    
    Returns:
        a tuple of loss_coeff, prob, fix_coeff
            loss_coeff is the weight of the Dirac delta function at allele
            frequency = 0 for this time point
            prob is the probability density at the given value of x
            fix_coeff is the weight of the Dirac delta function at allele
            frequency = 1 for this time point
    """
    tau = t/(2*N)
    
    # Stopping threshold. When a number being added to a sum is less than
    # this value, the infinite sum will be considered "converged."
    delta = 1e-20
    
    loss_tot = 0
    stop = False
    n = 0
    
    # http://docs.scipy.org/doc/scipy/reference/special.html
    
    # Coefficient of dirac delta function for allele loss
    while not stop:
        lambdacoeff = ((n+1)*(n+2))/2
        newsum = p * (2*n+3) * hyp2f1(-n, n+3, 2, p) *\
            math.exp(-lambdacoeff * tau)
        if newsum < delta:
            stop = True
        else:
            loss_tot += newsum
        n += 1
    loss_coeff = 1 - loss_tot
    loss_coeff *= (1-p)
    
    fix_tot = 0
    stop = False
    n = 0
    
    # Coefficient of dirac delta function for allele fixation
    while not stop:
        lambdacoeff = ((n+1)*(n+2))/2
        newsum = (1-p) * (2*n+3) * math.pow(-1, n+1) *\
            hyp2f1(-n, n+3, 2, p) * math.exp(-lambdacoeff * tau)
        if newsum < delta:
            stop = True
        else:
            fix_tot += newsum
        n += 1
    fix_coeff = 1 - fix_tot
    fix_coeff *= p
    
    # Main part of function (equivalent to Kimura's 1955 solution)
    
    tot = 0
    stop = False
    n = 0
    
    while not stop:
        lambdacoeff = ((n+1)*(n+2))/2
        newsum = (2*n+3) * (n+1) * (n+2) *\
            hyp2f1(-n, n+3, 2, p) * hyp2f1(-n, n+3, 2, x) *\
            math.exp(-lambdacoeff * tau)
        if newsum < delta:
            stop = True
        else:
            tot += newsum
        n += 1
    tot *= p*(1-p)
    
    return (loss_coeff, fix_coeff, tot)

def phi_diffusion_t(t, N, nbins=None):
    """
    Calculates probabilities of all possible allele frequencies at a given
    time point. Uses the above function (phi_diffusion) based on the diffusion
    approximation to the Wright-Fisher genetic drift problem. Since numbers
    can be weird in some cases, scales all results to ensure probabilities
    sum to 1.
    
    Parameters:
        t: the number of generations that genetic drift is allowed to take place
        N: the population size.
    
    Returns:
        a one-dimensional Numpy array of probabilities of various allele 
            frequencies, ranging from 0 (loss) to 2*N (fixation).
    """
    # Figure out how to bin data. This is necessary to make this computationally
    # tractable.
    if nbins is None or nbins > 2*N-1:
        nbins = 2*N-1
    
    # Determine width of bins in terms of number of chromosomes in the population,
    # so we can accurately integrate the PDF.
    binwidth = 1/nbins
    
    # Add bins for loss and fixation.
    nbins += 2
    
    # Define data structure to store results
    freq_vector = numpy.array([0] * nbins, dtype=numpy.float64)
    
    # These variables will store the weights of the Dirac delta functions
    # corresponding to allele loss and fixation.
    loss_coeff = 0
    fix_coeff = 0
    
    # Do calculation for every possible (discrete) allele frequency not 
    # corresponding to loss or fixation
    for bin_index in xrange(1, nbins-1):
    
        # The above function (Diffusion approximation) is a probability 
        # density function, so its results can only be interpreted across
        # ranges of values. Therefore, we will evaluate the PDF 0.5/(2*N) units
        # to the left and to the right of the desired frequency, and then
        # simulate integration by calculating the area under the curve via
        # the trapezoid rule. 
        
        bin_lower = bin_index * binwidth - 0.5 * binwidth
        bin_upper = bin_index * binwidth + 0.5 * binwidth
        
        # bin_lower = (freq-0.5)/(2*N)
        # bin_upper = (freq+0.5)/(2*N)
        
        loss_coeff, fix_coeff, prob_dens_lower = \
            phi_diffusion(bin_lower, t, N, 1/(2*N))
        loss_coeff, fix_coeff, prob_dens_upper = \
            phi_diffusion(bin_upper, t, N, 1/(2*N))
        
        if prob_dens_lower < prob_dens_upper:
            smaller = prob_dens_lower
            larger = prob_dens_upper
        else:
            smaller = prob_dens_upper
            larger = prob_dens_lower
        
        # Use trapezoid rule to get area under the curve and hence probability.
        p = (smaller * binwidth) + 0.5 * binwidth * (larger-smaller)
        
        # Store result.
        freq_vector[bin_index] = p
    
    # Use Dirac delta function weights as probabilities of loss and fixation
    freq_vector[0] = loss_coeff
    freq_vector[nbins-1] = fix_coeff
    
    # Build a "legend" where each index corresponds to an index in freq_vector
    # and denotes a population allele frequency.
    freq_legend = [0] * nbins
    freq_legend[nbins-1] = 1
    for bin_index in xrange(1, nbins-1):
        freq_legend[bin_index] = binwidth * bin_index
        
    return (freq_vector, freq_legend)

def phi_matrix_t(t, N):
    """
    Used just for double-checking/verifying the accuracy of the diffusion 
    approximation probabilities. Sets up the genetic drift problem as a 
    Markov chain transition probability matrix and, like phi_diffusion_t()
    above, returns a one-dimensional Numpy array of probabilities of each
    allele frequency.
    
    This relies on creating a huge (2N x 2N) matrix and performing t rounds
    of matrix multiplication on it. This should therefore not be attempted 
    for large values of t or N but instead should be used to evaluate the
    output of the diffusion approximation functions above.
    
    Parameters:
        t: the number of generations for which genetic drift is being simulated
        N: the population size
    
    Returns: 
        a one-dimensional Numpy array of probabilities of various allele 
            frequencies, ranging from 0 (loss) to 2*N (fixation).
    """
    # Define the matrix. This will represent transition probabilities in a 
    # Markov chain simulating diffusion. Each i,j entry in the matrix is the
    # probability of transitioning from i to j copies of an allele in a
    # population in one generation. Taking this matrix to the t power will
    # result in each matrix entry's representing the probability of going 
    # from i to j copies of an allele over the course of t generations.
    
    m = []
    for row in xrange(0, 2*N+1):
        m.append([0] * (2*N+1))
    m = numpy.matrix(m, dtype=numpy.float64)
    
    # Set "absorbing" conditions. One lost or fixed, an allele must remain
    # lost or fixed.
    for j in xrange(1, 2*N+1):
        m[0,j] = 0
    for j in xrange(0, 2*N):
        m[2*N,j] = 0
    m[0,0] = 1
    m[2*N, 2*N] = 1
    
    # Set other transition probabilities. This is a standard equation based 
    # on the Wright-Fisher model of drift and the binomial distribution. It 
    # can be found in the popular Hartl & Clark population genetics textbook.
    for i in xrange(1, 2*N):
        for j in xrange(0, 2*N+1):
            m[i,j] = binom(2*N, j) * math.pow((i/(2*N)), j) *\
                math.pow(((2*N-i)/(2*N)), (2*N-j))
    
    # Exponentiate matrix. This simulates letting an allele diffuse over the
    # course of t generations.
    m = numpy.linalg.matrix_power(m, t)
    
    # Now get the probabilities of transitioning from 1 to any other number
    # of copies of an allele in the population.
    freq_vector = numpy.array([0] * (2*N+1), dtype=numpy.float64)
    
    for freq in xrange(0, 2*N+1):
        freq_vector[freq] = m[1,freq]
    
    return freq_vector

def recomb_resample_prob(r, g, N, nthreads=1, matrix=False):
    """
    Given a recombination probability per site, a number of generations,
    and a population size, calculates the probability of sampling the same
    recombination event twice in the same individual.
    
    g represents a number of generations in a time period of interest: a 
    recombination event is allowed to have arisen at any generation in that
    time period and then drifted to any frequency in the current generation.
    
    Parameters:
        r: the recombination probability per site
        g: the number of generations in the time period of interest
        N: the population size
        nthreads: the number of threads to use to speed up processing
        matrix: whether or not to use the matrix/Markov chain formulation
            of the genetic drift problem. This is computationally intractable
            for large numbers and should therefore only be used for testing/
            comparison purposes.
    
    Returns:
        a single probability: the chance of resampling the same ancestral
            recombination event twice in the same (diploid) individual
    """
    # From http://www.pnas.org/content/102/22/7882.full
    # "Because the frequency of each new mutation is initially 1/(2N), 
    # the distribution of allele frequency at modern sites is equal to 
    # the number of new mutations entering the population in a generation, 
    # multiplied by the transient distribution given the time difference 
    # between that generation and the current population, summed across 
    # generations."
    
    # Note: can't resample the same event twice in generation 0, since the
    # recombination event has just happened and is thus at frequency 1/2N.
    
    resample_prob = 0
    
    t_end = g
    
    # Allow multiple threads to do the work. Threads will calculate the 
    # probability of resampling twice (in the current generation) a
    # recombination event that arose in a given generation between t = 0 
    # and t = g-1.
    
    pool = Pool(nthreads)
    threads = []
        
    for t_start in xrange(0, g):
        t_diff = t_end - t_start
        if matrix:
            thread = pool.apply_async(recomb_resample_prob_aux_matrix, [t_diff, N, r])
        else:
            thread = pool.apply_async(recomb_resample_prob_aux, [t_diff, N, r])
        threads.append(thread)
    
    # Allow threads to work.
    pool.close()
    pool.join()
    for thread in threads:
        resample_prob_gen = thread.get()
        resample_prob += resample_prob_gen
        
    return resample_prob

def recomb_resample_prob_aux(t_diff, N, r):
    """
    Auxilliary helper function for above. Called by subprocess module to perform
    calculation of probability of resampling a recombination event that arose
    t_diff generations before the current generation twice in an individual in
    the current generation.
    
    Parameters:
        t_diff: the number of generations between when the recombination event
            arose and the current generation
        N: the population size
        r: recombination probability per site
        
    Returns:
        the probability of resampling a recombination event that arose t_diff
            generations ago twice in a single individual in the current
            generation.
    """
    # Set number of bins to to make this computationally tractable.
    freq_probs, freq_legend = phi_diffusion_t(t_diff, N, nbins=500)
    
    # Probability of resampling has three components:
    # Recombination event arises in generation t_start (r)
    # Recombination event drifts to frequency f in generation t_end
    # Recombination event is sampled twice in an individual in 
    #   generation t_end (f/(2N))^2
    resample_prob = 0
    for bin_index in xrange(0, len(freq_probs)):
        freq = freq_legend[bin_index]
        resample_prob += (r * freq_probs[bin_index] * math.pow(freq, 2))
    
    return resample_prob

def recomb_resample_prob_aux_matrix(t_diff, N, r):
    """
    Same as above function, but uses matrix/Markov chain formulation of 
    the genetic drift problem to get an exact answer. Note: this is intractable
    for large numbers of generations and large population sizes, so it should
    only be used for testing/comparison purposes.
    
    Parameters:
        t_diff: the number of generations between when the recombination event
            arose and the current generation
        N: the population size
        r: recombination probability per site
        
    Returns:
        the probability of resampling a recombination event that arose t_diff
            generations ago twice in a single individual in the current
            generation.
    """
    freq_probs = phi_matrix_t(t_diff, N)
    # Probability of resampling has three components:
    # Recombination event arises in generation t_start (r)
    # Recombination event drifts to frequency f in generation t_end
    # Recombination event is sampled twice in an individual in 
    #   generation t_end (f/(2N))^2
    resample_prob = 0
    for bin_index in xrange(0, len(freq_probs)):
        freq = bin_index/(2*N)
        resample_prob += (r * freq_probs[bin_index] * math.pow(freq, 2))
    
    return resample_prob
    
def main(args):
    """Main method"""
    options = parse_args()
    
    N = options.pop_size
    r = options.recomb
    g = options.gens
    threads = options.threads
    
    print("{:.2e}".format(recomb_resample_prob(r, g, N, nthreads=threads)))
    
    
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))

