#! /usr/bin/env python3
import math
import numpy
from adlibs_dists import *
from pomegranate import *
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import scipy.stats as ss
import math
"""
adlibs_math.py

Math-related functions needed by different parts of AD-LIBS.
"""

# How much are distributions allowed to overlap each other?
overlap_thresh = 0.25

def exp_cdf(x, rate):
    """
    Cumulative distribution function for exponential distribution.
    
    Arguments:
        x -- (float) the value 
        rate -- (float) the exponential distribution parameter
    Returns:
        (float) the exponential CDF at x
    """
    return 1-math.exp(-rate*x)

def exp_quant(p, rate):
    """
    Quantile function for exponential distribution.
    
    Arguments:
        p -- (float) the value
        rate -- (float) the exponential distribution parameter
    Returns:
        (float) the quantile of the exponential distribution at x
    """
    if rate < 0 or rate == 1:
        return None
    return -math.log(1-p)/rate

def norm2lognorm(mu, sigma):
    """
    Converts normal distribution parameters to equivalent lognormal distribution
    parameters.
    
    Arguments:
        mu -- (float) the normal distribution mean
        sigma -- (float) the normal distribution standard deviation
    Returns:
        a tuple of (logmu, logsigma) where logmu is the lognormal equivalent of mu
            and logsigma is the lognormal equivalent of sigma
    """
    logmu = math.log(math.pow(mu,2)/math.sqrt(math.pow(sigma,2) + math.pow(mu,2)))
    logsigma = math.sqrt(math.log(math.pow(sigma,2)/math.pow(mu,2) + 1))
    return (logmu, logsigma)

def lognorm2norm(logmu, logsigma):
    """
    Converts lognormal distribution parameters to equivalent normal distribution
    parameters.
    
    Arguments:
        logmu -- (float) the lognormal distribution mean
        logsigma -- (float) the lognormal distribution sd
    Returns:
        a tuple of (mu, sigma) where mu is the normal equivalent of logmu and
            sigma is the normal equivalent of logsigma
    """
    mu = math.exp(logmu + math.pow(logsigma,2)/2)
    sigma = math.sqrt( math.exp( 2*logmu + math.pow(logsigma,2)*(math.exp(math.pow(logisgma,2)-1))))
    return (mu, sigma)

def norm_pdf(mu, sd, x):
    """
    PDF of the normal distribution.
    
    Arguments:
        mu -- (float) the normal distribution mean
        sd -- (float) the normal distribution standard deviation
        x -- (float) the value at which to evaluate.
    Returns:
        (float) the PDF of the normal distribution at x
    """
    return 1/(sd * math.sqrt(2*math.pi)) * math.exp(-(x-mu)**2/(2*(sd**2)))
    

def min_2_funs(fun1, args1, fun2, args2, x):
    """
    Suppose there are two functions, which each have their own parameters,
    but each is to be evaluated at a point x. Evaluates both functions at
    that point and returns the minimum of the two values.
    
    NOTE: each function must accept x as its last argument. If no additional
    arguments are required for the function, you must still provide an empty
    argument list.
    
    Arguments:
        fun1 -- the first function
        args1 -- the arguments to the first function (as a list)
        fun2 -- the second function
        args2 -- the arguments to the second function (as a list)
        x -- the last argument to be passed to both functions (should represent
            the point at which they're both being evaluated)
    Returns:
        the minimum of fun1 and fun2 evaluated at x
    """
    return min(( fun1(*(args1 + [x])), fun2(*(args2 + [x])) ))

def integrate_fun_trapezoid(fun, args, xmin, xmax):
    """
    Given a function whose last argument is the point at which to evaluate, uses
    the trapezoid rule to integrate that function between a minimum and maximum value.
    
    Last parameter of function being called must be x (point at which to 
    evaluate).
    
    Arguments:
        fun -- the function to integrate
        args -- the arguments to pass to the function (list)
        xmin -- (float) the minimum value for integration
        xmax -- (float) the maximum value for integration
    Returns:
        the approximate integral of fun from xmin to xmax.
    """
    # Use 100 steps.
    stepsize = (xmax-xmin)/100
    yprev = fun(*(args + [xmin]))
    x = xmin + stepsize
    area = 0
    while x <= xmax:
        y = fun(*(args + [x]))
        area += stepsize * min((y, yprev))
        area += 0.5 * stepsize * ( max((y, yprev)) - min((y, yprev)) )
        yprev = y
        x += stepsize
    return area
    
def get_dist_overlap(dist1, dist2, dist3):
    """
    Given three distributions (in Pomegranate module format),
    determines how much overlap there is between all three. This is done by
    defining a function that returns the minimum of two functions at
    any given point, then integrating that function from a minimum to maximum
    bound. This value is the percent overlap between two distributions.
    The overall percent overlap, then, is the sum of all three individual
    overlaps divided by three.
    
    NOTE: distributions should be provided in order of leftmost -> rightmost
    (first is BB, then AB, then AA).
    
    Arguments:
        dist1 -- a distribution object that has pdf and quantile functions
        dist2 -- a second distribution object with pdf and quantile functions
        dist3 -- a third distribution object with pdf and quantile functions
    
    Returns:
         a tuple of floats. The first is the percent overlap between dist1 and dist2, then
         the percent overlap between dist1 and dist3, then the percent overlap between
         dist2 and dist3.
    """
    
    # This came in handy here: http://stats.stackexchange.com/questions/12209/percentage-of-overlapping-regions-of-two-normal-distributions
    
    dist1_ov = integrate_fun_trapezoid(min_2_funs, \
        [dist1.pdf, [], dist2.pdf, []], \
        dist1.quantile(0.01), dist1.quantile(0.99))
    dist2_ov = integrate_fun_trapezoid(min_2_funs, \
        [dist2.pdf, [], dist3.pdf, []], \
        dist3.quantile(0.01), dist3.quantile(0.99))
    
    return((dist1_ov, dist2_ov))
    

def get_dist_params(pi_a1_orig, pi_a2_orig, pi_between_orig, window_size_orig, x_chr=False, debug=False, overlap_only=False):
    """
    Given pi within ancestral populations 1 and 2, and pi when one individual
    is sampled from each population, calculate mean and sd for all three 
    emission probability distributions.
    
    If the ancestral populations are similar to one another, it may be necessary
    to adjust the standard deviations of the distributions downward until the
    overlap between the distributions is minimized.
    
    Arguments:
        pi_a1_orig -- (float) the nucleotide diversity within ancestral population 1
        pi_a2_orig -- (float) the nucleotide diversity within ancestral population 2
        pi_between_orig -- (float) the nucleotide diversity when one haplotype is sampled
            from ancestral population 1 and another from ancestral population 2
        window_size_orig -- (int) the window size, in bp
        x_chr -- (boolean) -- does this chromosome/scaffold belong to a hemizygous sex
            chromosome?
        debug -- (boolean) -- should debug messages be print to the screen?
        overlap_only -- (boolean) -- instead of calculating the full parameters, do we
            just want to compute the overlap among emission probability distributions?
    Returns:
        a dict where keys are names of states (AA, AB, and BB) and values are dicts, 
            with keys mu and sd. Mu is keyed to the (float) mean of the distribution
            and sd is keyed to the (float) standard deviation of the distribution.
            
        Alternatively, if overlap_only is specified, this function returns the maximum
            overlap among any 2 of the 3 emission probability distributions (float).
    """
    window_size = window_size_orig
    
    if x_chr:
        # X chromosome effective population size is (3/4) the overall Ne
        pi_a1 = (3/4) * pi_a1_orig
        pi_a2 = (3/4) * pi_a2_orig
        pi_between = (3/4) * pi_between_orig
        #window_size = int(round((4/3) * window_size))
    else:
        pi_a1 = pi_a1_orig
        pi_a2 = pi_a2_orig
        pi_between = pi_between_orig
    
    ### Determine means of distributions
    
    het_pi_pos = 0.5*pi_a1 + 0.5*pi_between
    het_pi_neg = 0.5*pi_a2 + 0.5*pi_between
    
    het_ibs_mean = (1/window_size) * (1/het_pi_pos - 1/het_pi_neg)
    aa_ibs_mean = (1/window_size) * (1/pi_a1 - 1/pi_between)
    bb_ibs_mean = (1/window_size) * (1/pi_between - 1/pi_a2)
    
    ### Determine standard deviations of distributions
    
    pi1_sd = 1/(pi_a1 * math.sqrt(window_size * pi_a1-1))
    pi2_sd = 1/(pi_a2 * math.sqrt(window_size * pi_a2-1))
    pib_sd = 1/(pi_between * math.sqrt(window_size * pi_between-1))
    het_pos_sd = 1/(het_pi_pos * math.sqrt(window_size * het_pi_pos-1))
    het_neg_sd = 1/(het_pi_neg * math.sqrt(window_size * het_pi_neg-1))
    
    aa_sd = (1/window_size) * math.sqrt(math.pow(pi1_sd, 2) + math.pow(pib_sd, 2))
    bb_sd = (1/window_size) * math.sqrt(math.pow(pi2_sd, 2) + math.pow(pib_sd, 2))
    het_sd = (1/window_size) * math.sqrt(math.pow(het_pos_sd, 2) + math.pow(het_neg_sd, 2))
    #het_sd = (1/window_size) * math.sqrt((math.pow(het_pos_sd, 2) + math.pow(het_neg_sd, 2))/2)
    
    ### NOW figure everything out in log space.
    
    aa_ibs_mean_log = math.log(1/(pi_a1*window_size)) - math.log(1/(pi_between*window_size))
    bb_ibs_mean_log = math.log(1/(pi_between*window_size)) - math.log(1/(pi_a2*window_size))
    het_ibs_mean_log = math.log(1/(het_pi_pos*window_size)) - math.log(1/(het_pi_neg*window_size))
    
    mu, pi1_sd_log = norm2lognorm(1/(pi_a1*window_size), pi1_sd/window_size)
    mu, pi2_sd_log = norm2lognorm(1/(pi_a2*window_size), pi2_sd/window_size)
    mu, pib_sd_log = norm2lognorm(1/(pi_between*window_size), pib_sd/window_size)
    mu, het_pos_sd_log = norm2lognorm(1/(het_pi_pos*window_size), het_pos_sd/window_size)
    mu, het_neg_sd_log = norm2lognorm(1/(het_pi_neg*window_size), het_neg_sd/window_size)
    
    aa_sd_log = math.sqrt(math.pow(pi1_sd_log, 2) + math.pow(pib_sd_log, 2))
    bb_sd_log = math.sqrt(math.pow(pi2_sd_log, 2) + math.pow(pib_sd_log, 2))
    het_sd_log = math.sqrt(math.pow(het_pos_sd_log, 2) + math.pow(het_neg_sd_log, 2))
    
    # Here, we face the issue that the three distributions might overlap too much.
    # In the case that they overlap a lot, it's better to artificially change
    # the standard deviations to minimize overlap than to have them be too 
    # faithful to the true data.
    
    aa = ADLIBSNormalDistributionSkip(aa_ibs_mean_log, aa_sd_log, None)
    ab = ADLIBSNormalDistributionSkip(het_ibs_mean_log, het_sd_log, None)
    bb = ADLIBSNormalDistributionSkip(bb_ibs_mean_log, bb_sd_log, None)

    bb_ab_ov, ab_aa_ov = get_dist_overlap(bb, ab, aa)
    
    dist_overlap = max((bb_ab_ov, ab_aa_ov))
    
    if overlap_only:
        return dist_overlap
    
    # Determine the maximum amount of overlap between distributions to allow.
    global overlap_thresh
    
    while dist_overlap > overlap_thresh:
        if debug:
            print("# Reducing overlap from {}".format(dist_overlap), file=sys.stderr)
        
        # Reduce all standard deviations by half.
        bb_sd_log *= 0.5
        aa_sd_log *= 0.5
        het_sd_log *= 0.5
        
        aa = ADLIBSNormalDistributionSkip(aa_ibs_mean_log, aa_sd_log, float("Inf"))
        ab = ADLIBSNormalDistributionSkip(het_ibs_mean_log, het_sd_log, float("Inf"))
        bb = ADLIBSNormalDistributionSkip(bb_ibs_mean_log, bb_sd_log, float("Inf"))
        
        bb_ab_ov, ab_aa_ov = get_dist_overlap(bb, ab, aa)
        dist_overlap = max((bb_ab_ov, ab_aa_ov))
    
    if debug:
        print("# Overlap: {} {}".format(bb_ab_ov, ab_aa_ov))
    
    # Return a data structure containing means and standard deviations for
    # all three distributions.
    params = {'AA': {}, 'AB': {}, 'BB': {}}
    
    params['AA']['mu'] = aa_ibs_mean_log
    params['AA']['sd'] = aa_sd_log
    params['AB']['mu'] = het_ibs_mean_log
    params['AB']['sd'] = het_sd_log
    params['BB']['mu'] = bb_ibs_mean_log
    params['BB']['sd'] = bb_sd_log
    
    if debug:
        xmin = bb.quantile(0.001)
        xmax = aa.quantile(0.999)
        xpts = numpy.arange(xmin, xmax, 0.001)
        plt.xlabel('Score')
        plt.ylabel('Density')
        if not x_chr:
            plt.title("Emission probability distributions")
        else:
            plt.title("Emission probability distributions: X chromosome")
        
        aa_pdf = numpy.vectorize(aa.pdf, [numpy.float])
        ab_pdf = numpy.vectorize(ab.pdf, [numpy.float])
        bb_pdf = numpy.vectorize(bb.pdf, [numpy.float])
        
        plt.plot(xpts, bb_pdf(xpts), 'red')
        plt.plot(xpts, ab_pdf(xpts), 'purple')
        plt.plot(xpts, aa_pdf(xpts), 'blue')
        plt.show()
            
    return params
    

def get_trans_probs(r, g, p, resample_prob):
    """
    Computes transition probabilities between states.
    
    Arguments:
        r -- (float) per site, per generation recombination probability
        g -- (int) approximate number of generations since admixture (estimated
            beforehand)
        p -- (float) approximate percent ancestry the admixed population derives from
            ancestral population A (estimated beforehand)
        resample_prob -- (float) the probability of resampling the same ancestral
            recombination event twice in the same individual, in the generations since
            admixture. This is referred to as z in the paper.
    Returns:
        a dict where keys are strings of "fromstate_tostate" (i.e. aa_ab) and values
            are floating point approximate transition probabilities from fromstate to 
            tostate.
    """
    probs = {}
    
    probs['aa_ab'] = 2*g*r*(1-p)*(g*p*r - g*r + 1)
    probs['aa_bb'] = (1-p)*(-g**2 * p * r**2 + g**2 * r**2 + resample_prob)
    probs['ab_aa'] = 2*p*(g**2  * p * r**2 - g**2 * r**2 + g*r + resample_prob)
    probs['ab_bb'] = 2*(1-p)*(-g**2 * p * r**2 + g*r + resample_prob)
    probs['bb_aa'] = p*(g**2 * p * r**2 + resample_prob)
    probs['bb_ab'] = 2*g*p*r*(1 - g*p*r)
    
    return probs
    
