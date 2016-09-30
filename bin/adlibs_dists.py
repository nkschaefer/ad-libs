from __future__ import division, print_function
from collections import Counter
import math
import random
import pomegranate
import numpy
from scipy.special import erf
import scipy.stats
    
class ADLIBSLogNormalDistribution(pomegranate.Distribution):
    name="ADLIBSLogNormalDistribution"
    
    def __init__(self, mu, sigma, skip_score=None, mirror=False):
        # Convert mu and sigma into lognormal-friendly mu and sigma.
        self.init_from_normal(mu, sigma, mirror)
        
        # Allow a specially announced score to only be used for "skip" states.
        # If this score is ever seen, it's marked impossible.
        self.skip_score = None
        
    def init_from_normal(self, mu, sigma, mirror=False):
        
        # Determine if we need to shift the distribution left so it can handle negative
        # values.
        if mu - 3*sigma < 0:
            self.shift = 0-(mu-3*sigma)
            mu += self.shift
        else:
            self.shift = 0
        
        self.mu = mu
        self.sigma = sigma
        
        self.logmu = math.log(math.pow(mu,2)/math.sqrt(math.pow(sigma,2) + math.pow(mu,2)))
        self.logsigma = math.sqrt(math.log(math.pow(sigma,2)/math.pow(mu,2) + 1))
        
        # Calculate the mode (point of maximum value of the PDF).
        self.mode = math.exp(self.logmu - self.logsigma**2)
        self.median = math.exp(self.logmu)
        
        self.mirror = mirror
        
        # This is confusing as hell.
        # http://stackoverflow.com/questions/8870982/how-do-i-get-a-lognormal-distribution-in-python-with-mu-and-sigma
        
        self.dist = scipy.stats.lognorm(self.logsigma, scale=math.exp(self.logmu))
        
    def log_probability(self, symbol):
        """
        What's the probability of the given float under this distribution?
        """
        if symbol == self.skip_score:
            return float("-Inf")
        
        if self.mirror:
            x = 2*self.mode - symbol - self.shift
            #x = -(symbol - 2*self.mu) - self.shift
        else:
            x = symbol + self.shift
        
        cdf = self.dist.cdf(x)
        
        return log(cdf)
    
    def probability(self, symbol):
        return math.exp(self.log_probability(symbol))
    
    def pdf(self, x):
        """
        Probability density function.
        """
        if x == self.skip_score:
            return float("-Inf")
        
        if self.mirror:
            x = 2*self.mode - x - self.shift
        else:
            x = x + self.shift
        
        return self.dist.pdf(x)
    
    def cdf(self, x):
        """
        Cumulative distribution function.
        """
        if x == self.skip_score:
            return float('-Inf')
        if self.mirror:
            x = 2*self.mode - x- self.shift
        else:
            x = x + self.shift
        return self.dist.cdf(x)
         
    def quantile(self, quant):
        """
        Quantile function.
        """
        # Allow mirroring
        
        if self.mirror:
            # Need to get the value of 1-the given quantile first. This value
            # then needs to be shifted by the appropriate amount (2*mode - x).
            x = self.dist.ppf(1-quant)
            x = 2*self.mode - x
        else:
            x = self.dist.ppf(quant)
        
        x -= self.shift
        
        return x
        
        
    def sample(self):
        r = random.random()
        if self.mirror:
            # This will be mirrored but also left of 0. Need to shift up by
            # 2*mu, which will make it so the mean is the same as usual
            # but it's skewed left instead of right.
            draw = -self.dist.ppf(r) + 2*self.mode
            #draw = -self.dist.ppf(r) + 2*self.mu
        else:
            # Use percent point function (quantile function, inverse of CDF)
            draw = self.dist.ppf(r)
        
        if self.shift:
            # Move the random draw the right amount.
            draw -= self.shift
                
        return draw
    
    def from_sample(self, items, weights=None):
        # First calculate the mean and stuff normally. Then translate these into
        # lognormal parameters.
        items = numpy.array(items)
        mu = numpy.mean(items)
        sd = numpy.std(items)
        self.init_from_normal(mu, sd)

#ADLIBSLogNormalDistribution.register()

class ADLIBSFixedDistribution(pomegranate.Distribution):
    name="ADLIBSFixedDistribution"
    
    def __init__(self, score):
        self.score = score
    
    def log_probability(self, symbol):
        if symbol == self.score:
            return log(1)
        else:
            return float('-Inf')
    
    def sample(self):
        return self.score
    
    def from_sample(self, items, weights=None):
        # What to do here? I guess just set it equal to the first thing.
        self.score = items[0]

#ADLIBSFixedDistribution.register()

class ADLIBSNormalDistributionSkip(pomegranate.NormalDistribution):
    name="ADLIBSNormalDistributionSkip"
    
    def __init__(self, mean, std, skip_score):
        self.skip_score = skip_score
        self.mean = mean
        self.sd = std
        return super(ADLIBSNormalDistributionSkip, self).__init__(mean, std)
        
    def log_probability(self, symbol, epsilon=1E-9):
        if symbol == self.skip_score:
            return float('-Inf')
        else:
            return super(ADLIBSNormalDistributionSkip, self).log_probability(symbol)
    
    def probability(self, symbol, epsilon=1E-9):
        return math.exp(self.log_probability(symbol, epsilon))
        
    def quantile(self, quant):
        """
        Normal distribution quantile function.
        """
        return scipy.stats.norm.ppf(quant, self.mean, self.sd)
    
    def pdf(self, x):
        """
        Probability density function.
        """
        if x == self.skip_score:
            return float('-Inf')
        else:
            return scipy.stats.norm.pdf(x, self.mean, self.sd)
    
    def cdf(self, x):
        """
        Cumulative distribution function.
        """
        if x == self.skip_score:
            return float('-Inf')
        else:
            return scipy.stats.norm.cdf(x, self.mean, self.sd)
            
#ADLIBSNormalDistributionSkip.register()
            
class ADLIBSGumbelDistribution(pomegranate.Distribution):
    name="ADLIBSGumbelDistribution"
    
    def __init__(self, mu, sigma, skip_score=None, mirror=False):
        # mu and sigma are mean and sd of the data, parameters for a normal
        # distribution.
        
        # Here the mu provided is the mean of a normal distribution.
        # Mu in this distribution is actually the mode.
        self.mu = mu
        self.sigma = sigma
        self.beta = (sigma * math.sqrt(6))/math.pi
        # Convert from mean to mode.
        self.mu -= (self.beta * numpy.euler_gamma)
        self.mirror = mirror
        self.skip_score = skip_score
        
    def log_probability(self, symbol):
        """
        What's the probability of the given float under this distribution?
        """
        if symbol == self.skip_score:
            return float("-Inf")
        
        if self.mirror:
            # NOTE: mu is the mode here, so this is the mirror point.
            x = 2*self.mu - symbol
        else:
            x = symbol
        
        #cdf = math.exp(-math.exp(-(x-self.mu)/self.beta))
        #return log(cdf)
        try:
            return -math.exp((self.mu-x)/self.beta)
        except OverflowError:
            return float("-Inf")
            
    def probability(self, symbol):
        return math.exp(self.log_probability(symbol))
    
    def pdf(self, x):
        """
        Probability density function.
        """
        if x == self.skip_score:
            return float("-Inf")
        
        if self.mirror:
            # NOTE: mu is the mode here, so this is the mirror point.
            x = 2*self.mu - x
        
        z = (x-self.mu)/self.beta
        
        try:
            pdf = (1/self.beta) * math.exp(-(z + math.exp(-z)))
        except OverflowError:
            return float("-Inf")
            
        return pdf
    
    def cdf(self, x):
        """
        Cumulative distribution function.
        """
        if x == self.skip_score:
            return float('-Inf')
        
        if self.mirror:
            x = 2*self.mu - x
        
        try:
            cdf = math.exp(-math.exp(-(x - self.mu)/self.beta))
        except OverflowError:
            return float('-Inf')
        return cdf
        
    def ppf(self, quant):
        return self.mu - self.beta * math.log(-math.log(quant))
        
    def quantile(self, quant):
        """
        Quantile function.
        """
        # Allow mirroring
        
        if self.mirror:
            # Need to get the value of 1-the given quantile first. This value
            # then needs to be shifted by the appropriate amount (2*mode - x).
            x = self.ppf(1-quant)
            x = 2*self.mu - x
        else:
            x = self.ppf(quant)
        
        return x
        
    def sample(self):
        r = random.random()
        if self.mirror:
            # This will be mirrored but also left of 0. Need to shift up by
            # 2*mu, which will make it so the mean is the same as usual
            # but it's skewed left instead of right.
            draw = -self.ppf(r) + 2*self.mu
        else:
            # Use percent point function (quantile function, inverse of CDF)
            draw = self.ppf(r)
        
        print("{}\t{}".format(r, math.exp(self.log_probability(draw))))
        
        return draw
    
    def from_sample(self, items, weights=None):
        # First calculate the mean and stuff normally. Then translate these into
        # lognormal parameters.
        items = numpy.array(items)
        mu = numpy.mean(items)
        sigma = numpy.std(items)
        self.mu = mu
        self.beta = math.sqrt((6*math.sqrt(sigma))/(math.pi**2))

#ADLIBSGumbelDistribution.register()

