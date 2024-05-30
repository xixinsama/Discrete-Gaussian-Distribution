### Importing dependencies

# Estimators for moments
from scipy.stats import skew, kurtosis, moment
# Statistical (normality) tests
from scipy.stats import chisquare
# Numpy stuff
from numpy import cov, set_printoptions, diag, array, mean


# Math functions
from math import ceil, sqrt, exp, log
# Data management
from copy import deepcopy
from numpy import floor


# For debugging purposes
import sys
if sys.version_info >= (3, 4):
    from importlib import reload  # Python 3.4+ only.

# Tailcut rate
tau = 14
# Minimal size of a bucket for the chi-squared test (must be >= 5)
chi2_bucket = 10
# Minimal p-value
pmin = 0.001
# Print options
set_printoptions(precision=4)


def gaussian(x, mu, sigma):
    """
    Gaussian function of center mu and "standard deviation" sigma.
    """
    return exp(- ((x - mu) ** 2) / (2 * (sigma ** 2)))


def make_gaussian_pdt(mu, sigma):
    """
    Make the probability distribution table (PDT) of a discrete Gaussian.
    The output is a dictionary.
    """
    # The distribution is restricted to [-zmax, zmax).
    zmax = ceil(tau * sigma)
    pdt = dict()
    for z in range(int(floor(mu)) - zmax, int(ceil(mu)) + zmax):
        pdt[z] = gaussian(z, mu, sigma)
    gauss_sum = sum(pdt.values())
    for z in pdt:
        pdt[z] /= gauss_sum
    return pdt


class UnivariateSamples:
    """
    Class for computing statistics on univariate Gaussian samples.
    """

    def __init__(self, mu, sigma, list_samples):
        """
        Input:
        - the expected center mu of a discrete Gaussian over Z
        - the expected standard deviation sigma of a discrete Gaussian over Z
        - a list of samples defining an empiric distribution

        Output:
        - the means of the expected and empiric distributions
        - the standard deviations of the expected and empiric distributions
        - the skewness of the expected and empiric distributions
        - the kurtosis of the expected and empiric distributions
        - a chi-square test between the two distributions
        """
        zmax = ceil(tau * sigma)
        # Expected center standard variation.
        self.exp_mu = mu
        self.exp_sigma = sigma
        # Number of samples
        self.nsamples = len(list_samples)
        self.histogram = dict()
        self.outlier = 0
        # Initialize histogram
        for z in range(int(floor(mu)) - zmax, int(ceil(mu)) + zmax):
            self.histogram[z] = 0
        for z in list_samples:
            # Detect and count outliers (samples not in [-zmax, zmax))
            if z not in self.histogram:
                self.outlier += 1
            # Fill histogram according to the samples
            else:
                self.histogram[z] += 1
        # Empiric mean, variance, skewness, kurtosis and standard deviation
        self.mean = sum(list_samples) / self.nsamples
        self.variance = moment(list_samples, 2)
        self.skewness = skew(list_samples)
        self.kurtosis = kurtosis(list_samples)
        self.stdev = sqrt(self.variance)
        # Chi-square statistic and p-value
        self.chi2_stat, self.chi2_pvalue = self.chisquare()
        # Final assessment: the dataset is valid if:
        # - the chi-square p-value is higher than pmin
        # - there is no outlier
        self.is_valid = True
        self.is_valid &= (self.chi2_pvalue > pmin)
        self.is_valid &= (self.outlier == 0)


    def __repr__(self):
        """
        Print the sample statistics in a readable form.
        """
        rep = "\n"
        rep += "Testing a Gaussian sampler with center = {c} and sigma = {s}\n".format(c=self.exp_mu, s=self.exp_sigma)
        rep += "Number of samples: {nsamples}\n\n".format(nsamples=self.nsamples)
        rep += "Moments  |   Expected     Empiric\n"
        rep += "---------+----------------------\n"
        rep += "Mean:    |   {exp:.5f}      {emp:.5f}\n".format(exp=self.exp_mu, emp=self.mean)
        rep += "St. dev. |   {exp:.5f}      {emp:.5f}\n".format(exp=self.exp_sigma, emp=self.stdev)
        rep += "Skewness |   {exp:.5f}      {emp:.5f}\n".format(exp=0, emp=self.skewness)
        rep += "Kurtosis |   {exp:.5f}      {emp:.5f}\n".format(exp=0, emp=self.kurtosis)
        rep += "\n"
        rep += "Chi-2 statistic:   {stat}\n".format(stat=self.chi2_stat)
        rep += "Chi-2 p-value:     {pval}   (should be > {p})\n".format(pval=self.chi2_pvalue, p=pmin)
        rep += "\n"
        rep += "How many outliers? {o}".format(o=self.outlier)
        rep += "\n\n"
        rep += "Is the sample valid? {i}".format(i=self.is_valid)
        return rep

    def chisquare(self):
        """
        Run a chi-square test to compare the expected and empiric distributions
        """
        # We construct two histograms:
        # - the expected one (exp_histogram)
        # - the empirical one (histogram)
        histogram = deepcopy(self.histogram)
        # The chi-square test require buckets to have enough elements,
        # so we aggregate samples in the left and right tails in two buckets
        exp_histogram = make_gaussian_pdt(self.exp_mu, self.exp_sigma)
        obs = list(histogram.values())
        exp = list(exp_histogram.values())
        z = 0
        while(1):
            if (z >= len(exp) - 1):
                break
            while (z < len(exp) - 1) and (exp[z] < chi2_bucket / self.nsamples):
                obs[z + 1] += obs[z]
                exp[z + 1] += exp[z]
                obs.pop(z)
                exp.pop(z)
            z += 1
        obs[-2] += obs[-1]
        exp[-2] += exp[-1]
        obs.pop(-1)
        exp.pop(-1)
        exp = [round(prob * self.nsamples) for prob in exp]

        # 调整期望频率使其总和与观察频率的总和相等
        exp = [e * sum(obs) / sum(exp) for e in exp]

        diff = self.nsamples - sum(exp_histogram.values())
        exp_histogram[int(round(self.exp_mu))] += diff
        res = chisquare(obs, f_exp=exp)
        return res
    
