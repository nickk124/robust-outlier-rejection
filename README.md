# What is RCR?
RCR is advanced, but easy to use, outlier rejection.

The simplest form of outlier rejection is sigma clipping, where measurements that are more than a specified number of standard deviations from the mean are rejected from the sample. This number of standard deviations should not be chosen arbitrarily, but is a function of your sample’s size. A simple prescription for this was introduced by William Chauvenet in 1863. Sigma clipping plus this prescription, applied iteratively, is what we call traditional Chauvenet rejection.

However, both sigma clipping and traditional Chauvenet rejection make use of non-robust quantities: the mean and the standard deviation are both sensitive to the very outliers that they are being used to reject. This limits such techniques to samples with small contaminants or small contamination fractions.

Robust Chauvenet Rejection (RCR) instead first makes use of robust replacements for the mean, such as the median and the half-sample mode, and similar robust replacements that we have developed for the standard deviation.

RCR has been carefully calibrated, and extensively simulated (see [Maples et al. 2020](https://arxiv.org/abs/1807.05276)). It can be applied to samples with both large contaminants and large contaminant fractions (sometimes in excess of 90% contaminated).

# How do I use RCR?
We have boiled it down to two simple user choices:

**1.** Are your uncontaminated measurements distributed symmetrically, like a Gaussian (or mildy peaked or flat-topped), mildly asymmetrically, or neither?

**2.** Are the contaminants to your measurements high and low in equal proportions, all high or all low, or something in between?

RCR can be applied to weighted data, to functional data (e.g., x vs. y), and we have incorporated bulk rejection to decrease computation times with large samples.

Besides the source code here, try our online calculators ([single value](https://skynet.unc.edu/rcr/calculator/value) and [functional](https://skynet.unc.edu/rcr/calculator/functional)).

# Licensing and Citation

RCR is free to use for academic and non-commercial applications (see license in this repository). We only ask that you cite [Maples et al. 2020](https://arxiv.org/abs/1807.05276).

For commercial applications, or consultation, feel free to contact us.


There is no more fundamental act in science than measurement. There is no more fundamental problem in science than contaminated measurements. RCR is not a complete solution...but it is very close! We hope that you enjoy it.

Dan Reichart, Michael Maples, Nick Konz

Department of Physics and Astronomy

University of North Carolina at Chapel Hill
