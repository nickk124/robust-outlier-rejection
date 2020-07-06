.. rcr documentation master file, created by
   sphinx-quickstart on Fri Jun 19 09:51:42 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _index:

rcr
===

.. image:: https://img.shields.io/badge/GitHub-dfm%2Femcee-blue.svg?style=flat
    :target: https://github.com/nickk124/RCR
    :alt: Github Repository
.. image:: https://readthedocs.org/projects/rcr/badge/?version=latest
   :target: https://rcr.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://travis-ci.com/nickk124/RCR.svg?branch=master
    :target: https://travis-ci.com/nickk124/RCR
    :alt: Build Status
    
**RCR**, or Robust Chauvenet Rejection, is advanced, and easy to use, outlier rejection. Originally published 
in `Maples et al. 2018 <https://arxiv.org/abs/1807.05276>`_, this site will 
show you how to use RCR in Python. RCR can be applied to weighted data 
and used for model-fitting, and we have incorporated rejecting 
outliers in bulk to have the best of both computational efficiency and accuracy.

RCR has been carefully calibrated, and extensively simulated. It can be 
applied to samples with both large contaminants and large contaminant 
fractions (sometimes in excess of 90% contaminated). Finally, because RCR 
runs on a C++ backend, it is quite fast.

Introduction
------------
The simplest form of outlier rejection is sigma clipping, where measurements 
that are more than a specified number of standard deviations from the mean 
are rejected from the sample. This number of standard deviations should not 
be chosen arbitrarily, but is a function of your sample’s size. A simple 
prescription for this was introduced by William Chauvenet in 1863. 
Sigma clipping plus this prescription, applied iteratively, is what we call 
traditional Chauvenet rejection.

However, both sigma clipping and traditional Chauvenet rejection make use 
of non-robust quantities: the mean and the standard deviation are both 
sensitive to the very outliers that they are being used to reject. 
This limits such techniques to samples with small contaminants or small 
contamination fractions.

Robust Chauvenet Rejection (RCR) instead first makes use of robust 
replacements for the mean, such as the median and the half-sample mode, 
and similar robust replacements that we have developed for the standard 
deviation.

Basic Usage Example
-------------------
Here's a quick example of RCR in action: we have a dataset of :math:`N=1000` 
measurements, 85% of which are contaminants. The contaminants are sampled from
one side of a Gaussian/normal distribution with standard deviation :math:`\sigma=10`,
while the uncontaminated points are from a regular, symmetric Gaussian with :math:`\sigma=1`.
Both distributions are centered at :math:`\mu=0`.

The question is, how can we recover the :math:`\mu` and :math:`\sigma` of the underlying distribution, 
in the face of such heavy contamination? The example below shows how to do it with RCR.

.. literalinclude:: examples/intro/intro.py

Output:

.. image:: 
   _static/examples/intro/intro.*

For a more in-depth explanation of using RCR for this type of one-dimensional outlier rejection, see :ref:`singlevalue`.

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   guide/install
   guide/papers
   guide/faq
   modules


.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/singlevalue
   tutorials/functional