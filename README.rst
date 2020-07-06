Robust Chauvenet Outlier Rejection (RCR)
========================================
.. image:: https://readthedocs.org/projects/rcr/badge/?version=latest
   :target: https://rcr.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://travis-ci.com/nickk124/RCR.svg?branch=master
    :target: https://travis-ci.com/nickk124/RCR
    :alt: Build Status
    
.. image:: https://img.shields.io/badge/arXiv-1807.05276-orange.svg?style=flat
    :target: https://arxiv.org/abs/1807.05276
    :alt: arXiv Paper

What is RCR?
============
RCR is advanced, but easy to use, outlier rejection.

The simplest form of outlier rejection is sigma clipping, where measurements that are more than a specified number of standard deviations from the mean are rejected from the sample. This number of standard deviations should not be chosen arbitrarily, but is a function of your sample’s size. A simple prescription for this was introduced by William Chauvenet in 1863. Sigma clipping plus this prescription, applied iteratively, is what we call traditional Chauvenet rejection.

However, both sigma clipping and traditional Chauvenet rejection make use of non-robust quantities: the mean and the standard deviation are both sensitive to the very outliers that they are being used to reject. This limits such techniques to samples with small contaminants or small contamination fractions.

Robust Chauvenet Rejection (RCR) instead first makes use of robust replacements for the mean, such as the median and the half-sample mode, and similar robust replacements that we have developed for the standard deviation.

RCR has been carefully calibrated, and extensively simulated (see `Maples et al. 2018 <https://arxiv.org/abs/1807.05276>`_). It can be applied to samples with both large contaminants and large contaminant fractions (sometimes in excess of 90% contaminated).

Documentation
=============

The `documentation <rcr.readthedocs.io>`_ covers all of the RCR API, and provides thorough examples for using RCR in all of its forms.

Installation
============

Linux and macOS
---------------

RCR can be used most easily via Python, installed using ``python3 -m pip install rcr`` in the command line.

The C++ source code is also included here in ``/src``, with documentation in ``/docs/cpp_docs``.

Windows
-------

Before installing, you'll need to have **Microsoft Visual C++ 14.0**, found under the `Microsoft Visual C++ Build Tools <https://visualstudio.microsoft.com/downloads/>`_. If that doesn't work, you may need the latest `Windows SDK <https://developer.microsoft.com/en-us/windows/downloads/windows-10-sdk/>`_. (Both can be installed through the Visual Studio Installer.)

After that, run ``python3 -m pip install rcr`` in the command line.


Licensing and Citation
======================

RCR is free to use for academic and non-commercial applications (see license in this repository). We only ask that you cite `Maples et al. 2018 <https://arxiv.org/abs/1807.05276>`_ as:

.. code-block:: bib

   @article{maples2018robust,
       title={Robust Chauvenet Outlier Rejection},
       author={{Maples}, M.P. and {Reichart}, D.E. and {Konz}, N.C. and {Berger}, T.A. and {Trotter}, A.S. and {Martin}, J.R. and {Dutton}, D.A. and {Paggen}, M.L. and {Joyner}, R.E. and {Salemi}, C.P.},
       journal={The Astrophysical Journal Supplement Series},
       volume={238},
       number={1},
       pages={2},
       year={2018},
       publisher={IOP Publishing}
   }

For commercial applications, or consultation, feel free to contact us.

There is no more fundamental act in science than measurement. There is no more fundamental problem in science than contaminated measurements. RCR is not a complete solution...but it is very close! We hope that you enjoy it.

Nick Konz, Dan Reichart, Michael Maples

Department of Physics and Astronomy

University of North Carolina at Chapel Hill
