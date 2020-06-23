.. _singlevalue:

Rejecting 1D Outliers
=====================

This page gives a full tutorial for using RCR to detect and reject outliers
within one-dimensional datasets. Although this page avoids unnecessary statistical technicalities 
(see the paper **ADD LINK HERE**), a more bare-bones example is given on the main page :ref:`index`.

To begin, consider some dataset of :math:`N` measurements, made up of 1) samples from some contaminant 
distribution (outliers) and 2) samples from some underlying "true" uncontaminated distribution. 
RCR has various outlier rejection techniques (**add link to them here**) that have each been 
chosen to work best for different shapes of these distributions. The table below explains this.

.. _rejectiontechs:

Table of Rejection Techniques
-----------------------------

========================  ====================================  =============================
Best Rejection Technique  Uncontaminated ("true") Distribution  Contaminant Distribution
========================  ====================================  =============================
``SS_MEDIAN_DL``          Symmetric                             Two-Sided/Symmetric
``LS_MODE_68``            Symmetric                             One-Sided
``LS_MODE_DL``            Symmetric                             In-Between One- and Two-Sided
``ES_MODE_DL``            Mildly Asymmetric/Very low :math:`N`  (Any)                
========================  ====================================  =============================

*(Note that an uncontaminated distribution labeled as "symmetric" means approximately Gaussian/normal, 
mildly peaked, or mildly flat-topped, meaning an*
`exponential power distribution/generalized normal distribution <https://en.wikipedia.org/wiki/Generalized_normal_distribution>`_ 
*with positive and negative kurtosis, respectively.)*

For this tutorial, let's consider the case of both the uncontaminated and contaminated distributions being 
Gaussian/normal (so then, symmetric), both centered at :math:`\mu=0`. Being outliers, we'll give the contaminated 
distribution a standard deviation of :math:`\sigma=5`, and the uncontaminated distribution :math:`\sigma=1`. 
Referring to the table above, this means that the rejection technique that we'll need to
use is ``SS_MEDIAN_DL``, which we will show how to do shortly. Let's arbitrarily choose :math:`N = 500` datapoints total,
with a fraction of 50% contaminated. We can create the dataset in Python as follows:

.. code-block:: python

   import numpy as np

   np.random.seed(18318) # get consistent random results

   N = 500                  # total measurement count
   frac_contaminated = 0.5  # fraction of sample that will be contaminated

   # symmetric, uncontaminated distribution
   mu = 0 
   sigma_uncontaminated = 1
   uncontaminated_samples = np.random.normal(mu, sigma_uncontaminated, 
       int(N * (1 - frac_contaminated)))

   # symmetric, contaminated distribution
   sigma_contaminated = 5
   contaminated_samples = np.random.normal(mu, sigma_contaminated, 
        int(N * frac_contaminated))

   # combine to create overall dataset
   data = np.concatenate((uncontaminated_samples, contaminated_samples))
   np.random.shuffle(data)

To see what this dataset looks like, we'll plot it below (projected randomly along the :math:`y`-axis for added visibility).

.. code-block:: python

   plot data
   import matplotlib.pyplot as plt

   plt.figure(figsize=(8,5))
   ax = plt.subplot(111)

   ydata = np.random.uniform(0, 1, N) # project randomly into 2D for better visualization

   ax.plot(contaminated_samples, ydata[:int(N * frac_contaminated)], "k.", 
      label="Pre-RCR dataset", alpha=0.75, ms=4)
   ax.plot(uncontaminated_samples, ydata[int(N * frac_contaminated):], "k.", 
      alpha=0.75, ms=4)

   plt.xlim(-15, 15)
   plt.ylim(0, 1)
   plt.xlabel("data")
   plt.yticks([])

   box = ax.get_position()
   ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])
   ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

   plt.show()

Now, what do we do if we to estimate the :math:`\mu` and :math:`\sigma` of the underlying uncontaminated distribution?
Without RCR, we get:

.. code-block:: python

   # get results pre-RCR
   contaminated_mu = np.mean(data)
   contaminated_sigma = np.std(data)
   print(contaminated_mu, contaminated_sigma)

Output:

.. code-block:: python

    -0.3168378799621606 3.792535849537549

Unsurprisingly, the contaminants don't have a great effect on :math:`\mu`, as both the contaminants 
and the true distribution have the same :math:`\mu=0`. However, :math:`\sigma` is grossly
overestimated due to the contaminants, compared to the expected :math:`\sigma=1`.

So, how can we use RCR? After importing ``rcr`` (see :ref:`install`), we initialize the
``RCR`` object with the desired rejection technique; in our case ``SS_MEDIAN_DL``.
Next, we perform the outlier rejection (the, recommended, bulk rejection variant; see :ref:`bulk`)
using the `performBulkRejection()` method and the data (as well as optional weights for the data; see :ref:`weighting`), 
as follows:

.. code-block:: python

   # perform RCR
   import rcr

   # initialize RCR with rejection technique:
   # (chosen from shape of uncontaminated + contaminated distribution)
   r = rcr.RCR(rcr.SS_MEDIAN_DL)
   r.performBulkRejection(data) # perform outlier rejection

Next, we can obtain the results of RCR with the `result` member of ``RCR``. In our case, we're interested in the RCR-recovered
values for :math:`\mu` and :math:`\sigma` of the underlying uncontaminated distribution:

.. code-block:: python

   # View results post-RCR
   cleaned_mu = r.result.mu
   cleaned_sigma = r.result.stDev
   print(cleaned_mu, cleaned_sigma)

Output:

.. code-block:: python
   
   -0.1584668560834893 1.8260572902969874

Successfully, RCR managed to recover both a :math:`\mu` and :math:`\sigma` that are significantly 
closer to the true values of :math:`0` and :math:`1`, respectively, both by a factor of about 2.

We can also access the subsets of rejected and nonrejected datapoints of the dataset, as well as
the corresponding indices and flags thereof, from ``RCR.result``. For example, we can plot the
post-rejection dataset with:

.. code-block:: python

   # plot rejections
   cleaned_data = r.result.cleanY

   flags = r.result.flags 
   # list of booleans corresponding to the original dataset, 
   # true if the corresponding datapoint is not an outlier.

   cleaned_data_indices = r.result.indices 
   # indices of data in original dataset that are not outliers

   plt.figure(figsize=(8,5))
   ax = plt.subplot(111)
   ax.plot(data[cleaned_data_indices], ydata[cleaned_data_indices], "b.", 
        label="RCR-accepted points", alpha=0.75, ms=4)

   plt.xlim(-15, 15)
   plt.ylim(0, 1)
   plt.xlabel("data")
   plt.yticks([])

   box = ax.get_position()
   ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])
   ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

   plt.show()

Output:

.. image:: 
   ../_static/examples/singlevalue/singlevalue_postRCR.*

In the next section, we'll explore how we can apply weights to datapoints
to use with RCR.

.. _weighting:

Weighting Data
--------------

For both single-value/one-dimensional RCR, and the :math:`n`-dimensional
model-fitting/functional variant (see :ref:`functional`), numerical, non-negative weights can be
optionally provided for each of the datapoints. However, what does it really mean
to weight datapoints? If you have some datapoint :math:`y_n`, giving it a weight
of :math:`w_n=2` is simply analogous to counting it twice. Now, what's 
an example of where weighting can be useful?

Lets say that we'd like to perform RCR on the same dataset as above, except now
we somehow know that the true, uncontaminated datapoints should
be normally/Gaussian-distributed (again with :math:`\mu=0` and :math:`\sigma=1`) *a priori*.
We can use this prior knowledge to perform a sort of Bayesian outlier rejection,
by giving the datapoints weights that are proportional to the value of the
known normal probability density function. In Python, we can do this simply as:

.. code-block:: python

   from scipy.stats import norm

   # function to weight each datapoint according to the prior knowledge
   def weight_data(datapoint):
      return norm.pdf(datapoint, loc=mu, scale=sigma_uncontaminated)

   # create weights
   weights = weight_data(data)

Next we can perform RCR and view the results as usual, but now with providing the weights as the first argument
of ``performBulkRejection()``:

.. code-block:: python

   # perform RCR; same rejection technique
   r = rcr.RCR(rcr.SS_MEDIAN_DL)
   r.performBulkRejection(weights, data) # perform outlier rejection, now with weights

   # View results post-RCR
   cleaned_mu = r.result.mu
   cleaned_sigma = r.result.stDev
   print(cleaned_mu, cleaned_sigma)

Output:

.. code-block:: python

   -0.05519770432617514 0.7825197746126461

This is much closer to the expected values of :math:`\mu=0` and :math:`\sigma=1` than 
what we got with the unweighted/equally-weighted dataset above (this time actually,
:math:`\sigma` was slightly *under*-estimated).

We can then plot the cleaned dataset/non-rejected data as usual:

.. code-block:: python

   # plot rejections
   cleaned_data = r.result.cleanY
   cleaned_data_indices = r.result.indices

   plt.figure(figsize=(8,5))
   ax = plt.subplot(111)
   ax.plot(data[cleaned_data_indices], ydata[cleaned_data_indices], "b.", 
      label="RCR-accepted points,\nwith weights applied to data", alpha=0.75, ms=4)

   plt.xlim(-15, 15)
   plt.ylim(0, 1)
   plt.xlabel("data")
   plt.yticks([])

   box = ax.get_position()
   ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])
   ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

   plt.show()

Output:

.. image:: 
   ../_static/examples/singlevalue/singlevalue_postRCR_weight.*

As expected, the width of the cleaned dataset is noticeably smaller after applying weights.