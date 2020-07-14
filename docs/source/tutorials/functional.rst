.. _functional:

Rejecting Outliers While Model Fitting
======================================

Introduction
------------

In its most simple form, RCR is an excellent tool for
detecting and rejecting outliers within heavily contaminated one-dimensional 
datasets, as shown in :ref:`singlevalue`. However, this only scratches
the surface of RCR. In its more generalized form, RCR can also be used
to reject outliers within some :math:`n`-dimensional dataset 
while also *simultaneously* fitting a model to
that dataset. This section will explain how this can be done
fairly easily in practice, while avoiding going into unnessarily
technicalities. We recommend reading :ref:`singlevalue` before
tackling this section, as the following is essentially a generalization
of that section.

For the case of one-dimensional data (see :ref:`singlevalue`), 
RCR can be thought of as being used to reject outliers from some
dataset :math:`\left\{y_i\right\}_{i=1}^N`, distributed about
a single, parameterized "true" value :math:`y`. In this case, we often
wish to get a best estimate of :math:`y`, in order to properly characterize
the underlying measurement distribution; :math:`y` is just a one-dimensional 
*model* that we want to fit to the data, characterized by location
and scale parameters like :math:`\mu` and :math:`\sigma`.

If we generalize the dimensionality of this, we can imagine using 
RCR on measurements distributed about some :math:`n`-dimensional 
model *function* :math:`y(\vec{x}|\vec{\theta})`,
where :math:`\vec{x}` is an :math:`n`-dimensional vector of the model's 
independent variable(s), and :math:`\vec{\theta}` is an 
:math:`M`-dimensional vector of the model's parameters. In this case, 
we say that :math:`y(\vec{x}|\vec{\theta})` is an :math:`n`-dimensional,
:math:`M`-parameter model. For a more concrete example of this, consider
a simple linear model :math:`y = b + mx`. In this case,
:math:`n = 1`, and our parameter vector is just :math:`\vec{\theta} = (b ,m)`.

For this more-general case, what does our dataset look like? Each datapoint
will be some value of :math:`y(\vec{x}|\vec{\theta})` associated with
a value for :math:`\vec{x}`. As such, if we have :math:`N` datapoints in total,
indexing each by :math:`i`, our dataset can be written compactly 
as :math:`\left\{\left(\vec{x}_i, y_i\right)\right\}_{i=1}^N`
(be sure to avoid getting :math:`N` confused with :math:`n` here; 
the former is the number of datapoints that we're fitting the 
model to, while the latter is the dimensionality of the dataset/model).

Last but not least, before we get into the code, it's important to point out
that in order for RCR to fit any arbitrary model function to a dataset,
(partial) derivatives of the model function with respect to each model parameter must
be supplied (due to the specific algorithm that is used for fitting). For example,
consider a one-dimensional exponential model of the 
form :math:`y(\vec{x}|\vec{\theta})=be^{mx}`. If we choose to order the model
parameters as :math:`\vec{\theta} = (b ,m)`, then our model parameter derivatives are

.. math::
    \frac{\partial y(\vec{x}|\vec{\theta})}{\partial b} = e^{mx} \qquad \mathrm{ and } \qquad \frac{\partial y(\vec{x}|\vec{\theta})}{\partial m} = xbe^{mx} .


Implementation
--------------

Finally, we have everything that we need to use RCR for outlier rejection 
during model fitting. Although ``rcr`` supports any arbitrary
:math:`n`-dimensional nonlinear model function (as long as the model parameter 
derivatives are well-defined), for simplicity let's consider a simple linear model
:math:`y(\vec{x}|\vec{\theta})=b + mx`. The parameter partial derivatives
are then simply

.. math::
    \frac{\partial y(\vec{x}|\vec{\theta})}{\partial b} = 1 \qquad \mathrm{ and } \qquad \frac{\partial y(\vec{x}|\vec{\theta})}{\partial m} = x .

Before we start coding, it's important to consider the following:

.. note::

    Within ``rcr``, model functions and their derivatives must be defined
    exactly with arguments 1) ``x`` and
    2) ``params``, where ``x`` is the :math:`n`-dimensional list or numpy array 
    (or ``float``, in the case of :math:`n=1`) of independent variable(s), and 
    ``params`` is the :math:`M`-dimensional list/array of model parameters. 
    *Make sure to maintain consistent ordering of the model parameters vector
    throughout your code.*

Now, onto the code; let's start by defining our model function and its
parameter derivatives:

.. code-block:: python

    def linear(x, params): # model function
        return params[0] + x * params[1]

    def d_linear_1(x, params): # first model parameter derivative
        return 1

    def d_linear_2(x, params): # second model parameter derivative
        return x

Next, let's start creating our dataset. We'll have :math:`N=200` points
total, with 85% of the datapoints being outliers. Our "true" model that the datapoints
will be generated about will have parameters of :math:`b=0` and :math:`m=1`. 
In code, this is simply:

.. code-block:: python

    import numpy as np

    N = 200 # number of datapoints
    f = 0.85 # fraction of datapoints that are outliers

    params_true = [0, 1] # parameters of "true" model

We'll generate our datapoints in a certain range
of :math:`x` values about the "true" model line. For this example, we'll 
make uncontaminated datapoints that are Gaussian/normally distributed, with
standard deviation :math:`\sigma=1`, about the true model. In order to highlight
the power of RCR with dealing with especially difficult outliers, we'll generate
one-sided outliers/contaminants, sampled from the positive side of a 
Gaussian with :math:`\sigma=10`. In code, this will take the form of:

.. code-block:: python

    sigma_uncontaminated = 1 # standard deviations used to generate datapoints
    sigma_contaminated = 10

    # generate x-datapoints randomly in an interval
    x_range = (-10, 10)
    xdata_uncontaminated = np.random.uniform(
        x_range[0], x_range[1], int(N * (1 - f)))
    xdata_contaminated = np.random.uniform(
        x_range[0], x_range[1], int(N * f))


    # generate y-datapoints about the true model:
    # symmetric uncontaminated distribution
    ydata_uncontaminated = np.random.normal(
        loc=linear(xdata_uncontaminated, params_true),
        scale=sigma_uncontaminated
        )

    # one-sided contaminated distribution
    ydata_contaminated = linear(xdata_contaminated, params_true) + np.abs(
        np.random.normal(0, sigma_contaminated, int(N * f)))


    # combine dataset
    xdata = np.concatenate((xdata_contaminated, xdata_uncontaminated))
    ydata = np.concatenate((ydata_contaminated, ydata_uncontaminated))
   
Let's plot the dataset over the true, underlying model:

.. code-block:: python

    # plot dataset
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 5))
    ax = plt.subplot(111)

    ax.plot(xdata_contaminated, ydata_contaminated, "k.", 
        label="Pre-RCR dataset", alpha=0.75, ms=4)
    ax.plot(xdata_uncontaminated, ydata_uncontaminated, "k.", 
        alpha=0.75, ms=4)


    # plot model
    x_model = np.linspace(x_range[0], x_range[1], 1000)
    ax.plot(x_model, linear(x_model, params_true),
        "b--", label="True model", alpha=0.5, lw=2)

    plt.xlim(-10, 10)
    plt.ylim(-15, 25)
    plt.xlabel("$x$")
    plt.ylabel("$y$")

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

Output:

.. image:: 
   ../_static/examples/functional/preRCR.*

Clearly, these outliers are pretty nasty. This looks like a job for RCR. First,
we  need to supply an initial guess
for the model parameters, to give the fitting engine within RCR a starting place. 
Approaching this
dataset with no knowledge of what is or isn't an outlier, it would be hard
to tell what the true best fit should be; as such, we'll use an initial guess
that naively should work with the data, but is pretty far off of the true values of 
:math:`b=0` and :math:`m=1`; let's try :math:`b=5` and :math:`m=1.5`:

.. code-block:: python

    guess = [5, 1.5]

Next, we'll need to initialize the model, as an instance 
of the ``rcr.FunctionalForm`` class. The required arguments (in order)
to construct an instance of this class are 1) the model function,
2) the (:math:`n`-dimensional) :math:`x`-data, 3) the :math:`y`-data,
4) a list of the model parameter derivative functions, in order and 5)
the guess for the parameters. This is implemented as:

.. code-block:: python

    model = rcr.FunctionalForm(linear,
        xdata,
        ydata,
        [d_linear_1, d_linear_2],
        guess
    )

Now, we're finally ready to run RCR on the dataset/model.
Our uncontaminated distribution of data is
symmetric, while our contaminated distribution is
one-sided/completely asymmetric. Therefore, following
the :ref:`rejectiontechs`, the rejection technique
that will perform best on this dataset is ``LS_MODE_68``. Given this,
we'll perform RCR as usual, except now, we need to tell our instance of
the ``RCR`` class that we're fitting to our specific parametric model:

.. code-block:: python

    r = rcr.RCR(rcr.LS_MODE_68) # setting up for RCR with this rejection technique

    r.setParametricModel(model) 
    # tell RCR that we are model fitting, and give it the model of choice

    r.performBulkRejection(ydata) # perform RCR

That was only a few lines of code, but what actually happened here? Essentially,
(see :ref:`papers` for more details), RCR can iteratively reject outliers and
fit the model to the data at the same time. As such, we can access the same 
outlier-rejection results from ``r.result`` as in :ref:`singlevalue`, while also
having model-fitting results from our model, with the member ``model.result``:

.. code-block:: python

    best_fit_parameters = model.result.parameters # best fit parameters

    rejected_data = r.result.rejectedY # rejected and non-rejected data
    nonrejected_data = r.result.cleanY
    nonrejected_indices = r.result.indices 
    # indices from original dataset of nonrejected data

    print(best_fit_parameters)

Output:

.. code-block:: python

    [1.2367288755077883, 1.004037971689524]

Before we discuss this result, it's teaching to compare it to the
traditional method of ordinary least-squares fitting; we'll summarize 
this in a plot, as follows:

.. code-block:: python

    # plot results

    plt.figure(figsize=(8, 5))
    ax = plt.subplot(111)

    ax.plot(xdata_contaminated, ydata_contaminated, "k.", 
        label="Pre-RCR dataset", alpha=0.75, ms=4)
    ax.plot(xdata_uncontaminated, ydata_uncontaminated, "k.", 
        alpha=0.75, ms=4)

    ax.plot(xdata[nonrejected_indices], ydata[nonrejected_indices], "bo", 
        label="Post-RCR dataset", alpha=0.4, ms=4)

    # plot true model
    ax.plot(x_model, linear(x_model, params_true),
        "b--", label="True model", alpha=0.5, lw=2)

    # plot regular linear least squares best fit
    from scipy.stats import linregress

    slope_lsq, intercept_lsq, _, _, _ = linregress(xdata, ydata)

    ax.plot(x_model, linear(x_model, [intercept_lsq, slope_lsq]),
        "r-", label="Least-squares best fit", alpha=0.5, lw=2)

    # plot RCR-fitted model
    ax.plot(x_model, linear(x_model, best_fit_parameters),
        "g-", label="RCR best fit", alpha=0.5, lw=2)


    plt.xlim(-10, 10)
    plt.ylim(-15, 25)
    plt.xlabel("$x$")
    plt.ylabel("$y$")

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    print("Least-squares fit results:", intercept_lsq, slope_lsq)

    plt.show()

Output:

.. code-block:: python
   
   Least-squares fit results:
   7.202089027278407 1.0412871059773106

.. image:: 
   ../_static/examples/functional/postRCR.*

RCR gave us a best fit values of :math:`b=1.237` and :math:`m=1.004`, while
traditional linear least squares gave :math:`b=7.202` and :math:`m=1.041`.
The slope (true value of :math:`m=1`) was recovered very well in both cases, 
but this isn't super surprising, given that both the contaminated and uncontaminated
measurement distributions were generated without any scatter
along the :math:`x`-axis. However, due to the heavy scatter/contamination along the
:math:`y`-axis, the least-squares result for the intercept :math:`b` is, 
expectly, heavily biased by the outliers, very far off of the true value of 
:math:`b=1`. However, RCR was able to successfully reject many of the outliers,
while maintaining almost all of the uncontaminated distribution 
(shown in blue circles), giving a best fit :math:`b=1.237` 
that is significantly closer to the true value of :math:`b=1` than the least-squares
result. 

Overall, the RCR fit (green line) is clearly a much better fit 
(true best fit in blue dashed line) than the least squares best fit (red line).

.. _errorbars:

Data with Uncertainties and/or Weights
--------------------------------------

Realistically, many datasets will have measurements that have uncertainties, or *error bars*,
as practically all physical measurements cannot truly be made with exact precision.
In most model-fitting scenarios, only uncertainties in the *dependent* variable (:math:`y`) 
are considered, with any uncertainties in the independent variable(s) :math:`\vec{x}`
considered to be negligible (for a more generalized treatment, that includes such
:math:`\vec{x}`-uncertainties, as well as uncertainty in the dataset
that cannot solely be attributed to the data error bars, 
see e.g. `Konz 2020 <https://github.com/nickk124/seniorthesis/blob/master/konz_thesis_final.pdf>`_). In this case, which
we take for RCR, our dataset becomes 
:math:`\left\{\left(\vec{x}_i, y_i \pm \sigma_{y,i}\right)\right\}_{i=1}^N`, i.e.
our measurement error bars/uncertainties are :math:`\left\{\sigma_{y,i}\right\}_{i=1}^N`.

Just as in one-dimensional RCR, weights :math:`w_i`
can also be applied to model-fitting datasets (e.g. :ref:`weighting`). 
We note that the inclusion of error bars as described in the previous paragraph
is not mutual exclusive with such weighting; both weights and error bars can be used
in practice. 

To use a dataset with error bars and/or weights with model-fitting RCR, simply
use the optional arguments ``error_y`` and ``weights`` of the ``rcr.FunctionalForm()``
constructor, where the former is an ordered vector/list of measurement uncertainties 
:math:`\left\{\sigma_{y,i}\right\}_{i=1}^N`, and the latter is an ordered vector/list 
of measurement weights :math:`\left\{ w_i\right\}_{i=1}^N`. An example of this 
is given in the following section.

Model Parameter Uncertainties/Error Bars
----------------------------------------

In many cases, we often want not just best fit parameters for a model and dataset,
but also *uncertainties*, or "error bars" for these parameters. This is easily
available in ``rcr``, again via the ``model.result`` object, as 
``model.result.parameter_uncertainties``. However, before we go into a
worked code example, note the following:

.. note::

    In ``rcr``, best fit model parameter uncertainties can only be calculated if 
    error bars/uncertainties *and/or* weights were given for the dataset before fitting.

Now, let's try adding error bars to our linear dataset, same as above. First, 
we'll initialize the error bars, randomly, giving higher error, on average, to the contaminants:

.. code-block:: python

    error_y_uncontaminated = np.random.uniform(low=0.1, high=1, size=int(N * (1 - f)))
    error_y_contaminated = np.random.uniform(low=1, high=2, size=int(N * f))

    error_y = np.concatenate((error_y_contaminated, error_y_uncontaminated))

Next, let's initailize the model as before, except now using the optional keyword argument ``error_y``
to input the error bars. We then can perform RCR as usual.

.. code-block:: python

    # instantiate model
    model = rcr.FunctionalForm(linear,
        xdata,
        ydata,
        [d_linear_1, d_linear_2],
        guess,
        error_y=error_y
    )

    # initialize and perform RCR as usual
    r = rcr.RCR(rcr.LS_MODE_68) # setting up for RCR with this rejection technique
    r.setParametricModel(model) # tell RCR that we are model fitting
    r.performBulkRejection(ydata) # perform RCR

Let's check out the results:

.. code-block:: python

    # view results
    best_fit_parameters = model.result.parameters # best fit parameters
    best_fit_parameter_errors = model.result.parameter_uncertainties # and their uncertainties

    rejected_data = r.result.rejectedY # rejected and non-rejected data
    nonrejected_data = r.result.cleanY
    nonrejected_indices = r.result.indices

    print(best_fit_parameters)
    print(best_fit_parameter_errors)

Output:

.. code-block:: python

    [6.612942587028933, 0.9732622673909074]
    [1.6299290812536242, 0.3258511725157285]

So, our RCR-recovered best fit is :math:`b = 6.61 \pm 1.63` and :math:`m = 0.973 \pm 0.326`.
Unfortunately, this fit isn't nearly as good as when we didn't have measurement uncertainties.
But why? To see, let's plot the dataset alongside the fit:

.. code-block:: python

    # plot results

    plt.figure(figsize=(8, 5))
    ax = plt.subplot(111)

    ax.errorbar(xdata_contaminated, ydata_contaminated, yerr=error_y_contaminated, 
        fmt="k.", label="Pre-RCR dataset", alpha=0.75, ms=4)
    ax.errorbar(xdata_uncontaminated, ydata_uncontaminated, yerr=error_y_uncontaminated,
        fmt="k.", alpha=0.75, ms=4)

    ax.plot(xdata[nonrejected_indices], ydata[nonrejected_indices], "bo", 
        label="Post-RCR dataset", alpha=0.4, ms=4)

    # plot true model
    ax.plot(x_model, linear(x_model, params_true),
        "b--", label="True model", alpha=0.5, lw=2)

    # plot RCR-fitted model
    ax.plot(x_model, linear(x_model, best_fit_parameters),
        "g-", label="RCR best fit", alpha=0.5, lw=2)


    plt.xlim(-10, 10)
    plt.ylim(-15, 25)
    plt.xlabel("$x$")
    plt.ylabel("$y$")

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

Output:

.. image::
    ../_static/examples/functional/postRCR_erry.*

Adding error bars, or *intrinsic* uncertainties, to the measurements in
the dataset introduced even more overall uncertainty to the data, beyond just the *extrinsic* uncertainty,
or scatter/sample variance of the datapoints themselves. That, combined with the extremely high
contaminant fraction of 85%, made it so that RCR was unable to tell apart the contaminants from the non-outlier datapoints,
under-rejecting the outliers, as shown in the plot. As such, the final dataset that the
model was fit to included too many outliers, biasing the fitted line to have too high an intercept.
RCR would've worked better if either/both 1) there were smaller error bars
or 2) the fraction of contaminants was lower.

.. _priors:

Applying Prior Knowledge to Model Parameters (Advanced)
-------------------------------------------------------

Let's say that we want to fit some model to a dataset, and we know certain, *prior* information
about one of the parameters of the model, :math:`a`, in advance. 
From the point of view of `Bayesian inference <https://en.wikipedia.org/wiki/Bayes%27_theorem>`_,
this can be formalized by specifying the *prior probability distribution*, or *prior probability density function*
(PDF) of that parameter :math:`p(a)`. For example, let's say that for the linear dataset/model above, we know *a priori* 
that the intercept :math:`b` should be :math:`b=0`, with uncertainty of :math:`1`, 
i.e. :math:`b=0 \pm 1`. This translates to a *prior probability distribution* of a Gaussian
with mean :math:`\mu=0` and standard deviation :math:`\sigma=1`, i.e.

.. math::
    p(b) = \frac{1}{\sqrt{2\pi}}e^{-\frac{1}{2}b^2}.

However, let's say that we don't know anything in advance about the slope :math:`m`. In this case,
we say that the prior on :math:`m` is *uninformative*, i.e. all values are equally likely 
(again, this is before we even consider any data), which manifests mathematically as 

.. math::
    p(m) \propto 1.

In ``rcr``, prior probability distributions can be specified for any or all of the parameters of a model,
which will affect the rejection of outliers (essentially by modifying the rejection probability of certain measurements
according to the prior probabilities of all of the model parameter solutions that these measurements can contribute to). For simplicity and ease-of-use,
we've included two types of common priors within the library, as well as allowing for any sort of custom
prior PDF. These options are described in the table below.

.. _priorstypes:

Types of Model Parameter Priors in RCR
--------------------------------------

====================== =============================================================================
Prior Type             Parameters Needed To Specify        
====================== =============================================================================
``GAUSSIAN_PRIORS``    Means and standard deviations of some or all model parameters           
``CONSTRAINED_PRIORS`` Lower and/or upper bounds on some or all model parameters         
``MIXED_PRIORS``       Combination of some or all of the two above
``CUSTOM_PRIORS``      For some or all model parameters :math:`a_j`, custom prior PDF :math:`p(a_j)`
====================== =============================================================================

Now, how can we use these different types of priors in practice?

Gaussian Priors
^^^^^^^^^^^^^^^

Let's say that you want to apply Gaussian/normal prior probability distributions on some (or all) of your
model parameters. To do so, you'll first need to create a list, where each element of the list
corresponds to a model parameter, and is itself a list of 1) the mean of the Gaussian
for that parameter's prior and 2) the standard deviation of the same. If no Gaussian prior
is desired for a certain parameter, just give NaNs for those fields.

This is pretty dense, so we'll show a specific instance of this usage. Following the example within the introduction to
this section (:ref:`priors`), lets use the same linear model as before, and apply a Gaussian prior to the 
intercept :math:`b`, with mean :math:`\mu=0` and standard deviation :math:`\sigma=1`. We'll use no prior 
(uninformative) on the slope :math:`m`. From here, our list of parameters (not model parameters) 
that describe the Gaussian priors will be: 

.. code-block:: python

    gaussianParams = [
            [0,             1], # mu = 0, sigma = 1
            [float('nan'),  float('nan')] 
            # no prior on the slope parameter, so just use NaNs
        ]

Now, to introduce these priors before performing any fitting/RCR, we'll need to create an instance of 
the ``Priors`` class from ``rcr``, making sure to specify which type of prior we're implementing using
the correct object from the above table (in this case ``GAUSSIAN_PRIORS``). Here it is in code:

.. code-block:: python

    mypriors = rcr.Priors(rcr.GAUSSIAN_PRIORS, gaussianParams)

From here, RCR can be performed as usual, by 1) supplying the optional argument ``has_priors=True``
to the ``FunctionalForm`` constructor when initializing the model, and after that 2) 
initializing the ``priors`` attribute of your model with your ``Priors`` object:, e.g.:

.. code-block:: python

    model = rcr.FunctionalForm(linear,
        x,
        y,
        [linear_partial1, linear_partial2],
        guess,
        has_priors=True
    )

    model.priors = mypriors

From here RCR can be utilized with this model given the usual methods.

Constrained/Bounded Priors
^^^^^^^^^^^^^^^^^^^^^^^^^^

Another very common type of prior is to give hard constraints/bounds on 
certain model parameters. Following the same linear example, let's say
that we know that the slope :math:`m` of our model should be nonnegative 
(this type of prior is often for some physical reason), but we don't know anything about
the intercept :math:`b`. 

Similar to the usage of Gaussian priors, to implement this
we'll need to create a list where each element corresponds to a model parameter, and is
itself a list of 1) the lower bound and 2) the upper bounds that we want to give the corresponding parameter
if we only want to supply one (or neither) of the bounds, just use a NaN instead. 
Following our chosen example, this list can be coded as

.. code-block:: python

    paramBounds = [
            [float('nan'),  float('nan')] 
            [0,             float('nan')]  # constrain m > 0
        ]

Next, we need to instantiate an ``rcr.Priors`` object, in a similar
manner to the case of Gaussian priors (except now being sure to specify
``CONSTRAINED_PRIORS``):

.. code-block:: python

    mypriors = rcr.Priors(rcr.CONSTRAINED_PRIORS, paramBounds)

Finally, we'll need to initialize our model with the priors as in the end of the previous section 
(again with ``has_priors=True``),
and then we're good to go.

Both Gaussian and/or Constrained (Mixed) Priors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

What if we want to apply Gaussian priors to some model parameters, constrained priors to others,
or even a mix of both for certain parameters (e.g. force a parameter to be positive, while also
making it Gaussian-distributed)? To do this, simply create the lists that specify these priors---
``paramBounds`` and ``gaussianParams`` following the previous examples---and supply them both
to the constructor for your ``Priors`` object, making sure to specify the priors type as 
``MIXED_PRIORS``:

.. code-block:: python

    mypriors = rcr.Priors(rcr.MIXED_PRIORS, gaussianParams, paramBounds)

From here, RCR can be used as normal, after initializing our model (with ``has_priors=True``)
and supplying the model with the Priors object.

Custom Priors
^^^^^^^^^^^^^

In the most general case, RCR can work with any type of prior
probability distributions/density functions. To implement this,
you'll need a function :math:`\vec{p}(\vec{\theta})` that takes in
a vector of model parameters :math:`\vec{\theta}`, and returns a vector
of each parameter's prior probability density function evaluated 
given the corresponding parameter's value.

As an example, let's consider that for our linear model, we'd like to 1) place
an (unusual) prior on :math:`b`:

.. math::
    p(b) = e^{-|b|}\left|\cos^2b\right|,

and 2) constrain :math:`m` to be within the interval :math:`(0, 2]`. 
We can then implement :math:`\vec{p}(\vec{\theta})` as:

.. code-block:: python

    def prior_pdfs(model_parameters):
        pdfs = np.zeros(2) # vector of model parameter density function values
        b = model_parameters[0]
        pdfs[0] = np.exp(-np.abs(b)) * np.abs(np.cos(b)**2.)

        b = model_parameters[0]
        pdfs[1] = 1 if 0 < m <= 2 else 0 
        # p(m) = 0 if m is outside bounds of (0, 2]

        return pdfs

After such a :math:`\vec{p}(\vec{\theta})` is defined, we'll need to 
use it to instantiate an ``rcr.Priors`` object as usual, this time
declaring our type of priors as ``CUSTOM_PRIORS``:

.. code-block:: python

    mypriors = rcr.Priors(rcr.CUSTOM_PRIORS, prior_pdfs)

After creating our model (with ``has_priors=True``) and supplying it with
our Priors object ``mypriors``, RCR can then be used as usual.

.. _pivots:

Automatically Minimizing Correlation between Linear Model Parameters (Advanced)
-------------------------------------------------------------------------------

Let's again consider a linear model :math:`y = b + mx`. Usually
the fitted slope :math:`m` will be correlated with the fitted intercept :math:`b`.
Why is this? Consider redefining this model as :math:`y = b + m(x-x_p)`, with some
*pivot point* :math:`x_p`. Then, the intercept parameter is effectively :math:`b-mx_p`.
Therefore, given some fitted :math:`m`, the fitted intercept will be impacted by
:math:`m`, with the degree of this depending on choice of :math:`x_p`. 
As shown in `Trotter (2011) <https://cdr.lib.unc.edu/concern/dissertations/1544bq461>`_,
there exists some optimal :math:`x_p` that minimizes the correlation between :math:`b` and :math:`m`.

RCR has an algorithm for this; but first, why does this matter? One reason is
that if the linear parameters are uncorrelated with each other, then the uncertainty
in their estimation can be represented with simple, one-dimensional error bars. However,
if the two parameters *do* have some nontrivial correlation, then in order to properly present
their uncertainties, a full 2D correlation ellipse (`e.g. <http://www.math.wpi.edu/saspdf/insight/chap18.pdf>`_)
is needed, making things more complicated.

RCR's method for this correlation minimization works for any model that can be written as :math:`y = b + m(x-x_p)`.
So, for example, we could have some power-law model 

.. math::
    y(x|a_0, a_1) = a_0\left(\frac{x}{10^{x_p}}\right)^{a_1},

that can be actually be *linearized* (to be used with RCR's pivot point optimizing algorithm)
as follows:

.. math::
    y(x) &= a_0\left(\frac{x}{10^{x_p}}\right)^{a_1}

    \log_{10}y(x) &= \log_{10}\left[a_0\left(\frac{x}{10^{x_p}}\right)^{a_1}\right]

    &= \log_{10}a_0 + \log_{10}\left(\frac{x}{10^{x_p}}\right)^{a_1}

    &= \log_{10}a_0 + a_1\log_{10}\left(\frac{x}{10^{x_p}}\right)

    &= \log_{10}a_0 + a_1\left[\log_{10}x - \log_{10}10^{x_p}\right] 

    \log_{10}y(x) &\equiv \log_{10}a_0 + a_1\left[\log_{10}x - (\log_{10}x)_p\right]

So, the linearized version of this power law has intercept :math:`\log_{10}a_0`, slope :math:`a_1`,
pivot point :math:`(\log_{10}x)_p`, and the data transforms as :math:`\log_{10}y \rightarrow y` and :math:`\log_{10}x \rightarrow x`.
The formula for the pivot point is 

.. math::
    (\log_{10}x)_p = \frac{\sum\limits_iw_i(\log_{10} x_i)y^2(x_i)}{\sum\limits_i w_iy^2(x_i)}

(`Maples et al. 2018 <https://arxiv.org/abs/1807.05276>`_, Section 8.3.5); we'll need such a formula
for the pivot point of any model that we'd like to apply this procedure to. Keeping that in mind,
lets look at this in code.

To use the slope-intercept correlation minimization procedure with RCR, we'll need to define
this pivot point function. However, first read the following note:

.. note::
    Pivot point functions need to be defined with parameters of 1) `xdata`, 
    a list or array of the :math:`x`-data in the dataset, 2) `weights`, a list/array 
    of the weights of the dataset, 3) :math:`f`, the model function, and 4) :math:`params`,
    a list/array of model parameters. (To be made easier in a future patch)

Keeping this in mind, let's define our pivot point function for this power law model:

.. code-block:: python

    def get_pivot_powerlaw(xdata, weights, f, params):
        topsum = np.sum(weights * np.log10(xdata) * np.power(f(xdata, params), 2.))
        bottomsum = np.sum(weights * np.power(f(xdata, params), 2.))
        
        return topsum / bottomsum

Now, we need to define our model, making sure to use the pivot point. In order to use pivot points
within the function definition of a model (and its derivatives), we'll need to use the
static attribute ``pivot`` of the ``rcr.FunctionalForm`` class (or ``pivot_ND`` for the :math:`n`-dimensional case)
within these definitions. So, our power-law model and its parameter-derivatives can be defined as:

.. code-block:: python

    def powerlaw(x, params):
        a0 = params[0]
        a1 = params[1]
        return a0 * np.power(x / np.power(10., rcr.FunctionalForm.pivot), a1)

    def powerlaw_partial1(x, params):
        a1 = params[1]
        return np.power((x / np.power(10., rcr.FunctionalForm.pivot)), a1)

    def powerlaw_partial2(x, params):
        a0 = params[0]
        a1 = params[1]
        piv = rcr.FunctionalForm.pivot # renamed for brevity

        return a0 * np.power((x / np.power(10., piv)), a1) * np.log(x / np.power(10., piv))

Next, we can use this with RCR as normal, except now supplying additional arguments
of ``pivot_function`` and ``pivot_guess`` to the model constructor ``rcr.FunctionalForm``, 
where ``pivot_function`` is the function that returns
the pivot for model give ``xdata, weights, f, params``, and ``pivot_guess`` is a guess for the optimal
pivot point (for the iterative optimization algorithm), that should accompany the initial guess for the model parameters. 
For example, if our initial guess for the pivot point
is some ``pivot_guess = 1.5``, we could initialize our model as:

.. code-block:: python

    model = rcr.FunctionalForm(powerlaw,
        xdata,
        ydata,
        [powerlaw_partial1, powerlaw_partial2],
        guess,
        pivot_function=get_pivot_powerlaw,
        pivot_guess=pivot_guess
    )

From here, we can perform RCR as normal, and access the optimal value for the pivot points found by RCR
with ``model.result.pivot`` (or ``model.result.pivot_ND`` for the `n`-dimensional model case).

Finally, note that the support for :math:`n`-dimensional models (i.e. :math:`n` independent variables)
is still available when using this feature; in this case, your pivot point function
should return a list/array of :math:`n` pivot points.