.. _functional:

Rejecting Outliers While Model Fitting
======================================

Introduction
------------

In it's most simple form, RCR is an excellent tool for
detecting and rejecting outliers within heavily contaminated one-dimensional 
datasets, as shown in :ref:`singlevalue`. However, this only scratches
the surface of RCR. In it's more generalized form, RCR can also be used
to reject outliers within some :math:`n`-dimensional dataset 
while also *simultaneously* fitting a model to
that dataset. This section will explain how this can be done
fairly easily in practice, while avoiding going into unnessarily
technicalities. We recommend reading :ref:`singlevalue` before
tackling this section, as the following is essentially a generalization
on that section.

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

Next, let's create our dataset. We'll create a dataset with :math:`N=200` points
total, with 85% datapoints being outliers. Our "true" model that the datapoints
will be generated about will have parameters of :math:`b=0` and :math:`m=1`. 
In code, this is simply:

.. code-block:: python

    import numpy as np

    N = 200 # number of datapoints
    f = 0.85 # fraction of datapoints that are outliers

    params_true = [0, 1] # parameters of "true" model

Next, let's create the data. We'll generate our datapoints in a certain range
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
that naively should work with the data, but are pretty far off of the true values of 
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
(see the paper **HERE** for more details), RCR can iteratively reject outliers and
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

Model Parameter Uncertainties/Error Bars
----------------------------------------

In many cases, we often want not just best fit parameters for a model and dataset,
but also *uncertainties*, or "error bars" for these parameters. This is easily
available in ``rcr``, simply via the 

Weighting Data
--------------

Data with Uncertainties
-----------------------

Applying Prior Knowledge to Model Parameters
--------------------------------------------

Automatically Minimizing Correlation between Linear Model Parameters (Advanced)
-------------------------------------------------------------------------------

Why is it useful (think of Dan's examples...)