.. _functional:

Rejecting Outliers While Model Fitting
======================================

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
a simple, one-dimensional linear model :math:`y = b + mx`. In this case,
:math:`n = 1`, and our parameter vector is just :math:`\vec{\theta} = (b ,m)`.

For this more-general case, what does our dataset look like? Each datapoint
will be some value of :math:`y(\vec{x}|\vec{\theta})` associated with
a value for :math:`\vec{x}`. As such, if we have :math:`N` datapoints in total,
indexing each by :math:`i`, our dataset can be written compactly 
as :math:`\left\{\left(\vec{x}_i, y_i\right)\right\}_{i=1}^N`
(be sure to avoid getting :math:`N` confused with :math:`n` here; 
the former is the number of datapoints that we're fitting the 
model to, while the latter is the dimensionality of the dataset/model).


Automatically Minimizing Correlation between Linear Model Parameters
--------------------------------------------------------------------

Why is it useful (think of Dan's examples...)