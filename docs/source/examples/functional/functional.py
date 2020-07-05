import rcr
import numpy as np

np.random.seed(18318) # get consistent random results

# define model and model param derivatives

def linear(x, params):          # model function
    return params[0] + x * params[1]

def d_linear_1(x, params): # first model parameter derivative
    return 1

def d_linear_2(x, params): # second model parameter derivative
    return x


# build dataset
import numpy as np

N = 200 # number of datapoints
f = 0.85 # fraction of datapoints that are outliers

params_true = [0, 1] # parameters of "true" model

sigma_uncontaminated = 1
sigma_contaminated = 10

# generate x-datapoints randomly in an interval
x_range = (-10, 10)

xdata_uncontaminated = np.random.uniform(x_range[0], x_range[1], int(N * (1 - f)))
xdata_contaminated = np.random.uniform(x_range[0], x_range[1], int(N * f))


# generate y-datapoints about the true model

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

# plot dataset
import matplotlib.pyplot as plt

plt.figure(figsize=(8,5))
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

# plt.show()

# initialize model and guess for model parameters

guess = [5, 1.5]

model = rcr.FunctionalForm(linear,
    xdata,
    ydata,
    [d_linear_1, d_linear_2],
    guess
)


# initialize and perform RCR

r = rcr.RCR(rcr.LS_MODE_68) # setting up for RCR with this rejection technique

r.setParametricModel(model) # tell RCR that we are model fitting

r.performBulkRejection(ydata) # perform RCR


# view results
best_fit_parameters = model.result.parameters # best fit parameters

rejected_data = r.result.rejectedY # rejected and non-rejected data
nonrejected_data = r.result.cleanY
nonrejected_indices = r.result.indices

print(best_fit_parameters)


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

# plt.show()

# model param error bars

# create measurement uncertainties
error_y_uncontaminated = np.random.uniform(low=0.1, high=1, size=int(N * (1 - f)))

error_y_contaminated = np.random.uniform(low=1, high=2, size=int(N * f))

error_y = np.concatenate((error_y_contaminated, error_y_uncontaminated))


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


# view results
best_fit_parameters = model.result.parameters # best fit parameters
best_fit_parameter_errors = model.result.parameter_uncertainties # and their uncertainties

rejected_data = r.result.rejectedY # rejected and non-rejected data
nonrejected_data = r.result.cleanY
nonrejected_indices = r.result.indices

print(best_fit_parameters)
print(best_fit_parameter_errors)


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

# plt.show()



# priors 

# # Gaussian
# gaussianParams =    [[float('nan'), float('nan')], [1.0, 2.0] ]
# mypriors = rcr.Priors(rcr.MIXED_PRIORS, gaussianParams, boundedParams)

# model = rcr.FunctionalForm(linear,
#     x,
#     y,
#     [linear_partial1, linear_partial2],
#     guess,
#     tol=tolerance,
#     weights=w,
#     error_y=err_y,
#     has_priors=True
# )

# model.priors = mypriors

# #Bounded/constrained
# boundedParams =     [[0.0, float('nan')], [float('nan'), float('nan')] ]


# pivots

# define function that returns pivot point
def get_pivot_powerlaw(xdata, weights, f, params):
    topsum = np.sum(weights * np.log10(xdata) * np.power(f(xdata, params), 2.))
    bottomsum = np.sum(weights * np.power(f(xdata, params), 2.))
    
    return topsum / bottomsum


weights = np.ones(len(xdata))

xp = get_pivot_powerlaw(np.abs(xdata), weights, linear, guess)

print(xp)