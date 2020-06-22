import rcr
import numpy as np

np.random.seed(18318) # get consistent random results

# define model and model param derivatives

def linear(x, params):          # model function
    return params[0] + x * params[1]

def linear_partial1(x, params): # first model parameter derivative
    return 1

def linear_partial2(x, params): # second model parameter derivative
    return x


# build dataset
import numpy as np

N = 200 # number of datapoints
f = 0.85 # fraction of datapoints that are outliers

params_true = [0, 1] # parameters of "true" model

sigma_uncontaminated = 1
sigma_contaminated = 10

# generate x-datapoints randomly in an interval
xdata_uncontaminated = np.random.uniform(-10, 10, int(N * (1 - f)))
xdata_contaminated = np.random.uniform(-10, 10, int(N * f))

# generate y-datapoints
ydata_uncontaminated = np.random.normal(
    linear(xdata_uncontaminated, params_true),
    
    )

# plot dataset and underlying model
# (with outlier and nonoutlier dists?)