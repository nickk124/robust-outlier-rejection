import math as m
import rcr

def function_linear(x, params): # model function
    return params[0] + x * params[1]

def partial1_linear(x, params): # first model parameter derivative
    return 1

def partial2_linear(x, params): # second model parameter derivative
    return x


xdata = [1,2,3,4,5,6,7,8,9,10]
ydata = [-6,-5,-4,-3,-2,-1,25,1,2,3]
wdata = [1,1,1,1,1,1,0.2,1,1,1]

guess = [-7, 1]

model = rcr.FunctionalForm(
    function_linear,
    xdata,
    ydata,
    [partial1_linear, partial2_linear],
    guess,
    weights=wdata,
)

r = rcr.RCR(rcr.LS_MODE_68) # setting up for RCR with this rejection technique
r.setParametricModel(model) # tell RCR that we are model fitting
r.performBulkRejection(wdata, ydata) # perform RCR

best_fit_parameters = model.result.parameters # best fit parameters
rejected_data = r.result.rejectedY # rejected and non-rejected data
nonrejected_data = r.result.cleanY
nonrejected_indices = r.result.indices
# indices from original dataset of nonrejected data

print(best_fit_parameters, rejected_data, nonrejected_data, nonrejected_indices)