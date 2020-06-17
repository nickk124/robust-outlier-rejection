import rcr

def linear(x, params):
    return params[0] + x * params[1]

def linear_partial1(x, params):
    return 1

def linear_partial2(x, params):
    return x

def get_pivot(xdata, weights, f, params):
    sum = 0

    for x in xdata:
        sum += x
    return sum / len(xdata)

# if __name__ == "__main__":
y = [0.1, 0.2, 0, 0.3, -0.1, -0.2, -0.3, 11]
# w = [0.5, 1.0, 1.2, 3.0, 1.0, 0.6, 0.2, 1.0, 0.1, 1.0]
rcr_obj = rcr.RCR(rcr.ES_MODE_DL)
rcr_obj.performBulkRejection(y)

print(rcr_obj.result.flags)


tolerance = 0.01
guess = [2.1, 0.9]
y = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 20.0, 11.0]
x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
w = [0.5, 1.0, 1.2, 3.0, 1.0, 0.6, 0.2, 1.0, 0.1, 2.0]
err_y = [0.1 for xx in x]

gaussianParams =    [[float('nan'), float('nan')], [1.0, 2.0] ]
boundedParams =     [[0.0, float('nan')], [float('nan'), float('nan')] ]
priors = rcr.Priors(rcr.MIXED_PRIORS, gaussianParams, boundedParams)

# boundedParams =     [[0, 1], [-20, 20] ]
# priors = rcr.Priors(rcr.CONSTRAINED_PRIORS, boundedParams)

pivot_guess = 4.5

# // FunctionalForm(func_linear, x, y, partialsvector_lin, tolerance, guess, w, testPriors, getAvg, xb_guess);

model = rcr.FunctionalForm(linear,
    x,
    y,
    [linear_partial1, linear_partial2],
    guess,
    tol=tolerance,
    weights=w,
    error_y=err_y,
    has_priors=True,
    pivot_function=get_pivot,
    pivot_guess=pivot_guess
)

model.priors = priors

r = rcr.RCR(rcr.SS_MEDIAN_DL) #setting up for RCR with this rejection technique
r.setParametricModel(model)
r.performBulkRejection(w, y) # //running Bulk rejection RCR

final_parameters = model.parameters

print(final_parameters)

rejected_data = r.result.rejectedY
nonrejected_data = r.result.cleanY
flags = r.result.flags

print(rejected_data)
print(nonrejected_data)
print(flags)

pivot = model.result.pivot
print(pivot)