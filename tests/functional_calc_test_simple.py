import numpy as np
import rcr

test_server = False

if not test_server:

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

    r = rcr.RCR() # setting up for RCR with this rejection technique

    r.setRejectionTech(rcr.ES_MODE_DL)

    r.setParametricModel(model) # tell RCR that we are model fitting
    r.performBulkRejection(wdata, ydata) # perform RCR

    best_fit_parameters = model.result.parameters # best fit parameters
    rejected_data = r.result.rejectedY # rejected and non-rejected data
    nonrejected_data = r.result.cleanY
    nonrejected_indices = r.result.indices
    flags = r.result.flags
    # indices from original dataset of nonrejected data

    print(best_fit_parameters, rejected_data, nonrejected_data, nonrejected_indices, flags)

else:
    def function_linear(x, params): # model function
        return params[0] + x * params[1]

    def partial1_linear(x, params): # first model parameter derivative
        return 1

    def partial2_linear(x, params): # second model parameter derivative
        return x



    def function_quadratic(x,  params):
        a0 = params[0]
        a1 = params[1]
        a2 = params[2]

        return a0 + a1 * x + a2 * (x**2.0)


    def partial1_quadratic(x, params):

        return 1.0


    def partial2_quadratic(x, params):

        return x


    def partial3_quadratic(x, params):

        return np.power(x, 2.0)






    def function_cubic(x, params):
        a0 = params[0]
        a1 = params[1]
        a2 = params[2]
        a3 = params[3]

        return a0 + a1 * x + a2 * np.power(x, 2.0) + a3 * np.power(x, 3.0)


    def partial1_cubic(x, params):

        return 1.0


    def partial2_cubic(x, params):

        return x


    def partial3_cubic(x, params):

        return np.power(x, 2.0)


    def partial4_cubic(x, params):

        return np.power(x, 3.0)






    def function_powerlaw(x, params):
        a0 = params[0]
        a1 = params[1]

        return a0 * np.power(x, a1)


    def partial1_powerlaw(x, params):
        a1 = params[1]

        return np.power(x, a1)


    def partial2_powerlaw(x, params):
        a0 = params[0]
        a1 = params[1]

        return a0 * np.power(x, a1) * np.log(x)





    def function_exponential(x, params):
        a0 = params[0]
        a1 = params[1]

        return a0 * np.power(10, a1*x)


    def partial1_exponential(x, params):
        a1 = params[1]

        return np.power(10, a1*x)


    def partial2_exponential(x, params):
        a0 = params[0]
        a1 = params[1]

        return np.power(10, a1*x) * x * np.log(10)



    def function_logarithmic(x, params):
        a0 = params[0]

        return a0 * np.log10(x)

    def partial1_logarithmic(x, params):
        a0 = params[0]

        return np.log10(x)

    import rcr
    print('----------------------------------')

    mainData = request.vars.mainData
    xData = request.vars.xData
    guess = request.vars.guessData
    symDist = stringToBool(request.vars.symDist)
    equalCont = stringToBool(request.vars.equalCont)
    oneSidedCont = stringToBool(request.vars.oneSidedCont)
    inBetweenCont = stringToBool(request.vars.inBetweenCont)
    rejType = request.vars.rejectionType
    weighted = stringToBool(request.vars.weighted)
    functionType = request.vars.functionType
    # bar_guess = request.vars.barData
    priors = request.vars.priorsData
    hasPriors = stringToBool(request.vars.priorsCheck)
    print('request received')

    mainData = mainData.split("\n") #puts the input data into a list
    xData = xData.split("\n")
    guess = guess.split("\n")
    priorsData = priors.split("\n") #watch this line

    # bar_guess = bar_guess.rstrip()
    # bar_guess = bar_guess.lstrip() #removes any spaces on left or right
    # pivot_guess = float(bar_guess) #converts to float

    yValues = []
    xValues = []
    wValues = []
    guessValues = []

    lb_values = []
    ub_values = []
    mu_values = []
    dev_values = []

    hasLb = []
    hasUb = []
    hasMu = []
    hasDev = []

    for i in mainData:
        print(i)
        if len(i) == 0:
            continue
        i.replace("\t"," ") #replaces any potential tab delimiter with a space
        dataPoint = i.split() #list made from the row of inputData
        if weighted:
            yValues.append(float(dataPoint[0]))
            wValues.append(float(dataPoint[1]))
        else:
            yValues.append(float(dataPoint[0]))
            wValues.append(1.0)

    for i in xData: # doesn't account for 3D case
        if len(i) == 0:
            continue
        xValues.append(float(i))

    for i in guess:
        if len(i) == 0:
            continue
        guessValues.append(float(i))

    boundedParams = []
    gaussianParams = []
    myPriorType = None
    priors_obj = None

    if hasPriors: #creates list of each of the four prior params
        for i in priorsData:
            if len(i) == 0:
                continue
            i.replace("\t"," ") #replaces any potential tab delimiter with a space
            dataPoint = i.split() #list made from the row of inputData

            lb_s = dataPoint[0] #individual strings
            ub_s = dataPoint[1]
            mu_s = dataPoint[2]
            dev_s = dataPoint[3]

            boundedParams.append([formatPriorString(lb_s), formatPriorString(ub_s)])
            gaussianParams.append([formatPriorString(mu_s), formatPriorString(dev_s)])

            if lb_s == "x" or lb_s == "X" :
                hasLb.append(0)
            else:
                hasLb.append(1)

            if ub_s == "x" or ub_s == "X":
                hasUb.append(0)
            else:
                hasUb.append(1)

            if mu_s == "x" or mu_s == "X":
                hasMu.append(0)
            else:
                hasMu.append(1)

            if dev_s == "x" or dev_s == "X":
                hasDev.append(0)
            else:
                hasDev.append(1)

    CountArr = [0, 0, 0, 0]
    hasArray = [hasLb, hasUb, hasMu, hasDev]

    ct = 0
    for list in hasArray:
        for i in list:
            if i==1:
                CountArr[ct] += 1
        ct += 1

    if (CountArr[0] > 0 or CountArr[1] > 0) and (CountArr[2] == 0 and CountArr[3] == 0):
        myPriorType = rcr.CONSTRAINED_PRIORS
        priors_obj = rcr.Priors(myPriorType, boundedParams)

    elif (CountArr[0] == 0 and CountArr[1] == 0) and (CountArr[2] > 0 and CountArr[3] > 0):
        myPriorType = rcr.GAUSSIAN_PRIORS
        priors_obj = rcr.Priors(myPriorType, gaussianParams)

    elif (CountArr[0] > 0 or CountArr[1] > 0) and (CountArr[2] > 0 and CountArr[3] > 0):
        myPriorType = rcr.MIXED_PRIORS
        priors_obj = rcr.Priors(myPriorType, gaussianParams, boundedParams)
    else:
        priors_obj = rcr.Priors()

    #create int that defines the func type
    model_function = None
    model_partials = None
    # model_pivotfunc = get_pivot_default


    if functionType == 'Linear':
        model_function = function_linear
        model_partials = [partial1_linear, partial2_linear]

    elif functionType == 'Quadratic':
        model_function = function_quadratic
        model_partials = [partial1_quadratic, partial2_quadratic, partial3_quadratic]

    elif functionType == 'Cubic':
        model_function = function_cubic
        model_partials = [partial1_cubic, partial2_cubic, partial3_cubic, partial4_cubic]

    elif functionType == 'PowerLaw':
        model_function = function_powerlaw
        model_partials = [partial1_powerlaw, partial2_powerlaw]
        # model_pivotfunc = get_pivot_powerlaw

    elif functionType == 'Exponential':
        model_function = function_exponential
        model_partials = [partial1_exponential, partial2_exponential]
        # model_pivotfunc = get_pivot_exponential

    elif functionType == 'Logarithmic':
        model_function = function_logarithmic
        model_partials = [partial1_logarithmic]
        # model_pivotfunc = get_pivot_logarithmic

    else:
        # return
        pass

    #Rejection tech determination:
    rejTechInt = int(rejType)

    r  =  rcr.RCR()

    rejTech = ''
    if symDist:
        if equalCont:
            r.setRejectionTech(rcr.SS_MEDIAN_DL)
            rejTech = 'ss_med_dl'
        elif oneSidedCont:
            r.setRejectionTech(rcr.LS_MODE_68)
            rejTech = 'ls_mode_68'
        elif inBetweenCont:
            r.setRejectionTech(rcr.LS_MODE_DL)
            rejTech = 'ls_mode_dl'
    else:
        r.setRejectionTech(rcr.ES_MODE_DL)
        rejTech = 'es_mode_dl'

    
    print(model_function, 
        xValues, 
        yValues, 
        model_partials, 
        guessValues, 
        wValues, 
        hasPriors
    )

    model = rcr.FunctionalForm(
        model_function, 
        xValues, 
        yValues, 
        model_partials, 
        guessValues, 
        weights=wValues, 
        has_priors=hasPriors
        )
    # model = rcr.FunctionalForm(model_function, xValues, yValues, model_partials, guessValues, weights=wValues, has_priors=hasPriors, pivot_function=model_pivotfunc, pivot_guess=pivot_guess)
    model.priors = priors_obj

    r.setParametricModel(model)
    r.performBulkRejection(wValues, yValues)

    parameterslist = [str(round_to_n(param, 5)) for param in model.result.parameters] #strings
    flags = r.result.flags


    #add pivot value to end of parameterslist
    # parameterslist.append(str(round_to_n(model.result.pivot, 5)))

    #creates lists of the nonrejected and the rejected data

    nonRejectedList = [] #formatted nonweighted as ['x1 y1', 'x2 y2', 'x3 y3'...] strings or weighted as ['x1 y1 w1', 'x2 y2 w2', ...]
    rejectedList = []
    allDataList = [] #formatted nonweighted as ['x1 y1 f1', 'x2 y2 f2', 'x3 y3 f3'...] strings or weighted as ['x1 y1 w1 f1', 'x2 y2 w2 f2', ...]

    for k in range(len(xValues)):
        if flags[k]:
            nonRejectedList.append(str(xValues[k]) + '\t' + str(yValues[k]))
            allDataList.append(str(xValues[k]) + '\t' + str(yValues[k]) + '\t' + 'True')
        else:
            rejectedList.append(str(xValues[k]) + '\t' + str(yValues[k]))
            allDataList.append(str(xValues[k]) + '\t' + str(yValues[k]) + '\t' + 'False')





    #finally, formats the data to be sent via json

    nonRejectedData = '\n'.join(nonRejectedList)
    rejectedData = '\n'.join(rejectedList)
    allData = '\n'.join(allDataList)
    finalParametersData = '\n'.join(parameterslist)

    print('=-=-=-=-=-=')
    jsonData = dict(allData=allData, nonRejectedData=nonRejectedData, rejectedData=rejectedData, finalParametersData=finalParametersData)
    