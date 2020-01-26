import smtplib
import math
import numpy as np
from gluon.serializers import json

def validateInputData(input):
    input = input.split("\n");

    for i in input:
        print i
        try:
            float(i)
        except ValueError:
            raise Http(400, 'Must provide data')
    print 'test'

def stringToBool(s):
	return s == "true"

def round_to_n(x, n):
    if x == 0:
        return x
    else:
        return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def formatPoints(yVal, wVal=None, flags=None):
    s = ''
    print len(yVal)
    for i in range(len(yVal)):
        s += str(yVal[i])
        if wVal is not None:
            s += '\t'
            s += str(wVal[i])
        if flags is not None:
            s += '\t'
            s += str(flags[i])
        s += '\n'

    return s

def calculateStDev(delta, mu, data):
    sum = 0.0
    for i in range(len(data)):
        sum += (data[i] - mu)**2

    return math.sqrt(sum / (len(data) - delta))

def calculateWeightedStDev(delta, mu, weights, data):
    top = 0.0
    wSum = 0.0
    wSqSum = 0.0
    for i in range(len(data)):
        w = weights[i]
        top += w * ((data[i] - mu)**2)
        wSum += w
        wSqSum += w**2

    return math.sqrt(top / (wSum - delta*wSqSum / wSum));

def characterizeData():
    inputData = request.vars.inputData
    weighted = False
    inputData = inputData.split("\n")
    yValues = []
    wValues = []
    for i in inputData:
        if len(i) == 0:
            continue
        i.replace("\t"," ")
        dataPoint = i.split()
        if weighted:
            yValues.append(float(dataPoint[0]))
            wValues.append(float(dataPoint[1]))
        else:
            yValues.append(float(dataPoint[0]))

    ySum = 0.0
    if weighted:
        wSum = 0.0
        for i in range(len(yValues)):
            ySum += yValues[i] * wValues[i]
            wSum += wValues[i]
    else:
        for i in range(len(yValues)):
            ySum += yValues[i]
        wSum = len(yValues)
    mu = ySum / wSum

    yValuesAbove = []
    yValuesBelow = []
    wValuesAbove = []
    wValuesBelow = []

    for i in range(len(yValues)):
        val = yValues[i]
        if weighted:
            wVal = wValues[i]
        else:
            wVal = 1.0
        if val > mu:
            yValuesAbove.append(val)
            wValuesAbove.append(wVal)
        elif val < mu:
            yValuesBelow.append(val)
            wValuesBelow.append(wVal)
        else:
            yValuesAbove.append(val)
            yValuesBelow.append(val)
            wValuesAbove.append(.5*wVal)
            wValuesBelow.append(.5*wVal)


    if weighted:
        stDevTotal = calculateWeightedStDev(1.0, mu, wValues, yValues)
    else:
        stDevTotal = calculateStDev(1.0, mu, yValues)

    stDevAbove = calculateWeightedStDev(.5, mu, wValuesAbove, yValuesAbove)
    stDevBelow = calculateWeightedStDev(.5, mu, wValuesBelow, yValuesBelow)

    jsonData = dict(mu=float('%.6g' % mu), stDevBelow=float('%.6g' % stDevBelow), stDevAbove=float('%.6g' % stDevAbove), stDevTotal=float('%.6g' % stDevTotal))
    #from gluon.serializers import json
    return json(jsonData)

def callRCR():
    import RCRLib as R
    print '----------------------------------'
    #if not validateInputData(request.vars.inputData): raise Http(400, 'Must provide data')
    #validateInputData(request.vars.inputData)

    inputData = request.vars.inputData #data received
    symDist = stringToBool(request.vars.symDist) #necessary bools received
    equalCont = stringToBool(request.vars.equalCont)
    oneSidedCont = stringToBool(request.vars.oneSidedCont)
    inBetweenCont = stringToBool(request.vars.inBetweenCont)
    weighted = stringToBool(request.vars.weighted)

    print inputData
    inputData = inputData.split("\n")

    yValues = []
    wValues = []

    for i in inputData:
        if len(i) == 0:
            continue
        i.replace("\t"," ")
        dataPoint = i.split()
        if weighted:
            yValues.append(float(dataPoint[0]))
            wValues.append(float(dataPoint[1]))
        else:
            yValues.append(float(dataPoint[0]))
    r  =  R.RCR()
    rejTech = ''
    if symDist:
        if equalCont:
            r.setRejectionTech(R.SS_MEDIAN_DL)
            rejTech = 'ss_med_dl'
        elif oneSidedCont:
            r.setRejectionTech(R.LS_MODE_68)
            rejTech = 'ls_mode_68'
        elif inBetweenCont:
            r.setRejectionTech(R.LS_MODE_DL)
            rejTech = 'ls_mode_dl'
    else:
        r.setRejectionTech(R.ES_MODE_DL)
        rejTech = 'es_mode_dl'

    yValues = R.DubVec(yValues)

    print '=-=-=-=-=-='
    if weighted:
        wValues = R.DubVec(wValues)
        r.performBulkRejection(wValues, yValues)
        nonRejectedData = formatPoints(r.result.cleanY, wVal=r.result.cleanW)
        rejectedData = formatPoints(r.result.rejectedY, wVal=r.result.rejectedW)
        allData = formatPoints(r.result.originalY, wVal=r.result.originalW, flags=r.result.flags)
    else:
        r.performBulkRejection(yValues)
        nonRejectedData = formatPoints(r.result.cleanY)
        rejectedData = formatPoints(r.result.rejectedY)
        allData = formatPoints(r.result.originalY, flags=r.result.flags)
        print allData

    print '=-=-=-=-=-='
    jsonData = dict(mu=float('%.6g' % r.result.mu), sigma=float('%.6g' % r.result.sigma), stDev=float('%.6g' % r.result.stDev),
                    sigmaBelow=float('%.6g' % r.result.sigmaBelow), sigmaAbove=float('%.6g' % r.result.sigmaAbove),
                    stDevBelow=float('%.6g' % r.result.stDevBelow), stDevAbove=float('%.6g' % r.result.stDevAbove), stDevTotal=float('%.6g' % r.result.stDevTotal),
                    allData=allData, nonRejectedData=nonRejectedData, rejectedData=rejectedData, rejTech=rejTech);
    #from gluon.serializers import json
    return json(jsonData)

def callFunctionalRCR():
    #==============================
    import RCRwebutils as rs
    print '----------------------------------'

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
    bar_guess = request.vars.barData
    priors = request.vars.priorsData
    hasPriors = stringToBool(request.vars.priorsCheck)
    print 'request received'

    mainData = mainData.split("\n") #puts the input data into a list
    xData = xData.split("\n")
    guess = guess.split("\n")
    priorsData = priors.split("\n") #watch this line

    bar_guess = bar_guess.rstrip()
    bar_guess = bar_guess.lstrip() #removes any spaces on left or right
    bar_guess = float(bar_guess) #converts to float

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
        print i
        if len(i) == 0:
            continue
        i.replace("\t"," ") #replaces any potential tab delimiter with a space
        dataPoint = i.split() #list made from the row of inputData
        if weighted:
            yValues.append(float(dataPoint[0]))
            wValues.append(float(dataPoint[1]))
        else:
            yValues.append(float(dataPoint[0]))

    for i in xData: # doesn't account for 3D case
        if len(i) == 0:
            continue
        xValues.append(float(i))

    for i in guess:
        if len(i) == 0:
            continue
        guessValues.append(float(i))

    if hasPriors: #creates list of each of the four prior params
        for i in priorsData:
            if len(i) == 0:
                continue
            i.replace("\t"," ") #replaces any potential tab delimiter with a space
            dataPoint = i.split() #list made from the row of inputData

            iList = []

            lb_s = dataPoint[0] #individual strings
            ub_s = dataPoint[1]
            mu_s = dataPoint[2]
            dev_s = dataPoint[3]

            if lb_s == "x" or lb_s == "X" :
                hasLb.append(0)
                lb_values.append(0)
            else:
                hasLb.append(1)
                lb_values.append(float(lb_s))

            if ub_s == "x" or ub_s == "X":
                hasUb.append(0)
                ub_values.append(0)
            else:
                hasUb.append(1)
                ub_values.append(float(ub_s))

            if mu_s == "x" or mu_s == "X":
                hasMu.append(0)
                mu_values.append(0)
            else:
                hasMu.append(1)
                mu_values.append(float(mu_s))

            if dev_s == "x" or dev_s == "X":
                hasDev.append(0)
                dev_values.append(0)
            else:
                hasDev.append(1)
                dev_values.append(float(dev_s))

    onlyBounded = False
    onlyGaussian = False
    both = False

    whichPrior = 0

    CountArr = [0, 0, 0, 0]
    hasArray = [hasLb, hasUb, hasMu, hasDev]

    ct = 0
    for list in hasArray:
        for i in list:
            if i==1:
                CountArr[ct] += 1
        ct += 1

    if (CountArr[0] > 0 or CountArr[1] > 0) and (CountArr[2] == 0 and CountArr[3] == 0):
        onlyBounded = True
        whichPrior = 2
    elif (CountArr[0] == 0 and CountArr[1] == 0) and (CountArr[2] > 0 and CountArr[3] > 0):
        onlyGaussian = True
        whichPrior = 1
    elif (CountArr[0] > 0 or CountArr[1] > 0) and (CountArr[2] > 0 and CountArr[3] > 0):
        both = True
        whichPrior = 3

    #create int that defines the func type
    funcTypeInt = 0

    if functionType == 'Linear':
        funcTypeInt = 1
    elif functionType == 'Quadratic':
        funcTypeInt = 2
    elif functionType == 'Cubic':
        funcTypeInt = 3
    elif functionType == 'PowerLaw':
        funcTypeInt = 4
    elif functionType == 'Exponential':
        funcTypeInt = 5
    elif functionType == 'Logarithmic':
        funcTypeInt = 6
    else:
        return

    #Rejection tech determination:
    rejTechInt = int(rejType)

    #converts the data into c++ std::vectors
    dataSize = len(xValues) #number of datapoints
    paramCount = len(guessValues) #number of parameters in the function

    xarray = rs.DoubleVector(dataSize)
    yarray = rs.DoubleVector(dataSize)
    guessarray = rs.DoubleVector(paramCount)

    warray = rs.DoubleVector(dataSize) #if nonweighted, will just be an empty array.

    priorsParamsArray = rs.DoubleVector(paramCount * 4) #only gaussian or only bounded will just ust the first half of this array
    hasPriorsArray = rs.IntVector(paramCount * 4)

    for i in range(0,dataSize):
        xarray[i] = xValues[i]
        yarray[i] = yValues[i]

    for j in range(0,paramCount):
        guessarray[j] = guessValues[j]

    if weighted: #creates weight vector if needed
        for k in range(0,dataSize):
            warray[k] = wValues[k]

    if onlyBounded:
        for j in range(0,paramCount):
            priorsParamsArray[2*j] = lb_values[j]
            priorsParamsArray[2*j+1] = ub_values[j]

            hasPriorsArray[2*j] = hasLb[j]
            hasPriorsArray[2*j+1] = hasUb[j]

    elif onlyGaussian:
        for j in range(0,paramCount):
            priorsParamsArray[2*j] = mu_values[j]
            priorsParamsArray[2*j+1] = dev_values[j]

            hasPriorsArray[2*j] = hasMu[j]
            hasPriorsArray[2*j+1] = hasDev[j]

    elif both:
        for j in range(0,paramCount):
            priorsParamsArray[2*j] = mu_values[j]
            priorsParamsArray[2*j+1] = dev_values[j]

            hasPriorsArray[2*j] = hasMu[j]
            hasPriorsArray[2*j+1] = hasDev[j]

        for j in range(0,paramCount):
            priorsParamsArray[paramCount*2 + 2*j] = lb_values[j]
            priorsParamsArray[paramCount*2 + 2*j+1] = ub_values[j]

            hasPriorsArray[paramCount*2 + 2*j] = hasLb[j]
            hasPriorsArray[paramCount*2 + 2*j+1] = hasUb[j]


    """
    xarray_np = np.array(xValues, dtype=np.float64)
    yarray_np = np.array(yValues, dtype=np.float64)
    guessarray_np = np.array(guessValues, dtype=np.float64)

    warray_np = np.array(wValues, dtype=np.float64)
    """
    """
    testx = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], dtype=np.float64)
    testy = np.array([2.0, 3.0, 400.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0], dtype=np.float64)
    testguess = np.array([2.1, 0.9], dtype=np.float64)
    testfType = 1;
    testdataSize = 10;
    testrejTechNo = 1;
    """
    # Getting the unique thread name:

    # Call the callable RCR!
    #order: (x, y, sigy, guess, fType, dataSize, rejTechNo)

    result = rs.DoubleVector()

    print 'performing RCR'

    if weighted:
        result = rs.requestHandlerWeighted(xarray, yarray, guessarray, warray, funcTypeInt, dataSize, rejTechInt, whichPrior, priorsParamsArray, hasPriorsArray, bar_guess)
    else:
        result = rs.requestHandlerUnWeighted(xarray, yarray, guessarray, funcTypeInt, dataSize, rejTechInt, whichPrior, priorsParamsArray, hasPriorsArray, bar_guess)

    #there will now be a text file of data in the current directory, with the final calculated function parameters, and binary numbers specifying which data was rejected, and which wasn't.


    parameterslist = [] #strings
    flagslist = [] #strings

    for i in range(paramCount):
        parameterslist.append(str(result[i])) #converts to string for usage via json
    for j in range(dataSize):
        flagslist.append(str(int(result[paramCount + 1 + j]))) #double -> int -> string

    #add bar value to end of parameterslist
    parameterslist.append(str(result[-1]))

    #creates lists of the nonrejected and the rejected data

    nonRejectedList = [] #formatted nonweighted as ['x1 y1', 'x2 y2', 'x3 y3'...] strings or weighted as ['x1 y1 w1', 'x2 y2 w2', ...]
    rejectedList = []
    allDataList = [] #formatted nonweighted as ['x1 y1 f1', 'x2 y2 f2', 'x3 y3 f3'...] strings or weighted as ['x1 y1 w1 f1', 'x2 y2 w2 f2', ...]

    if (weighted == False):
        for k in range(dataSize):
            if flagslist[k] == '1': #flagged true
                nonRejectedList.append(str(xValues[k]) + '\t' + str(yValues[k]))
                allDataList.append(str(xValues[k]) + '\t' + str(yValues[k]) + '\t' + 'True')
            elif flagslist[k] == '0': #flagged false
                rejectedList.append(str(xValues[k]) + '\t' + str(yValues[k]))
                allDataList.append(str(xValues[k]) + '\t' + str(yValues[k]) + '\t' + 'False')
    elif (weighted == True):
        for k in range(dataSize):
            if flagslist[k] == '1': #flagged true
                nonRejectedList.append(str(xValues[k]) + '\t' + str(yValues[k]) + '\t' + str(wValues[k]))
                allDataList.append(str(xValues[k]) + '\t' + str(yValues[k]) + '\t' + str(wValues[k]) + '\t' + 'True')
            elif flagslist[k] == '0': #flagged false
                rejectedList.append(str(xValues[k]) + '\t' + str(yValues[k]) + '\t' + str(wValues[k]))
                allDataList.append(str(xValues[k]) + '\t' + str(yValues[k]) + '\t' + str(wValues[k]) + '\t' + 'False')





    #finally, formats the data to be sent via json

    nonRejectedData = '\n'.join(nonRejectedList)
    rejectedData = '\n'.join(rejectedList)
    allData = '\n'.join(allDataList)
    finalParametersData = '\n'.join(parameterslist)

    print '=-=-=-=-=-='
    jsonData = dict(allData=allData, nonRejectedData=nonRejectedData, rejectedData=rejectedData, finalParametersData = finalParametersData)
    #from gluon.serializers import json
    return json(jsonData)

def callTRK():
    import TRKwebpageutils as TRK
    print '----------------------------------'
    mainData = request.vars.mainData
    guess = request.vars.guessData
    weighted = stringToBool(request.vars.weighted)
    functionType = request.vars.functionType
    priors = request.vars.priorsData
    hasPriors = stringToBool(request.vars.hasPriors)
    findPivot = stringToBool(request.vars.findPivot)
    optimizeScale = stringToBool(request.vars.optimizeScale)
    findUncertainties = stringToBool(request.vars.findUncertainties)
    scale = request.vars.scale
    print 'request received'

    scale = float(scale)

    #data and guess request strings to python 

    mainData = mainData.split("\n")
    guess = guess.split("\n")
    priorsData = priors.split("\n")

    M = len(guess) - 2
    N = len(mainData)

    x = []
    sx = []
    y = []
    sy = []
    w = []

    allparamsguess = []

    for i in mainData:
        print i
        if len(i) == 0:
            continue
        i.replace("\t"," ") #replaces any potential tab delimiter with a space
        dataPoint = i.split() #list made from the row of inputData

        x.append(float(dataPoint[0]))
        sx.append(float(dataPoint[1]))
        y.append(float(dataPoint[2]))
        sy.append(float(dataPoint[3]))

        if weighted:
            w.append(float(dataPoint[4]))
        else:
            w.append(1.0)

    for i in guess:
        if len(i) == 0:
            continue
        allparamsguess.append(float(i))

    #priors request strings to python 

    lb_values = []
    ub_values = []
    mu_values = []
    dev_values = []

    hasLb = []
    hasUb = []
    hasMu = []
    hasDev = []

    if hasPriors: #creates list of each of the four prior params
        for i in priorsData:
            if len(i) == 0:
                continue
            i.replace("\t"," ") #replaces any potential tab delimiter with a space
            dataPoint = i.split() #list made from the row of inputData

            iList = []

            lb_s = dataPoint[0] #individual strings
            ub_s = dataPoint[1]
            mu_s = dataPoint[2]
            dev_s = dataPoint[3]

            if lb_s == "x" or lb_s == "X" :
                hasLb.append(0)
                lb_values.append(0)
            else:
                hasLb.append(1)
                lb_values.append(float(lb_s))

            if ub_s == "x" or ub_s == "X":
                hasUb.append(0)
                ub_values.append(0)
            else:
                hasUb.append(1)
                ub_values.append(float(ub_s))

            if mu_s == "x" or mu_s == "X":
                hasMu.append(0)
                mu_values.append(0)
            else:
                hasMu.append(1)
                mu_values.append(float(mu_s))

            if dev_s == "x" or dev_s == "X":
                hasDev.append(0)
                dev_values.append(0)
            else:
                hasDev.append(1)
                dev_values.append(float(dev_s))

    onlyBounded = False
    onlyGaussian = False
    both = False

    whichPrior = 0

    CountArr = [0, 0, 0, 0]
    hasArray = [hasLb, hasUb, hasMu, hasDev]

    ct = 0
    for list in hasArray:
        for i in list:
            if i==1:
                CountArr[ct] += 1
        ct += 1

    if (CountArr[0] > 0 or CountArr[1] > 0) and (CountArr[2] == 0 and CountArr[3] == 0):
        onlyBounded = True
        whichPrior = 2
    elif (CountArr[0] == 0 and CountArr[1] == 0) and (CountArr[2] > 0 and CountArr[3] > 0):
        onlyGaussian = True
        whichPrior = 1
    elif (CountArr[0] > 0 or CountArr[1] > 0) and (CountArr[2] > 0 and CountArr[3] > 0):
        both = True
        whichPrior = 3

    #create int that defines the func type
    funcTypeInt = 0

    if functionType == 'Linear':
        funcTypeInt = 1
    elif functionType == 'Quadratic':
        funcTypeInt = 2
    elif functionType == 'Cubic':
        funcTypeInt = 3
    elif functionType == 'PowerLaw':
        funcTypeInt = 4
    elif functionType == 'Exponential':
        funcTypeInt = 5
    elif functionType == 'Logarithmic':
        funcTypeInt = 6
    else:
        return

    xVec = TRK.DoubleVector(N)
    sxVec = TRK.DoubleVector(N)
    yVec = TRK.DoubleVector(N)
    syVec = TRK.DoubleVector(N)
    wVec = TRK.DoubleVector(N)

    priorsParamsVec = TRK.DoubleVector(M * 4) #only gaussian or only bounded will just ust the first half of this array
    hasPriorsVec = TRK.IntVector(M * 4)

    allparamsguessVec = TRK.DoubleVector(M+2)

    for i in range(N):
        xVec[i] = x[i]
        sxVec[i] = sx[i]
        yVec[i] = y[i]
        syVec[i] = sy[i]
        wVec[i] = w[i]

    for j in range(M+2):
        allparamsguessVec[j] = allparamsguess[j]

    if onlyBounded:
        for j in range(0,paramCount):
            priorsParamsVec[2*j] = lb_values[j]
            priorsParamsVec[2*j+1] = ub_values[j]

            hasPriorsVec[2*j] = hasLb[j]
            hasPriorsVec[2*j+1] = hasUb[j]

    elif onlyGaussian:
        for j in range(0,paramCount):
            priorsParamsVec[2*j] = mu_values[j]
            priorsParamsVec[2*j+1] = dev_values[j]

            hasPriorsVec[2*j] = hasMu[j]
            hasPriorsVec[2*j+1] = hasDev[j]

    elif both:
        for j in range(0,paramCount):
            priorsParamsVec[2*j] = mu_values[j]
            priorsParamsVec[2*j+1] = dev_values[j]

            hasPriorsVec[2*j] = hasMu[j]
            hasPriorsVec[2*j+1] = hasDev[j]

        for j in range(0,paramCount):
            priorsParamsVec[paramCount*2 + 2*j] = lb_values[j]
            priorsParamsVec[paramCount*2 + 2*j+1] = ub_values[j]

            hasPriorsVec[paramCount*2 + 2*j] = hasLb[j]
            hasPriorsVec[paramCount*2 + 2*j+1] = hasUb[j]

    result = TRK.DoubleVector()        
    print 'performing TRK'

    priorsCheckInt = 0
    pivotCheckInt = 0
    opScaleInt = 0
    doMCMCInt = 0

    if hasPriors:
        priorsCheckInt = 1
    if findPivot:
        pivotCheckInt = 1
    if optimizeScale:
        opScaleInt = 1
    if findUncertainties:
        doMCMCInt = 1

    result = TRK.requestHandler(funcTypeInt, xVec, yVec, wVec, sxVec, syVec, allparamsguessVec, N, pivotCheckInt, priorsCheckInt, priorsParamsVec, hasPriorsVec, opScaleInt, doMCMCInt, scale)
    #result = {best fit params, slop, - 1 2 3, + 1 2 3 sigmas, s0, a, b, pivot, bincount1, bincount2 ... , hist1, edges1, hist2, edges2 ...

    
    """
    needed outputs:
    scalesData = "s0\na\nb"
    pivotData = str(float(pivot))
    paramHistogramData = "hist11 hist12 ... hist1n\tedges11 edges12 ... edges1n+1\nhist21 hist 22 ..." = "[hist1]\t[edges1]\n[hist2]..."
    bestFit_123SigmasData = "-1 -2 -3\t1 2 3\n-1 -2 -3\t1 2 3\n..." where +/- 1, 2, 3 are +/- 1, 2, 3 standard deviations of each successive param
    finalParametersData = "a0\na1\na2 ... \nslopx\nslopy"

    """

    #best fit params

    finalParametersData = []

    for j in range(M + 2):
        finalParametersData.append(result[j])

    finalParametersData = [str(i) for i in finalParametersData]
    finalParametersData = "\n".join(finalParametersData)

    #best fit uncertainties

    bestFit_123SigmasData = []

    for j in range(M):
        minusSigmas = []
        plusSigmas = []
        for i in range(3):
            minusSigmas.append(str(round_to_n(float(result[M + 2 + i + j*6]), 2)))
        for i in range(3,6):
            plusSigmas.append(str(round_to_n(float(result[M + 2 + i + j*6]), 2)))
        minusSigmas = " ".join(minusSigmas)
        plusSigmas = " ".join(plusSigmas)


        bestFit_123SigmasData.append(minusSigmas + "\t" + plusSigmas)
        
    bestFit_123SigmasData = "\n".join(bestFit_123SigmasData)


    #scale data
    scales = [result[M + 2 + 6 * M], result[M + 2 + 6 * M + 1], result[M + 2 + 6 * M + 2]]
    scalesData = "\n".join([str(s) for s in scales])

    #pivot data

    pivotData = str(result[M + 2 + 6 * M + 3])

    #histogram data

    bincounts = []

    for i in range(M):
        bincounts.append(int(result[M + 2 + 6 * M + 4 + i]))

    totalBins = sum(bincounts)
    firstHistIndex = M + 2 + 6 * M + 4 + M

    paramHistogramData = []

    for i in range(M):
        hist = []
        edge = []
        bincount = bincounts[i]

        for k in range(bincount):
            hist.append(result[firstHistIndex + 2*sum(bincounts[0:i]) + i + k])
            edge.append(result[firstHistIndex + 2*sum(bincounts[0:i]) + i + k + bincount])
        
        edge.append(result[firstHistIndex + 2*sum(bincounts[0:i]) + i + bincount + bincount])

        hist = [str(i) for i in hist]
        edge = [str(i) for i in edge]

        hist = " ".join(hist)
        edge = " ".join(edge)

        paramHistogramData.append(hist + "\t" + edge)


    paramHistogramData = "\n".join(paramHistogramData)

    print '=-=-=-=-=-='
    jsonData = dict(finalParametersData = finalParametersData, bestFit_123SigmasData = bestFit_123SigmasData, paramHistogramData = paramHistogramData, pivotData = pivotData, scales = scalesData, )
    from gluon.serializers import json
    return json(jsonData)


def value():
    return dict(message=T("Hello World"))

def functional():
    return dict(message=T("Hello World"))

def trk():
    return dict(message=T("Hello World"))

def index():
    return dict(message=T("Hello World"))

def downloads():
    return dict(message=T("Hello World"))

def contact():
    return dict(message=T("Hello World"))
