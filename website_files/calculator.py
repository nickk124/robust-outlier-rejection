import RCRLib as R
import RCRSWIGFull as rs
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
    print '----------------------------------'
    #if not validateInputData(request.vars.inputData): raise Http(400, 'Must provide data')
    #validateInputData(request.vars.inputData)

    inputData = request.vars.inputData #data received
    symDist = stringToBool(request.vars.symDist) #necessary bools received
    equalCont = stringToBool(request.vars.equalCont)
    oneSidedCont = stringToBool(request.vars.oneSidedCont)
    inBetweenCont = stringToBool(request.vars.inBetweenCont)
    weighted = stringToBool(request.vars.weighted)

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
    print '----------------------------------'

    mainData = request.vars.mainData
    xData = request.vars.xData
    guess = request.vars.guessData
    symDist = stringToBool(request.vars.symDist)
    equalCont = stringToBool(request.vars.equalCont)
    oneSidedCont = stringToBool(request.vars.oneSidedCont)
    inBetweenCont = stringToBool(request.vars.inBetweenCont)
    weighted = stringToBool(request.vars.weighted)
    functionType = request.vars.functionType
    print 'request received'

    mainData = mainData.split("\n") #puts the input data into a list
    xData = xData.split("\n")
    guess = guess.split("\n")

    yValues = []
    xValues = []
    sigma_yValues = []
    wValues = []
    guessValues = []

    for i in mainData:
        if len(i) == 0:
            continue
        i.replace("\t"," ") #replaces any potential tab delimiter with a space
        dataPoint = i.split() #list made from the row of inputData
        if weighted:
            yValues.append(float(dataPoint[0]))
            sigma_yValues.append(float(dataPoint[1]))
            wValues.append(float(dataPoint[2]))
        else:
            yValues.append(float(dataPoint[0]))
            sigma_yValues.append(float(dataPoint[1]))

    for i in xData: # doesn't account for 3D case
        if len(i) == 0:
            continue
        xValues.append(float(i))

    for i in guess:
        if len(i) == 0:
            continue
        guessValues.append(float(i))

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
    rejTechInt = 0

    if symDist:
        if equalCont:
            rejTechInt = 1
        elif oneSidedCont:
            rejTechInt = 2
        elif inBetweenCont:
            rejTechInt = 3
    else:
        rejTechInt = 4

    #converts the data into c++ std::vectors
    dataSize = len(xValues) #number of datapoints
    paramCount = len(guessValues) #number of parameters in the function

    xarray = rs.DoubleVector(dataSize)
    yarray = rs.DoubleVector(dataSize)
    sigmayarray = rs.DoubleVector(dataSize)
    guessarray = rs.DoubleVector(paramCount)

    warray = rs.DoubleVector(dataSize) #if nonweighted, will just be an empty array.

    for i in range(0,dataSize):
        xarray[i] = xValues[i]
        yarray[i] = yValues[i]
        sigmayarray[i] = sigma_yValues[i]

    for j in range(0,paramCount):
        guessarray[j] = guessValues[j]

    if weighted: #creates weight vector if needed
        for k in range(0,dataSize):
            warray[k] = wValues[k]

    """
    xarray_np = np.array(xValues, dtype=np.float64)
    yarray_np = np.array(yValues, dtype=np.float64)
    sigmayarray_np = np.array(sigma_yValues, dtype=np.float64)
    guessarray_np = np.array(guessValues, dtype=np.float64)

    warray_np = np.array(wValues, dtype=np.float64)
    """
    """
    testx = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], dtype=np.float64)
    testy = np.array([2.0, 3.0, 400.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0], dtype=np.float64)
    testsigmay = np.array([0.1, 0.2, 0.4, 0.2, 0.3, 0.5, 0.10, 0.11, 0.6, 0.2], dtype=np.float64)
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
        result = rs.requestHandlerWeighted(xarray, yarray, sigmayarray, guessarray, warray, funcTypeInt, dataSize, rejTechInt)
    else:
        result = rs.requestHandlerUnWeighted(xarray, yarray, sigmayarray, guessarray, funcTypeInt, dataSize, rejTechInt)

    #there will now be a text file of data in the current directory, with the final calculated function parameters, and binary numbers specifying which data was rejected, and which wasn't.


    parameterslist = [] #strings
    flagslist = [] #strings

    for i in range(paramCount):
        parameterslist.append(str(result[i])) #converts to string for usage via json
    for j in range(dataSize):
        flagslist.append(str(int(result[paramCount + 1 + j]))) #double -> int -> string

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
    jsonData = dict(allData=allData, nonRejectedData=nonRejectedData, rejectedData=rejectedData, finalParametersData = finalParametersData);
    #from gluon.serializers import json
    return json(jsonData)

def value():
    return dict(message=T("Hello World"))

def functional():
    return dict(message=T("Hello World"))

def index():
    return dict(message=T("Hello World"))

def downloads():
    return dict(message=T("Hello World"))

def contact():
    return dict(message=T("Hello World"))
