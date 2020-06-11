//#include "stdafx.h"
#include "RCR.h"

namespace RCRLib {

    //Non-member Functions
    bool RCR::unityTablesLoaded = false;
    std::vector<double> RCR::ESUnity;
    std::vector<double> RCR::SSUnity;
    std::vector<double> RCR::LSUnity;
    std::vector<double> RCR::ESDLUnityCF;
    std::vector<double> RCR::LSDLUnityCF;
    std::vector<double> RCR::LS68UnityCF;
    std::vector<double> RCR::SSDLUnityCF;
    std::vector<std::vector<double> > RCR::SSConstants;


    //manipulation tools:
    void swap(int a, int b, std::vector<double> &y)
    {
        double tmp;
        tmp = y[a];
        y[a] = y[b];
        y[b] = tmp;
    }
    void swap(int a, int b, std::vector<int> &y)
    {
        int tmp;
        tmp = y[a];
        y[a] = y[b];
        y[b] = tmp;
    }
    void QS(int left, int right, std::vector<double> &y)
    {
        int i = left, j = right;
        double pivot = y[(left + right) / 2];

        while (i <= j)
        {
            while (y[i] < pivot)
            {
                i++;
            }
            while (y[j] > pivot)
            {
                j--;
            }
            if (i <= j)
            {
                swap(i, j, y);
                i++;
                j--;
            }
        }

        if (left < j)
        {
            QS(left, j, y);
        }
        if (i < right)
        {
            QS(i, right, y);
        }
    }
    void QS(int left, int right, std::vector<double> &w, std::vector<double> &y)
    {
        int i = left, j = right;
        double pivot = y[(left + right) / 2];

        while (i <= j)
        {
            while (y[i] < pivot)
            {
                i++;
            }
            while (y[j] > pivot)
            {
                j--;
            }
            if (i <= j)
            {
                swap(i, j, y);
                swap(i, j, w);
                i++;
                j--;
            }
        }

        if (left < j)
        {
            QS(left, j, w, y);
        }
        if (i < right)
        {
            QS(i, right, w, y);
        }
    }
    void QS(int left, int right, std::vector<int> &w, std::vector<double> &y)
    {
        int i = left, j = right;
        double pivot = y[(left + right) / 2];

        while (i <= j)
        {
            while (y[i] < pivot)
            {
                i++;
            }
            while (y[j] > pivot)
            {
                j--;
            }
            if (i <= j)
            {
                swap(i, j, y);
                swap(i, j, w);
                i++;
                j--;
            }
        }

        if (left < j)
        {
            QS(left, j, w, y);
        }
        if (i < right)
        {
            QS(i, right, w, y);
        }
    }
    void sort(std::vector<double> &y)
    {
        QS(0, (int) y.size() - 1, y);
    }
    void sort(std::vector<double> &w, std::vector<double> &y)
    {
        QS(0, (int) y.size() - 1, w, y);
    }
    void sort(std::vector<int> &w, std::vector<double> &y)
    {
        QS(0, (int) y.size() - 1, w, y);
    }


    //FN & CF Models:
    double getCFRatio(std::vector<double> &w)
    {
        int size = (int) w.size();
        double mean = 0, stDev = 0;
        for (int i = 0; i < size; i++)
        {
            mean += w[i];
        }
        mean = mean / size;
        for (int i = 0; i < size; i++)
        {
            stDev += pow(w[i] - mean, 2);
        }
        stDev = sqrt(stDev / (size - 1));
        return stDev / mean;
    }
    double getFNRatio(std::vector<double> &x, std::vector<double> &w)
    {
        int counter = 0;
        double mean = 0, stDev = 0;
        while (x[counter] < 1.0)
        {
            mean += w[counter];
            counter++;
        }
        mean = mean / counter;
        for (int i = 0; i < counter; i++)
        {
            stDev += pow(w[i] - mean, 2);
        }
        stDev = sqrt(stDev / (counter - 1));
        return stDev / mean;
    }


    //statistics tools
    double inverf(double x)
    {
        //x = x*x;
        //return sqrt(2.0) * sqrt(-2.0 / (PI*inverfMult) - log(1 - x) / 2.0 + sqrt(pow((2.0 / (PI*inverfMult) + log(1 - x) / 2.0), 2) - (1.0 / inverfMult)*log(1 - x)));

        x = log(1 - x*x);
        return squareRootOf2 * sqrt(-2.0 / (PI*inverfMult) - x * .5 + sqrt(pow((2.0 / (PI*inverfMult) + x * .5), 2) - x / inverfMult));

    }
    double erfcCustom(double x)
    {
        x = x / sqrt(2.0);
        return 1.0 / pow((1 + x*(.0705230784 + x*(.0422820123 + x*(.0092705272 + x*(.0001520143 + x*(.0002765672 + .0000430638*x)))))), 16);
    }
    double getDiff(double mu, double datum)
    {
        return std::abs(datum - mu);
    }
    //bool isEqual(double x, double y, double maxRelativeError = .000000001, double maxAbsoluteError = DBL_MIN)
    bool isEqual(double x, double y, double maxRelativeError = .00000001, double maxAbsoluteError = DBL_MIN)// .000001; .0000001;.00000001
    {
        if (std::abs(x - y) < maxAbsoluteError)
        {
            return true;
        }
        double relativeError = (std::abs(y) > std::abs(x) ? std::abs((x - y) / y) : std::abs((x - y) / x));
        if (relativeError <= maxRelativeError)
        {
            return true;
        }
        return false;
    }
    bool distinctValuesCheck(std::vector<bool> &flags, std::vector<double>& y)
    {
        int i = 0, size = (int) y.size();
        double a, b;

        while (!flags[i])
        {
            i++;
        }

        a = y[i];
        b = a;
        while (i < size)
        {
            if (!isEqual(y[i], a) && flags[i])
            {
                b = y[i];
                break;
            }
            i++;
        }

        while (i < size)
        {
            if (!isEqual(y[i], a) && !isEqual(y[i], b) && flags[i])
                //if (y[i] != a && y[i] != b && flags[i])
            {
                return true;
            }
            i++;
        }
        return false;

    }
    bool distinctValuesCheck(int paramCount, std::vector<bool> &flags, std::vector<double>& y)
    {
        bool distinctValue;
        std::vector<double> distincts;
        for (size_t index = 0; index < y.size(); index++)
        {
            if (flags[index])
            {
                distinctValue = true;
                for (size_t j = 0; j < distincts.size(); j++)
                {
                    if (isEqual(distincts[j], y[index]))
                    {
                        distinctValue = false;
                    }
                }
                if (distinctValue)
                {
                    distincts.push_back(y[index]);
                }
            }
            if (distincts.size() > paramCount)
            {
                return true;
            }
        }
        return false;
    }

    //algorithm tools:
    void setTrueVec(std::vector<bool> &flags, std::vector<int> &indices, std::vector<double> &w, std::vector<double> &y, std::vector<double> &trueW, std::vector<double> &trueY)
    {
        int trueCount = 0, currentIndex;
        std::vector<int> indicesVec;
        std::vector<double> trueWVec, trueYVec;
        for (size_t i = 0; i < flags.size(); i++)
        {
            if (flags[i])
            {
                trueCount += 1;
            }
        }
        trueWVec.resize(trueCount);
        trueYVec.resize(trueCount);
        indicesVec.resize(trueCount);
        currentIndex = 0;
        for (size_t i = 0; i < flags.size(); i++)
        {
            if (flags[i])
            {
                trueWVec[currentIndex] = (w[i]);
                trueYVec[currentIndex] = (y[i]);
                indicesVec[currentIndex] = (int) i;

                currentIndex += 1;
            }
        }
        trueY = trueYVec;
        trueW = trueWVec;
        indices = indicesVec;
    }
    void setTrueVec(std::vector<bool> &flags, std::vector<int> &indices, std::vector<double> &y, std::vector<double> &trueY)
    {
        int trueCount = 0, currentIndex;
        std::vector<int> indicesVec;
        std::vector<double> trueYVec;
        for (size_t i = 0; i < flags.size(); i++)
        {
            if (flags[i])
            {
                trueCount += 1;
            }
        }
        trueYVec.resize(trueCount);
        indicesVec.resize(trueCount);
        currentIndex = 0;
        for (size_t i = 0; i < flags.size(); i++)
        {
            if (flags[i])
            {
                trueYVec[currentIndex] = (y[i]);
                indicesVec[currentIndex] = (int) i;
                currentIndex += 1;
            }
        }
        trueY = trueYVec;
        indices = indicesVec;

    }
    int countAmountLessThanOne(std::vector<double> &x)
    {
        if (x.size() == 1)
        {
            if (x[0] >= 1)
            {
                return 0;
            }
            else
            {
                return 1;
            }
        }
        int amountUnderOne = 0;
        while (x[amountUnderOne] < 1.0)
        {
            amountUnderOne++;
        }
        return amountUnderOne;
    }
    int binarySearch(bool searchUp, int minimumIndex, double toFind, std::vector<double> &toSearch)
    {
        int low, high, midPoint = 0, lowIn = -1, highIn = -1;
        if (searchUp)
        {
            low = minimumIndex;
            high = (int) toSearch.size();
        }
        else
        {
            low = 0;
            high = minimumIndex;
        }
        while (low != lowIn || high != highIn)
        {
            lowIn = low;
            highIn = high;
            midPoint = (int)(low + (high - low) / 2.0);

            if (isEqual(toFind, toSearch[midPoint]))
            {
                low = midPoint;
                high = midPoint;
                
            }
            else if (toFind > toSearch[midPoint])
            {
                low = midPoint;
            }
            else if (toFind < toSearch[midPoint])
            {
                high = midPoint;
            }
            
        }
        if (searchUp)
        {
            return low;
        }
        else
        {
            return high;
        }
    }
    double min(double a, double b)
    {
        return (a < b ? a : b);
    }
    double max(double a, double b)
    {
        return (a > b ? a : b);
    }
    std::vector<double> getXVec(int size, std::vector<double> &w)
    {
        double wSum = w[0], sSum = .682689 * wSum;

        std::vector<double> x;
        x.resize(size);

        for (int i = 1; i < size; i++)
        {
            wSum += w[i];
        }

        wSum = 1.0 / wSum;

        x[0] = inverf(sSum * wSum);

        for (int i = 1; i < size; i++)
        {
            sSum += .317311 * w[i - 1] + .682689 * w[i];
            x[i] = inverf(sSum * wSum);
        }

        return x;
    }
    std::vector<double> getXVec(int size)
    {
        double sSum = .682689;
        std::vector<double> x;
        x.resize(size);
        x[0] = inverf(sSum / size);
        for (int i = 1; i < size; i++)
        {
            sSum += 1;
            x[i] = inverf(sSum / size);
        }
        return x;
    }


    //mu calculations:
    double getMean(int trueCount, std::vector<double> &y)
    {
        double top = 0, bottom = trueCount;
        for (int i = 0; i < trueCount; i++)
        {
            top += y[i];
        }
        return top / bottom;
    }
    double getMean(int trueCount, std::vector<double> &w, std::vector<double> &y)
    {
        double top = 0, bottom = 0;
        for (int i = 0; i < trueCount; i++)
        {
            top += w[i] * y[i];
            bottom += w[i];
        }

        return top / bottom;
    }
    double getMedian(int trueCount, std::vector<double> &w, std::vector<double> &y)
    {
        size_t sumCounter = 0;
        double median = 0, totalSum = 0, runningSum = 0;
        for (int i = 0; i < trueCount; i++)
        {
            totalSum += w[i];
        }
        if (trueCount > 1)
        {
            runningSum = w[sumCounter] * .5;
            while (runningSum < .5*totalSum)
            {
                sumCounter++;
                runningSum += w[sumCounter - 1] * .5 + w[sumCounter] * .5;
            }
            if (sumCounter == 0)
            {
                median = y[0];
                std::cout << median << std::endl;
            }
            else
            {
                median = y[sumCounter - 1] + (.5*totalSum - (runningSum - (w[sumCounter - 1] * .5 + w[sumCounter] * .5))) / (w[sumCounter - 1] * .5 + w[sumCounter] * .5)*(y[sumCounter] - y[sumCounter - 1]);
                std::cout << median << std::endl;
            }
        }
        else
        {
            median = y[0];
            std::cout << median << std::endl;
        }
        return median;
    }
    double getMedian(std::vector<double> &y)
    {
        int high = (int)(floor(y.size() / 2));
        int low = high - 1;
        double runningSum = 0, median = 0;
        double totalSum = y.size();
        if (y.size() > 1)
        {
            if (y.size() % 2 == 0)
            {
                runningSum = y.size() / 2.0 + .5;
            }
            else
            {
                runningSum = y.size() / 2.0;
            }
            median = y[low] + (.5*totalSum - runningSum + 1.0)* (y[high] - y[low]);
        }

        else
        {
            median = y[0];
        }
        return median;

    }
    double getMode(int trueCount, std::vector<double> &w, std::vector<double> &y)
    {
        int k, lowerLimit = 0, upperLimit = trueCount - 1, lowerLimitIn = -1, upperLimitIn = -1, finalLower = -1, finalUpper = -1, size;
        double halfWeightSum = 0, sSum, total, minDist = 999999;
        std::vector<double> sVec;

        while (lowerLimit != lowerLimitIn || upperLimit != upperLimitIn)
        {
            //std::cout<< lowerLimit << "\t" << upperLimit << "\n";
            lowerLimitIn = lowerLimit;
            upperLimitIn = upperLimit;
            size = upperLimit - lowerLimit + 1;
            minDist = 999999;
            halfWeightSum = 0;
            for (int i = lowerLimit; i < upperLimit + 1; i++)
            {
                halfWeightSum += w[i];
            }
            halfWeightSum *= .5;

            sVec.resize(size, 0.0);
            sSum = .5 * w[lowerLimit];
            sVec[0] = sSum;
            for (int i = lowerLimit + 1; i < lowerLimit + size; i++)
            {
                sSum += w[i - 1] * .5 + w[i] * .5;
                sVec[i - lowerLimit] = sSum;
            }

            for (size_t i = 0; i < sVec.size(); i++)
            {
                if ((sVec[i] < halfWeightSum) || isEqual(sVec[i],halfWeightSum))
                {
                    total = sVec[i] + halfWeightSum;
                    k = (int) i; // was 0
                    while (k < sVec.size() && ((sVec[k] < total) || isEqual(sVec[k], total)))
                    {
                        k++;
                    }
                    k--;
                    total = std::abs(y[k + lowerLimit] - y[i + lowerLimit]);

                    
                    if (isEqual(total,minDist))
                    {
                        finalLower = (int)(min(finalLower, i + lowerLimit));
                        finalUpper = (int)(max(finalUpper, k + lowerLimit));
                    }
                    else if (total < minDist)
                    {
                        minDist = total;
                        finalLower = (int) i + lowerLimit;
                        finalUpper = k + lowerLimit;
                    }
                }
                if ((sVec[i] > halfWeightSum) || isEqual(sVec[i], halfWeightSum))
                {
                    total = sVec[i] - halfWeightSum;
                    k = (int) i; // was svec.size() - 1
                    while (k > -1 && ((sVec[k] > total) || isEqual(sVec[k], total)))
                    {
                        k--;
                    }
                    k++;
                    total = std::abs(y[i + lowerLimit] - y[k + lowerLimit]);

                    
                    if (isEqual(total,minDist))
                    {
                        finalLower = (int)(min(finalLower, k + lowerLimit));
                        finalUpper = (int)(max(finalUpper, i + lowerLimit));
                    }
                    else if (total < minDist)
                    {
                        minDist = total;
                        finalLower = k + lowerLimit;
                        finalUpper = (int) i + lowerLimit;
                    }
                }
            }

            lowerLimit = finalLower;
            upperLimit = finalUpper;

            sVec.clear();
        }

        std::vector<double> newValues(y.begin() + lowerLimit, y.begin() + upperLimit + 1);
        std::vector<double> newWeights(w.begin() + lowerLimit, w.begin() + upperLimit + 1);
        return getMedian((int) newWeights.size(), newWeights, newValues);
    }
    double getMode(int trueCount, std::vector<double> &y)
    {
        int k, lowerLimit = 0, upperLimit = trueCount - 1, lowerLimitIn = -1, upperLimitIn = -1, finalLower = -1, finalUpper = -1, size;
        double halfWeightSum = 0, sSum, total, minDist = 999999;
        std::vector<double> sVec;
        while (lowerLimit != lowerLimitIn || upperLimit != upperLimitIn)
        {
            lowerLimitIn = lowerLimit;
            upperLimitIn = upperLimit;
            size = upperLimit - lowerLimit + 1;
            minDist = 999999;
            halfWeightSum = 0;
            halfWeightSum = size;

            halfWeightSum *= 0.5;
            sVec.resize(size, 0.0);
            sSum = .5;
            sVec[0] = sSum;
            for (int i = lowerLimit + 1; i < lowerLimit + size; i++)
            {
                sSum += 1;
                sVec[i - lowerLimit] = sSum;
            }
            for (size_t i = 0; i < sVec.size(); i++)
            {
                if ((sVec[i] < halfWeightSum) || isEqual(sVec[i],halfWeightSum))
                {
                    total = sVec[i] + halfWeightSum;
                    /*k = 0;
                    while (k < sVec.size() && sVec[k] <= total)
                    {
                    k++;
                    }
                    k--;*/
                    k = binarySearch(true, (int) i, total, sVec);
                    total = std::abs(y[k + lowerLimit] - y[i + lowerLimit]);


                    if (isEqual(total,minDist))
                    {
                        finalLower = (int)(min(finalLower, i + lowerLimit));
                        finalUpper = (int)(max(finalUpper, k + lowerLimit));
                    }
                    else if (total < minDist)
                    {
                        minDist = total;
                        finalLower = (int) i + lowerLimit;
                        finalUpper = k + lowerLimit;
                    }
                }
                if ((sVec[i] > halfWeightSum) || isEqual(sVec[i], halfWeightSum))
                {
                    total = sVec[i] - halfWeightSum;
                    /*k = sVec.size() - 1;
                    while (k > -1 && sVec[k] >= total)
                    {
                    k--;
                    }
                    k++;*/
                    k = binarySearch(false, (int) i, total, sVec);

                    total = std::abs(y[i + lowerLimit] - y[k + lowerLimit]);

                    if (isEqual(total, minDist))
                    {
                        finalLower = (int)(min(finalLower, k + lowerLimit));
                        finalUpper = (int)(max(finalUpper, i + lowerLimit));
                    }
                    else if (total < minDist)
                    {
                        minDist = total;
                        finalLower = k + lowerLimit;
                        finalUpper = (int) i + lowerLimit;
                    }

                }
            }
            lowerLimit = finalLower;
            upperLimit = finalUpper;
            sVec.clear();
        }
        std::vector<double> newValues(y.begin() + lowerLimit, y.begin() + upperLimit + 1);
        return getMedian(newValues);
    }


    //sigma calculations:
    double getOriginFixedRegressionLine(int start, int end, std::vector<double> &w, std::vector<double> &x, std::vector<double> &y)
    {
        double wxProd, prodX = 0, prodY = 0;
        for (int i = start; i < end; i++)
        {
            wxProd = w[i] * x[i];
            prodX += wxProd * x[i];
            prodY += wxProd * y[i];
        }
        return prodY / prodX;
    }
    double getOriginFixedRegressionLine(int start, int end, std::vector<double> &x, std::vector<double> &y)
    {
        double prodSum = 0, xSquares = 0;
        for (int i = start; i < end; i++)
        {
            prodSum += x[i] * y[i];
            xSquares += x[i] * x[i];
        }
        return prodSum / xSquares;
    }
    double mFinder(int low, int high, int lastXUnderOne, int increment, std::vector<double> &w, std::vector<double> &x, std::vector<double> &y)
    {
        bool stop = false;
        int bestM = -1, mLow = -1, mHigh = -1;
        double a, b, c, d, e, f, error, xAtM, xAtIxAtMDiff, wAtI, wxProd, tau, sigma, factor, minError = DBL_MAX;
        while (!stop)
        {
            for (int m = low; m < high; m += increment)
            {
                error = 0;
                a = 0;
                b = 0;
                c = 0;
                d = 0;
                e = 0;
                f = 0;
                xAtM = x[m];
                for (int i = 0; i <= m; i++)
                {
                    wxProd = w[i] * x[i];
                    a += wxProd * x[i];
                    e += wxProd * y[i];
                }
                for (int i = m + 1; i < lastXUnderOne + 1; i++)
                {
                    xAtIxAtMDiff = (x[i] - xAtM);
                    wAtI = w[i];
                    a += xAtM * xAtM * wAtI;
                    b += xAtM * wAtI * xAtIxAtMDiff;
                    d += wAtI * xAtIxAtMDiff * xAtIxAtMDiff;
                    e += xAtM * wAtI * y[i];
                    f += wAtI * y[i] * xAtIxAtMDiff;
                }
                c = b;
                if (a != 0 && d - c*b != 0)
                {
                    tau = (f - e*c / a) / (d - c*b / a);
                    sigma = (e - tau*b) / a;
                    for (int i = 0; i <= m; i++)
                    {
                        factor = (sigma * x[i] - y[i]);
                        error += w[i] * factor*factor;
                    }
                    for (int i = m + 1; i < lastXUnderOne + 1; i++)
                    {
                        factor = sigma*xAtM + tau*(x[i] - xAtM) - y[i];
                        error += w[i] * factor*factor;
                    }
                    if (error < minError)
                    {
                        minError = error;
                        bestM = m;
                        mLow = (int)(max(bestM - increment - 1, 1));
                        mHigh = (int)(min(bestM + increment + 1, lastXUnderOne));
                    }
                }
            }
            if (increment > 1)
            {
                increment = (int)(max(floor((double)(mHigh - mLow) / 6.36), 1.0));
            }
            else
            {
                stop = true;
            }
            low = mLow;
            high = mHigh;
            minError = DBL_MAX;
        }
        if (low == high)
        {
            return low;
        }
        return bestM;
    }
    double mFinder(int low, int high, int lastXUnderOne, int increment, std::vector<double> &x, std::vector<double> &y)
    {
        bool stop = false;
        int bestM = -1, mLow = -1, mHigh = -1;
        double a, b, c, d, e, f, error, xAtM, xAtIxAtMDiff, tau, sigma, factor, minError = 999999;
        while (!stop)
        {
            for (int m = low; m < high; m += increment)
            {
                error = 0;
                a = 0;
                b = 0;
                c = 0;
                d = 0;
                e = 0;
                f = 0;
                xAtM = x[m];
                for (int i = 0; i <= m; i++)
                {
                    a += x[i] * x[i];
                    e += x[i] * y[i];
                }
                for (int i = m + 1; i < lastXUnderOne + 1; i++)
                {
                    xAtIxAtMDiff = (x[i] - xAtM);
                    a += xAtM * xAtM;
                    b += xAtM * xAtIxAtMDiff;
                    d += xAtIxAtMDiff * xAtIxAtMDiff;
                    e += xAtM * y[i];
                    f += y[i] * xAtIxAtMDiff;
                }
                c = b;
                if (a != 0 && d - c*b != 0)
                {
                    tau = (f - e*c / a) / (d - c*b / a);
                    sigma = (e - tau*b) / a;
                    for (int i = 0; i <= m; i++)
                    {
                        factor = (sigma * x[i] - y[i]);
                        error += factor*factor;
                    }
                    for (int i = m + 1; i < lastXUnderOne + 1; i++)
                    {
                        factor = sigma*xAtM + tau*(x[i] - xAtM) - y[i];
                        error += factor*factor;
                    }
                    if (error < minError)
                    {
                        minError = error;
                        bestM = m;
                        mLow = (int)(max(bestM - increment - 1, 1));
                        mHigh = (int)(min(bestM + increment + 1, lastXUnderOne));
                    }
                }
            }
            if (increment > 1)
            {
                increment = (int)(max(floor((double)(mHigh - mLow) / 6.36), 1.0));
            }
            else
            {
                stop = true;
            }
            low = mLow;
            high = mHigh;
            minError = 999999;
        }
        if (low == high)
        {
            return low;
        }
        return bestM;
    }
    double getStDev(double delta, std::vector<double> &w, std::vector<double> &y)
    {
        int size = (int) w.size();
        double top = 0, wSum = 0, wSumSq = 0, weight;
        for (int i = 0; i < size; i++)
        {
            weight = w[i];
            top += weight * y[i] * y[i];
            wSum += weight;
            wSumSq += weight*weight;
        }
        return sqrt(top / (wSum - delta*wSumSq / wSum));
    }
    double getStDev(double delta, std::vector<double> &y)
    {
        int size = (int) y.size();
        double top = 0;
        for (int i = 0; i < size; i++)
        {
            top += y[i] * y[i];
        }
        return sqrt(top / (size - delta));
    }
    double get68th(std::vector<double> &w, std::vector<double> &y)
    {
        size_t sumCounter = 0;
        double stDev = 0, totalSum = 0, runningSum; //, temp = 0, weightTemp = 0;
        sort(w, y);
        for (size_t i = 0; i < y.size(); i++)
        {
            totalSum += w[i];
        }
        if (y.size() > 1)
        {
            runningSum = w[sumCounter] * .682689;
            while (runningSum < .682689*totalSum)
            {
                sumCounter++;
                runningSum += w[sumCounter - 1] * .317311 + w[sumCounter] * .682689;
            }
            if (sumCounter == 0)
            {
                stDev = y[0];
            }
            else
            {
                stDev = y[sumCounter - 1] + (.682689*totalSum - (runningSum - (w[sumCounter - 1] * .317311 + w[sumCounter] * .682689))) / (w[sumCounter - 1] * .317311 + w[sumCounter] * .682689)*(y[sumCounter] - y[sumCounter - 1]);
            }
        }
        else
        {
            stDev = y[0];
        }

        return stDev;
    }
    double get68th(std::vector<double> &y)
    {
        size_t sumCounter = 0;
        double stDev = 0, totalSum = 0, runningSum; //, temp = 0, weightTemp = 0;
        sort(y);
        for (size_t i = 0; i < y.size(); i++)
        {
            totalSum += 1.0;
        }
        if (y.size() > 1)
        {
            runningSum = 1.0 * .682689;
            while (runningSum < .682689*totalSum)
            {
                sumCounter++;
                runningSum += 1.0 * .317311 + 1.0 * .682689;
            }
            if (sumCounter == 0)
            {
                stDev = y[0];
            }
            else
            {
                stDev = y[sumCounter - 1] + (.682689*totalSum - (runningSum - (1.0 * .317311 + 1.0 * .682689))) / (1.0 * .317311 + 1.0 * .682689)*(y[sumCounter] - y[sumCounter - 1]);
            }
        }
        else
        {
            stDev = y[0];
        }

        return stDev;

    }
    double fitSL(std::vector<double> &w, std::vector<double> &x, std::vector<double> &y)
    {
        return getOriginFixedRegressionLine(0, countAmountLessThanOne(x), w, x, y);
    }
    double fitSL(std::vector<double> &x, std::vector<double> &y)
    {
        return getOriginFixedRegressionLine(0, countAmountLessThanOne(x), x, y);
    }


    //constructors
    RCR::RCR()
    {
        this->rejectionTech = LS_MODE_DL;
        this->delta = 1.0;
        alignTechniques();
        //FunctionalForm f;
        //this->parametricModel = &f;
        loadUnityTables();
    }
    RCR::RCR(RejectionTechs rejectionTech)
    {
        this->rejectionTech = rejectionTech;
        this->delta = 1.0;
        alignTechniques();
        loadUnityTables();
        this->muType = VALUE;
        //setparametricModel(this->mF);
        //this->parametricModel = &mF;
        //revertparametricModel();
    }


    std::vector<int> getHalfSampleBounds(int lowerLimit, int upperLimit, std::vector<double> w, std::vector<double> y)
    {
        int	size = upperLimit - lowerLimit + 1, k, finalLower, finalUpper;
        double minDist = 999999, halfWeightSum = 0, sSum, total;
        std::vector<double> sVec;
        halfWeightSum = 0;
        for (int i = lowerLimit; i < upperLimit + 1; i++)
        {
            halfWeightSum += w[i];
        }
        halfWeightSum *= .5;
        sVec.resize(size, 0.0);
        sSum = .5 * w[lowerLimit];
        sVec[0] = sSum;
        for (int i = lowerLimit + 1; i < lowerLimit + size; i++)
        {
            sSum += w[i - 1] * .5 + w[i] * .5;
            sVec[i - lowerLimit] = sSum;
        }
        for (int i = 0; i < sVec.size(); i++)
        {
            if (sVec[i] < halfWeightSum || isEqual(sVec[i],halfWeightSum))
            {
                total = sVec[i] + halfWeightSum;
                k = binarySearch(true, i, total, sVec);
                total = std::abs(y[k + lowerLimit] - y[i + lowerLimit]);

                if (isEqual(total,minDist))
                {
                    finalLower = std::min(finalLower, i + lowerLimit);
                    finalUpper = std::max(finalUpper, k + lowerLimit);
                }
                else if (total < minDist)
                {
                    minDist = total;
                    finalLower = i + lowerLimit;
                    finalUpper = k + lowerLimit;
                }
            }
            if (sVec[i] > halfWeightSum || isEqual(sVec[i], halfWeightSum))
            {
                total = sVec[i] - halfWeightSum;
                k = binarySearch(false, i, total, sVec);
                total = std::abs(y[i + lowerLimit] - y[k + lowerLimit]);
                if (isEqual(total, minDist))
                {
                    finalLower = std::min(finalLower, k + lowerLimit);
                    finalUpper = std::max(finalUpper, i + lowerLimit);
                }
                else if (total < minDist)
                {
                    minDist = total;
                    finalLower = k + lowerLimit;
                    finalUpper = i + lowerLimit;
                }
            }
        }
        std::vector<int> toRet;
        toRet.push_back(finalLower);
        toRet.push_back(finalUpper);
        return toRet;
    }
    std::vector<double> get2DMode(FunctionalForm &f)
    {
        std::vector<int> xBounds, yBounds;
        std::vector<double> sortedXHold, sortedYHold, sortedWXHold, sortedWYHold;
        std::vector<double> sortedX, sortedY, sortedWX, sortedWY;
        std::vector<double> x = f.parameterSpace[0], y = f.parameterSpace[1], xW = f.weightSpace[0], yW = f.weightSpace[1];
        sortedX = x;
        sortedY = y;
        sortedWX = xW;
        sortedWY = yW;
        int xLowIndex = 0, xHighIndex = (int) sortedX.size() - 1, yLowIndex = 0, yHighIndex = (int) sortedY.size() - 1, xLowIndexIn = -1, xHighIndexIn = -1, yLowIndexIn = -1, yHighIndexIn = -1;
        double xLowValue, xHighValue, yLowValue, yHighValue;
        bool nullSet = false;
        while (!nullSet && (xLowIndex != xLowIndexIn || xHighIndex != xHighIndexIn || yLowIndex != yLowIndexIn || yHighIndex != yHighIndexIn))
        {

            xLowIndexIn = xLowIndex;
            xHighIndexIn = xHighIndex;
            yLowIndexIn = yLowIndex;
            yHighIndexIn = yHighIndex;

            sort(sortedWX, sortedX);
            xBounds = getHalfSampleBounds(0, (int) sortedX.size() - 1, sortedWX, sortedX);
            xLowIndex = xBounds[0];
            xHighIndex = xBounds[1];
            xLowValue = sortedX[xLowIndex];
            xHighValue = sortedX[xHighIndex];


            sort(sortedWY, sortedY);
            yBounds = getHalfSampleBounds(0, (int) sortedY.size() - 1, sortedWY, sortedY);
            yLowIndex = yBounds[0];
            yHighIndex = yBounds[1];
            yLowValue = sortedY[yLowIndex];
            yHighValue = sortedY[yHighIndex];

            //std::cout << xLowIndex << "\t" << xHighIndex << "\t" << yLowIndex << "\t" << yHighIndex << "\n";
            //std::cout << xLowValue << "\t" << xHighValue << "\t" << yLowValue << "\t" << yHighValue << "\n";

            sortedXHold = sortedX;
            sortedYHold = sortedY;
            sortedWXHold = sortedWX;
            sortedWYHold = sortedWY;

            sortedX.clear();
            sortedY.clear();
            sortedWX.clear();
            sortedWY.clear();
            for (size_t i = 0; i < x.size(); i++)
            {
                if (x[i] >= xLowValue && x[i] <= xHighValue && y[i] >= yLowValue && y[i] <= yHighValue)
                {
                    sortedX.push_back(x[i]);
                    sortedY.push_back(y[i]);
                    sortedWX.push_back(xW[i]);
                    sortedWY.push_back(yW[i]);
                }
            }
            nullSet = false;
            if (sortedX.size() == 0 || sortedY.size() == 0)
            {
                nullSet = true;
                sortedX = sortedXHold;
                sortedY = sortedYHold;
                sortedWX = sortedWXHold;
                sortedWY = sortedWYHold;
            }
        }
        sort(sortedWX, sortedX);
        sort(sortedWY, sortedY);
        std::vector<double> toRet;
        double m = getMedian((int) sortedX.size(), sortedWX, sortedX), b = getMedian((int) sortedY.size(), sortedWY, sortedY);
        //std::cout << "MODE        a: " << m << "\t" << "b: " << b << "\n";
        toRet.push_back(m);
        toRet.push_back(b);

        return toRet;

    }
    std::vector<double> get2DMedian(FunctionalForm &f)
    {
        std::vector<double> sortedX, sortedY, sortedWX, sortedWY, x = f.parameterSpace[0], y = f.parameterSpace[1], xW = f.weightSpace[0], yW = f.weightSpace[1];

        //adding the extra exception parameter values from buildModelSpace:
        /*
        for (int i = 0; i < f.extraParameterSpace[0].size(); i++)
        {
            x.push_back(f.extraParameterSpace[0][i]);
            y.push_back(f.extraParameterSpace[1][i]);

            xW.push_back(f.extraWeightSpace[0][i]);
            yW.push_back(f.extraWeightSpace[1][i]);
        }
        */
        sortedX = x;
        sortedY = y;

        sortedWX = xW;
        sortedWY = yW;

        sort(sortedWX, sortedX);
        sort(sortedWY, sortedY);

        std::vector<double> toRet;
        double m = getMedian((int) sortedX.size(), sortedWX, sortedX), b = getMedian((int) sortedY.size(), sortedWY, sortedY);
        //std::cout << "MEDIAN        a: " << m << "\t" << "b: " << b << "\n";
        toRet.push_back(m);
        toRet.push_back(b);

        f.meanstartingpoint = toRet; //used for the guess for the generalized mean calculation

        return toRet;
    }
    std::vector<double> get3DMode(FunctionalForm &f)
    {
        std::vector<int> xBounds, yBounds, zBounds;
        std::vector<double> sortedXHold, sortedYHold, sortedZHold, sortedXWHold, sortedYWHold, sortedZWHold;
        std::vector<double> sortedX, sortedY, sortedZ;
        std::vector<double> sortedXW, sortedYW, sortedZW;
        std::vector<double> x = f.parameterSpace[0], y = f.parameterSpace[1], z = f.parameterSpace[2];
        std::vector<double> xw = f.weightSpace[0], yw = f.weightSpace[1], zw = f.weightSpace[2];

        sortedX = x;
        sortedY = y;
        sortedZ = z;
        sortedXW = xw;
        sortedYW = yw;
        sortedZW = zw;
        int xLowIndex = 0, xHighIndex = (int) sortedX.size() - 1, yLowIndex = 0, yHighIndex = (int) sortedY.size() - 1, zLowIndex = 0, zHighIndex = (int) sortedZ.size() - 1;
        int xLowIndexIn = -1, xHighIndexIn = -1, yLowIndexIn = -1, yHighIndexIn = -1, zLowIndexIn = -1, zHighIndexIn = -1;
        double xLowValue, xHighValue, yLowValue, yHighValue, zLowValue, zHighValue;
        bool nullSet = false;
        while (!nullSet && (xLowIndex != xLowIndexIn || xHighIndex != xHighIndexIn || yLowIndex != yLowIndexIn || yHighIndex != yHighIndexIn || zLowIndex != zLowIndexIn || zHighIndex != zHighIndexIn))
        {

            xLowIndexIn = xLowIndex;
            xHighIndexIn = xHighIndex;
            yLowIndexIn = yLowIndex;
            yHighIndexIn = yHighIndex;
            zLowIndexIn = zLowIndex;
            zHighIndexIn = zHighIndex;

            sort(sortedXW, sortedX);
            xBounds = getHalfSampleBounds(0, (int) sortedX.size() - 1, sortedXW, sortedX);
            xLowIndex = xBounds[0];
            xHighIndex = xBounds[1];
            xLowValue = sortedX[xLowIndex];
            xHighValue = sortedX[xHighIndex];


            sort(sortedYW, sortedY);
            yBounds = getHalfSampleBounds(0, (int) sortedY.size() - 1, sortedYW, sortedY);
            yLowIndex = yBounds[0];
            yHighIndex = yBounds[1];
            yLowValue = sortedY[yLowIndex];
            yHighValue = sortedY[yHighIndex];

            sort(sortedZW, sortedZ);
            zBounds = getHalfSampleBounds(0, (int) sortedZ.size() - 1, sortedZW, sortedZ);
            zLowIndex = zBounds[0];
            zHighIndex = zBounds[1];
            zLowValue = sortedZ[zLowIndex];
            zHighValue = sortedZ[zHighIndex];

            //std::cout << xLowIndex << "\t" << xHighIndex << "\t" << yLowIndex << "\t" << yHighIndex << "\n";
            //std::cout << xLowValue << "\t" << xHighValue << "\t" << yLowValue << "\t" << yHighValue << "\n";

            sortedXHold = sortedX;
            sortedYHold = sortedY;
            sortedZHold = sortedZ;

            sortedXWHold = sortedXW;
            sortedYWHold = sortedYW;
            sortedZWHold = sortedZW;

            sortedX.clear();
            sortedY.clear();
            sortedZ.clear();

            sortedXW.clear();
            sortedYW.clear();
            sortedZW.clear();

            for (size_t i = 0; i < x.size(); i++)
            {
                if (x[i] >= xLowValue && x[i] <= xHighValue && y[i] >= yLowValue && y[i] <= yHighValue && z[i] >= zLowValue && z[i] <= zHighValue)
                {
                    sortedX.push_back(x[i]);
                    sortedY.push_back(y[i]);
                    sortedZ.push_back(z[i]);

                    sortedXW.push_back(xw[i]);
                    sortedYW.push_back(yw[i]);
                    sortedZW.push_back(zw[i]);

                }
            }
            nullSet = false;
            if (sortedX.size() == 0 || sortedY.size() == 0 || sortedZ.size() == 0)
            {
                nullSet = true;
                sortedX = sortedXHold;
                sortedY = sortedYHold;
                sortedZ = sortedZHold;

                sortedXW = sortedXWHold;
                sortedYW = sortedYWHold;
                sortedZW = sortedZWHold;
            }
        }
        sort(sortedXW, sortedX);
        sort(sortedYW, sortedY);
        sort(sortedZW, sortedZ);
        std::vector<double> toRet;
        double b = getMedian((int) sortedX.size(), sortedXW, sortedX), m = getMedian((int) sortedY.size(), sortedYW, sortedY), xNot = getMedian((int) sortedZ.size(), sortedZW, sortedZ);
        //std::cout << "MODE        A: " << b << "\t" << "B: " << m << "\tC: " << xNot << "\n";
        toRet.push_back(b);
        toRet.push_back(m);
        toRet.push_back(xNot);

        
        return toRet;

    }
    std::vector<double> get3DMedian(FunctionalForm &f)
    {
        std::vector<double> sortedX, sortedY, sortedZ, x = f.parameterSpace[0], y = f.parameterSpace[1], z = f.parameterSpace[2];
        std::vector<double> sortedXW, sortedYW, sortedZW, xw = f.weightSpace[0], yw = f.weightSpace[1], zw = f.weightSpace[2];

        for (int i = 0; i < f.extraParameterSpace[0].size(); i++)
        {
            x.push_back(f.extraParameterSpace[0][i]);
            y.push_back(f.extraParameterSpace[1][i]);
            z.push_back(f.extraParameterSpace[2][i]);

            xw.push_back(f.extraWeightSpace[0][i]);
            yw.push_back(f.extraWeightSpace[1][i]);
            zw.push_back(f.extraWeightSpace[2][i]);
        }


        sortedX = x;
        sortedY = y;
        sortedZ = z;

        sortedXW = xw;
        sortedYW = yw;
        sortedZW = zw;

        sort(sortedXW, sortedX);
        sort(sortedYW, sortedY);
        sort(sortedZW, sortedZ);

        std::vector<double> toRet;
        double b = getMedian((int) sortedX.size(), sortedXW, sortedX), m = getMedian((int) sortedY.size(), sortedYW, sortedY), xNot = getMedian((int) sortedZ.size(), sortedZW, sortedZ);
        //std::cout << "MEDIAN        A: " << b << "\t" << "B: " << m << "\tC: " << xNot << "\n";
        toRet.push_back(b);
        toRet.push_back(m);
        toRet.push_back(xNot);

        f.meanstartingpoint = toRet; //used for the guess for the generalized mean calculation

        return toRet;
    }
    std::vector<double> getNDMedian(FunctionalForm &f)
    {
        std::vector<double> yHold, sortedYHold, wHold, sortedWHold, nDimMedian;

        nDimMedian.resize(f.parameterSpace.size());

        for (size_t i = 0; i < f.parameterSpace.size(); i++)
        {
            yHold = f.parameterSpace[i];
            wHold = f.weightSpace[i];

            for (int j = 0; j < f.extraParameterSpace[0].size(); j++) //adds the extra exception parameter data
            {
                yHold.push_back(f.extraParameterSpace[i][j]);

                wHold.push_back(f.extraWeightSpace[i][j]);
            }

            sortedYHold = yHold;
            sortedWHold = wHold;

            sort(sortedWHold, sortedYHold);

            nDimMedian[i] = getMedian((int) sortedYHold.size(), sortedWHold, sortedYHold);
        }

        f.meanstartingpoint = nDimMedian; //used for the guess for the generalized mean calculation

        return nDimMedian;
    }
    std::vector<double> getNDMode(FunctionalForm &f)
    {
        std::vector<int> indexBounds, indexBoundIn;
        std::vector<double> w, y, sortedW, sortedY, sortedYHold, sortedWHold, halfSampleSubsetOld, halfSampleSubsetNew;
        std::vector<std::vector<double> > inBoundsParameterSpace, inBoundsWeightSpace, sortedYHoldVec, sortedWHoldVec, sortedYVec, sortedWVec;
        struct Bound
        {
            int lowIndex, highIndex;
            double lowValue, highValue;
        };

        Bound paramBound, paramBoundIn;
        std::vector<Bound> halfSampleBounds, halfSampleBoundsIn;

        inBoundsParameterSpace = f.parameterSpace;
        inBoundsWeightSpace = f.weightSpace;

        for (int j = 0; j < f.parameterSpace.size(); j++) {
            sortedYVec.push_back(f.parameterSpace[j]);
            sortedWVec.push_back(f.weightSpace[j]);

            sortedY = sortedYVec[j];
            sortedW = sortedWVec[j];

            //indexBounds = getHalfSampleBounds(0, sortedY.size() - 1, sortedW, sortedY);
            paramBound.lowIndex = 0;
            paramBound.highIndex = (int) sortedY.size() - 1;

            //paramBound.lowValue = sortedY[indexBounds[0]];
            //paramBound.highValue = sortedY[indexBounds[1]];

            paramBoundIn.lowIndex = -1;
            paramBoundIn.highIndex = -1;

            

            halfSampleBounds.push_back(paramBound);

            halfSampleBoundsIn.push_back(paramBoundIn);
        }

        bool nullSet = false;

        while (true)
        {
            int checker = 0;
            for (int j = 0; j < f.parameterSpace.size(); j++) {
                if (halfSampleBounds[j].lowIndex != halfSampleBoundsIn[j].lowIndex) {
                    checker += 1;
                }
            }
            if (nullSet) { //checks both conditions for the continuation of the loop
                break;
            }
            if (checker < 1) {
                break;
            }

            for (int j = 0; j < f.parameterSpace.size(); j++)
            {
                halfSampleBoundsIn[j].lowIndex = halfSampleBounds[j].lowIndex;
                halfSampleBoundsIn[j].highIndex = halfSampleBounds[j].highIndex;

                sortedW = sortedWVec[j];
                sortedY = sortedYVec[j];

                sort(sortedW, sortedY);

                indexBounds = getHalfSampleBounds(0, (int) sortedY.size() - 1, sortedW, sortedY);
                paramBound.lowIndex = indexBounds[0];
                paramBound.highIndex = indexBounds[1];
                paramBound.lowValue = sortedY[indexBounds[0]];
                paramBound.highValue = sortedY[indexBounds[1]];

                halfSampleBounds[j] = paramBound;

                sortedYHold = sortedY;
                sortedWHold = sortedW;

                sortedYHoldVec.push_back(sortedYHold);
                sortedWHoldVec.push_back(sortedWHold);

                sortedY.clear();
                sortedW.clear();

            }
            for (size_t i = 0; i < f.parameterSpace[0].size(); i++) //generalizes the code on line 1227
            {
                int check = 0;
                for (int j = 0; j < f.parameterSpace.size(); j++) { // j indexes each parameter
                    if (f.parameterSpace[j][i] >= halfSampleBounds[j].lowValue && f.parameterSpace[j][i] <= halfSampleBounds[j].highValue) {
                        check += 1;
                    }
                }
                if (check == f.parameterSpace.size()) {
                    for (int j = 0; j < f.parameterSpace.size(); j++) {
                        sortedYVec[j].push_back(f.parameterSpace[j][i]);
                        sortedWVec[j].push_back(f.weightSpace[j][i]);
                    }
                }
            }
            nullSet = false;

            int sortedcheck = 0;
            for (int j = 0; j < f.parameterSpace.size(); j++) {
                if (sortedYVec[j].size() == 0) {
                    sortedcheck += 1;
                }
            }
            if (sortedcheck >= 1) {
                nullSet = true;

                for (int j = 0; j < f.parameterSpace.size(); j++) {
                    sortedYVec[j] = sortedYHoldVec[j];
                    sortedWVec[j] = sortedWHoldVec[j];
                }
            }
        }
        
        for (int j = 0; j < f.parameterSpace.size(); j++)
        {
            sort(sortedWVec[j], sortedYVec[j]);
        }
        std::vector <double> toRet;

        for (int j = 0; j < f.parameterSpace.size(); j++)
        {
            toRet.push_back(getMedian((int) sortedYVec[j].size(), sortedWVec[j], sortedYVec[j]));
        }

        return toRet;
    }
    void RCR::setInitialModel(std::vector<double> &y)
    {
        bool found;
        int indexFound = -1;
        std::vector<double> model;
        std::vector<std::vector<double> > models;
        this->result.flags.resize(y.size(), true);
        this->parametricModel->setTrueVec(result.flags, y);
        found = false;
        while (!found)
        {
            this->parametricModel->buildModelSpace();
            this->delta = parametricModel->parameterSpace.size();
            model = get2DMode(*parametricModel);

            this->parametricModel->setModel(model);
            found = false;
            for (size_t i = 0; i < models.size(); i++)
            {
                if (std::abs(models[i][0] - model[0]) < .0001 && std::abs(models[i][1] - model[1]) < .0001)// && std::abs(models[i][2] - model[2]) < .0001)
                {
                    found = true;
                    indexFound = (int) i;
                }
            }
            if (!found)
            {
                models.push_back(model);
            }
        }
        double mSum = 0, bSum = 0; //, xSum = 0;
        for (size_t i = indexFound; i < models.size(); i++)
        {
            bSum += models[i][0];
            mSum += models[i][1];
            //xSum += models[i][2];
        }

        bSum = bSum / (double)(models.size() - indexFound);
        mSum = mSum / (double)(models.size() - indexFound);
        //xSum = xSum / (double)(models.size() - indexFound);

        std::vector<double> finalModel;
        finalModel.push_back(bSum);
        finalModel.push_back(mSum);
        //finalModel.push_back(xSum);
        this->parametricModel->setModel(finalModel);
    }

    //standard procedures:
    void RCR::alignTechniques()
    {
        switch (this->rejectionTech)
        {
        case SS_MEDIAN_DL:
            this->muTech = MEDIAN;
            this->sigmaTech = DOUBLE_LINE;
            this->sigmaChoice = SINGLE;
            break;
        case LS_MODE_68:
            this->muTech = MODE;
            this->sigmaTech = SIXTY_EIGHTH_PERCENTILE;
            this->sigmaChoice = LOWER;
            break;
        case LS_MODE_DL:
            this->muTech = MODE;
            this->sigmaTech = DOUBLE_LINE;
            this->sigmaChoice = LOWER;
            break;
        case ES_MODE_DL:
            this->muTech = MODE;
            this->sigmaTech = DOUBLE_LINE;
            this->sigmaChoice = EACH;
        }
    }
    void RCR::loadUnityTables()
    {
        if (!unityTablesLoaded)
        {
            std::vector<double> dubFiller;
            ESUnity.resize(1001, 0.0);
            SSUnity.resize(1001, 0.0);
            LSUnity.resize(1001, 0.0);
            ESDLUnityCF.resize(101, 0.0);
            LSDLUnityCF.resize(101, 0.0);
            LS68UnityCF.resize(101, 0.0);
            SSDLUnityCF.resize(101, 0.0);
            SSConstants.resize(2, dubFiller);
            SSConstants[0].resize(8, 0.0);
            SSConstants[1].resize(8, 0.0);
            SSConstants[0][4] = .202399;
            SSConstants[1][4] = .464231;
            SSConstants[0][5] = -.29158;
            SSConstants[1][5] = .26031;
            SSConstants[0][6] = -.03321;
            SSConstants[1][6] = .363776;
            SSConstants[0][7] = -.18181;
            SSConstants[1][7] = .454703;
            ESUnity[0] = 0;
            SSUnity[0] = 0;
            LSUnity[0] = 0;
            ESUnity[1] = 0;
            SSUnity[1] = 0;
            LSUnity[1] = 0;
            ESUnity[2] = 0;
            SSUnity[2] = 0;
            LSUnity[2] = 0;
            ESUnity[3] = 0;
            SSUnity[3] = 0;
            LSUnity[3] = 0;
            ESUnity[4] = 0;
            SSUnity[4] = 2.15681;
            LSUnity[4] = 0;
            ESUnity[5] = 8.1666;
            SSUnity[5] = 8.81255;
            LSUnity[5] = 36.8534;
            ESUnity[6] = 5.0808;
            SSUnity[6] = 2.72685;
            LSUnity[6] = 35.9209;
            ESUnity[7] = 6.60453;
            SSUnity[7] = 3.45589;
            LSUnity[7] = 34.563;
            ESUnity[8] = 6.57429;
            SSUnity[8] = 2.32749;
            LSUnity[8] = 33.7065;
            ESUnity[9] = 5.86414;
            SSUnity[9] = 2.90262;
            LSUnity[9] = 27.7158;
            ESUnity[10] = 4.75354;
            SSUnity[10] = 2.41888;
            LSUnity[10] = 19.4552;
            ESUnity[11] = 4.31763;
            SSUnity[11] = 2.76084;
            LSUnity[11] = 15.9112;
            ESUnity[12] = 3.95886;
            SSUnity[12] = 2.37034;
            LSUnity[12] = 12.9506;
            ESUnity[13] = 3.67176;
            SSUnity[13] = 2.47951;
            LSUnity[13] = 11.4666;
            ESUnity[14] = 3.61546;
            SSUnity[14] = 2.3207;
            LSUnity[14] = 10.9081;
            ESUnity[15] = 3.27558;
            SSUnity[15] = 2.41864;
            LSUnity[15] = 9.57823;
            ESUnity[16] = 3.1039;
            SSUnity[16] = 2.21943;
            LSUnity[16] = 8.10273;
            ESUnity[17] = 2.9816;
            SSUnity[17] = 2.33924;
            LSUnity[17] = 7.26361;
            ESUnity[18] = 2.7951;
            SSUnity[18] = 2.18108;
            LSUnity[18] = 6.34613;
            ESUnity[19] = 2.71315;
            SSUnity[19] = 2.25062;
            LSUnity[19] = 5.77615;
            ESUnity[20] = 2.65111;
            SSUnity[20] = 2.14112;
            LSUnity[20] = 5.31425;
            ESUnity[21] = 2.60358;
            SSUnity[21] = 2.1984;
            LSUnity[21] = 4.88903;
            ESUnity[22] = 2.55717;
            SSUnity[22] = 2.09772;
            LSUnity[22] = 4.54831;
            ESUnity[23] = 2.48317;
            SSUnity[23] = 2.13145;
            LSUnity[23] = 4.4066;
            ESUnity[24] = 2.44219;
            SSUnity[24] = 2.08836;
            LSUnity[24] = 4.14193;
            ESUnity[25] = 2.42666;
            SSUnity[25] = 2.12187;
            LSUnity[25] = 3.98687;
            ESUnity[26] = 2.34112;
            SSUnity[26] = 2.03728;
            LSUnity[26] = 3.88145;
            ESUnity[27] = 2.34151;
            SSUnity[27] = 2.08679;
            LSUnity[27] = 3.78224;
            ESUnity[28] = 2.3212;
            SSUnity[28] = 2.05492;
            LSUnity[28] = 3.53612;
            ESUnity[29] = 2.27048;
            SSUnity[29] = 2.05858;
            LSUnity[29] = 3.51816;
            ESUnity[30] = 2.2469;
            SSUnity[30] = 2.04114;
            LSUnity[30] = 3.41364;
            ESUnity[31] = 2.21856;
            SSUnity[31] = 2.06094;
            LSUnity[31] = 3.31962;
            ESUnity[32] = 2.18775;
            SSUnity[32] = 2.02109;
            LSUnity[32] = 3.23311;
            ESUnity[33] = 2.17785;
            SSUnity[33] = 2.02455;
            LSUnity[33] = 3.23204;
            ESUnity[34] = 2.18495;
            SSUnity[34] = 1.99134;
            LSUnity[34] = 3.1276;
            ESUnity[35] = 2.14245;
            SSUnity[35] = 2.01884;
            LSUnity[35] = 3.09802;
            ESUnity[36] = 2.13495;
            SSUnity[36] = 2.01627;
            LSUnity[36] = 3.0242;
            ESUnity[37] = 2.12995;
            SSUnity[37] = 2.0241;
            LSUnity[37] = 3.01162;
            ESUnity[38] = 2.12228;
            SSUnity[38] = 1.96923;
            LSUnity[38] = 2.91956;
            ESUnity[39] = 2.06979;
            SSUnity[39] = 2.01425;
            LSUnity[39] = 2.92877;
            ESUnity[40] = 2.08234;
            SSUnity[40] = 1.9681;
            LSUnity[40] = 2.87747;
            ESUnity[41] = 2.06751;
            SSUnity[41] = 1.98454;
            LSUnity[41] = 2.82572;
            ESUnity[42] = 2.05053;
            SSUnity[42] = 1.97773;
            LSUnity[42] = 2.84595;
            ESUnity[43] = 2.04679;
            SSUnity[43] = 1.99193;
            LSUnity[43] = 2.82653;
            ESUnity[44] = 2.04511;
            SSUnity[44] = 1.95293;
            LSUnity[44] = 2.76119;
            ESUnity[45] = 2.03733;
            SSUnity[45] = 1.96469;
            LSUnity[45] = 2.73985;
            ESUnity[46] = 2.00788;
            SSUnity[46] = 1.97371;
            LSUnity[46] = 2.7262;
            ESUnity[47] = 1.99989;
            SSUnity[47] = 1.99114;
            LSUnity[47] = 2.70399;
            ESUnity[48] = 1.98854;
            SSUnity[48] = 1.93276;
            LSUnity[48] = 2.68546;
            ESUnity[49] = 1.99807;
            SSUnity[49] = 1.96951;
            LSUnity[49] = 2.66641;
            ESUnity[50] = 1.97196;
            SSUnity[50] = 1.9455;
            LSUnity[50] = 2.66605;
            ESUnity[51] = 1.96665;
            SSUnity[51] = 1.96897;
            LSUnity[51] = 2.62692;
            ESUnity[52] = 1.96873;
            SSUnity[52] = 1.96819;
            LSUnity[52] = 2.59297;
            ESUnity[53] = 1.96407;
            SSUnity[53] = 1.94352;
            LSUnity[53] = 2.58996;
            ESUnity[54] = 1.97134;
            SSUnity[54] = 1.9391;
            LSUnity[54] = 2.58574;
            ESUnity[55] = 1.97476;
            SSUnity[55] = 1.97577;
            LSUnity[55] = 2.56806;
            ESUnity[56] = 1.9479;
            SSUnity[56] = 1.94596;
            LSUnity[56] = 2.52424;
            ESUnity[57] = 1.93348;
            SSUnity[57] = 1.94278;
            LSUnity[57] = 2.54642;
            ESUnity[58] = 1.93729;
            SSUnity[58] = 1.95162;
            LSUnity[58] = 2.51156;
            ESUnity[59] = 1.94159;
            SSUnity[59] = 1.94947;
            LSUnity[59] = 2.52289;
            ESUnity[60] = 1.93376;
            SSUnity[60] = 1.9354;
            LSUnity[60] = 2.51645;
            ESUnity[61] = 1.93016;
            SSUnity[61] = 1.95432;
            LSUnity[61] = 2.50626;
            ESUnity[62] = 1.92467;
            SSUnity[62] = 1.94189;
            LSUnity[62] = 2.47096;
            ESUnity[63] = 1.90928;
            SSUnity[63] = 1.93058;
            LSUnity[63] = 2.47922;
            ESUnity[64] = 1.91914;
            SSUnity[64] = 1.92598;
            LSUnity[64] = 2.48032;
            ESUnity[65] = 1.90782;
            SSUnity[65] = 1.93135;
            LSUnity[65] = 2.46965;
            ESUnity[66] = 1.89641;
            SSUnity[66] = 1.92705;
            LSUnity[66] = 2.45449;
            ESUnity[67] = 1.91532;
            SSUnity[67] = 1.94832;
            LSUnity[67] = 2.47265;
            ESUnity[68] = 1.90886;
            SSUnity[68] = 1.94477;
            LSUnity[68] = 2.46335;
            ESUnity[69] = 1.92046;
            SSUnity[69] = 1.93304;
            LSUnity[69] = 2.46095;
            ESUnity[70] = 1.89849;
            SSUnity[70] = 1.91891;
            LSUnity[70] = 2.46207;
            ESUnity[71] = 1.90122;
            SSUnity[71] = 1.92816;
            LSUnity[71] = 2.44411;
            ESUnity[72] = 1.89933;
            SSUnity[72] = 1.93071;
            LSUnity[72] = 2.44908;
            ESUnity[73] = 1.89731;
            SSUnity[73] = 1.92636;
            LSUnity[73] = 2.44387;
            ESUnity[74] = 1.89772;
            SSUnity[74] = 1.92686;
            LSUnity[74] = 2.42582;
            ESUnity[75] = 1.88835;
            SSUnity[75] = 1.94204;
            LSUnity[75] = 2.43031;
            ESUnity[76] = 1.87918;
            SSUnity[76] = 1.91158;
            LSUnity[76] = 2.41448;
            ESUnity[77] = 1.87955;
            SSUnity[77] = 1.93024;
            LSUnity[77] = 2.40721;
            ESUnity[78] = 1.88165;
            SSUnity[78] = 1.93606;
            LSUnity[78] = 2.4291;
            ESUnity[79] = 1.88726;
            SSUnity[79] = 1.94116;
            LSUnity[79] = 2.3857;
            ESUnity[80] = 1.86683;
            SSUnity[80] = 1.94119;
            LSUnity[80] = 2.41122;
            ESUnity[81] = 1.87656;
            SSUnity[81] = 1.92857;
            LSUnity[81] = 2.38928;
            ESUnity[82] = 1.88028;
            SSUnity[82] = 1.92096;
            LSUnity[82] = 2.3845;
            ESUnity[83] = 1.88149;
            SSUnity[83] = 1.94413;
            LSUnity[83] = 2.37153;
            ESUnity[84] = 1.87248;
            SSUnity[84] = 1.90939;
            LSUnity[84] = 2.39368;
            ESUnity[85] = 1.88399;
            SSUnity[85] = 1.90963;
            LSUnity[85] = 2.39749;
            ESUnity[86] = 1.86493;
            SSUnity[86] = 1.91973;
            LSUnity[86] = 2.36566;
            ESUnity[87] = 1.87471;
            SSUnity[87] = 1.92188;
            LSUnity[87] = 2.36452;
            ESUnity[88] = 1.86693;
            SSUnity[88] = 1.92554;
            LSUnity[88] = 2.36156;
            ESUnity[89] = 1.86201;
            SSUnity[89] = 1.92088;
            LSUnity[89] = 2.3567;
            ESUnity[90] = 1.86646;
            SSUnity[90] = 1.91027;
            LSUnity[90] = 2.37369;
            ESUnity[91] = 1.86772;
            SSUnity[91] = 1.93532;
            LSUnity[91] = 2.36823;
            ESUnity[92] = 1.86229;
            SSUnity[92] = 1.92246;
            LSUnity[92] = 2.38324;
            ESUnity[93] = 1.85577;
            SSUnity[93] = 1.92766;
            LSUnity[93] = 2.36444;
            ESUnity[94] = 1.87048;
            SSUnity[94] = 1.90794;
            LSUnity[94] = 2.35227;
            ESUnity[95] = 1.86354;
            SSUnity[95] = 1.91497;
            LSUnity[95] = 2.35479;
            ESUnity[96] = 1.84507;
            SSUnity[96] = 1.92607;
            LSUnity[96] = 2.36963;
            ESUnity[97] = 1.8591;
            SSUnity[97] = 1.93512;
            LSUnity[97] = 2.34846;
            ESUnity[98] = 1.85134;
            SSUnity[98] = 1.91289;
            LSUnity[98] = 2.36358;
            ESUnity[99] = 1.8582;
            SSUnity[99] = 1.94698;
            LSUnity[99] = 2.36588;
            ESUnity[100] = 1.83742;
            SSUnity[100] = 1.92398;
            LSUnity[100] = 2.33237;
            ESUnity[101] = 1.85959;
            SSUnity[101] = 1.93819;
            LSUnity[101] = 2.31503;
            ESUnity[102] = 1.84882;
            SSUnity[102] = 1.91674;
            LSUnity[102] = 2.32805;
            ESUnity[103] = 1.85987;
            SSUnity[103] = 1.93119;
            LSUnity[103] = 2.33006;
            ESUnity[104] = 1.86173;
            SSUnity[104] = 1.93532;
            LSUnity[104] = 2.30914;
            ESUnity[105] = 1.83974;
            SSUnity[105] = 1.91059;
            LSUnity[105] = 2.33301;
            ESUnity[106] = 1.86365;
            SSUnity[106] = 1.91698;
            LSUnity[106] = 2.32591;
            ESUnity[107] = 1.84821;
            SSUnity[107] = 1.91589;
            LSUnity[107] = 2.32715;
            ESUnity[108] = 1.85502;
            SSUnity[108] = 1.92052;
            LSUnity[108] = 2.33602;
            ESUnity[109] = 1.85445;
            SSUnity[109] = 1.91982;
            LSUnity[109] = 2.3315;
            ESUnity[110] = 1.85968;
            SSUnity[110] = 1.90968;
            LSUnity[110] = 2.32701;
            ESUnity[111] = 1.85108;
            SSUnity[111] = 1.90836;
            LSUnity[111] = 2.31637;
            ESUnity[112] = 1.85655;
            SSUnity[112] = 1.91082;
            LSUnity[112] = 2.31514;
            ESUnity[113] = 1.85065;
            SSUnity[113] = 1.91337;
            LSUnity[113] = 2.30695;
            ESUnity[114] = 1.85595;
            SSUnity[114] = 1.8941;
            LSUnity[114] = 2.32584;
            ESUnity[115] = 1.85469;
            SSUnity[115] = 1.92452;
            LSUnity[115] = 2.33675;
            ESUnity[116] = 1.87006;
            SSUnity[116] = 1.92157;
            LSUnity[116] = 2.33606;
            ESUnity[117] = 1.87031;
            SSUnity[117] = 1.90225;
            LSUnity[117] = 2.32676;
            ESUnity[118] = 1.87044;
            SSUnity[118] = 1.90111;
            LSUnity[118] = 2.3061;
            ESUnity[119] = 1.83924;
            SSUnity[119] = 1.89789;
            LSUnity[119] = 2.31024;
            ESUnity[120] = 1.8517;
            SSUnity[120] = 1.90852;
            LSUnity[120] = 2.32392;
            ESUnity[121] = 1.84286;
            SSUnity[121] = 1.91614;
            LSUnity[121] = 2.30527;
            ESUnity[122] = 1.84867;
            SSUnity[122] = 1.90442;
            LSUnity[122] = 2.32664;
            ESUnity[123] = 1.86211;
            SSUnity[123] = 1.90563;
            LSUnity[123] = 2.32271;
            ESUnity[124] = 1.85371;
            SSUnity[124] = 1.89635;
            LSUnity[124] = 2.3182;
            ESUnity[125] = 1.86309;
            SSUnity[125] = 1.91467;
            LSUnity[125] = 2.30015;
            ESUnity[126] = 1.85043;
            SSUnity[126] = 1.89847;
            LSUnity[126] = 2.30045;
            ESUnity[127] = 1.85959;
            SSUnity[127] = 1.91982;
            LSUnity[127] = 2.30848;
            ESUnity[128] = 1.85356;
            SSUnity[128] = 1.90748;
            LSUnity[128] = 2.33988;
            ESUnity[129] = 1.86067;
            SSUnity[129] = 1.93264;
            LSUnity[129] = 2.30429;
            ESUnity[130] = 1.86836;
            SSUnity[130] = 1.89417;
            LSUnity[130] = 2.30381;
            ESUnity[131] = 1.86463;
            SSUnity[131] = 1.91169;
            LSUnity[131] = 2.30484;
            ESUnity[132] = 1.86499;
            SSUnity[132] = 1.89229;
            LSUnity[132] = 2.30502;
            ESUnity[133] = 1.86168;
            SSUnity[133] = 1.92045;
            LSUnity[133] = 2.30216;
            ESUnity[134] = 1.85547;
            SSUnity[134] = 1.9015;
            LSUnity[134] = 2.29704;
            ESUnity[135] = 1.85889;
            SSUnity[135] = 1.90308;
            LSUnity[135] = 2.32589;
            ESUnity[136] = 1.85769;
            SSUnity[136] = 1.91328;
            LSUnity[136] = 2.30058;
            ESUnity[137] = 1.86705;
            SSUnity[137] = 1.90585;
            LSUnity[137] = 2.32206;
            ESUnity[138] = 1.86334;
            SSUnity[138] = 1.90799;
            LSUnity[138] = 2.32104;
            ESUnity[139] = 1.84789;
            SSUnity[139] = 1.90027;
            LSUnity[139] = 2.31013;
            ESUnity[140] = 1.86541;
            SSUnity[140] = 1.90334;
            LSUnity[140] = 2.30736;
            ESUnity[141] = 1.84999;
            SSUnity[141] = 1.90647;
            LSUnity[141] = 2.29894;
            ESUnity[142] = 1.87105;
            SSUnity[142] = 1.90834;
            LSUnity[142] = 2.31059;
            ESUnity[143] = 1.86007;
            SSUnity[143] = 1.90034;
            LSUnity[143] = 2.30877;
            ESUnity[144] = 1.87412;
            SSUnity[144] = 1.91302;
            LSUnity[144] = 2.31979;
            ESUnity[145] = 1.85801;
            SSUnity[145] = 1.91276;
            LSUnity[145] = 2.31853;
            ESUnity[146] = 1.86571;
            SSUnity[146] = 1.9122;
            LSUnity[146] = 2.30831;
            ESUnity[147] = 1.86972;
            SSUnity[147] = 1.91081;
            LSUnity[147] = 2.31743;
            ESUnity[148] = 1.87379;
            SSUnity[148] = 1.88134;
            LSUnity[148] = 2.31835;
            ESUnity[149] = 1.8744;
            SSUnity[149] = 1.902;
            LSUnity[149] = 2.31956;
            ESUnity[150] = 1.8554;
            SSUnity[150] = 1.91332;
            LSUnity[150] = 2.29383;
            ESUnity[151] = 1.86018;
            SSUnity[151] = 1.90927;
            LSUnity[151] = 2.3098;
            ESUnity[152] = 1.8826;
            SSUnity[152] = 1.92343;
            LSUnity[152] = 2.32044;
            ESUnity[153] = 1.84798;
            SSUnity[153] = 1.91227;
            LSUnity[153] = 2.31076;
            ESUnity[154] = 1.88937;
            SSUnity[154] = 1.91807;
            LSUnity[154] = 2.32562;
            ESUnity[155] = 1.89387;
            SSUnity[155] = 1.90075;
            LSUnity[155] = 2.31578;
            ESUnity[156] = 1.86612;
            SSUnity[156] = 1.90756;
            LSUnity[156] = 2.31019;
            ESUnity[157] = 1.87548;
            SSUnity[157] = 1.91373;
            LSUnity[157] = 2.31885;
            ESUnity[158] = 1.88121;
            SSUnity[158] = 1.90331;
            LSUnity[158] = 2.31229;
            ESUnity[159] = 1.88415;
            SSUnity[159] = 1.90487;
            LSUnity[159] = 2.31798;
            ESUnity[160] = 1.87015;
            SSUnity[160] = 1.90561;
            LSUnity[160] = 2.31418;
            ESUnity[161] = 1.8861;
            SSUnity[161] = 1.91607;
            LSUnity[161] = 2.31711;
            ESUnity[162] = 1.88231;
            SSUnity[162] = 1.92025;
            LSUnity[162] = 2.29978;
            ESUnity[163] = 1.86763;
            SSUnity[163] = 1.91531;
            LSUnity[163] = 2.30664;
            ESUnity[164] = 1.87352;
            SSUnity[164] = 1.91193;
            LSUnity[164] = 2.32066;
            ESUnity[165] = 1.87648;
            SSUnity[165] = 1.89599;
            LSUnity[165] = 2.32196;
            ESUnity[166] = 1.86567;
            SSUnity[166] = 1.91862;
            LSUnity[166] = 2.31512;
            ESUnity[167] = 1.87065;
            SSUnity[167] = 1.90261;
            LSUnity[167] = 2.30798;
            ESUnity[168] = 1.86588;
            SSUnity[168] = 1.90792;
            LSUnity[168] = 2.32331;
            ESUnity[169] = 1.89753;
            SSUnity[169] = 1.90866;
            LSUnity[169] = 2.30176;
            ESUnity[170] = 1.87486;
            SSUnity[170] = 1.89406;
            LSUnity[170] = 2.29932;
            ESUnity[171] = 1.87875;
            SSUnity[171] = 1.91595;
            LSUnity[171] = 2.30248;
            ESUnity[172] = 1.85801;
            SSUnity[172] = 1.89274;
            LSUnity[172] = 2.30617;
            ESUnity[173] = 1.8802;
            SSUnity[173] = 1.89967;
            LSUnity[173] = 2.31901;
            ESUnity[174] = 1.87334;
            SSUnity[174] = 1.90876;
            LSUnity[174] = 2.32557;
            ESUnity[175] = 1.88844;
            SSUnity[175] = 1.89717;
            LSUnity[175] = 2.31063;
            ESUnity[176] = 1.89476;
            SSUnity[176] = 1.90859;
            LSUnity[176] = 2.31588;
            ESUnity[177] = 1.90011;
            SSUnity[177] = 1.91391;
            LSUnity[177] = 2.3159;
            ESUnity[178] = 1.8868;
            SSUnity[178] = 1.90559;
            LSUnity[178] = 2.31772;
            ESUnity[179] = 1.90092;
            SSUnity[179] = 1.90095;
            LSUnity[179] = 2.31842;
            ESUnity[180] = 1.88328;
            SSUnity[180] = 1.90033;
            LSUnity[180] = 2.31526;
            ESUnity[181] = 1.89275;
            SSUnity[181] = 1.90289;
            LSUnity[181] = 2.3276;
            ESUnity[182] = 1.87378;
            SSUnity[182] = 1.90295;
            LSUnity[182] = 2.32912;
            ESUnity[183] = 1.89406;
            SSUnity[183] = 1.90713;
            LSUnity[183] = 2.34188;
            ESUnity[184] = 1.90745;
            SSUnity[184] = 1.90025;
            LSUnity[184] = 2.27696;
            ESUnity[185] = 1.89517;
            SSUnity[185] = 1.88309;
            LSUnity[185] = 2.30697;
            ESUnity[186] = 1.89659;
            SSUnity[186] = 1.90272;
            LSUnity[186] = 2.30451;
            ESUnity[187] = 1.89442;
            SSUnity[187] = 1.91291;
            LSUnity[187] = 2.32676;
            ESUnity[188] = 1.89757;
            SSUnity[188] = 1.90952;
            LSUnity[188] = 2.30071;
            ESUnity[189] = 1.91347;
            SSUnity[189] = 1.8967;
            LSUnity[189] = 2.31046;
            ESUnity[190] = 1.873;
            SSUnity[190] = 1.9081;
            LSUnity[190] = 2.31215;
            ESUnity[191] = 1.91664;
            SSUnity[191] = 1.90577;
            LSUnity[191] = 2.32906;
            ESUnity[192] = 1.90175;
            SSUnity[192] = 1.91587;
            LSUnity[192] = 2.32624;
            ESUnity[193] = 1.90809;
            SSUnity[193] = 1.88008;
            LSUnity[193] = 2.31233;
            ESUnity[194] = 1.90104;
            SSUnity[194] = 1.89939;
            LSUnity[194] = 2.31996;
            ESUnity[195] = 1.89359;
            SSUnity[195] = 1.90212;
            LSUnity[195] = 2.31474;
            ESUnity[196] = 1.91446;
            SSUnity[196] = 1.93413;
            LSUnity[196] = 2.33245;
            ESUnity[197] = 1.90136;
            SSUnity[197] = 1.89113;
            LSUnity[197] = 2.32916;
            ESUnity[198] = 1.89328;
            SSUnity[198] = 1.92183;
            LSUnity[198] = 2.32184;
            ESUnity[199] = 1.90306;
            SSUnity[199] = 1.88261;
            LSUnity[199] = 2.32315;
            ESUnity[200] = 1.91074;
            SSUnity[200] = 1.91105;
            LSUnity[200] = 2.33396;
            ESUnity[201] = 1.92347;
            SSUnity[201] = 1.91124;
            LSUnity[201] = 2.33788;
            ESUnity[202] = 1.88847;
            SSUnity[202] = 1.91656;
            LSUnity[202] = 2.32957;
            ESUnity[203] = 1.90926;
            SSUnity[203] = 1.91666;
            LSUnity[203] = 2.32283;
            ESUnity[204] = 1.90457;
            SSUnity[204] = 1.89339;
            LSUnity[204] = 2.32555;
            ESUnity[205] = 1.90546;
            SSUnity[205] = 1.90816;
            LSUnity[205] = 2.324;
            ESUnity[206] = 1.92029;
            SSUnity[206] = 1.90318;
            LSUnity[206] = 2.31826;
            ESUnity[207] = 1.90872;
            SSUnity[207] = 1.90783;
            LSUnity[207] = 2.31881;
            ESUnity[208] = 1.91051;
            SSUnity[208] = 1.89488;
            LSUnity[208] = 2.32999;
            ESUnity[209] = 1.89849;
            SSUnity[209] = 1.91292;
            LSUnity[209] = 2.31921;
            ESUnity[210] = 1.90966;
            SSUnity[210] = 1.91516;
            LSUnity[210] = 2.32627;
            ESUnity[211] = 1.90013;
            SSUnity[211] = 1.9127;
            LSUnity[211] = 2.32733;
            ESUnity[212] = 1.91262;
            SSUnity[212] = 1.89191;
            LSUnity[212] = 2.32873;
            ESUnity[213] = 1.89846;
            SSUnity[213] = 1.88983;
            LSUnity[213] = 2.33075;
            ESUnity[214] = 1.91245;
            SSUnity[214] = 1.90066;
            LSUnity[214] = 2.33536;
            ESUnity[215] = 1.91623;
            SSUnity[215] = 1.90265;
            LSUnity[215] = 2.32293;
            ESUnity[216] = 1.92408;
            SSUnity[216] = 1.90542;
            LSUnity[216] = 2.31342;
            ESUnity[217] = 1.92469;
            SSUnity[217] = 1.9038;
            LSUnity[217] = 2.32611;
            ESUnity[218] = 1.91827;
            SSUnity[218] = 1.8986;
            LSUnity[218] = 2.33335;
            ESUnity[219] = 1.91272;
            SSUnity[219] = 1.89256;
            LSUnity[219] = 2.35003;
            ESUnity[220] = 1.92312;
            SSUnity[220] = 1.88351;
            LSUnity[220] = 2.33859;
            ESUnity[221] = 1.90166;
            SSUnity[221] = 1.8852;
            LSUnity[221] = 2.32431;
            ESUnity[222] = 1.93004;
            SSUnity[222] = 1.92514;
            LSUnity[222] = 2.32131;
            ESUnity[223] = 1.92804;
            SSUnity[223] = 1.88772;
            LSUnity[223] = 2.36006;
            ESUnity[224] = 1.91557;
            SSUnity[224] = 1.9043;
            LSUnity[224] = 2.32564;
            ESUnity[225] = 1.92108;
            SSUnity[225] = 1.9022;
            LSUnity[225] = 2.34198;
            ESUnity[226] = 1.9281;
            SSUnity[226] = 1.89706;
            LSUnity[226] = 2.31628;
            ESUnity[227] = 1.91955;
            SSUnity[227] = 1.89796;
            LSUnity[227] = 2.33134;
            ESUnity[228] = 1.91109;
            SSUnity[228] = 1.89817;
            LSUnity[228] = 2.35282;
            ESUnity[229] = 1.90594;
            SSUnity[229] = 1.88593;
            LSUnity[229] = 2.3304;
            ESUnity[230] = 1.94078;
            SSUnity[230] = 1.9075;
            LSUnity[230] = 2.33832;
            ESUnity[231] = 1.9506;
            SSUnity[231] = 1.88818;
            LSUnity[231] = 2.3339;
            ESUnity[232] = 1.92283;
            SSUnity[232] = 1.89838;
            LSUnity[232] = 2.32174;
            ESUnity[233] = 1.92585;
            SSUnity[233] = 1.89365;
            LSUnity[233] = 2.35165;
            ESUnity[234] = 1.92885;
            SSUnity[234] = 1.89466;
            LSUnity[234] = 2.34869;
            ESUnity[235] = 1.92967;
            SSUnity[235] = 1.90708;
            LSUnity[235] = 2.35611;
            ESUnity[236] = 1.94285;
            SSUnity[236] = 1.90425;
            LSUnity[236] = 2.3312;
            ESUnity[237] = 1.92227;
            SSUnity[237] = 1.91575;
            LSUnity[237] = 2.32226;
            ESUnity[238] = 1.93435;
            SSUnity[238] = 1.89347;
            LSUnity[238] = 2.33333;
            ESUnity[239] = 1.93533;
            SSUnity[239] = 1.89518;
            LSUnity[239] = 2.34529;
            ESUnity[240] = 1.92872;
            SSUnity[240] = 1.90846;
            LSUnity[240] = 2.33104;
            ESUnity[241] = 1.93116;
            SSUnity[241] = 1.91872;
            LSUnity[241] = 2.33178;
            ESUnity[242] = 1.95184;
            SSUnity[242] = 1.90166;
            LSUnity[242] = 2.34451;
            ESUnity[243] = 1.95929;
            SSUnity[243] = 1.89268;
            LSUnity[243] = 2.35776;
            ESUnity[244] = 1.93613;
            SSUnity[244] = 1.90013;
            LSUnity[244] = 2.33515;
            ESUnity[245] = 1.93753;
            SSUnity[245] = 1.90088;
            LSUnity[245] = 2.36907;
            ESUnity[246] = 1.95813;
            SSUnity[246] = 1.90508;
            LSUnity[246] = 2.36208;
            ESUnity[247] = 1.94259;
            SSUnity[247] = 1.90721;
            LSUnity[247] = 2.3577;
            ESUnity[248] = 1.94857;
            SSUnity[248] = 1.92669;
            LSUnity[248] = 2.35587;
            ESUnity[249] = 1.93827;
            SSUnity[249] = 1.90787;
            LSUnity[249] = 2.37169;
            ESUnity[250] = 1.9403;
            SSUnity[250] = 1.87833;
            LSUnity[250] = 2.35331;
            ESUnity[251] = 1.95039;
            SSUnity[251] = 1.90566;
            LSUnity[251] = 2.34057;
            ESUnity[252] = 1.9599;
            SSUnity[252] = 1.88696;
            LSUnity[252] = 2.36204;
            ESUnity[253] = 1.94687;
            SSUnity[253] = 1.8908;
            LSUnity[253] = 2.32979;
            ESUnity[254] = 1.94243;
            SSUnity[254] = 1.8995;
            LSUnity[254] = 2.35296;
            ESUnity[255] = 1.93983;
            SSUnity[255] = 1.89255;
            LSUnity[255] = 2.34675;
            ESUnity[256] = 1.95263;
            SSUnity[256] = 1.88287;
            LSUnity[256] = 2.35479;
            ESUnity[257] = 1.94929;
            SSUnity[257] = 1.89618;
            LSUnity[257] = 2.34909;
            ESUnity[258] = 1.96734;
            SSUnity[258] = 1.89197;
            LSUnity[258] = 2.36341;
            ESUnity[259] = 1.95456;
            SSUnity[259] = 1.89878;
            LSUnity[259] = 2.36852;
            ESUnity[260] = 1.96577;
            SSUnity[260] = 1.88857;
            LSUnity[260] = 2.35247;
            ESUnity[261] = 1.94643;
            SSUnity[261] = 1.88001;
            LSUnity[261] = 2.35652;
            ESUnity[262] = 1.95825;
            SSUnity[262] = 1.88883;
            LSUnity[262] = 2.34908;
            ESUnity[263] = 1.95878;
            SSUnity[263] = 1.88522;
            LSUnity[263] = 2.34402;
            ESUnity[264] = 1.95404;
            SSUnity[264] = 1.91297;
            LSUnity[264] = 2.34519;
            ESUnity[265] = 1.95153;
            SSUnity[265] = 1.89276;
            LSUnity[265] = 2.35936;
            ESUnity[266] = 1.94547;
            SSUnity[266] = 1.89648;
            LSUnity[266] = 2.34917;
            ESUnity[267] = 1.96893;
            SSUnity[267] = 1.89694;
            LSUnity[267] = 2.35581;
            ESUnity[268] = 1.9594;
            SSUnity[268] = 1.9099;
            LSUnity[268] = 2.37457;
            ESUnity[269] = 1.97205;
            SSUnity[269] = 1.89939;
            LSUnity[269] = 2.36935;
            ESUnity[270] = 1.99737;
            SSUnity[270] = 1.90364;
            LSUnity[270] = 2.37624;
            ESUnity[271] = 1.96382;
            SSUnity[271] = 1.91509;
            LSUnity[271] = 2.35;
            ESUnity[272] = 1.976;
            SSUnity[272] = 1.88928;
            LSUnity[272] = 2.36405;
            ESUnity[273] = 1.98593;
            SSUnity[273] = 1.90324;
            LSUnity[273] = 2.37997;
            ESUnity[274] = 1.9721;
            SSUnity[274] = 1.89895;
            LSUnity[274] = 2.36499;
            ESUnity[275] = 1.97738;
            SSUnity[275] = 1.90954;
            LSUnity[275] = 2.36181;
            ESUnity[276] = 1.98635;
            SSUnity[276] = 1.89269;
            LSUnity[276] = 2.37889;
            ESUnity[277] = 1.96346;
            SSUnity[277] = 1.88783;
            LSUnity[277] = 2.37869;
            ESUnity[278] = 1.97599;
            SSUnity[278] = 1.89777;
            LSUnity[278] = 2.37991;
            ESUnity[279] = 1.97878;
            SSUnity[279] = 1.90334;
            LSUnity[279] = 2.37276;
            ESUnity[280] = 1.99636;
            SSUnity[280] = 1.88776;
            LSUnity[280] = 2.35524;
            ESUnity[281] = 1.97553;
            SSUnity[281] = 1.87254;
            LSUnity[281] = 2.38919;
            ESUnity[282] = 1.98378;
            SSUnity[282] = 1.90645;
            LSUnity[282] = 2.34321;
            ESUnity[283] = 1.96246;
            SSUnity[283] = 1.90498;
            LSUnity[283] = 2.36569;
            ESUnity[284] = 1.97884;
            SSUnity[284] = 1.90964;
            LSUnity[284] = 2.35995;
            ESUnity[285] = 1.9695;
            SSUnity[285] = 1.89432;
            LSUnity[285] = 2.38254;
            ESUnity[286] = 1.98307;
            SSUnity[286] = 1.90426;
            LSUnity[286] = 2.37014;
            ESUnity[287] = 1.97067;
            SSUnity[287] = 1.91295;
            LSUnity[287] = 2.34551;
            ESUnity[288] = 1.99561;
            SSUnity[288] = 1.91371;
            LSUnity[288] = 2.38429;
            ESUnity[289] = 1.99559;
            SSUnity[289] = 1.90281;
            LSUnity[289] = 2.39698;
            ESUnity[290] = 1.97476;
            SSUnity[290] = 1.88529;
            LSUnity[290] = 2.37262;
            ESUnity[291] = 1.98626;
            SSUnity[291] = 1.89394;
            LSUnity[291] = 2.36594;
            ESUnity[292] = 1.9757;
            SSUnity[292] = 1.91211;
            LSUnity[292] = 2.36637;
            ESUnity[293] = 1.98518;
            SSUnity[293] = 1.89553;
            LSUnity[293] = 2.38661;
            ESUnity[294] = 1.9889;
            SSUnity[294] = 1.88646;
            LSUnity[294] = 2.38222;
            ESUnity[295] = 1.98579;
            SSUnity[295] = 1.89605;
            LSUnity[295] = 2.35873;
            ESUnity[296] = 1.9907;
            SSUnity[296] = 1.89684;
            LSUnity[296] = 2.37462;
            ESUnity[297] = 1.98393;
            SSUnity[297] = 1.89165;
            LSUnity[297] = 2.38352;
            ESUnity[298] = 1.98949;
            SSUnity[298] = 1.89019;
            LSUnity[298] = 2.37973;
            ESUnity[299] = 1.97859;
            SSUnity[299] = 1.90702;
            LSUnity[299] = 2.39807;
            ESUnity[300] = 1.98883;
            SSUnity[300] = 1.90584;
            LSUnity[300] = 2.37303;
            ESUnity[301] = 1.99809;
            SSUnity[301] = 1.8982;
            LSUnity[301] = 2.38523;
            ESUnity[302] = 2.01698;
            SSUnity[302] = 1.89138;
            LSUnity[302] = 2.39735;
            ESUnity[303] = 1.97457;
            SSUnity[303] = 1.88391;
            LSUnity[303] = 2.39825;
            ESUnity[304] = 2.00691;
            SSUnity[304] = 1.91162;
            LSUnity[304] = 2.38887;
            ESUnity[305] = 1.99878;
            SSUnity[305] = 1.89593;
            LSUnity[305] = 2.39099;
            ESUnity[306] = 1.98738;
            SSUnity[306] = 1.87554;
            LSUnity[306] = 2.36962;
            ESUnity[307] = 1.98681;
            SSUnity[307] = 1.89157;
            LSUnity[307] = 2.39097;
            ESUnity[308] = 2.02094;
            SSUnity[308] = 1.89109;
            LSUnity[308] = 2.38772;
            ESUnity[309] = 2.00732;
            SSUnity[309] = 1.91159;
            LSUnity[309] = 2.37725;
            ESUnity[310] = 2.01242;
            SSUnity[310] = 1.88954;
            LSUnity[310] = 2.39564;
            ESUnity[311] = 2.00143;
            SSUnity[311] = 1.89924;
            LSUnity[311] = 2.38023;
            ESUnity[312] = 2.01719;
            SSUnity[312] = 1.90355;
            LSUnity[312] = 2.3893;
            ESUnity[313] = 2.00147;
            SSUnity[313] = 1.88566;
            LSUnity[313] = 2.39098;
            ESUnity[314] = 2.01433;
            SSUnity[314] = 1.91358;
            LSUnity[314] = 2.39579;
            ESUnity[315] = 1.98978;
            SSUnity[315] = 1.89022;
            LSUnity[315] = 2.36855;
            ESUnity[316] = 2.01423;
            SSUnity[316] = 1.9111;
            LSUnity[316] = 2.37588;
            ESUnity[317] = 2.02388;
            SSUnity[317] = 1.90806;
            LSUnity[317] = 2.41048;
            ESUnity[318] = 1.99763;
            SSUnity[318] = 1.90755;
            LSUnity[318] = 2.39376;
            ESUnity[319] = 2.01322;
            SSUnity[319] = 1.90232;
            LSUnity[319] = 2.38372;
            ESUnity[320] = 2.01779;
            SSUnity[320] = 1.91273;
            LSUnity[320] = 2.40223;
            ESUnity[321] = 2.00585;
            SSUnity[321] = 1.8973;
            LSUnity[321] = 2.37983;
            ESUnity[322] = 2.0262;
            SSUnity[322] = 1.9001;
            LSUnity[322] = 2.38342;
            ESUnity[323] = 2.01612;
            SSUnity[323] = 1.89809;
            LSUnity[323] = 2.37161;
            ESUnity[324] = 2.015;
            SSUnity[324] = 1.87097;
            LSUnity[324] = 2.39749;
            ESUnity[325] = 2.01545;
            SSUnity[325] = 1.906;
            LSUnity[325] = 2.40345;
            ESUnity[326] = 2.02878;
            SSUnity[326] = 1.90316;
            LSUnity[326] = 2.41587;
            ESUnity[327] = 1.9953;
            SSUnity[327] = 1.89894;
            LSUnity[327] = 2.36511;
            ESUnity[328] = 2.02431;
            SSUnity[328] = 1.88788;
            LSUnity[328] = 2.37592;
            ESUnity[329] = 2.00385;
            SSUnity[329] = 1.89129;
            LSUnity[329] = 2.39139;
            ESUnity[330] = 2.01857;
            SSUnity[330] = 1.91002;
            LSUnity[330] = 2.41094;
            ESUnity[331] = 2.01121;
            SSUnity[331] = 1.89125;
            LSUnity[331] = 2.39567;
            ESUnity[332] = 2.00012;
            SSUnity[332] = 1.88322;
            LSUnity[332] = 2.40709;
            ESUnity[333] = 2.04458;
            SSUnity[333] = 1.88575;
            LSUnity[333] = 2.41242;
            ESUnity[334] = 2.01686;
            SSUnity[334] = 1.91677;
            LSUnity[334] = 2.38928;
            ESUnity[335] = 2.03361;
            SSUnity[335] = 1.89967;
            LSUnity[335] = 2.41823;
            ESUnity[336] = 2.03484;
            SSUnity[336] = 1.88543;
            LSUnity[336] = 2.39035;
            ESUnity[337] = 2.01569;
            SSUnity[337] = 1.90811;
            LSUnity[337] = 2.40312;
            ESUnity[338] = 2.02847;
            SSUnity[338] = 1.90059;
            LSUnity[338] = 2.41389;
            ESUnity[339] = 2.03367;
            SSUnity[339] = 1.89444;
            LSUnity[339] = 2.40169;
            ESUnity[340] = 2.01199;
            SSUnity[340] = 1.89978;
            LSUnity[340] = 2.38431;
            ESUnity[341] = 2.03147;
            SSUnity[341] = 1.89799;
            LSUnity[341] = 2.41435;
            ESUnity[342] = 2.02107;
            SSUnity[342] = 1.90427;
            LSUnity[342] = 2.39744;
            ESUnity[343] = 2.02272;
            SSUnity[343] = 1.90272;
            LSUnity[343] = 2.43255;
            ESUnity[344] = 2.02693;
            SSUnity[344] = 1.9061;
            LSUnity[344] = 2.40659;
            ESUnity[345] = 2.02361;
            SSUnity[345] = 1.90819;
            LSUnity[345] = 2.42163;
            ESUnity[346] = 2.05638;
            SSUnity[346] = 1.90026;
            LSUnity[346] = 2.42772;
            ESUnity[347] = 2.05819;
            SSUnity[347] = 1.89927;
            LSUnity[347] = 2.39879;
            ESUnity[348] = 2.03968;
            SSUnity[348] = 1.90479;
            LSUnity[348] = 2.41731;
            ESUnity[349] = 2.06116;
            SSUnity[349] = 1.88755;
            LSUnity[349] = 2.40878;
            ESUnity[350] = 2.03919;
            SSUnity[350] = 1.88131;
            LSUnity[350] = 2.40211;
            ESUnity[351] = 2.057;
            SSUnity[351] = 1.89867;
            LSUnity[351] = 2.43587;
            ESUnity[352] = 2.04157;
            SSUnity[352] = 1.90644;
            LSUnity[352] = 2.38383;
            ESUnity[353] = 2.03594;
            SSUnity[353] = 1.88228;
            LSUnity[353] = 2.40603;
            ESUnity[354] = 2.04629;
            SSUnity[354] = 1.88071;
            LSUnity[354] = 2.43596;
            ESUnity[355] = 2.04504;
            SSUnity[355] = 1.90087;
            LSUnity[355] = 2.40904;
            ESUnity[356] = 2.05291;
            SSUnity[356] = 1.90146;
            LSUnity[356] = 2.42474;
            ESUnity[357] = 2.0337;
            SSUnity[357] = 1.87144;
            LSUnity[357] = 2.4032;
            ESUnity[358] = 2.03813;
            SSUnity[358] = 1.91352;
            LSUnity[358] = 2.41303;
            ESUnity[359] = 2.05083;
            SSUnity[359] = 1.9141;
            LSUnity[359] = 2.41279;
            ESUnity[360] = 2.05647;
            SSUnity[360] = 1.89858;
            LSUnity[360] = 2.44204;
            ESUnity[361] = 2.04183;
            SSUnity[361] = 1.89862;
            LSUnity[361] = 2.43031;
            ESUnity[362] = 2.06155;
            SSUnity[362] = 1.90519;
            LSUnity[362] = 2.40928;
            ESUnity[363] = 2.03937;
            SSUnity[363] = 1.90058;
            LSUnity[363] = 2.41103;
            ESUnity[364] = 2.03977;
            SSUnity[364] = 1.89918;
            LSUnity[364] = 2.41396;
            ESUnity[365] = 2.0622;
            SSUnity[365] = 1.90559;
            LSUnity[365] = 2.4162;
            ESUnity[366] = 2.05104;
            SSUnity[366] = 1.89944;
            LSUnity[366] = 2.41199;
            ESUnity[367] = 2.03419;
            SSUnity[367] = 1.90062;
            LSUnity[367] = 2.41941;
            ESUnity[368] = 2.03653;
            SSUnity[368] = 1.89251;
            LSUnity[368] = 2.39934;
            ESUnity[369] = 2.04335;
            SSUnity[369] = 1.90507;
            LSUnity[369] = 2.4284;
            ESUnity[370] = 2.05126;
            SSUnity[370] = 1.88644;
            LSUnity[370] = 2.43677;
            ESUnity[371] = 2.05407;
            SSUnity[371] = 1.89033;
            LSUnity[371] = 2.41108;
            ESUnity[372] = 2.04052;
            SSUnity[372] = 1.88505;
            LSUnity[372] = 2.42056;
            ESUnity[373] = 2.0543;
            SSUnity[373] = 1.88696;
            LSUnity[373] = 2.43162;
            ESUnity[374] = 2.06232;
            SSUnity[374] = 1.87862;
            LSUnity[374] = 2.42774;
            ESUnity[375] = 2.06837;
            SSUnity[375] = 1.90548;
            LSUnity[375] = 2.42649;
            ESUnity[376] = 2.05954;
            SSUnity[376] = 1.89442;
            LSUnity[376] = 2.4378;
            ESUnity[377] = 2.07001;
            SSUnity[377] = 1.90166;
            LSUnity[377] = 2.43639;
            ESUnity[378] = 2.06376;
            SSUnity[378] = 1.89409;
            LSUnity[378] = 2.44567;
            ESUnity[379] = 2.05613;
            SSUnity[379] = 1.89856;
            LSUnity[379] = 2.42983;
            ESUnity[380] = 2.08693;
            SSUnity[380] = 1.88847;
            LSUnity[380] = 2.42878;
            ESUnity[381] = 2.06129;
            SSUnity[381] = 1.90425;
            LSUnity[381] = 2.43716;
            ESUnity[382] = 2.06122;
            SSUnity[382] = 1.90832;
            LSUnity[382] = 2.443;
            ESUnity[383] = 2.07765;
            SSUnity[383] = 1.88451;
            LSUnity[383] = 2.42191;
            ESUnity[384] = 2.06185;
            SSUnity[384] = 1.9019;
            LSUnity[384] = 2.43273;
            ESUnity[385] = 2.07341;
            SSUnity[385] = 1.90881;
            LSUnity[385] = 2.42527;
            ESUnity[386] = 2.09579;
            SSUnity[386] = 1.91162;
            LSUnity[386] = 2.44494;
            ESUnity[387] = 2.07121;
            SSUnity[387] = 1.90154;
            LSUnity[387] = 2.43479;
            ESUnity[388] = 2.06034;
            SSUnity[388] = 1.90307;
            LSUnity[388] = 2.4229;
            ESUnity[389] = 2.05612;
            SSUnity[389] = 1.91407;
            LSUnity[389] = 2.44183;
            ESUnity[390] = 2.07538;
            SSUnity[390] = 1.88398;
            LSUnity[390] = 2.42446;
            ESUnity[391] = 2.09243;
            SSUnity[391] = 1.91621;
            LSUnity[391] = 2.46586;
            ESUnity[392] = 2.07411;
            SSUnity[392] = 1.88462;
            LSUnity[392] = 2.43181;
            ESUnity[393] = 2.0629;
            SSUnity[393] = 1.90465;
            LSUnity[393] = 2.42846;
            ESUnity[394] = 2.07861;
            SSUnity[394] = 1.87257;
            LSUnity[394] = 2.46485;
            ESUnity[395] = 2.08307;
            SSUnity[395] = 1.89214;
            LSUnity[395] = 2.43435;
            ESUnity[396] = 2.07181;
            SSUnity[396] = 1.90098;
            LSUnity[396] = 2.44537;
            ESUnity[397] = 2.06602;
            SSUnity[397] = 1.90786;
            LSUnity[397] = 2.42507;
            ESUnity[398] = 2.068;
            SSUnity[398] = 1.88021;
            LSUnity[398] = 2.44751;
            ESUnity[399] = 2.09641;
            SSUnity[399] = 1.90179;
            LSUnity[399] = 2.45449;
            ESUnity[400] = 2.07932;
            SSUnity[400] = 1.90588;
            LSUnity[400] = 2.43665;
            ESUnity[401] = 2.08298;
            SSUnity[401] = 1.88598;
            LSUnity[401] = 2.46273;
            ESUnity[402] = 2.08073;
            SSUnity[402] = 1.88188;
            LSUnity[402] = 2.43859;
            ESUnity[403] = 2.08581;
            SSUnity[403] = 1.9037;
            LSUnity[403] = 2.46295;
            ESUnity[404] = 2.0851;
            SSUnity[404] = 1.894;
            LSUnity[404] = 2.47624;
            ESUnity[405] = 2.10295;
            SSUnity[405] = 1.90562;
            LSUnity[405] = 2.42944;
            ESUnity[406] = 2.08136;
            SSUnity[406] = 1.9023;
            LSUnity[406] = 2.46189;
            ESUnity[407] = 2.07008;
            SSUnity[407] = 1.90257;
            LSUnity[407] = 2.44163;
            ESUnity[408] = 2.09438;
            SSUnity[408] = 1.91506;
            LSUnity[408] = 2.41992;
            ESUnity[409] = 2.09176;
            SSUnity[409] = 1.88724;
            LSUnity[409] = 2.43553;
            ESUnity[410] = 2.08998;
            SSUnity[410] = 1.89661;
            LSUnity[410] = 2.4875;
            ESUnity[411] = 2.08718;
            SSUnity[411] = 1.89036;
            LSUnity[411] = 2.43687;
            ESUnity[412] = 2.10746;
            SSUnity[412] = 1.91581;
            LSUnity[412] = 2.46288;
            ESUnity[413] = 2.07256;
            SSUnity[413] = 1.89118;
            LSUnity[413] = 2.4557;
            ESUnity[414] = 2.08746;
            SSUnity[414] = 1.89605;
            LSUnity[414] = 2.42297;
            ESUnity[415] = 2.0944;
            SSUnity[415] = 1.88842;
            LSUnity[415] = 2.45249;
            ESUnity[416] = 2.10528;
            SSUnity[416] = 1.8849;
            LSUnity[416] = 2.45016;
            ESUnity[417] = 2.08079;
            SSUnity[417] = 1.89656;
            LSUnity[417] = 2.49956;
            ESUnity[418] = 2.08737;
            SSUnity[418] = 1.89162;
            LSUnity[418] = 2.43977;
            ESUnity[419] = 2.10175;
            SSUnity[419] = 1.88283;
            LSUnity[419] = 2.45717;
            ESUnity[420] = 2.09705;
            SSUnity[420] = 1.89988;
            LSUnity[420] = 2.46192;
            ESUnity[421] = 2.08389;
            SSUnity[421] = 1.89194;
            LSUnity[421] = 2.47636;
            ESUnity[422] = 2.10981;
            SSUnity[422] = 1.89482;
            LSUnity[422] = 2.46011;
            ESUnity[423] = 2.11035;
            SSUnity[423] = 1.9133;
            LSUnity[423] = 2.46693;
            ESUnity[424] = 2.11369;
            SSUnity[424] = 1.89462;
            LSUnity[424] = 2.46462;
            ESUnity[425] = 2.08858;
            SSUnity[425] = 1.89691;
            LSUnity[425] = 2.46019;
            ESUnity[426] = 2.10032;
            SSUnity[426] = 1.89646;
            LSUnity[426] = 2.4574;
            ESUnity[427] = 2.10094;
            SSUnity[427] = 1.90788;
            LSUnity[427] = 2.45756;
            ESUnity[428] = 2.10825;
            SSUnity[428] = 1.91017;
            LSUnity[428] = 2.45052;
            ESUnity[429] = 2.11005;
            SSUnity[429] = 1.87967;
            LSUnity[429] = 2.44069;
            ESUnity[430] = 2.10437;
            SSUnity[430] = 1.90209;
            LSUnity[430] = 2.47748;
            ESUnity[431] = 2.11447;
            SSUnity[431] = 1.90721;
            LSUnity[431] = 2.45291;
            ESUnity[432] = 2.12606;
            SSUnity[432] = 1.91069;
            LSUnity[432] = 2.45148;
            ESUnity[433] = 2.11678;
            SSUnity[433] = 1.89462;
            LSUnity[433] = 2.4435;
            ESUnity[434] = 2.1413;
            SSUnity[434] = 1.88581;
            LSUnity[434] = 2.4392;
            ESUnity[435] = 2.09585;
            SSUnity[435] = 1.88022;
            LSUnity[435] = 2.46629;
            ESUnity[436] = 2.12634;
            SSUnity[436] = 1.88962;
            LSUnity[436] = 2.44054;
            ESUnity[437] = 2.11713;
            SSUnity[437] = 1.89843;
            LSUnity[437] = 2.47332;
            ESUnity[438] = 2.10869;
            SSUnity[438] = 1.89884;
            LSUnity[438] = 2.45765;
            ESUnity[439] = 2.12638;
            SSUnity[439] = 1.88286;
            LSUnity[439] = 2.48058;
            ESUnity[440] = 2.12402;
            SSUnity[440] = 1.90381;
            LSUnity[440] = 2.47736;
            ESUnity[441] = 2.12494;
            SSUnity[441] = 1.89093;
            LSUnity[441] = 2.47677;
            ESUnity[442] = 2.10982;
            SSUnity[442] = 1.89283;
            LSUnity[442] = 2.47976;
            ESUnity[443] = 2.11702;
            SSUnity[443] = 1.92003;
            LSUnity[443] = 2.45027;
            ESUnity[444] = 2.13528;
            SSUnity[444] = 1.89169;
            LSUnity[444] = 2.47531;
            ESUnity[445] = 2.12185;
            SSUnity[445] = 1.91138;
            LSUnity[445] = 2.44251;
            ESUnity[446] = 2.1013;
            SSUnity[446] = 1.89636;
            LSUnity[446] = 2.49847;
            ESUnity[447] = 2.12729;
            SSUnity[447] = 1.91967;
            LSUnity[447] = 2.4591;
            ESUnity[448] = 2.11141;
            SSUnity[448] = 1.8953;
            LSUnity[448] = 2.45593;
            ESUnity[449] = 2.11924;
            SSUnity[449] = 1.90448;
            LSUnity[449] = 2.49194;
            ESUnity[450] = 2.11256;
            SSUnity[450] = 1.90181;
            LSUnity[450] = 2.46848;
            ESUnity[451] = 2.11446;
            SSUnity[451] = 1.90356;
            LSUnity[451] = 2.46689;
            ESUnity[452] = 2.12816;
            SSUnity[452] = 1.89381;
            LSUnity[452] = 2.44947;
            ESUnity[453] = 2.12969;
            SSUnity[453] = 1.90262;
            LSUnity[453] = 2.46825;
            ESUnity[454] = 2.11749;
            SSUnity[454] = 1.87854;
            LSUnity[454] = 2.47907;
            ESUnity[455] = 2.14418;
            SSUnity[455] = 1.89806;
            LSUnity[455] = 2.46224;
            ESUnity[456] = 2.10116;
            SSUnity[456] = 1.92212;
            LSUnity[456] = 2.46655;
            ESUnity[457] = 2.13384;
            SSUnity[457] = 1.90101;
            LSUnity[457] = 2.48657;
            ESUnity[458] = 2.13099;
            SSUnity[458] = 1.88927;
            LSUnity[458] = 2.52373;
            ESUnity[459] = 2.14597;
            SSUnity[459] = 1.88755;
            LSUnity[459] = 2.48287;
            ESUnity[460] = 2.13887;
            SSUnity[460] = 1.89137;
            LSUnity[460] = 2.46254;
            ESUnity[461] = 2.13683;
            SSUnity[461] = 1.89204;
            LSUnity[461] = 2.49295;
            ESUnity[462] = 2.14822;
            SSUnity[462] = 1.90466;
            LSUnity[462] = 2.48062;
            ESUnity[463] = 2.15896;
            SSUnity[463] = 1.89143;
            LSUnity[463] = 2.46829;
            ESUnity[464] = 2.13273;
            SSUnity[464] = 1.90659;
            LSUnity[464] = 2.48751;
            ESUnity[465] = 2.13917;
            SSUnity[465] = 1.90207;
            LSUnity[465] = 2.49493;
            ESUnity[466] = 2.1528;
            SSUnity[466] = 1.9059;
            LSUnity[466] = 2.47246;
            ESUnity[467] = 2.13894;
            SSUnity[467] = 1.90186;
            LSUnity[467] = 2.47503;
            ESUnity[468] = 2.14315;
            SSUnity[468] = 1.89992;
            LSUnity[468] = 2.49402;
            ESUnity[469] = 2.12441;
            SSUnity[469] = 1.90816;
            LSUnity[469] = 2.47679;
            ESUnity[470] = 2.13479;
            SSUnity[470] = 1.89709;
            LSUnity[470] = 2.47079;
            ESUnity[471] = 2.13607;
            SSUnity[471] = 1.90251;
            LSUnity[471] = 2.46039;
            ESUnity[472] = 2.14126;
            SSUnity[472] = 1.8906;
            LSUnity[472] = 2.47582;
            ESUnity[473] = 2.14049;
            SSUnity[473] = 1.8974;
            LSUnity[473] = 2.49556;
            ESUnity[474] = 2.16284;
            SSUnity[474] = 1.90416;
            LSUnity[474] = 2.4796;
            ESUnity[475] = 2.13623;
            SSUnity[475] = 1.89726;
            LSUnity[475] = 2.47735;
            ESUnity[476] = 2.15132;
            SSUnity[476] = 1.88361;
            LSUnity[476] = 2.49906;
            ESUnity[477] = 2.16334;
            SSUnity[477] = 1.89477;
            LSUnity[477] = 2.49884;
            ESUnity[478] = 2.15033;
            SSUnity[478] = 1.87703;
            LSUnity[478] = 2.49199;
            ESUnity[479] = 2.14009;
            SSUnity[479] = 1.89657;
            LSUnity[479] = 2.50527;
            ESUnity[480] = 2.14827;
            SSUnity[480] = 1.9032;
            LSUnity[480] = 2.4915;
            ESUnity[481] = 2.14371;
            SSUnity[481] = 1.8782;
            LSUnity[481] = 2.48734;
            ESUnity[482] = 2.14474;
            SSUnity[482] = 1.90761;
            LSUnity[482] = 2.49249;
            ESUnity[483] = 2.147;
            SSUnity[483] = 1.90883;
            LSUnity[483] = 2.4917;
            ESUnity[484] = 2.15871;
            SSUnity[484] = 1.89922;
            LSUnity[484] = 2.46335;
            ESUnity[485] = 2.13863;
            SSUnity[485] = 1.88931;
            LSUnity[485] = 2.48479;
            ESUnity[486] = 2.15408;
            SSUnity[486] = 1.90294;
            LSUnity[486] = 2.4946;
            ESUnity[487] = 2.15092;
            SSUnity[487] = 1.91728;
            LSUnity[487] = 2.49776;
            ESUnity[488] = 2.13819;
            SSUnity[488] = 1.8975;
            LSUnity[488] = 2.48784;
            ESUnity[489] = 2.18943;
            SSUnity[489] = 1.90712;
            LSUnity[489] = 2.46877;
            ESUnity[490] = 2.16572;
            SSUnity[490] = 1.90092;
            LSUnity[490] = 2.48966;
            ESUnity[491] = 2.15846;
            SSUnity[491] = 1.90051;
            LSUnity[491] = 2.49595;
            ESUnity[492] = 2.14524;
            SSUnity[492] = 1.90849;
            LSUnity[492] = 2.47949;
            ESUnity[493] = 2.17515;
            SSUnity[493] = 1.88904;
            LSUnity[493] = 2.51807;
            ESUnity[494] = 2.15688;
            SSUnity[494] = 1.90892;
            LSUnity[494] = 2.49209;
            ESUnity[495] = 2.1765;
            SSUnity[495] = 1.88063;
            LSUnity[495] = 2.50142;
            ESUnity[496] = 2.1592;
            SSUnity[496] = 1.88996;
            LSUnity[496] = 2.51039;
            ESUnity[497] = 2.1673;
            SSUnity[497] = 1.90751;
            LSUnity[497] = 2.51679;
            ESUnity[498] = 2.15585;
            SSUnity[498] = 1.86961;
            LSUnity[498] = 2.51214;
            ESUnity[499] = 2.17656;
            SSUnity[499] = 1.90242;
            LSUnity[499] = 2.50556;
            ESUnity[500] = 2.16368;
            SSUnity[500] = 1.90487;
            LSUnity[500] = 2.50435;
            ESUnity[501] = 2.17599;
            SSUnity[501] = 1.89017;
            LSUnity[501] = 2.48686;
            ESUnity[502] = 2.17565;
            SSUnity[502] = 1.89402;
            LSUnity[502] = 2.48369;
            ESUnity[503] = 2.16999;
            SSUnity[503] = 1.90003;
            LSUnity[503] = 2.49709;
            ESUnity[504] = 2.17862;
            SSUnity[504] = 1.88839;
            LSUnity[504] = 2.49323;
            ESUnity[505] = 2.17052;
            SSUnity[505] = 1.89566;
            LSUnity[505] = 2.50155;
            ESUnity[506] = 2.19717;
            SSUnity[506] = 1.90426;
            LSUnity[506] = 2.50268;
            ESUnity[507] = 2.16832;
            SSUnity[507] = 1.90072;
            LSUnity[507] = 2.50462;
            ESUnity[508] = 2.16011;
            SSUnity[508] = 1.89577;
            LSUnity[508] = 2.50216;
            ESUnity[509] = 2.18991;
            SSUnity[509] = 1.88399;
            LSUnity[509] = 2.50117;
            ESUnity[510] = 2.17384;
            SSUnity[510] = 1.89207;
            LSUnity[510] = 2.52367;
            ESUnity[511] = 2.16112;
            SSUnity[511] = 1.90252;
            LSUnity[511] = 2.51911;
            ESUnity[512] = 2.17405;
            SSUnity[512] = 1.88967;
            LSUnity[512] = 2.52693;
            ESUnity[513] = 2.1751;
            SSUnity[513] = 1.90138;
            LSUnity[513] = 2.50703;
            ESUnity[514] = 2.19438;
            SSUnity[514] = 1.91235;
            LSUnity[514] = 2.50557;
            ESUnity[515] = 2.16823;
            SSUnity[515] = 1.90589;
            LSUnity[515] = 2.51375;
            ESUnity[516] = 2.16627;
            SSUnity[516] = 1.89749;
            LSUnity[516] = 2.52358;
            ESUnity[517] = 2.16182;
            SSUnity[517] = 1.91055;
            LSUnity[517] = 2.5186;
            ESUnity[518] = 2.18889;
            SSUnity[518] = 1.88218;
            LSUnity[518] = 2.492;
            ESUnity[519] = 2.158;
            SSUnity[519] = 1.90668;
            LSUnity[519] = 2.53665;
            ESUnity[520] = 2.18684;
            SSUnity[520] = 1.90062;
            LSUnity[520] = 2.50627;
            ESUnity[521] = 2.19179;
            SSUnity[521] = 1.88645;
            LSUnity[521] = 2.50134;
            ESUnity[522] = 2.20725;
            SSUnity[522] = 1.91537;
            LSUnity[522] = 2.52079;
            ESUnity[523] = 2.1885;
            SSUnity[523] = 1.88989;
            LSUnity[523] = 2.52392;
            ESUnity[524] = 2.18414;
            SSUnity[524] = 1.87472;
            LSUnity[524] = 2.50891;
            ESUnity[525] = 2.16727;
            SSUnity[525] = 1.90244;
            LSUnity[525] = 2.52027;
            ESUnity[526] = 2.19537;
            SSUnity[526] = 1.91055;
            LSUnity[526] = 2.51075;
            ESUnity[527] = 2.19736;
            SSUnity[527] = 1.90391;
            LSUnity[527] = 2.52028;
            ESUnity[528] = 2.16647;
            SSUnity[528] = 1.8828;
            LSUnity[528] = 2.53086;
            ESUnity[529] = 2.19236;
            SSUnity[529] = 1.89432;
            LSUnity[529] = 2.52228;
            ESUnity[530] = 2.17266;
            SSUnity[530] = 1.90331;
            LSUnity[530] = 2.50621;
            ESUnity[531] = 2.2053;
            SSUnity[531] = 1.88782;
            LSUnity[531] = 2.52003;
            ESUnity[532] = 2.18514;
            SSUnity[532] = 1.88651;
            LSUnity[532] = 2.55193;
            ESUnity[533] = 2.20062;
            SSUnity[533] = 1.91428;
            LSUnity[533] = 2.52878;
            ESUnity[534] = 2.20196;
            SSUnity[534] = 1.87854;
            LSUnity[534] = 2.50474;
            ESUnity[535] = 2.20121;
            SSUnity[535] = 1.88879;
            LSUnity[535] = 2.50677;
            ESUnity[536] = 2.20077;
            SSUnity[536] = 1.90204;
            LSUnity[536] = 2.52377;
            ESUnity[537] = 2.18712;
            SSUnity[537] = 1.88806;
            LSUnity[537] = 2.52075;
            ESUnity[538] = 2.1855;
            SSUnity[538] = 1.9132;
            LSUnity[538] = 2.49925;
            ESUnity[539] = 2.203;
            SSUnity[539] = 1.87505;
            LSUnity[539] = 2.52755;
            ESUnity[540] = 2.19592;
            SSUnity[540] = 1.89064;
            LSUnity[540] = 2.52391;
            ESUnity[541] = 2.19697;
            SSUnity[541] = 1.89937;
            LSUnity[541] = 2.51727;
            ESUnity[542] = 2.19238;
            SSUnity[542] = 1.90317;
            LSUnity[542] = 2.53087;
            ESUnity[543] = 2.21358;
            SSUnity[543] = 1.87959;
            LSUnity[543] = 2.54322;
            ESUnity[544] = 2.20131;
            SSUnity[544] = 1.87987;
            LSUnity[544] = 2.53557;
            ESUnity[545] = 2.20271;
            SSUnity[545] = 1.91334;
            LSUnity[545] = 2.52718;
            ESUnity[546] = 2.2132;
            SSUnity[546] = 1.89087;
            LSUnity[546] = 2.53461;
            ESUnity[547] = 2.20322;
            SSUnity[547] = 1.90313;
            LSUnity[547] = 2.52847;
            ESUnity[548] = 2.18998;
            SSUnity[548] = 1.88701;
            LSUnity[548] = 2.53709;
            ESUnity[549] = 2.18898;
            SSUnity[549] = 1.88566;
            LSUnity[549] = 2.50712;
            ESUnity[550] = 2.2032;
            SSUnity[550] = 1.91023;
            LSUnity[550] = 2.53503;
            ESUnity[551] = 2.21105;
            SSUnity[551] = 1.89054;
            LSUnity[551] = 2.52062;
            ESUnity[552] = 2.20624;
            SSUnity[552] = 1.8945;
            LSUnity[552] = 2.53159;
            ESUnity[553] = 2.21042;
            SSUnity[553] = 1.89948;
            LSUnity[553] = 2.53295;
            ESUnity[554] = 2.2051;
            SSUnity[554] = 1.89636;
            LSUnity[554] = 2.50698;
            ESUnity[555] = 2.21223;
            SSUnity[555] = 1.88011;
            LSUnity[555] = 2.5267;
            ESUnity[556] = 2.25237;
            SSUnity[556] = 1.8943;
            LSUnity[556] = 2.51318;
            ESUnity[557] = 2.21964;
            SSUnity[557] = 1.89894;
            LSUnity[557] = 2.53213;
            ESUnity[558] = 2.20576;
            SSUnity[558] = 1.89616;
            LSUnity[558] = 2.51027;
            ESUnity[559] = 2.19588;
            SSUnity[559] = 1.89537;
            LSUnity[559] = 2.54457;
            ESUnity[560] = 2.20689;
            SSUnity[560] = 1.88978;
            LSUnity[560] = 2.54982;
            ESUnity[561] = 2.20957;
            SSUnity[561] = 1.88551;
            LSUnity[561] = 2.53728;
            ESUnity[562] = 2.20479;
            SSUnity[562] = 1.90539;
            LSUnity[562] = 2.53525;
            ESUnity[563] = 2.21863;
            SSUnity[563] = 1.89958;
            LSUnity[563] = 2.53594;
            ESUnity[564] = 2.20926;
            SSUnity[564] = 1.88902;
            LSUnity[564] = 2.5311;
            ESUnity[565] = 2.19959;
            SSUnity[565] = 1.8978;
            LSUnity[565] = 2.51435;
            ESUnity[566] = 2.22527;
            SSUnity[566] = 1.88383;
            LSUnity[566] = 2.54877;
            ESUnity[567] = 2.23337;
            SSUnity[567] = 1.89353;
            LSUnity[567] = 2.53915;
            ESUnity[568] = 2.23426;
            SSUnity[568] = 1.87327;
            LSUnity[568] = 2.52995;
            ESUnity[569] = 2.21943;
            SSUnity[569] = 1.89261;
            LSUnity[569] = 2.52399;
            ESUnity[570] = 2.21984;
            SSUnity[570] = 1.89132;
            LSUnity[570] = 2.54886;
            ESUnity[571] = 2.23413;
            SSUnity[571] = 1.91017;
            LSUnity[571] = 2.55376;
            ESUnity[572] = 2.20928;
            SSUnity[572] = 1.8932;
            LSUnity[572] = 2.56368;
            ESUnity[573] = 2.21516;
            SSUnity[573] = 1.90353;
            LSUnity[573] = 2.5365;
            ESUnity[574] = 2.21589;
            SSUnity[574] = 1.89165;
            LSUnity[574] = 2.52975;
            ESUnity[575] = 2.21113;
            SSUnity[575] = 1.90393;
            LSUnity[575] = 2.55837;
            ESUnity[576] = 2.23569;
            SSUnity[576] = 1.89115;
            LSUnity[576] = 2.56132;
            ESUnity[577] = 2.23778;
            SSUnity[577] = 1.89779;
            LSUnity[577] = 2.54107;
            ESUnity[578] = 2.22916;
            SSUnity[578] = 1.90125;
            LSUnity[578] = 2.53877;
            ESUnity[579] = 2.23161;
            SSUnity[579] = 1.89697;
            LSUnity[579] = 2.53377;
            ESUnity[580] = 2.23586;
            SSUnity[580] = 1.88323;
            LSUnity[580] = 2.5597;
            ESUnity[581] = 2.2355;
            SSUnity[581] = 1.8885;
            LSUnity[581] = 2.52129;
            ESUnity[582] = 2.21636;
            SSUnity[582] = 1.89575;
            LSUnity[582] = 2.5522;
            ESUnity[583] = 2.21727;
            SSUnity[583] = 1.91027;
            LSUnity[583] = 2.52495;
            ESUnity[584] = 2.23536;
            SSUnity[584] = 1.89907;
            LSUnity[584] = 2.53839;
            ESUnity[585] = 2.24593;
            SSUnity[585] = 1.88814;
            LSUnity[585] = 2.574;
            ESUnity[586] = 2.23659;
            SSUnity[586] = 1.90582;
            LSUnity[586] = 2.52807;
            ESUnity[587] = 2.24643;
            SSUnity[587] = 1.88155;
            LSUnity[587] = 2.55233;
            ESUnity[588] = 2.22351;
            SSUnity[588] = 1.90385;
            LSUnity[588] = 2.55734;
            ESUnity[589] = 2.23788;
            SSUnity[589] = 1.88722;
            LSUnity[589] = 2.55201;
            ESUnity[590] = 2.22869;
            SSUnity[590] = 1.89629;
            LSUnity[590] = 2.56279;
            ESUnity[591] = 2.23329;
            SSUnity[591] = 1.89471;
            LSUnity[591] = 2.54147;
            ESUnity[592] = 2.23923;
            SSUnity[592] = 1.90323;
            LSUnity[592] = 2.54816;
            ESUnity[593] = 2.23153;
            SSUnity[593] = 1.91022;
            LSUnity[593] = 2.5406;
            ESUnity[594] = 2.22083;
            SSUnity[594] = 1.87994;
            LSUnity[594] = 2.54299;
            ESUnity[595] = 2.23567;
            SSUnity[595] = 1.88495;
            LSUnity[595] = 2.53876;
            ESUnity[596] = 2.24503;
            SSUnity[596] = 1.90931;
            LSUnity[596] = 2.54489;
            ESUnity[597] = 2.23027;
            SSUnity[597] = 1.89433;
            LSUnity[597] = 2.54608;
            ESUnity[598] = 2.25075;
            SSUnity[598] = 1.89826;
            LSUnity[598] = 2.52993;
            ESUnity[599] = 2.2466;
            SSUnity[599] = 1.85809;
            LSUnity[599] = 2.5548;
            ESUnity[600] = 2.2261;
            SSUnity[600] = 1.89566;
            LSUnity[600] = 2.54196;
            ESUnity[601] = 2.26533;
            SSUnity[601] = 1.90934;
            LSUnity[601] = 2.55706;
            ESUnity[602] = 2.23839;
            SSUnity[602] = 1.89164;
            LSUnity[602] = 2.56675;
            ESUnity[603] = 2.24523;
            SSUnity[603] = 1.88081;
            LSUnity[603] = 2.55882;
            ESUnity[604] = 2.25686;
            SSUnity[604] = 1.8992;
            LSUnity[604] = 2.54778;
            ESUnity[605] = 2.25003;
            SSUnity[605] = 1.87815;
            LSUnity[605] = 2.55636;
            ESUnity[606] = 2.26011;
            SSUnity[606] = 1.89481;
            LSUnity[606] = 2.56739;
            ESUnity[607] = 2.23715;
            SSUnity[607] = 1.90496;
            LSUnity[607] = 2.55791;
            ESUnity[608] = 2.24906;
            SSUnity[608] = 1.87239;
            LSUnity[608] = 2.54213;
            ESUnity[609] = 2.2572;
            SSUnity[609] = 1.89234;
            LSUnity[609] = 2.57212;
            ESUnity[610] = 2.24812;
            SSUnity[610] = 1.89378;
            LSUnity[610] = 2.5587;
            ESUnity[611] = 2.26422;
            SSUnity[611] = 1.89824;
            LSUnity[611] = 2.54781;
            ESUnity[612] = 2.25018;
            SSUnity[612] = 1.87987;
            LSUnity[612] = 2.55629;
            ESUnity[613] = 2.23922;
            SSUnity[613] = 1.88727;
            LSUnity[613] = 2.55217;
            ESUnity[614] = 2.26886;
            SSUnity[614] = 1.88735;
            LSUnity[614] = 2.56069;
            ESUnity[615] = 2.249;
            SSUnity[615] = 1.8909;
            LSUnity[615] = 2.57395;
            ESUnity[616] = 2.25936;
            SSUnity[616] = 1.88136;
            LSUnity[616] = 2.55695;
            ESUnity[617] = 2.27496;
            SSUnity[617] = 1.88493;
            LSUnity[617] = 2.58945;
            ESUnity[618] = 2.27758;
            SSUnity[618] = 1.89821;
            LSUnity[618] = 2.57694;
            ESUnity[619] = 2.24577;
            SSUnity[619] = 1.89608;
            LSUnity[619] = 2.57801;
            ESUnity[620] = 2.254;
            SSUnity[620] = 1.90322;
            LSUnity[620] = 2.56366;
            ESUnity[621] = 2.25192;
            SSUnity[621] = 1.89722;
            LSUnity[621] = 2.55821;
            ESUnity[622] = 2.26225;
            SSUnity[622] = 1.90471;
            LSUnity[622] = 2.55877;
            ESUnity[623] = 2.26882;
            SSUnity[623] = 1.90275;
            LSUnity[623] = 2.56047;
            ESUnity[624] = 2.25137;
            SSUnity[624] = 1.90514;
            LSUnity[624] = 2.58122;
            ESUnity[625] = 2.25483;
            SSUnity[625] = 1.89977;
            LSUnity[625] = 2.54434;
            ESUnity[626] = 2.29103;
            SSUnity[626] = 1.89659;
            LSUnity[626] = 2.57663;
            ESUnity[627] = 2.26272;
            SSUnity[627] = 1.8967;
            LSUnity[627] = 2.56465;
            ESUnity[628] = 2.26391;
            SSUnity[628] = 1.89102;
            LSUnity[628] = 2.56682;
            ESUnity[629] = 2.2784;
            SSUnity[629] = 1.89847;
            LSUnity[629] = 2.5836;
            ESUnity[630] = 2.2738;
            SSUnity[630] = 1.87468;
            LSUnity[630] = 2.5914;
            ESUnity[631] = 2.26043;
            SSUnity[631] = 1.89986;
            LSUnity[631] = 2.58477;
            ESUnity[632] = 2.25974;
            SSUnity[632] = 1.89101;
            LSUnity[632] = 2.56708;
            ESUnity[633] = 2.25427;
            SSUnity[633] = 1.90189;
            LSUnity[633] = 2.58162;
            ESUnity[634] = 2.26472;
            SSUnity[634] = 1.89183;
            LSUnity[634] = 2.58168;
            ESUnity[635] = 2.26555;
            SSUnity[635] = 1.88167;
            LSUnity[635] = 2.54164;
            ESUnity[636] = 2.27494;
            SSUnity[636] = 1.88848;
            LSUnity[636] = 2.60683;
            ESUnity[637] = 2.28194;
            SSUnity[637] = 1.8998;
            LSUnity[637] = 2.5715;
            ESUnity[638] = 2.28255;
            SSUnity[638] = 1.90078;
            LSUnity[638] = 2.55917;
            ESUnity[639] = 2.26372;
            SSUnity[639] = 1.88937;
            LSUnity[639] = 2.5675;
            ESUnity[640] = 2.27794;
            SSUnity[640] = 1.89129;
            LSUnity[640] = 2.58735;
            ESUnity[641] = 2.2836;
            SSUnity[641] = 1.89438;
            LSUnity[641] = 2.57018;
            ESUnity[642] = 2.28825;
            SSUnity[642] = 1.88236;
            LSUnity[642] = 2.58657;
            ESUnity[643] = 2.29874;
            SSUnity[643] = 1.88634;
            LSUnity[643] = 2.60211;
            ESUnity[644] = 2.28286;
            SSUnity[644] = 1.90459;
            LSUnity[644] = 2.55943;
            ESUnity[645] = 2.27046;
            SSUnity[645] = 1.88899;
            LSUnity[645] = 2.58734;
            ESUnity[646] = 2.29208;
            SSUnity[646] = 1.9121;
            LSUnity[646] = 2.57386;
            ESUnity[647] = 2.26638;
            SSUnity[647] = 1.89335;
            LSUnity[647] = 2.54582;
            ESUnity[648] = 2.27392;
            SSUnity[648] = 1.88067;
            LSUnity[648] = 2.55171;
            ESUnity[649] = 2.28982;
            SSUnity[649] = 1.89845;
            LSUnity[649] = 2.59369;
            ESUnity[650] = 2.27952;
            SSUnity[650] = 1.90379;
            LSUnity[650] = 2.58087;
            ESUnity[651] = 2.29714;
            SSUnity[651] = 1.89761;
            LSUnity[651] = 2.56467;
            ESUnity[652] = 2.27305;
            SSUnity[652] = 1.88957;
            LSUnity[652] = 2.58683;
            ESUnity[653] = 2.28039;
            SSUnity[653] = 1.89529;
            LSUnity[653] = 2.57931;
            ESUnity[654] = 2.29303;
            SSUnity[654] = 1.89919;
            LSUnity[654] = 2.56709;
            ESUnity[655] = 2.29785;
            SSUnity[655] = 1.89856;
            LSUnity[655] = 2.59857;
            ESUnity[656] = 2.28801;
            SSUnity[656] = 1.87728;
            LSUnity[656] = 2.59497;
            ESUnity[657] = 2.274;
            SSUnity[657] = 1.90058;
            LSUnity[657] = 2.60257;
            ESUnity[658] = 2.29605;
            SSUnity[658] = 1.8912;
            LSUnity[658] = 2.59976;
            ESUnity[659] = 2.29291;
            SSUnity[659] = 1.89529;
            LSUnity[659] = 2.58045;
            ESUnity[660] = 2.28823;
            SSUnity[660] = 1.90655;
            LSUnity[660] = 2.60188;
            ESUnity[661] = 2.27533;
            SSUnity[661] = 1.89635;
            LSUnity[661] = 2.57599;
            ESUnity[662] = 2.30257;
            SSUnity[662] = 1.89044;
            LSUnity[662] = 2.58117;
            ESUnity[663] = 2.2948;
            SSUnity[663] = 1.90135;
            LSUnity[663] = 2.57755;
            ESUnity[664] = 2.29129;
            SSUnity[664] = 1.89168;
            LSUnity[664] = 2.5852;
            ESUnity[665] = 2.29115;
            SSUnity[665] = 1.88689;
            LSUnity[665] = 2.55517;
            ESUnity[666] = 2.29554;
            SSUnity[666] = 1.88565;
            LSUnity[666] = 2.58127;
            ESUnity[667] = 2.30859;
            SSUnity[667] = 1.88905;
            LSUnity[667] = 2.56754;
            ESUnity[668] = 2.29487;
            SSUnity[668] = 1.88723;
            LSUnity[668] = 2.59936;
            ESUnity[669] = 2.29137;
            SSUnity[669] = 1.87604;
            LSUnity[669] = 2.54267;
            ESUnity[670] = 2.29391;
            SSUnity[670] = 1.91439;
            LSUnity[670] = 2.57909;
            ESUnity[671] = 2.29261;
            SSUnity[671] = 1.90133;
            LSUnity[671] = 2.58412;
            ESUnity[672] = 2.30659;
            SSUnity[672] = 1.90071;
            LSUnity[672] = 2.61833;
            ESUnity[673] = 2.27934;
            SSUnity[673] = 1.90261;
            LSUnity[673] = 2.60542;
            ESUnity[674] = 2.30991;
            SSUnity[674] = 1.89697;
            LSUnity[674] = 2.60639;
            ESUnity[675] = 2.28966;
            SSUnity[675] = 1.89686;
            LSUnity[675] = 2.60918;
            ESUnity[676] = 2.3074;
            SSUnity[676] = 1.90015;
            LSUnity[676] = 2.56979;
            ESUnity[677] = 2.29334;
            SSUnity[677] = 1.90973;
            LSUnity[677] = 2.60049;
            ESUnity[678] = 2.30214;
            SSUnity[678] = 1.89383;
            LSUnity[678] = 2.55136;
            ESUnity[679] = 2.30023;
            SSUnity[679] = 1.89142;
            LSUnity[679] = 2.58865;
            ESUnity[680] = 2.29928;
            SSUnity[680] = 1.8829;
            LSUnity[680] = 2.58476;
            ESUnity[681] = 2.30966;
            SSUnity[681] = 1.87547;
            LSUnity[681] = 2.56743;
            ESUnity[682] = 2.29758;
            SSUnity[682] = 1.88567;
            LSUnity[682] = 2.59835;
            ESUnity[683] = 2.30252;
            SSUnity[683] = 1.87654;
            LSUnity[683] = 2.58465;
            ESUnity[684] = 2.3262;
            SSUnity[684] = 1.90044;
            LSUnity[684] = 2.6009;
            ESUnity[685] = 2.31343;
            SSUnity[685] = 1.89668;
            LSUnity[685] = 2.60052;
            ESUnity[686] = 2.31017;
            SSUnity[686] = 1.88688;
            LSUnity[686] = 2.60446;
            ESUnity[687] = 2.30876;
            SSUnity[687] = 1.8951;
            LSUnity[687] = 2.61713;
            ESUnity[688] = 2.30309;
            SSUnity[688] = 1.90837;
            LSUnity[688] = 2.60886;
            ESUnity[689] = 2.30082;
            SSUnity[689] = 1.90282;
            LSUnity[689] = 2.58914;
            ESUnity[690] = 2.30767;
            SSUnity[690] = 1.8893;
            LSUnity[690] = 2.58108;
            ESUnity[691] = 2.30642;
            SSUnity[691] = 1.88972;
            LSUnity[691] = 2.61292;
            ESUnity[692] = 2.31165;
            SSUnity[692] = 1.89732;
            LSUnity[692] = 2.5944;
            ESUnity[693] = 2.33401;
            SSUnity[693] = 1.89699;
            LSUnity[693] = 2.58891;
            ESUnity[694] = 2.31506;
            SSUnity[694] = 1.89371;
            LSUnity[694] = 2.58533;
            ESUnity[695] = 2.30913;
            SSUnity[695] = 1.91153;
            LSUnity[695] = 2.61508;
            ESUnity[696] = 2.30838;
            SSUnity[696] = 1.90202;
            LSUnity[696] = 2.59526;
            ESUnity[697] = 2.31412;
            SSUnity[697] = 1.89278;
            LSUnity[697] = 2.59377;
            ESUnity[698] = 2.31942;
            SSUnity[698] = 1.88933;
            LSUnity[698] = 2.59972;
            ESUnity[699] = 2.32585;
            SSUnity[699] = 1.88005;
            LSUnity[699] = 2.61878;
            ESUnity[700] = 2.31608;
            SSUnity[700] = 1.9068;
            LSUnity[700] = 2.60688;
            ESUnity[701] = 2.34261;
            SSUnity[701] = 1.89558;
            LSUnity[701] = 2.6091;
            ESUnity[702] = 2.32858;
            SSUnity[702] = 1.89649;
            LSUnity[702] = 2.59368;
            ESUnity[703] = 2.33181;
            SSUnity[703] = 1.89931;
            LSUnity[703] = 2.61317;
            ESUnity[704] = 2.31739;
            SSUnity[704] = 1.89903;
            LSUnity[704] = 2.59413;
            ESUnity[705] = 2.32328;
            SSUnity[705] = 1.89768;
            LSUnity[705] = 2.60147;
            ESUnity[706] = 2.32405;
            SSUnity[706] = 1.8823;
            LSUnity[706] = 2.60645;
            ESUnity[707] = 2.34106;
            SSUnity[707] = 1.88875;
            LSUnity[707] = 2.61759;
            ESUnity[708] = 2.33121;
            SSUnity[708] = 1.89121;
            LSUnity[708] = 2.60378;
            ESUnity[709] = 2.35545;
            SSUnity[709] = 1.90444;
            LSUnity[709] = 2.6127;
            ESUnity[710] = 2.34668;
            SSUnity[710] = 1.90323;
            LSUnity[710] = 2.60983;
            ESUnity[711] = 2.33791;
            SSUnity[711] = 1.87558;
            LSUnity[711] = 2.61437;
            ESUnity[712] = 2.36319;
            SSUnity[712] = 1.90504;
            LSUnity[712] = 2.60514;
            ESUnity[713] = 2.33368;
            SSUnity[713] = 1.90253;
            LSUnity[713] = 2.62028;
            ESUnity[714] = 2.32436;
            SSUnity[714] = 1.90053;
            LSUnity[714] = 2.58728;
            ESUnity[715] = 2.33373;
            SSUnity[715] = 1.89214;
            LSUnity[715] = 2.60484;
            ESUnity[716] = 2.31906;
            SSUnity[716] = 1.90118;
            LSUnity[716] = 2.5839;
            ESUnity[717] = 2.32447;
            SSUnity[717] = 1.88344;
            LSUnity[717] = 2.61168;
            ESUnity[718] = 2.33932;
            SSUnity[718] = 1.89658;
            LSUnity[718] = 2.59679;
            ESUnity[719] = 2.3463;
            SSUnity[719] = 1.89599;
            LSUnity[719] = 2.63831;
            ESUnity[720] = 2.31751;
            SSUnity[720] = 1.89615;
            LSUnity[720] = 2.60361;
            ESUnity[721] = 2.32482;
            SSUnity[721] = 1.88639;
            LSUnity[721] = 2.59929;
            ESUnity[722] = 2.34584;
            SSUnity[722] = 1.88376;
            LSUnity[722] = 2.62249;
            ESUnity[723] = 2.31332;
            SSUnity[723] = 1.88869;
            LSUnity[723] = 2.6155;
            ESUnity[724] = 2.33577;
            SSUnity[724] = 1.89978;
            LSUnity[724] = 2.62211;
            ESUnity[725] = 2.32814;
            SSUnity[725] = 1.90987;
            LSUnity[725] = 2.63329;
            ESUnity[726] = 2.33068;
            SSUnity[726] = 1.88091;
            LSUnity[726] = 2.62473;
            ESUnity[727] = 2.33119;
            SSUnity[727] = 1.91121;
            LSUnity[727] = 2.62245;
            ESUnity[728] = 2.33357;
            SSUnity[728] = 1.88613;
            LSUnity[728] = 2.5998;
            ESUnity[729] = 2.32783;
            SSUnity[729] = 1.88337;
            LSUnity[729] = 2.61497;
            ESUnity[730] = 2.33489;
            SSUnity[730] = 1.89822;
            LSUnity[730] = 2.62168;
            ESUnity[731] = 2.33765;
            SSUnity[731] = 1.88458;
            LSUnity[731] = 2.62665;
            ESUnity[732] = 2.32657;
            SSUnity[732] = 1.91133;
            LSUnity[732] = 2.62334;
            ESUnity[733] = 2.34187;
            SSUnity[733] = 1.88807;
            LSUnity[733] = 2.61483;
            ESUnity[734] = 2.34125;
            SSUnity[734] = 1.89796;
            LSUnity[734] = 2.60848;
            ESUnity[735] = 2.33287;
            SSUnity[735] = 1.89759;
            LSUnity[735] = 2.63601;
            ESUnity[736] = 2.36289;
            SSUnity[736] = 1.90655;
            LSUnity[736] = 2.63043;
            ESUnity[737] = 2.37368;
            SSUnity[737] = 1.88328;
            LSUnity[737] = 2.63654;
            ESUnity[738] = 2.33286;
            SSUnity[738] = 1.89338;
            LSUnity[738] = 2.64424;
            ESUnity[739] = 2.36672;
            SSUnity[739] = 1.90343;
            LSUnity[739] = 2.60307;
            ESUnity[740] = 2.32112;
            SSUnity[740] = 1.89968;
            LSUnity[740] = 2.63131;
            ESUnity[741] = 2.35406;
            SSUnity[741] = 1.8861;
            LSUnity[741] = 2.61724;
            ESUnity[742] = 2.37977;
            SSUnity[742] = 1.90912;
            LSUnity[742] = 2.63227;
            ESUnity[743] = 2.3474;
            SSUnity[743] = 1.89127;
            LSUnity[743] = 2.64063;
            ESUnity[744] = 2.34621;
            SSUnity[744] = 1.87304;
            LSUnity[744] = 2.64034;
            ESUnity[745] = 2.35117;
            SSUnity[745] = 1.88516;
            LSUnity[745] = 2.64152;
            ESUnity[746] = 2.33558;
            SSUnity[746] = 1.88347;
            LSUnity[746] = 2.61469;
            ESUnity[747] = 2.36639;
            SSUnity[747] = 1.91319;
            LSUnity[747] = 2.62656;
            ESUnity[748] = 2.343;
            SSUnity[748] = 1.89017;
            LSUnity[748] = 2.6106;
            ESUnity[749] = 2.33953;
            SSUnity[749] = 1.8873;
            LSUnity[749] = 2.60489;
            ESUnity[750] = 2.36846;
            SSUnity[750] = 1.90613;
            LSUnity[750] = 2.60702;
            ESUnity[751] = 2.36019;
            SSUnity[751] = 1.89508;
            LSUnity[751] = 2.62082;
            ESUnity[752] = 2.36664;
            SSUnity[752] = 1.87199;
            LSUnity[752] = 2.61918;
            ESUnity[753] = 2.3685;
            SSUnity[753] = 1.90935;
            LSUnity[753] = 2.63803;
            ESUnity[754] = 2.34536;
            SSUnity[754] = 1.90179;
            LSUnity[754] = 2.64378;
            ESUnity[755] = 2.38779;
            SSUnity[755] = 1.89347;
            LSUnity[755] = 2.63804;
            ESUnity[756] = 2.34615;
            SSUnity[756] = 1.88066;
            LSUnity[756] = 2.60679;
            ESUnity[757] = 2.35924;
            SSUnity[757] = 1.88283;
            LSUnity[757] = 2.63425;
            ESUnity[758] = 2.35619;
            SSUnity[758] = 1.91498;
            LSUnity[758] = 2.65165;
            ESUnity[759] = 2.35674;
            SSUnity[759] = 1.9015;
            LSUnity[759] = 2.66037;
            ESUnity[760] = 2.38124;
            SSUnity[760] = 1.886;
            LSUnity[760] = 2.61499;
            ESUnity[761] = 2.35764;
            SSUnity[761] = 1.90733;
            LSUnity[761] = 2.64311;
            ESUnity[762] = 2.3527;
            SSUnity[762] = 1.89976;
            LSUnity[762] = 2.64223;
            ESUnity[763] = 2.37836;
            SSUnity[763] = 1.88605;
            LSUnity[763] = 2.62837;
            ESUnity[764] = 2.35709;
            SSUnity[764] = 1.89132;
            LSUnity[764] = 2.63828;
            ESUnity[765] = 2.37681;
            SSUnity[765] = 1.8883;
            LSUnity[765] = 2.63047;
            ESUnity[766] = 2.3599;
            SSUnity[766] = 1.90426;
            LSUnity[766] = 2.65818;
            ESUnity[767] = 2.35993;
            SSUnity[767] = 1.88447;
            LSUnity[767] = 2.63292;
            ESUnity[768] = 2.37643;
            SSUnity[768] = 1.89635;
            LSUnity[768] = 2.65188;
            ESUnity[769] = 2.34256;
            SSUnity[769] = 1.88897;
            LSUnity[769] = 2.64348;
            ESUnity[770] = 2.37494;
            SSUnity[770] = 1.90566;
            LSUnity[770] = 2.63508;
            ESUnity[771] = 2.36611;
            SSUnity[771] = 1.89833;
            LSUnity[771] = 2.62616;
            ESUnity[772] = 2.3678;
            SSUnity[772] = 1.90698;
            LSUnity[772] = 2.64721;
            ESUnity[773] = 2.38232;
            SSUnity[773] = 1.89831;
            LSUnity[773] = 2.64219;
            ESUnity[774] = 2.37103;
            SSUnity[774] = 1.88439;
            LSUnity[774] = 2.63744;
            ESUnity[775] = 2.38294;
            SSUnity[775] = 1.8969;
            LSUnity[775] = 2.63433;
            ESUnity[776] = 2.36588;
            SSUnity[776] = 1.89293;
            LSUnity[776] = 2.66162;
            ESUnity[777] = 2.36488;
            SSUnity[777] = 1.88194;
            LSUnity[777] = 2.62189;
            ESUnity[778] = 2.39838;
            SSUnity[778] = 1.87484;
            LSUnity[778] = 2.64546;
            ESUnity[779] = 2.37459;
            SSUnity[779] = 1.89074;
            LSUnity[779] = 2.66013;
            ESUnity[780] = 2.39982;
            SSUnity[780] = 1.90319;
            LSUnity[780] = 2.63064;
            ESUnity[781] = 2.37237;
            SSUnity[781] = 1.91112;
            LSUnity[781] = 2.64863;
            ESUnity[782] = 2.39043;
            SSUnity[782] = 1.88657;
            LSUnity[782] = 2.64538;
            ESUnity[783] = 2.37584;
            SSUnity[783] = 1.9033;
            LSUnity[783] = 2.63274;
            ESUnity[784] = 2.36484;
            SSUnity[784] = 1.89549;
            LSUnity[784] = 2.61465;
            ESUnity[785] = 2.36998;
            SSUnity[785] = 1.88448;
            LSUnity[785] = 2.63777;
            ESUnity[786] = 2.39262;
            SSUnity[786] = 1.90199;
            LSUnity[786] = 2.65034;
            ESUnity[787] = 2.39616;
            SSUnity[787] = 1.89799;
            LSUnity[787] = 2.63136;
            ESUnity[788] = 2.38953;
            SSUnity[788] = 1.87417;
            LSUnity[788] = 2.66236;
            ESUnity[789] = 2.38486;
            SSUnity[789] = 1.88978;
            LSUnity[789] = 2.63136;
            ESUnity[790] = 2.37668;
            SSUnity[790] = 1.8832;
            LSUnity[790] = 2.65056;
            ESUnity[791] = 2.40776;
            SSUnity[791] = 1.88876;
            LSUnity[791] = 2.65583;
            ESUnity[792] = 2.39609;
            SSUnity[792] = 1.89214;
            LSUnity[792] = 2.64835;
            ESUnity[793] = 2.38562;
            SSUnity[793] = 1.88684;
            LSUnity[793] = 2.64503;
            ESUnity[794] = 2.39695;
            SSUnity[794] = 1.90323;
            LSUnity[794] = 2.6594;
            ESUnity[795] = 2.38721;
            SSUnity[795] = 1.9013;
            LSUnity[795] = 2.63609;
            ESUnity[796] = 2.36337;
            SSUnity[796] = 1.89762;
            LSUnity[796] = 2.64269;
            ESUnity[797] = 2.4008;
            SSUnity[797] = 1.87905;
            LSUnity[797] = 2.647;
            ESUnity[798] = 2.37085;
            SSUnity[798] = 1.88232;
            LSUnity[798] = 2.65979;
            ESUnity[799] = 2.40936;
            SSUnity[799] = 1.90718;
            LSUnity[799] = 2.65428;
            ESUnity[800] = 2.37824;
            SSUnity[800] = 1.87753;
            LSUnity[800] = 2.65347;
            ESUnity[801] = 2.39494;
            SSUnity[801] = 1.9065;
            LSUnity[801] = 2.64373;
            ESUnity[802] = 2.37432;
            SSUnity[802] = 1.90041;
            LSUnity[802] = 2.65751;
            ESUnity[803] = 2.39886;
            SSUnity[803] = 1.88791;
            LSUnity[803] = 2.63515;
            ESUnity[804] = 2.40034;
            SSUnity[804] = 1.9134;
            LSUnity[804] = 2.65027;
            ESUnity[805] = 2.38783;
            SSUnity[805] = 1.8941;
            LSUnity[805] = 2.65868;
            ESUnity[806] = 2.37574;
            SSUnity[806] = 1.89439;
            LSUnity[806] = 2.65317;
            ESUnity[807] = 2.3829;
            SSUnity[807] = 1.90251;
            LSUnity[807] = 2.67415;
            ESUnity[808] = 2.39571;
            SSUnity[808] = 1.89448;
            LSUnity[808] = 2.65978;
            ESUnity[809] = 2.37724;
            SSUnity[809] = 1.90285;
            LSUnity[809] = 2.66258;
            ESUnity[810] = 2.38703;
            SSUnity[810] = 1.89068;
            LSUnity[810] = 2.66796;
            ESUnity[811] = 2.38019;
            SSUnity[811] = 1.88837;
            LSUnity[811] = 2.68042;
            ESUnity[812] = 2.38043;
            SSUnity[812] = 1.89769;
            LSUnity[812] = 2.66031;
            ESUnity[813] = 2.41066;
            SSUnity[813] = 1.88539;
            LSUnity[813] = 2.63989;
            ESUnity[814] = 2.37786;
            SSUnity[814] = 1.88133;
            LSUnity[814] = 2.66617;
            ESUnity[815] = 2.41555;
            SSUnity[815] = 1.90799;
            LSUnity[815] = 2.65777;
            ESUnity[816] = 2.38219;
            SSUnity[816] = 1.89069;
            LSUnity[816] = 2.65101;
            ESUnity[817] = 2.3991;
            SSUnity[817] = 1.89388;
            LSUnity[817] = 2.67553;
            ESUnity[818] = 2.39931;
            SSUnity[818] = 1.89132;
            LSUnity[818] = 2.66463;
            ESUnity[819] = 2.41541;
            SSUnity[819] = 1.90184;
            LSUnity[819] = 2.65496;
            ESUnity[820] = 2.40905;
            SSUnity[820] = 1.87038;
            LSUnity[820] = 2.65991;
            ESUnity[821] = 2.38504;
            SSUnity[821] = 1.90038;
            LSUnity[821] = 2.67406;
            ESUnity[822] = 2.37377;
            SSUnity[822] = 1.90289;
            LSUnity[822] = 2.66285;
            ESUnity[823] = 2.39863;
            SSUnity[823] = 1.89678;
            LSUnity[823] = 2.66322;
            ESUnity[824] = 2.39556;
            SSUnity[824] = 1.88616;
            LSUnity[824] = 2.66159;
            ESUnity[825] = 2.38955;
            SSUnity[825] = 1.90152;
            LSUnity[825] = 2.65272;
            ESUnity[826] = 2.40792;
            SSUnity[826] = 1.89762;
            LSUnity[826] = 2.66115;
            ESUnity[827] = 2.40224;
            SSUnity[827] = 1.89782;
            LSUnity[827] = 2.66574;
            ESUnity[828] = 2.41852;
            SSUnity[828] = 1.9074;
            LSUnity[828] = 2.65239;
            ESUnity[829] = 2.40423;
            SSUnity[829] = 1.89764;
            LSUnity[829] = 2.67235;
            ESUnity[830] = 2.42243;
            SSUnity[830] = 1.89548;
            LSUnity[830] = 2.6834;
            ESUnity[831] = 2.3974;
            SSUnity[831] = 1.88731;
            LSUnity[831] = 2.68162;
            ESUnity[832] = 2.42301;
            SSUnity[832] = 1.90852;
            LSUnity[832] = 2.66852;
            ESUnity[833] = 2.41379;
            SSUnity[833] = 1.90634;
            LSUnity[833] = 2.64011;
            ESUnity[834] = 2.40933;
            SSUnity[834] = 1.88099;
            LSUnity[834] = 2.65396;
            ESUnity[835] = 2.43269;
            SSUnity[835] = 1.91536;
            LSUnity[835] = 2.67284;
            ESUnity[836] = 2.40928;
            SSUnity[836] = 1.89966;
            LSUnity[836] = 2.65112;
            ESUnity[837] = 2.42571;
            SSUnity[837] = 1.87809;
            LSUnity[837] = 2.6639;
            ESUnity[838] = 2.40922;
            SSUnity[838] = 1.88994;
            LSUnity[838] = 2.66532;
            ESUnity[839] = 2.38589;
            SSUnity[839] = 1.88381;
            LSUnity[839] = 2.66533;
            ESUnity[840] = 2.43387;
            SSUnity[840] = 1.89081;
            LSUnity[840] = 2.67687;
            ESUnity[841] = 2.43566;
            SSUnity[841] = 1.88562;
            LSUnity[841] = 2.66272;
            ESUnity[842] = 2.40253;
            SSUnity[842] = 1.89853;
            LSUnity[842] = 2.67128;
            ESUnity[843] = 2.40858;
            SSUnity[843] = 1.90482;
            LSUnity[843] = 2.67501;
            ESUnity[844] = 2.43298;
            SSUnity[844] = 1.88468;
            LSUnity[844] = 2.67628;
            ESUnity[845] = 2.3988;
            SSUnity[845] = 1.90035;
            LSUnity[845] = 2.66219;
            ESUnity[846] = 2.41652;
            SSUnity[846] = 1.88939;
            LSUnity[846] = 2.6629;
            ESUnity[847] = 2.42128;
            SSUnity[847] = 1.89295;
            LSUnity[847] = 2.68378;
            ESUnity[848] = 2.41883;
            SSUnity[848] = 1.90462;
            LSUnity[848] = 2.68541;
            ESUnity[849] = 2.41381;
            SSUnity[849] = 1.89681;
            LSUnity[849] = 2.65607;
            ESUnity[850] = 2.42592;
            SSUnity[850] = 1.88103;
            LSUnity[850] = 2.68323;
            ESUnity[851] = 2.43009;
            SSUnity[851] = 1.88996;
            LSUnity[851] = 2.65984;
            ESUnity[852] = 2.41947;
            SSUnity[852] = 1.89812;
            LSUnity[852] = 2.69041;
            ESUnity[853] = 2.42634;
            SSUnity[853] = 1.90114;
            LSUnity[853] = 2.67527;
            ESUnity[854] = 2.43983;
            SSUnity[854] = 1.88874;
            LSUnity[854] = 2.6632;
            ESUnity[855] = 2.41319;
            SSUnity[855] = 1.90127;
            LSUnity[855] = 2.68841;
            ESUnity[856] = 2.4213;
            SSUnity[856] = 1.87943;
            LSUnity[856] = 2.67395;
            ESUnity[857] = 2.43972;
            SSUnity[857] = 1.88883;
            LSUnity[857] = 2.65894;
            ESUnity[858] = 2.40919;
            SSUnity[858] = 1.89715;
            LSUnity[858] = 2.65999;
            ESUnity[859] = 2.43407;
            SSUnity[859] = 1.90068;
            LSUnity[859] = 2.66519;
            ESUnity[860] = 2.43571;
            SSUnity[860] = 1.8923;
            LSUnity[860] = 2.6833;
            ESUnity[861] = 2.43616;
            SSUnity[861] = 1.88333;
            LSUnity[861] = 2.67351;
            ESUnity[862] = 2.42489;
            SSUnity[862] = 1.91323;
            LSUnity[862] = 2.68194;
            ESUnity[863] = 2.41693;
            SSUnity[863] = 1.89676;
            LSUnity[863] = 2.69945;
            ESUnity[864] = 2.42638;
            SSUnity[864] = 1.88366;
            LSUnity[864] = 2.66126;
            ESUnity[865] = 2.44588;
            SSUnity[865] = 1.89779;
            LSUnity[865] = 2.66342;
            ESUnity[866] = 2.40138;
            SSUnity[866] = 1.91161;
            LSUnity[866] = 2.65907;
            ESUnity[867] = 2.44067;
            SSUnity[867] = 1.88421;
            LSUnity[867] = 2.66376;
            ESUnity[868] = 2.41251;
            SSUnity[868] = 1.89144;
            LSUnity[868] = 2.66585;
            ESUnity[869] = 2.4253;
            SSUnity[869] = 1.90206;
            LSUnity[869] = 2.68072;
            ESUnity[870] = 2.44413;
            SSUnity[870] = 1.88643;
            LSUnity[870] = 2.69919;
            ESUnity[871] = 2.45559;
            SSUnity[871] = 1.90517;
            LSUnity[871] = 2.69388;
            ESUnity[872] = 2.4249;
            SSUnity[872] = 1.90919;
            LSUnity[872] = 2.67666;
            ESUnity[873] = 2.44772;
            SSUnity[873] = 1.89141;
            LSUnity[873] = 2.6774;
            ESUnity[874] = 2.42393;
            SSUnity[874] = 1.88577;
            LSUnity[874] = 2.67711;
            ESUnity[875] = 2.42194;
            SSUnity[875] = 1.89182;
            LSUnity[875] = 2.67085;
            ESUnity[876] = 2.42421;
            SSUnity[876] = 1.885;
            LSUnity[876] = 2.65603;
            ESUnity[877] = 2.43479;
            SSUnity[877] = 1.90218;
            LSUnity[877] = 2.67537;
            ESUnity[878] = 2.44285;
            SSUnity[878] = 1.89239;
            LSUnity[878] = 2.66262;
            ESUnity[879] = 2.43805;
            SSUnity[879] = 1.90625;
            LSUnity[879] = 2.7217;
            ESUnity[880] = 2.43972;
            SSUnity[880] = 1.88104;
            LSUnity[880] = 2.6979;
            ESUnity[881] = 2.44475;
            SSUnity[881] = 1.89626;
            LSUnity[881] = 2.6986;
            ESUnity[882] = 2.44152;
            SSUnity[882] = 1.89945;
            LSUnity[882] = 2.67352;
            ESUnity[883] = 2.43006;
            SSUnity[883] = 1.88308;
            LSUnity[883] = 2.67068;
            ESUnity[884] = 2.42564;
            SSUnity[884] = 1.90484;
            LSUnity[884] = 2.69057;
            ESUnity[885] = 2.4388;
            SSUnity[885] = 1.89541;
            LSUnity[885] = 2.67911;
            ESUnity[886] = 2.46016;
            SSUnity[886] = 1.87879;
            LSUnity[886] = 2.66769;
            ESUnity[887] = 2.44581;
            SSUnity[887] = 1.89445;
            LSUnity[887] = 2.67372;
            ESUnity[888] = 2.45389;
            SSUnity[888] = 1.89553;
            LSUnity[888] = 2.68415;
            ESUnity[889] = 2.44211;
            SSUnity[889] = 1.88854;
            LSUnity[889] = 2.70737;
            ESUnity[890] = 2.44729;
            SSUnity[890] = 1.88014;
            LSUnity[890] = 2.67506;
            ESUnity[891] = 2.42274;
            SSUnity[891] = 1.91178;
            LSUnity[891] = 2.67741;
            ESUnity[892] = 2.46739;
            SSUnity[892] = 1.91681;
            LSUnity[892] = 2.68076;
            ESUnity[893] = 2.4538;
            SSUnity[893] = 1.90214;
            LSUnity[893] = 2.68482;
            ESUnity[894] = 2.43008;
            SSUnity[894] = 1.89259;
            LSUnity[894] = 2.6812;
            ESUnity[895] = 2.43388;
            SSUnity[895] = 1.89686;
            LSUnity[895] = 2.6784;
            ESUnity[896] = 2.44605;
            SSUnity[896] = 1.88722;
            LSUnity[896] = 2.67889;
            ESUnity[897] = 2.4527;
            SSUnity[897] = 1.89506;
            LSUnity[897] = 2.69316;
            ESUnity[898] = 2.4542;
            SSUnity[898] = 1.88619;
            LSUnity[898] = 2.69867;
            ESUnity[899] = 2.45912;
            SSUnity[899] = 1.89281;
            LSUnity[899] = 2.72472;
            ESUnity[900] = 2.43721;
            SSUnity[900] = 1.89551;
            LSUnity[900] = 2.69523;
            ESUnity[901] = 2.44794;
            SSUnity[901] = 1.8855;
            LSUnity[901] = 2.68489;
            ESUnity[902] = 2.44753;
            SSUnity[902] = 1.88574;
            LSUnity[902] = 2.67129;
            ESUnity[903] = 2.45321;
            SSUnity[903] = 1.90563;
            LSUnity[903] = 2.68269;
            ESUnity[904] = 2.4611;
            SSUnity[904] = 1.90001;
            LSUnity[904] = 2.67858;
            ESUnity[905] = 2.45027;
            SSUnity[905] = 1.88354;
            LSUnity[905] = 2.69434;
            ESUnity[906] = 2.42977;
            SSUnity[906] = 1.88374;
            LSUnity[906] = 2.69306;
            ESUnity[907] = 2.45646;
            SSUnity[907] = 1.90892;
            LSUnity[907] = 2.70113;
            ESUnity[908] = 2.43959;
            SSUnity[908] = 1.87765;
            LSUnity[908] = 2.71862;
            ESUnity[909] = 2.44658;
            SSUnity[909] = 1.89739;
            LSUnity[909] = 2.71305;
            ESUnity[910] = 2.46308;
            SSUnity[910] = 1.90702;
            LSUnity[910] = 2.68519;
            ESUnity[911] = 2.4674;
            SSUnity[911] = 1.89368;
            LSUnity[911] = 2.70542;
            ESUnity[912] = 2.44387;
            SSUnity[912] = 1.89146;
            LSUnity[912] = 2.69408;
            ESUnity[913] = 2.47152;
            SSUnity[913] = 1.89516;
            LSUnity[913] = 2.68624;
            ESUnity[914] = 2.479;
            SSUnity[914] = 1.8783;
            LSUnity[914] = 2.68845;
            ESUnity[915] = 2.45405;
            SSUnity[915] = 1.89885;
            LSUnity[915] = 2.7162;
            ESUnity[916] = 2.4617;
            SSUnity[916] = 1.89899;
            LSUnity[916] = 2.71963;
            ESUnity[917] = 2.48285;
            SSUnity[917] = 1.88313;
            LSUnity[917] = 2.70849;
            ESUnity[918] = 2.46095;
            SSUnity[918] = 1.90027;
            LSUnity[918] = 2.7113;
            ESUnity[919] = 2.45023;
            SSUnity[919] = 1.90857;
            LSUnity[919] = 2.70849;
            ESUnity[920] = 2.45919;
            SSUnity[920] = 1.89404;
            LSUnity[920] = 2.69117;
            ESUnity[921] = 2.4687;
            SSUnity[921] = 1.90807;
            LSUnity[921] = 2.70271;
            ESUnity[922] = 2.47058;
            SSUnity[922] = 1.90308;
            LSUnity[922] = 2.69867;
            ESUnity[923] = 2.45913;
            SSUnity[923] = 1.87233;
            LSUnity[923] = 2.68649;
            ESUnity[924] = 2.48053;
            SSUnity[924] = 1.88823;
            LSUnity[924] = 2.70945;
            ESUnity[925] = 2.47406;
            SSUnity[925] = 1.8987;
            LSUnity[925] = 2.70773;
            ESUnity[926] = 2.45633;
            SSUnity[926] = 1.87895;
            LSUnity[926] = 2.71815;
            ESUnity[927] = 2.47714;
            SSUnity[927] = 1.88939;
            LSUnity[927] = 2.68871;
            ESUnity[928] = 2.48744;
            SSUnity[928] = 1.8966;
            LSUnity[928] = 2.70891;
            ESUnity[929] = 2.46831;
            SSUnity[929] = 1.90515;
            LSUnity[929] = 2.69272;
            ESUnity[930] = 2.49996;
            SSUnity[930] = 1.8941;
            LSUnity[930] = 2.69751;
            ESUnity[931] = 2.48066;
            SSUnity[931] = 1.8914;
            LSUnity[931] = 2.69163;
            ESUnity[932] = 2.46329;
            SSUnity[932] = 1.88204;
            LSUnity[932] = 2.70256;
            ESUnity[933] = 2.47937;
            SSUnity[933] = 1.902;
            LSUnity[933] = 2.69113;
            ESUnity[934] = 2.49151;
            SSUnity[934] = 1.91323;
            LSUnity[934] = 2.72412;
            ESUnity[935] = 2.4485;
            SSUnity[935] = 1.89711;
            LSUnity[935] = 2.71537;
            ESUnity[936] = 2.44884;
            SSUnity[936] = 1.87751;
            LSUnity[936] = 2.706;
            ESUnity[937] = 2.4903;
            SSUnity[937] = 1.89697;
            LSUnity[937] = 2.70608;
            ESUnity[938] = 2.47449;
            SSUnity[938] = 1.89793;
            LSUnity[938] = 2.69244;
            ESUnity[939] = 2.49261;
            SSUnity[939] = 1.88654;
            LSUnity[939] = 2.71601;
            ESUnity[940] = 2.52088;
            SSUnity[940] = 1.88197;
            LSUnity[940] = 2.69698;
            ESUnity[941] = 2.47321;
            SSUnity[941] = 1.89874;
            LSUnity[941] = 2.69159;
            ESUnity[942] = 2.50311;
            SSUnity[942] = 1.89499;
            LSUnity[942] = 2.69689;
            ESUnity[943] = 2.47982;
            SSUnity[943] = 1.88469;
            LSUnity[943] = 2.70018;
            ESUnity[944] = 2.48593;
            SSUnity[944] = 1.88786;
            LSUnity[944] = 2.72669;
            ESUnity[945] = 2.48404;
            SSUnity[945] = 1.88063;
            LSUnity[945] = 2.70069;
            ESUnity[946] = 2.48011;
            SSUnity[946] = 1.89704;
            LSUnity[946] = 2.69244;
            ESUnity[947] = 2.46479;
            SSUnity[947] = 1.89433;
            LSUnity[947] = 2.7193;
            ESUnity[948] = 2.46686;
            SSUnity[948] = 1.89479;
            LSUnity[948] = 2.70408;
            ESUnity[949] = 2.47759;
            SSUnity[949] = 1.89705;
            LSUnity[949] = 2.74286;
            ESUnity[950] = 2.49516;
            SSUnity[950] = 1.8853;
            LSUnity[950] = 2.71983;
            ESUnity[951] = 2.46586;
            SSUnity[951] = 1.90603;
            LSUnity[951] = 2.714;
            ESUnity[952] = 2.48583;
            SSUnity[952] = 1.89423;
            LSUnity[952] = 2.70431;
            ESUnity[953] = 2.47507;
            SSUnity[953] = 1.91276;
            LSUnity[953] = 2.73146;
            ESUnity[954] = 2.49261;
            SSUnity[954] = 1.88896;
            LSUnity[954] = 2.71164;
            ESUnity[955] = 2.51366;
            SSUnity[955] = 1.89314;
            LSUnity[955] = 2.71262;
            ESUnity[956] = 2.49466;
            SSUnity[956] = 1.89921;
            LSUnity[956] = 2.73303;
            ESUnity[957] = 2.48851;
            SSUnity[957] = 1.88988;
            LSUnity[957] = 2.70071;
            ESUnity[958] = 2.48084;
            SSUnity[958] = 1.90855;
            LSUnity[958] = 2.76026;
            ESUnity[959] = 2.47907;
            SSUnity[959] = 1.8792;
            LSUnity[959] = 2.74374;
            ESUnity[960] = 2.49216;
            SSUnity[960] = 1.89072;
            LSUnity[960] = 2.71454;
            ESUnity[961] = 2.49586;
            SSUnity[961] = 1.88432;
            LSUnity[961] = 2.71946;
            ESUnity[962] = 2.48606;
            SSUnity[962] = 1.87707;
            LSUnity[962] = 2.7078;
            ESUnity[963] = 2.48748;
            SSUnity[963] = 1.91261;
            LSUnity[963] = 2.72343;
            ESUnity[964] = 2.50659;
            SSUnity[964] = 1.90356;
            LSUnity[964] = 2.69729;
            ESUnity[965] = 2.50226;
            SSUnity[965] = 1.88729;
            LSUnity[965] = 2.71403;
            ESUnity[966] = 2.49907;
            SSUnity[966] = 1.8989;
            LSUnity[966] = 2.73635;
            ESUnity[967] = 2.47668;
            SSUnity[967] = 1.88147;
            LSUnity[967] = 2.71698;
            ESUnity[968] = 2.49408;
            SSUnity[968] = 1.88846;
            LSUnity[968] = 2.7057;
            ESUnity[969] = 2.4747;
            SSUnity[969] = 1.90401;
            LSUnity[969] = 2.71981;
            ESUnity[970] = 2.49685;
            SSUnity[970] = 1.90417;
            LSUnity[970] = 2.70349;
            ESUnity[971] = 2.47902;
            SSUnity[971] = 1.88635;
            LSUnity[971] = 2.71896;
            ESUnity[972] = 2.48556;
            SSUnity[972] = 1.88752;
            LSUnity[972] = 2.7107;
            ESUnity[973] = 2.52355;
            SSUnity[973] = 1.90154;
            LSUnity[973] = 2.73126;
            ESUnity[974] = 2.51863;
            SSUnity[974] = 1.87506;
            LSUnity[974] = 2.72404;
            ESUnity[975] = 2.50303;
            SSUnity[975] = 1.8958;
            LSUnity[975] = 2.73664;
            ESUnity[976] = 2.49677;
            SSUnity[976] = 1.89085;
            LSUnity[976] = 2.72766;
            ESUnity[977] = 2.49893;
            SSUnity[977] = 1.88669;
            LSUnity[977] = 2.71329;
            ESUnity[978] = 2.495;
            SSUnity[978] = 1.90684;
            LSUnity[978] = 2.72124;
            ESUnity[979] = 2.48216;
            SSUnity[979] = 1.87381;
            LSUnity[979] = 2.70917;
            ESUnity[980] = 2.49511;
            SSUnity[980] = 1.88418;
            LSUnity[980] = 2.71533;
            ESUnity[981] = 2.49817;
            SSUnity[981] = 1.8935;
            LSUnity[981] = 2.7305;
            ESUnity[982] = 2.5022;
            SSUnity[982] = 1.89989;
            LSUnity[982] = 2.73146;
            ESUnity[983] = 2.50799;
            SSUnity[983] = 1.87953;
            LSUnity[983] = 2.74711;
            ESUnity[984] = 2.48675;
            SSUnity[984] = 1.89611;
            LSUnity[984] = 2.7491;
            ESUnity[985] = 2.51814;
            SSUnity[985] = 1.88679;
            LSUnity[985] = 2.74117;
            ESUnity[986] = 2.50713;
            SSUnity[986] = 1.89886;
            LSUnity[986] = 2.73801;
            ESUnity[987] = 2.52177;
            SSUnity[987] = 1.8935;
            LSUnity[987] = 2.73024;
            ESUnity[988] = 2.49259;
            SSUnity[988] = 1.89282;
            LSUnity[988] = 2.73127;
            ESUnity[989] = 2.51174;
            SSUnity[989] = 1.9088;
            LSUnity[989] = 2.73464;
            ESUnity[990] = 2.50533;
            SSUnity[990] = 1.91984;
            LSUnity[990] = 2.70061;
            ESUnity[991] = 2.50638;
            SSUnity[991] = 1.88528;
            LSUnity[991] = 2.73093;
            ESUnity[992] = 2.5212;
            SSUnity[992] = 1.90889;
            LSUnity[992] = 2.7378;
            ESUnity[993] = 2.52129;
            SSUnity[993] = 1.89211;
            LSUnity[993] = 2.74728;
            ESUnity[994] = 2.51768;
            SSUnity[994] = 1.88465;
            LSUnity[994] = 2.73877;
            ESUnity[995] = 2.51824;
            SSUnity[995] = 1.9064;
            LSUnity[995] = 2.72774;
            ESUnity[996] = 2.49155;
            SSUnity[996] = 1.88655;
            LSUnity[996] = 2.73405;
            ESUnity[997] = 2.50844;
            SSUnity[997] = 1.88969;
            LSUnity[997] = 2.72133;
            ESUnity[998] = 2.50625;
            SSUnity[998] = 1.89699;
            LSUnity[998] = 2.72068;
            ESUnity[999] = 2.48858;
            SSUnity[999] = 1.90052;
            LSUnity[999] = 2.73524;
            ESUnity[1000] = 2.50689;
            SSUnity[1000] = 1.87937;
            LSUnity[1000] = 2.71521;

            ESDLUnityCF[0] = 1;
            LSDLUnityCF[0] = 1;
            LS68UnityCF[0] = 1;
            SSDLUnityCF[0] = 1;
            ESDLUnityCF[1] = 1;
            LSDLUnityCF[1] = 1;
            LS68UnityCF[1] = 1;
            SSDLUnityCF[1] = 1;
            ESDLUnityCF[2] = 0.799351;
            LSDLUnityCF[2] = 0.799351;
            LS68UnityCF[2] = 0.799351;
            SSDLUnityCF[2] = 0.799351;
            ESDLUnityCF[3] = 0.944325;
            LSDLUnityCF[3] = 0.285521;
            LS68UnityCF[3] = 0.285521;
            SSDLUnityCF[3] = 0.521918;
            ESDLUnityCF[4] = 0.854403;
            LSDLUnityCF[4] = 0.340779;
            LS68UnityCF[4] = 0.340779;
            SSDLUnityCF[4] = 0.577339;
            ESDLUnityCF[5] = 0.785926;
            LSDLUnityCF[5] = 0.41601;
            LS68UnityCF[5] = 0.415995;
            SSDLUnityCF[5] = 0.528513;
            ESDLUnityCF[6] = 0.736956;
            LSDLUnityCF[6] = 0.447745;
            LS68UnityCF[6] = 0.447563;
            SSDLUnityCF[6] = 0.635396;
            ESDLUnityCF[7] = 0.722211;
            LSDLUnityCF[7] = 0.478807;
            LS68UnityCF[7] = 0.478511;
            SSDLUnityCF[7] = 0.659658;
            ESDLUnityCF[8] = 0.720667;
            LSDLUnityCF[8] = 0.489568;
            LS68UnityCF[8] = 0.489028;
            SSDLUnityCF[8] = 0.675483;
            ESDLUnityCF[9] = 0.725184;
            LSDLUnityCF[9] = 0.519816;
            LS68UnityCF[9] = 0.519352;
            SSDLUnityCF[9] = 0.696538;
            ESDLUnityCF[10] = 0.75123;
            LSDLUnityCF[10] = 0.530939;
            LS68UnityCF[10] = 0.530078;
            SSDLUnityCF[10] = 0.73329;
            ESDLUnityCF[11] = 0.755033;
            LSDLUnityCF[11] = 0.558881;
            LS68UnityCF[11] = 0.558642;
            SSDLUnityCF[11] = 0.727713;
            ESDLUnityCF[12] = 0.765896;
            LSDLUnityCF[12] = 0.545599;
            LS68UnityCF[12] = 0.566581;
            SSDLUnityCF[12] = 0.757271;
            ESDLUnityCF[13] = 0.780944;
            LSDLUnityCF[13] = 0.555364;
            LS68UnityCF[13] = 0.591561;
            SSDLUnityCF[13] = 0.768152;
            ESDLUnityCF[14] = 0.785441;
            LSDLUnityCF[14] = 0.548858;
            LS68UnityCF[14] = 0.592315;
            SSDLUnityCF[14] = 0.775964;
            ESDLUnityCF[15] = 0.787935;
            LSDLUnityCF[15] = 0.568839;
            LS68UnityCF[15] = 0.614014;
            SSDLUnityCF[15] = 0.783458;
            ESDLUnityCF[16] = 0.797322;
            LSDLUnityCF[16] = 0.573853;
            LS68UnityCF[16] = 0.621671;
            SSDLUnityCF[16] = 0.80077;
            ESDLUnityCF[17] = 0.802557;
            LSDLUnityCF[17] = 0.589571;
            LS68UnityCF[17] = 0.644098;
            SSDLUnityCF[17] = 0.797932;
            ESDLUnityCF[18] = 0.810293;
            LSDLUnityCF[18] = 0.59025;
            LS68UnityCF[18] = 0.640643;
            SSDLUnityCF[18] = 0.810314;
            ESDLUnityCF[19] = 0.812691;
            LSDLUnityCF[19] = 0.612567;
            LS68UnityCF[19] = 0.658897;
            SSDLUnityCF[19] = 0.81869;
            ESDLUnityCF[20] = 0.817311;
            LSDLUnityCF[20] = 0.615455;
            LS68UnityCF[20] = 0.659947;
            SSDLUnityCF[20] = 0.820111;
            ESDLUnityCF[21] = 0.820638;
            LSDLUnityCF[21] = 0.626054;
            LS68UnityCF[21] = 0.674892;
            SSDLUnityCF[21] = 0.829854;
            ESDLUnityCF[22] = 0.828034;
            LSDLUnityCF[22] = 0.62771;
            LS68UnityCF[22] = 0.670956;
            SSDLUnityCF[22] = 0.834026;
            ESDLUnityCF[23] = 0.828629;
            LSDLUnityCF[23] = 0.643118;
            LS68UnityCF[23] = 0.691888;
            SSDLUnityCF[23] = 0.841326;
            ESDLUnityCF[24] = 0.830826;
            LSDLUnityCF[24] = 0.646556;
            LS68UnityCF[24] = 0.692636;
            SSDLUnityCF[24] = 0.843973;
            ESDLUnityCF[25] = 0.838871;
            LSDLUnityCF[25] = 0.658871;
            LS68UnityCF[25] = 0.707093;
            SSDLUnityCF[25] = 0.847922;
            ESDLUnityCF[26] = 0.841844;
            LSDLUnityCF[26] = 0.66007;
            LS68UnityCF[26] = 0.700091;
            SSDLUnityCF[26] = 0.858153;
            ESDLUnityCF[27] = 0.839754;
            LSDLUnityCF[27] = 0.668774;
            LS68UnityCF[27] = 0.71671;
            SSDLUnityCF[27] = 0.856903;
            ESDLUnityCF[28] = 0.848066;
            LSDLUnityCF[28] = 0.671144;
            LS68UnityCF[28] = 0.717005;
            SSDLUnityCF[28] = 0.858529;
            ESDLUnityCF[29] = 0.848166;
            LSDLUnityCF[29] = 0.682491;
            LS68UnityCF[29] = 0.729161;
            SSDLUnityCF[29] = 0.868727;
            ESDLUnityCF[30] = 0.85345;
            LSDLUnityCF[30] = 0.679353;
            LS68UnityCF[30] = 0.720492;
            SSDLUnityCF[30] = 0.862837;
            ESDLUnityCF[31] = 0.854075;
            LSDLUnityCF[31] = 0.692286;
            LS68UnityCF[31] = 0.734685;
            SSDLUnityCF[31] = 0.871493;
            ESDLUnityCF[32] = 0.858788;
            LSDLUnityCF[32] = 0.693122;
            LS68UnityCF[32] = 0.739639;
            SSDLUnityCF[32] = 0.876271;
            ESDLUnityCF[33] = 0.860696;
            LSDLUnityCF[33] = 0.699773;
            LS68UnityCF[33] = 0.744988;
            SSDLUnityCF[33] = 0.875507;
            ESDLUnityCF[34] = 0.863621;
            LSDLUnityCF[34] = 0.699438;
            LS68UnityCF[34] = 0.742101;
            SSDLUnityCF[34] = 0.880854;
            ESDLUnityCF[35] = 0.862731;
            LSDLUnityCF[35] = 0.712464;
            LS68UnityCF[35] = 0.755613;
            SSDLUnityCF[35] = 0.878536;
            ESDLUnityCF[36] = 0.868193;
            LSDLUnityCF[36] = 0.711153;
            LS68UnityCF[36] = 0.754763;
            SSDLUnityCF[36] = 0.885307;
            ESDLUnityCF[37] = 0.863195;
            LSDLUnityCF[37] = 0.713465;
            LS68UnityCF[37] = 0.761168;
            SSDLUnityCF[37] = 0.887318;
            ESDLUnityCF[38] = 0.874256;
            LSDLUnityCF[38] = 0.717557;
            LS68UnityCF[38] = 0.76062;
            SSDLUnityCF[38] = 0.890606;
            ESDLUnityCF[39] = 0.872882;
            LSDLUnityCF[39] = 0.726389;
            LS68UnityCF[39] = 0.768433;
            SSDLUnityCF[39] = 0.890264;
            ESDLUnityCF[40] = 0.877161;
            LSDLUnityCF[40] = 0.726362;
            LS68UnityCF[40] = 0.767037;
            SSDLUnityCF[40] = 0.893034;
            ESDLUnityCF[41] = 0.873557;
            LSDLUnityCF[41] = 0.732999;
            LS68UnityCF[41] = 0.780292;
            SSDLUnityCF[41] = 0.898365;
            ESDLUnityCF[42] = 0.879086;
            LSDLUnityCF[42] = 0.726632;
            LS68UnityCF[42] = 0.772455;
            SSDLUnityCF[42] = 0.896602;
            ESDLUnityCF[43] = 0.880324;
            LSDLUnityCF[43] = 0.735104;
            LS68UnityCF[43] = 0.778726;
            SSDLUnityCF[43] = 0.899384;
            ESDLUnityCF[44] = 0.882034;
            LSDLUnityCF[44] = 0.736027;
            LS68UnityCF[44] = 0.780526;
            SSDLUnityCF[44] = 0.901048;
            ESDLUnityCF[45] = 0.881067;
            LSDLUnityCF[45] = 0.742498;
            LS68UnityCF[45] = 0.787753;
            SSDLUnityCF[45] = 0.904264;
            ESDLUnityCF[46] = 0.88604;
            LSDLUnityCF[46] = 0.744048;
            LS68UnityCF[46] = 0.787528;
            SSDLUnityCF[46] = 0.902966;
            ESDLUnityCF[47] = 0.883767;
            LSDLUnityCF[47] = 0.75067;
            LS68UnityCF[47] = 0.790481;
            SSDLUnityCF[47] = 0.906581;
            ESDLUnityCF[48] = 0.888839;
            LSDLUnityCF[48] = 0.745452;
            LS68UnityCF[48] = 0.791204;
            SSDLUnityCF[48] = 0.908435;
            ESDLUnityCF[49] = 0.887088;
            LSDLUnityCF[49] = 0.754027;
            LS68UnityCF[49] = 0.799592;
            SSDLUnityCF[49] = 0.907805;
            ESDLUnityCF[50] = 0.893795;
            LSDLUnityCF[50] = 0.75414;
            LS68UnityCF[50] = 0.790551;
            SSDLUnityCF[50] = 0.91133;
            ESDLUnityCF[51] = 0.890033;
            LSDLUnityCF[51] = 0.761648;
            LS68UnityCF[51] = 0.80264;
            SSDLUnityCF[51] = 0.913216;
            ESDLUnityCF[52] = 0.89339;
            LSDLUnityCF[52] = 0.756361;
            LS68UnityCF[52] = 0.805621;
            SSDLUnityCF[52] = 0.912063;
            ESDLUnityCF[53] = 0.892609;
            LSDLUnityCF[53] = 0.768596;
            LS68UnityCF[53] = 0.806967;
            SSDLUnityCF[53] = 0.913785;
            ESDLUnityCF[54] = 0.896503;
            LSDLUnityCF[54] = 0.762763;
            LS68UnityCF[54] = 0.803342;
            SSDLUnityCF[54] = 0.916749;
            ESDLUnityCF[55] = 0.896514;
            LSDLUnityCF[55] = 0.767644;
            LS68UnityCF[55] = 0.811288;
            SSDLUnityCF[55] = 0.915789;
            ESDLUnityCF[56] = 0.897336;
            LSDLUnityCF[56] = 0.771832;
            LS68UnityCF[56] = 0.810548;
            SSDLUnityCF[56] = 0.918895;
            ESDLUnityCF[57] = 0.897637;
            LSDLUnityCF[57] = 0.770402;
            LS68UnityCF[57] = 0.816211;
            SSDLUnityCF[57] = 0.919865;
            ESDLUnityCF[58] = 0.900154;
            LSDLUnityCF[58] = 0.77238;
            LS68UnityCF[58] = 0.810334;
            SSDLUnityCF[58] = 0.920497;
            ESDLUnityCF[59] = 0.899507;
            LSDLUnityCF[59] = 0.780689;
            LS68UnityCF[59] = 0.820914;
            SSDLUnityCF[59] = 0.921825;
            ESDLUnityCF[60] = 0.90176;
            LSDLUnityCF[60] = 0.773244;
            LS68UnityCF[60] = 0.818009;
            SSDLUnityCF[60] = 0.92345;
            ESDLUnityCF[61] = 0.903267;
            LSDLUnityCF[61] = 0.780089;
            LS68UnityCF[61] = 0.824183;
            SSDLUnityCF[61] = 0.92432;
            ESDLUnityCF[62] = 0.908935;
            LSDLUnityCF[62] = 0.781802;
            LS68UnityCF[62] = 0.818457;
            SSDLUnityCF[62] = 0.923952;
            ESDLUnityCF[63] = 0.901649;
            LSDLUnityCF[63] = 0.784195;
            LS68UnityCF[63] = 0.830039;
            SSDLUnityCF[63] = 0.925951;
            ESDLUnityCF[64] = 0.906029;
            LSDLUnityCF[64] = 0.781847;
            LS68UnityCF[64] = 0.824181;
            SSDLUnityCF[64] = 0.926259;
            ESDLUnityCF[65] = 0.905089;
            LSDLUnityCF[65] = 0.788101;
            LS68UnityCF[65] = 0.832788;
            SSDLUnityCF[65] = 0.927963;
            ESDLUnityCF[66] = 0.914035;
            LSDLUnityCF[66] = 0.783922;
            LS68UnityCF[66] = 0.823924;
            SSDLUnityCF[66] = 0.928288;
            ESDLUnityCF[67] = 0.906141;
            LSDLUnityCF[67] = 0.793554;
            LS68UnityCF[67] = 0.833908;
            SSDLUnityCF[67] = 0.930179;
            ESDLUnityCF[68] = 0.90988;
            LSDLUnityCF[68] = 0.794663;
            LS68UnityCF[68] = 0.831748;
            SSDLUnityCF[68] = 0.930933;
            ESDLUnityCF[69] = 0.909883;
            LSDLUnityCF[69] = 0.794824;
            LS68UnityCF[69] = 0.837215;
            SSDLUnityCF[69] = 0.925191;
            ESDLUnityCF[70] = 0.910747;
            LSDLUnityCF[70] = 0.797726;
            LS68UnityCF[70] = 0.83172;
            SSDLUnityCF[70] = 0.934905;
            ESDLUnityCF[71] = 0.91246;
            LSDLUnityCF[71] = 0.796223;
            LS68UnityCF[71] = 0.841664;
            SSDLUnityCF[71] = 0.932677;
            ESDLUnityCF[72] = 0.911705;
            LSDLUnityCF[72] = 0.797372;
            LS68UnityCF[72] = 0.837954;
            SSDLUnityCF[72] = 0.933801;
            ESDLUnityCF[73] = 0.9112;
            LSDLUnityCF[73] = 0.802599;
            LS68UnityCF[73] = 0.84414;
            SSDLUnityCF[73] = 0.934811;
            ESDLUnityCF[74] = 0.914824;
            LSDLUnityCF[74] = 0.799013;
            LS68UnityCF[74] = 0.836922;
            SSDLUnityCF[74] = 0.934255;
            ESDLUnityCF[75] = 0.915461;
            LSDLUnityCF[75] = 0.803278;
            LS68UnityCF[75] = 0.845512;
            SSDLUnityCF[75] = 0.936249;
            ESDLUnityCF[76] = 0.914989;
            LSDLUnityCF[76] = 0.802559;
            LS68UnityCF[76] = 0.84334;
            SSDLUnityCF[76] = 0.935739;
            ESDLUnityCF[77] = 0.914235;
            LSDLUnityCF[77] = 0.807167;
            LS68UnityCF[77] = 0.849;
            SSDLUnityCF[77] = 0.931335;
            ESDLUnityCF[78] = 0.91751;
            LSDLUnityCF[78] = 0.808534;
            LS68UnityCF[78] = 0.843154;
            SSDLUnityCF[78] = 0.93934;
            ESDLUnityCF[79] = 0.917076;
            LSDLUnityCF[79] = 0.808257;
            LS68UnityCF[79] = 0.851267;
            SSDLUnityCF[79] = 0.939084;
            ESDLUnityCF[80] = 0.917417;
            LSDLUnityCF[80] = 0.807184;
            LS68UnityCF[80] = 0.849053;
            SSDLUnityCF[80] = 0.939209;
            ESDLUnityCF[81] = 0.918592;
            LSDLUnityCF[81] = 0.812206;
            LS68UnityCF[81] = 0.853751;
            SSDLUnityCF[81] = 0.939624;
            ESDLUnityCF[82] = 0.920662;
            LSDLUnityCF[82] = 0.812232;
            LS68UnityCF[82] = 0.846056;
            SSDLUnityCF[82] = 0.94091;
            ESDLUnityCF[83] = 0.92628;
            LSDLUnityCF[83] = 0.814208;
            LS68UnityCF[83] = 0.855709;
            SSDLUnityCF[83] = 0.940528;
            ESDLUnityCF[84] = 0.917178;
            LSDLUnityCF[84] = 0.813478;
            LS68UnityCF[84] = 0.852672;
            SSDLUnityCF[84] = 0.940834;
            ESDLUnityCF[85] = 0.920872;
            LSDLUnityCF[85] = 0.817792;
            LS68UnityCF[85] = 0.859149;
            SSDLUnityCF[85] = 0.942504;
            ESDLUnityCF[86] = 0.923064;
            LSDLUnityCF[86] = 0.814158;
            LS68UnityCF[86] = 0.852491;
            SSDLUnityCF[86] = 0.942799;
            ESDLUnityCF[87] = 0.922202;
            LSDLUnityCF[87] = 0.819564;
            LS68UnityCF[87] = 0.85993;
            SSDLUnityCF[87] = 0.943895;
            ESDLUnityCF[88] = 0.927969;
            LSDLUnityCF[88] = 0.81924;
            LS68UnityCF[88] = 0.857938;
            SSDLUnityCF[88] = 0.943777;
            ESDLUnityCF[89] = 0.924312;
            LSDLUnityCF[89] = 0.82205;
            LS68UnityCF[89] = 0.863711;
            SSDLUnityCF[89] = 0.945166;
            ESDLUnityCF[90] = 0.923029;
            LSDLUnityCF[90] = 0.818275;
            LS68UnityCF[90] = 0.854912;
            SSDLUnityCF[90] = 0.944398;
            ESDLUnityCF[91] = 0.929101;
            LSDLUnityCF[91] = 0.82352;
            LS68UnityCF[91] = 0.865029;
            SSDLUnityCF[91] = 0.945379;
            ESDLUnityCF[92] = 0.927667;
            LSDLUnityCF[92] = 0.823616;
            LS68UnityCF[92] = 0.861066;
            SSDLUnityCF[92] = 0.946653;
            ESDLUnityCF[93] = 0.927667;
            LSDLUnityCF[93] = 0.825935;
            LS68UnityCF[93] = 0.866604;
            SSDLUnityCF[93] = 0.94639;
            ESDLUnityCF[94] = 0.928384;
            LSDLUnityCF[94] = 0.823225;
            LS68UnityCF[94] = 0.861614;
            SSDLUnityCF[94] = 0.946253;
            ESDLUnityCF[95] = 0.924162;
            LSDLUnityCF[95] = 0.828422;
            LS68UnityCF[95] = 0.868421;
            SSDLUnityCF[95] = 0.947046;
            ESDLUnityCF[96] = 0.927004;
            LSDLUnityCF[96] = 0.825587;
            LS68UnityCF[96] = 0.865454;
            SSDLUnityCF[96] = 0.948222;
            ESDLUnityCF[97] = 0.931712;
            LSDLUnityCF[97] = 0.831217;
            LS68UnityCF[97] = 0.869022;
            SSDLUnityCF[97] = 0.944239;
            ESDLUnityCF[98] = 0.932939;
            LSDLUnityCF[98] = 0.826585;
            LS68UnityCF[98] = 0.863025;
            SSDLUnityCF[98] = 0.949949;
            ESDLUnityCF[99] = 0.931967;
            LSDLUnityCF[99] = 0.83475;
            LS68UnityCF[99] = 0.870568;
            SSDLUnityCF[99] = 0.948752;
            ESDLUnityCF[100] = 0.928198;
            LSDLUnityCF[100] = 0.832441;
            LS68UnityCF[100] = 0.868542;
            SSDLUnityCF[100] = 0.950152;
            unityTablesLoaded = true;
        }
    }
    void RCR::setRejectionTech(RejectionTechs rejectionTech)
    {
        this->rejectionTech = rejectionTech;
        alignTechniques();
    }
    void RCR::setParametricModel(FunctionalForm &parametricModel)
    {
        this->parametricModel = &parametricModel;
        this->muType = PARAMETRIC;
    }
    void RCR::setNonParametricModel(NonParametric &nonParametricModel)
    {
        this->nonParametricModel = &nonParametricModel;
        this->muType = NONPARAMETRIC;
    }
    void RCR::setMuType(MuTypes muType)
    {
        this->muType = muType;
    }
    void RCR::performRejection(std::vector<double> &y)
    {
        this->result.flags.clear();
        this->result.indices.clear();
        this->result.flags.resize(y.size(), true);
        handleRCRLoopSelect(false, y);
        this->muTech = MEDIAN;
        this->sigmaTech = SIXTY_EIGHTH_PERCENTILE;
        handleRCRLoopSelect(false, y);
        this->muTech = MEAN;
        this->sigmaTech = STANDARD_DEVIATION;
        handleRCRLoopSelect(false, y);
        setTrueVec(result.flags, result.indices, y, result.cleanY);
        alignTechniques();
    }
    void RCR::performRejection(std::vector<double> &w, std::vector<double> &y)
    {
        this->result.flags.clear();
        this->result.indices.clear();
        this->result.flags.resize(y.size(), true);
        handleRCRLoopSelect(false, w, y);
        this->muTech = MEDIAN;
        this->sigmaTech = SIXTY_EIGHTH_PERCENTILE;
        handleRCRLoopSelect(false, w, y);
        this->muTech = MEAN;
        this->sigmaTech = STANDARD_DEVIATION;
        handleRCRLoopSelect(false, w, y);
        setTrueVec(result.flags, result.indices, w, y, result.cleanW, result.cleanY);
        alignTechniques();
    }
    void RCR::performBulkRejection(std::vector<double> &y)
    {
        this->result.flags.clear();
        this->result.flags.resize(y.size(), true);
        handleRCRLoopSelect(true, y);
        handleRCRLoopSelect(false, y);
        this->muTech = MEDIAN;
        this->sigmaTech = SIXTY_EIGHTH_PERCENTILE;
        handleRCRLoopSelect(false, y);
        this->muTech = MEAN;
        this->sigmaTech = STANDARD_DEVIATION;
        handleRCRLoopSelect(false, y);
        setTrueVec(result.flags, result.indices, y, result.cleanY);
        alignTechniques();
        setFinalVectors(y);
    }
    void RCR::performBulkRejection(std::vector<double> &w, std::vector<double> &y)
    {
        this->result.flags.clear();
        this->result.flags.resize(y.size(), true);
        handleRCRLoopSelect(true, w, y);
        handleRCRLoopSelect(false, w, y);
        this->muTech = MEDIAN;
        this->sigmaTech = SIXTY_EIGHTH_PERCENTILE;
        handleRCRLoopSelect(false, w, y);
        this->muTech = MEAN;
        this->sigmaTech = STANDARD_DEVIATION;
        handleRCRLoopSelect(false, w, y);
        setTrueVec(result.flags, result.indices, w, y, result.cleanW, result.cleanY);
        alignTechniques();
        setFinalVectors(w, y);
    }
    void RCR::setFinalVectors(std::vector<double> &w, std::vector<double> &y)
    {
        std::vector<bool> flagsHold = this->result.flags;
        std::vector<double> rejectedY, rejectedW;
        this->result.originalW = w;
        this->result.originalY = y;
        for (size_t i = 0; i < w.size(); i++)
        {
            if (!flagsHold[i])
            {
                rejectedW.push_back(w[i]);
                rejectedY.push_back(y[i]);
            }
        }
        this->result.rejectedW = rejectedW;
        this->result.rejectedY = rejectedY;

        std::vector<double> cleanYHold = this->result.cleanY, cleanYAboveMeanHold, cleanYBelowMeanHold;
        std::vector<double> cleanWHold = this->result.cleanW, cleanWAboveMeanHold, cleanWBelowMeanHold;
        double meanHold = getMean((int) cleanYHold.size(), cleanWHold, cleanYHold);
        for (size_t i = 0; i < cleanYHold.size(); i++)
        {
            cleanYHold[i] -= meanHold;
            if (isEqual(cleanYHold[i], 0.0))
            {
                cleanYAboveMeanHold.push_back(cleanYHold[i]);
                cleanYBelowMeanHold.push_back(cleanYHold[i]);
                cleanWBelowMeanHold.push_back(0.5*cleanWHold[i]);
                cleanWAboveMeanHold.push_back(0.5*cleanWHold[i]);
            }
            else if (cleanYHold[i] > 0.0)
            {
                cleanYAboveMeanHold.push_back(cleanYHold[i]);
                cleanWAboveMeanHold.push_back(cleanWHold[i]);
            }
            else
            {
                cleanYBelowMeanHold.push_back(cleanYHold[i]);
                cleanWBelowMeanHold.push_back(cleanWHold[i]);
            }
        }
        //double stDevTotalHold = getStDev(1.0, cleanWHold, cleanYHold);
        //double stDevBelowHold = getStDev(.5, cleanWBelowMeanHold, cleanYBelowMeanHold);
        //double stDevAboveHold = getStDev(.5, cleanWAboveMeanHold, cleanYAboveMeanHold);

        double stDevTotalHold = getStDev(delta, cleanWHold, cleanYHold);
        double stDevBelowHold = getStDev(delta/2.0, cleanWBelowMeanHold, cleanYBelowMeanHold);
        double stDevAboveHold = getStDev(delta/2.0, cleanWAboveMeanHold, cleanYAboveMeanHold);

        this->result.stDevAbove = stDevAboveHold;
        this->result.stDevBelow = stDevBelowHold;
        this->result.stDevTotal = stDevTotalHold;
    }
    void RCR::setFinalVectors(std::vector<double> &y)
    {
        std::vector<bool> flagsHold = this->result.flags;
        std::vector<double> rejectedY;
        this->result.originalY = y;
        for (size_t i = 0; i < y.size(); i++)
        {
            if (!flagsHold[i])
            {
                rejectedY.push_back(y[i]);
            }
        }
        this->result.rejectedY = rejectedY;
        std::vector<double> cleanYHold = this->result.cleanY, cleanYAboveMeanHold, cleanYBelowMeanHold;
        std::vector<double> cleanWAboveMeanHold, cleanWBelowMeanHold;
        double meanHold = getMean((int) cleanYHold.size(), cleanYHold);
        for (size_t i = 0; i < cleanYHold.size(); i++)
        {
            cleanYHold[i] -= meanHold;
            if (isEqual(cleanYHold[i], 0.0))
            {
                cleanYAboveMeanHold.push_back(cleanYHold[i]);
                cleanYBelowMeanHold.push_back(cleanYHold[i]);
                cleanWBelowMeanHold.push_back(0.5);
                cleanWAboveMeanHold.push_back(0.5);
            }
            else if (cleanYHold[i] > 0.0)
            {
                cleanYAboveMeanHold.push_back(cleanYHold[i]);
                cleanWAboveMeanHold.push_back(1.0);
            }
            else
            {
                cleanYBelowMeanHold.push_back(cleanYHold[i]);
                cleanWBelowMeanHold.push_back(1.0);
            }
        }
        //double stDevTotalHold = getStDev(1.0, cleanYHold);
        //double stDevBelowHold = getStDev(.5, cleanWBelowMeanHold, cleanYBelowMeanHold);
        //double stDevAboveHold = getStDev(.5, cleanWAboveMeanHold, cleanYAboveMeanHold);

        double stDevTotalHold = getStDev(delta, cleanYHold);
        double stDevBelowHold = getStDev(delta/2.0, cleanWBelowMeanHold, cleanYBelowMeanHold);
        double stDevAboveHold = getStDev(delta/2.0, cleanWAboveMeanHold, cleanYAboveMeanHold);

        this->result.stDevAbove = stDevAboveHold;
        this->result.stDevBelow = stDevBelowHold;
        this->result.stDevTotal = stDevTotalHold;

    }



    //FN & CF Models:
    double RCR::getLower68CF(int n, std::vector<double> &w)
    {
        double x = getCFRatio(w), y1, cf, a1, b1, logx, logn = log10(n);
        logx = log10(x);

        if (n < 101)
        {
            y1 = 1.0 / RCR::LS68UnityCF[n];
        }
        else
        {
            y1 = 1.0 / (1.0 - 2.3525*pow(n, -.627));
        }
        if (x == 0)
        {
            return y1;
        }
        if (n == 2)
        {
            b1 = 1.08149771508934;
            a1 = -0.398375456223868;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else if (n == 3)
        {
            b1 = 1.06994341034958;
            a1 = -1.14618991625901;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else if (n == 4)
        {
            cf = y1;
        }
        else if (n == 5)
        {
            b1 = 0.45972894692595;
            a1 = -1.11957357644441;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else if (5 < n && n < 101)
        {
            b1 = -1.4528*pow(logn, 3) + 5.3519*pow(logn, 2) - 5.33*logn + 2.2902 + pow(-1, n) * 0.1879*pow(logn, 0.9521);
            a1 = -1.1937*pow(logn, 4) + 6.5268*pow(logn, 3) - 13.308*pow(logn, 2) + 11.432*logn - 4.4769;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else
        {
            b1 = 1.4154 + pow(-1, n) * 0.363528;
            a1 = -0.5408 * logn - 0.6482;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        return cf;
    }
    double RCR::getLower68CF(int n)
    {
        if (n < 101)
        {
            return 1.0 / RCR::LS68UnityCF[n];
        }
        else
        {
            return 1.0 / (1.0 - 2.3525*pow(n, -.627));
        }
    }
    double RCR::getLowerDLCF(int n, std::vector<double> &w)
    {
        double x = getCFRatio(w), y1, cf, a1, b1, logx, logn = log10(n);
        logx = log10(x);

        if (n < 101)
        {
            y1 = 1.0 / RCR::LSDLUnityCF[n];
        }
        else
        {
            y1 = 1.0 / (1.0 - 3.3245*pow(n, -.65));
        }
        if (x == 0)
        {
            return y1;
        }
        if (n == 2)
        {
            b1 = 1.08149771508934;
            a1 = -0.398375456223868;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else if (n == 3)
        {
            b1 = 1.51433669748424;
            a1 = -1.10939332689999;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else if (n == 4)
        {
            cf = y1;
        }
        else if (n == 5)
        {
            b1 = 0.339404852185332;
            a1 = -1.1445790996528;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else if (5 < n && n < 21)
        {
            b1 = 43.179*pow(logn, 6) - 331.85*pow(logn, 5) + 968.25*pow(logn, 4) - 1399.1*pow(logn, 3) + 1070.7*pow(logn, 2) - 415.81*logn + 65.002 + pow(-1, n) * 0.1365 * pow(logn, 2.4716);
            a1 = -0.2683*pow(logn, 4) + 1.9174*pow(logn, 3) - 5.062*pow(logn, 2) + 5.452*logn - 2.9999;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else if (20 < n && n < 101)
        {
            b1 = 1.5144*logn - 0.0448 + pow(-1, n) * 0.1365 * pow(logn, 2.4716);
            a1 = -0.2683*pow(logn, 4) + 1.9174*pow(logn, 3) - 5.062*pow(logn, 2) + 5.452*logn - 2.9999;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else
        {
            b1 = 2.988077 + pow(-1, n) * 0.753032;
            a1 = -0.4282*logn - 0.4412;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        return cf;
    }
    double RCR::getLowerDLCF(int n)
    {
        if (n < 101)
        {
            return 1.0 / RCR::LSDLUnityCF[n];
        }
        else
        {
            return 1.0 / (1.0 - 3.3245*pow(n, -.65));
        }
    }
    double RCR::getSingleDLCF(int n, std::vector<double> &w)
    {
        double x = getCFRatio(w), y1, cf, a1, b1, logx, logn = log10(n);
        logx = log10(x);

        if (n < 101)
        {
            y1 = 1.0 / RCR::SSDLUnityCF[n];
        }
        else
        {
            y1 = 1.0 / (1.0 - 3.578*pow(n, -.942));
        }
        if (x == 0)
        {
            return y1;
        }
        if (n == 2)
        {
            b1 = 0.273907084639124;
            a1 = -3.15279135630884;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else if (n == 3)
        {
            b1 = 0.448654915529039;
            a1 = -1.19134294551807;

            cf = y1 / pow(10, pow(10, (a1 + b1*logx)));
        }
        else if (n == 4)
        {
            b1 = 3.38253309393705;
            a1 = -1.05087405984868;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else if (n == 5)
        {
            b1 = 0.118507989164207;
            a1 = -1.41453721585464;

            cf = y1 / pow(10, pow(10, (a1 + b1*logx)));
        }
        else
        {
            b1 = 0.1196*logn + 4.5073;
            a1 = -0.7914*logn + 0.0243;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }

        return cf;
    }
    double RCR::getSingleDLCF(int n)
    {
        if (n < 101)
        {
            return 1.0 / RCR::SSDLUnityCF[n];
        }
        else
        {
            return 1.0 / (1.0 - 3.578*pow(n, -.942));
        }
    }
    double RCR::getEachDLCF(int n, std::vector<double> &w)
    {
        double x = getCFRatio(w), y1, cf, a1, b1, logx, logn = log10(n);
        logx = log10(x);

        if (n < 101)
        {
            y1 = 1.0 / RCR::ESDLUnityCF[n];
        }
        else
        {
            y1 = 1.0 / (1.0 - 3.1666*pow(n, -.833));
        }
        if (x == 0)
        {
            return y1;
        }
        if (n == 2)
        {
            b1 = 0.733632602759432;
            a1 = -2.59506757852134;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else if (n == 3)
        {
            b1 = 0.816017988131836;
            a1 = -0.854637214866955;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else if (n == 4)
        {
            b1 = 1.16048439814909;
            a1 = -0.954253810365265;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else if (4 < n && n < 20)
        {
            b1 = 4.0458*pow(logn, 2) - 6.4354*logn + 2.7667;
            a1 = -1.3993*pow(logn, 3) + 6.5746*pow(logn, 2) - 9.8844*logn + 2.8572;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else if (19 < n && n < 101)
        {
            b1 = 1.7394*logn - 1.0435;
            a1 = -1.3993*pow(logn, 3) + 6.5746*pow(logn, 2) - 9.8844*logn + 2.8572;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }
        else
        {
            b1 = 1.4123*logn - 0.3893;
            a1 = -0.5989*logn - 0.6097;

            cf = y1 * pow(10, pow(10, (a1 + b1*logx)));
        }

        return cf;
    }
    double RCR::getEachDLCF(int n)
    {
        if (n < 101)
        {
            return 1.0 / RCR::ESDLUnityCF[n];
        }
        else
        {
            return 1.0 / (1.0 - 3.1666*pow(n, -.833));
        }
    }
    double RCR::getSingleFN(int n, std::vector<double> &x, std::vector<double> &w)
    {
        double ratio = getFNRatio(x, w), f_n, a1, y1, b1, logn = log10(n), logx;
        logx = log10(ratio);

        if (n < 1001)
        {
            y1 = RCR::SSUnity[n];
        }
        else
        {
            y1 = getSingleFN(n, x);
        }
        if (ratio == 0)
        {
            return y1;
        }
        if (3 < n && n < 8)
        {
            a1 = RCR::SSConstants[0][n];
            b1 = RCR::SSConstants[1][n];
            f_n = y1 * pow(10, pow(10, a1 + b1 * logx));
        }
        else if (7 < n && n < 1001)
        {

            b1 = -0.3556*pow(logn, 6) + 3.7036*pow(logn, 5) - 14.932*pow(logn, 4) + 29.176*pow(logn, 3) - 28.81*pow(logn, 2) + 14.397*logn - 2.6451;
            a1 = 0.2313*pow(logn, 6) - 3.02*pow(logn, 5) + 15.997*pow(logn, 4) - 43.713*pow(logn, 3) + 64.629*pow(logn, 2) - 49.976*logn + 15.484 + pow(-1, n)*0.1513*pow(n, -0.471);
            f_n = y1 * pow(10, pow(10, a1 + b1 * logx));

        }
        else if (n > 1000)
        {
            f_n = y1;
        }
        else
        {
            f_n = -999999;
        }
        return f_n;
    }
    double RCR::getSingleFN(int n, std::vector<double> &x)
    {
        if (n < 1001)
        {
            return RCR::SSUnity[n];
        }
        else
        {
            return 39.2519*pow(n, -.7969)  + 1.8688;
        }
    }
    double RCR::getLowerFN(int n, std::vector<double> &x, std::vector<double> &w)
    {
        double ratio = getFNRatio(x, w), f_n, a1, y1, b1, a2, b2, logn = log10(n), logx;
        logx = log10(ratio);

        if (n < 1001)
        {
            y1 = RCR::LSUnity[n];
        }
        else
        {
            y1 = getLowerFN(n, x);
        }
        if (ratio == 0)
        {
            return y1;
        }
        if (n == 5)
        {
            y1 = 36.8534;
            b1 = -0.300348560626506;
            a1 = -0.0828207791627729;
            f_n = y1 * pow(10, pow(10, a1 + b1*logx));
        }
        else if (n == 6)
        {
            b1 = -0.244262599183892;
            a1 = -0.267502535686502;
            f_n = y1 * pow(10, pow(10, a1 + b1*logx));
        }
        else if (n == 7)
        {
            b1 = -0.409677351330214;
            a1 = -0.558845927943435;
            f_n = y1 * pow(10, pow(10, a1 + b1*logx));
        }
        else if (n == 8)
        {
            b1 = -0.488354948027081;
            a1 = -0.889342857411619;
            f_n = y1 * pow(10, pow(10, a1 + b1*logx));
        }
        else if (8 < n && n < 1001)
        {
            b1 = 0.1462*pow(logn, 3.0) - 4.2139*pow(logn, 2.0) + 14.366*logn - 10.658;
            a1 = -0.541*pow(logn, 5.0) + 4.6943*pow(logn, 4.0) - 15.407*pow(logn, 3.0) + 21.875*pow(logn, 2.0) - 11.211*logn - 0.3798;
            b2 = 26.945*pow(logn, 3.0) - 221.42*pow(logn, 2.0) + 606.91*logn - 553.89;
            a2 = 18.149*pow(logn, 3.0) - 149.27*pow(logn, 2.0) + 410.15*logn - 378.47;

            if (n <= 264 || (n < 568 && (a1 + b1*logx > a2 + b2*logx)))
            {
                f_n = y1 * pow(10, pow(10, a1 + b1* logx));
            }
            else
            {
                f_n = y1 / pow(10, pow(10, a2 + b2 * logx));
            }
        }
        else
        {
            b1 = 0.0424*logn + 1.4479;
            a1 = 0.3861*logn - 2.5852;
            f_n = y1 / pow(10, pow(10, a1 + b1 * logx));
        }

        return f_n;
    }
    double RCR::getLowerFN(int n, std::vector<double> &x)
    {
        if (n < 1001)
        {
            return RCR::LSUnity[n];
        }
        else
        {
            return pow(1.3399, pow(n, .1765));
        }
    }
    double RCR::getEachFN(int n, std::vector<double> &x, std::vector<double> &w)
    {
        double ratio = getFNRatio(x, w), f_n, a1, y1, b1, a2, b2, logn = log10(n), logx;
        logx = log10(ratio);

        if (n < 1001)
        {
            y1 = RCR::ESUnity[n];
        }
        else
        {
            y1 = getEachFN(n, x);
        }
        if (ratio == 0)
        {
            return y1;
        }
        if (n == 5)
        {
            b1 = 2.13417275654528;
            a1 = -0.466431459550531;
            f_n = y1 * (pow(10.0, pow(10.0, a1 + b1*logx)));
        }
        else if (n == 6)
        {
            b1 = 1.0196951775215;
            a1 = -0.312373723591738;
            f_n = y1 * (pow(10.0, pow(10.0, a1 + b1*logx)));
        }
        else if (n == 7)
        {
            b1 = 0.579724747519776;
            a1 = -0.750190463040497;
            f_n = y1 * (pow(10.0, pow(10.0, a1 + b1*logx)));
        }
        else if (7 < n && n < 1001)
        {
            b1 = 5.8718*pow(logn, 4.0) - 47.049*pow(logn, 3.0) + 131.12*pow(logn, 2.0) - 150.24*logn + 61.727;
            a1 = 3.1767*pow(logn, 6.0) - 34.561*pow(logn, 5.0) + 152.16*pow(logn, 4.0) - 347.96*pow(logn, 3.0) + 435.59*pow(logn, 2.0) - 282.57*logn + 73.696;
            b2 = -2.7584*pow(logn, 2.0) + 17.078*logn - 24.602;
            a2 = -1.8953*pow(logn, 2.0) + 11.745*logn - 19.36;

            if (n < 191 || (n < 306 && a1 + b1*logx > a2 + b2*logx))
            {
                f_n = y1 * (pow(10.0, pow(10.0, a1 + b1*logx)));
            }
            else
            {
                f_n = y1 / (pow(10.0, pow(10.0, a2 + b2*logx)));
            }
        }
        else
        {
            b1 = 1.8064;
            a1 = -1.1827;
            f_n = y1 / (pow(10.0, pow(10.0, a1 + b1*logx)));
        }

        return f_n;
    }
    double RCR::getEachFN(int n, std::vector<double> &x)
    {
        if (n < 1001)
        {
            return RCR::ESUnity[n];
        }
        else
        {

            return pow(1.2591, pow(n, .2052));
        }
    }


    //statistic tools:
    double RCR::nCorrect(int n, std::vector<double> &w)
    {
        switch (this->rejectionTech)
        {
        case SS_MEDIAN_DL:
            return getSingleDLCF(n, w);
        case LS_MODE_68:
            return getLower68CF(n, w);
        case LS_MODE_DL:
            return getLowerDLCF(n, w);
        case ES_MODE_DL:
            return getEachDLCF(n, w);
        }
        return -1;
    }
    double RCR::nCorrect(int n)
    {
        switch (this->rejectionTech)
        {
        case SS_MEDIAN_DL:
            return getSingleDLCF(n);
        case LS_MODE_68:
            return getLower68CF(n);
        case LS_MODE_DL:
            return getLowerDLCF(n);
        case ES_MODE_DL:
            return getEachDLCF(n);
        }
        return -1;
    }
    double RCR::getFN(int n, std::vector<double> &x, std::vector<double> &w)
    {
        switch (this->sigmaChoice)
        {
        case SINGLE:
            return getSingleFN(n, x, w);
        case LOWER:
            return getLowerFN(n, x, w);
        case EACH:
            return getEachFN(n, x, w);
        }
        return -1;
    }
    double RCR::getFN(int n, std::vector<double> &x)
    {
        switch (this->sigmaChoice)
        {
        case SINGLE:
            return getSingleFN(n, x);
        case LOWER:
            return getLowerFN(n, x);
        case EACH:
            return getEachFN(n, x);
        }
        return -1;
    }


    //rejectors:
    bool RCR::reject(int paramCount, int trueCount, int &index, double max, std::vector<bool> &flags, std::vector<double> &y)
    {
        if (distinctValuesCheck(paramCount, flags, y) && trueCount*erfcCustom(max) < .5)
        {
            flags[index] = false;
            //	rejectedBy[index] = muTech;
            return false;
        }
        return true;
    }
    bool RCR::reject(int trueCount, int &index, double max, std::vector<bool> &flags, std::vector<double> &y)
    {
        if (distinctValuesCheck(flags, y) && trueCount*erfcCustom(max) < .5)
        {
            flags[index] = false;
            //	rejectedBy[index] = muTech;
            return false;
        }
        return true;
    }
    bool RCR::bulkReject(int paramCount, std::vector<bool> &flags, std::vector<int> &index, std::vector<double> &max, std::vector<double> &y)
    {
        bool noPointsRejected = true; //, threeDistinct = true;
        int size = (int) max.size();
        int i = size - 1;
        while (distinctValuesCheck(paramCount, flags, y) && size*erfcCustom(max[i]) < .5)
        {
            flags[index[i]] = false;
            noPointsRejected = false;
            //	rejectedBy[index[i]] = 5;
            i--;
        }
        return noPointsRejected;
    }
    bool RCR::bulkReject(std::vector<bool> &flags, std::vector<int> &index, std::vector<double> &max, std::vector<double> &y)
    {
        bool noPointsRejected = true; //, threeDistinct = true;
        int size = (int) max.size();
        int i = size - 1;
        while (distinctValuesCheck(flags, y) && size*erfcCustom(max[i]) < .5)
        {
            flags[index[i]] = false;
            noPointsRejected = false;
            //	rejectedBy[index[i]] = 5;
            i--;
        }
        return noPointsRejected;
    }


    //sigma calculations
    double RCR::fitDL(double counter, std::vector<double> &w, std::vector<double> &x, std::vector<double> &y)
    {
        int amountXUnderOne = 0;
        double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, doubleLineError = 0, singleLineError = 0, singleLineFit, xAtM, xAtIxAtMDiff, wxProd, sigma, tau, factor, deltaChiSquared, errorAtI;
        amountXUnderOne = countAmountLessThanOne(x);
        singleLineFit = getOriginFixedRegressionLine(0, amountXUnderOne, w, x, y);
        if (x.size() < 4)
        {
            return singleLineFit;
        }
        int m = (int)(mFinder(1, amountXUnderOne - 1, amountXUnderOne - 1, (int)(max((double)(y.size()) / 6.36, 1.0)), w, x, y));
        xAtM = x[m];
        for (int i = 0; i <= m; i++)
        {
            wxProd = w[i] * x[i];
            a += wxProd * x[i];
            e += wxProd * y[i];
        }
        for (int i = m + 1; i < amountXUnderOne; i++)
        {
            xAtIxAtMDiff = (x[i] - xAtM);
            a += xAtM * xAtM * w[i];
            b += xAtM * w[i] * xAtIxAtMDiff;
            d += w[i] * xAtIxAtMDiff * xAtIxAtMDiff;
            e += xAtM*w[i] * y[i];
            f += w[i] * y[i] * xAtIxAtMDiff;
        }

        c = b;

        tau = (f - e*c / a) / (d - c*b / a);
        sigma = (e - tau*b) / a;
        for (int i = 0; i <= m; i++)
        {
            factor = (sigma*x[i] - y[i]);
            doubleLineError += w[i] * factor*factor;
        }
        for (int i = m + 1; i < amountXUnderOne; i++)
        {
            factor = sigma*xAtM + tau*(x[i] - xAtM) - y[i];
            doubleLineError += w[i] * factor*factor;
        }

        for (int i = 0; i < amountXUnderOne; i++)
        {
            errorAtI = (singleLineFit * x[i] - y[i]);
            singleLineError += w[i] * errorAtI * errorAtI;
        }
        deltaChiSquared = (singleLineError - doubleLineError) / doubleLineError;

        if (sigma < 0)
        {
            sigma = .0000000001;
        }
        if (deltaChiSquared < getFN((int)counter, x, w))
        {
            sigma = singleLineFit;
        }
        return sigma;
    }
    double RCR::fitDL(double counter, std::vector<double> &x, std::vector<double> &y)
    {
        int amountXUnderOne = 0;
        double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, doubleLineError = 0, singleLineError = 0, singleLineFit, xAtM, xAtIxAtMDiff, sigma, tau, factor, deltaChiSquared, errorAtI;
        amountXUnderOne = countAmountLessThanOne(x);
        singleLineFit = getOriginFixedRegressionLine(0, amountXUnderOne, x, y);
        if (x.size() < 4)
        {
            return singleLineFit;
        }
        int m = (int)(mFinder(1, amountXUnderOne - 1, amountXUnderOne - 1, (int)(max((double)(y.size()) / 6.36, 1.0)), x, y));
        xAtM = x[m];
        for (int i = 0; i <= m; i++)
        {
            a += x[i] * x[i];
            e += x[i] * y[i];
        }
        for (int i = m + 1; i < amountXUnderOne; i++)
        {
            xAtIxAtMDiff = (x[i] - xAtM);
            a += xAtM * xAtM;
            b += xAtM * xAtIxAtMDiff;
            d += xAtIxAtMDiff * xAtIxAtMDiff;
            e += xAtM * y[i];
            f += y[i] * xAtIxAtMDiff;
        }

        c = b;

        tau = (f - e*c / a) / (d - c*b / a);
        sigma = (e - tau*b) / a;
        for (int i = 0; i <= m; i++)
        {
            factor = (sigma*x[i] - y[i]);
            doubleLineError += factor*factor;
        }
        for (int i = m + 1; i < amountXUnderOne; i++)
        {
            factor = sigma*xAtM + tau*(x[i] - xAtM) - y[i];
            doubleLineError += factor*factor;
        }

        for (int i = 0; i < amountXUnderOne; i++)
        {
            errorAtI = (singleLineFit * x[i] - y[i]);
            singleLineError += errorAtI * errorAtI;
        }
        deltaChiSquared = (singleLineError - doubleLineError) / doubleLineError;

        if (sigma < 0)
        {
            sigma = .0000000001;
        }
        if (deltaChiSquared < getFN((int)(counter), x))
        {
            sigma = singleLineFit;
        }
        return sigma;
    }


    //technique handlers:
    double RCR::handleMuTechSelect(std::vector<double> &w, std::vector<double> &y)
    {
        switch (muTech)
        {
        case MEAN:
            return getMean((int) w.size(), w, y);
        case MEDIAN:
            sort(w, y);
            return getMedian((int) w.size(), w, y);
        case MODE:
            sort(w, y);
            return 	getMode((int) w.size(), w, y);
        }
        return -1;
    }
    double RCR::handleMuTechSelect(int trueCount, std::vector<double> &y)
    {
        switch (muTech)
        {
        case MEAN:
            return getMean(trueCount, y);
        case MEDIAN:
            sort(y);
            return getMedian(y);
        case MODE:
            sort(y);
            return getMode(trueCount, y);
        }
        return -1;
    }
    std::vector<double> RCR::handleMuTechSelect()
    {
        std::vector<double> line, errors;
        switch (muTech)
        {
        case MEAN:
            line = parametricModel->regression();
            break;
        case MEDIAN:
            if (parametricModel->parameterSpace.size() == 2)
            {
                line = get2DMedian(*parametricModel);
            }
            else if (parametricModel->parameterSpace.size() == 3)
            {
                line = get3DMedian(*parametricModel);
            }
            else
            {
                line = getNDMedian(*parametricModel);
            }
            break;
        case MODE:
            if (parametricModel->parameterSpace.size() == 2)
            {
                line = get2DMode(*parametricModel);
            }
            else if (parametricModel->parameterSpace.size() == 3)
            {
                line = get3DMode(*parametricModel);
            }
            else  // ND case
            {
                line = getNDMode(*parametricModel);
            }
            break;
        }
        //parametricModel->setModel(line);
        parametricModel->parameters = line; //most recent calculation of parameters ("line" vector) is stored)

        if (parametricModel->NDcheck == true) {
            errors = parametricModel->getErrors_ND(line);
        }
        else if (parametricModel->NDcheck == false) {
            errors = parametricModel->getErrors(line);
        }
        return errors;
    }

    double RCR::handleSigmaTechSelect(int counter, std::vector<double> &w, std::vector<double> &y)
    {
        std::vector<double> x = getXVec((int) y.size(), w);
        int xBelow1Count = countAmountLessThanOne(x);
//        double nCorrection = nCorrect(counter, w);
        switch (sigmaTech)
        {
        case STANDARD_DEVIATION:
            //return (sigmaChoice == SINGLE) ? getStDev(1.0, w, y) : getStDev(0.5, w, y);
            return (sigmaChoice == SINGLE) ? getStDev(delta, w, y) : getStDev(delta/2.0, w, y);
        case SIXTY_EIGHTH_PERCENTILE:
            return get68th(w, y);
        case SINGLE_LINE:
            return (xBelow1Count > 1) ? fitSL(w, x, y) : get68th(w, y);
        case DOUBLE_LINE:
            return (xBelow1Count > 2) ? fitDL(counter, w, x, y) : ((xBelow1Count > 1) ? fitSL(w, x, y) : get68th(w, y));
        }
        return -1;
    }
    double RCR::handleSigmaTechSelect(int counter, std::vector<double> &y)
    {
        std::vector<double> x = getXVec((int) y.size());
        int xBelow1Count = countAmountLessThanOne(x);
//        double nCorrection = nCorrect(counter);
        switch (sigmaTech)
        {
        case STANDARD_DEVIATION:
            //return (sigmaChoice == SINGLE) ? getStDev(2.0, y) : getStDev(1.0, y);
            return (sigmaChoice == SINGLE) ? getStDev(delta, y) : getStDev(delta/2.0, y);
        case SIXTY_EIGHTH_PERCENTILE:
            return get68th(y);
        case SINGLE_LINE:
            return (xBelow1Count > 1) ? fitSL(x, y) : get68th(y);
        case DOUBLE_LINE:
            return (xBelow1Count > 2) ? fitDL(counter, x, y) : ((xBelow1Count > 1) ? fitSL(x, y) : get68th(y));
        }
        return -1;
    }
    double RCR::handleBulkSigmaTechSelect(int counter, std::vector<double> &w, std::vector<double> &y)
    {
        std::vector<double> x = getXVec((int) y.size(), w);
        int xBelowOne = countAmountLessThanOne(x);
        if (xBelowOne > 2)
        {
            return max(fitDL(counter, w, x, y), fitSL(w, x, y));
        }
        else if (xBelowOne > 1)
        {
            return fitSL(w, x, y);
        }
        else
        {
            return get68th(w, y);
        }
    }
    double RCR::handleBulkSigmaTechSelect(int counter, std::vector<double> &y)
    {
        std::vector<double> x = getXVec((int) y.size());
        int xBelowOne = countAmountLessThanOne(x);
        if (xBelowOne > 2)
        {
            return max(fitDL(counter, x, y), fitSL(x, y));
        }
        else if (xBelowOne > 1)
        {
            return fitSL(x, y);
        }
        else
        {
            return get68th(y);
        }
    }


    //iterative chauvenet loops:
    void RCR::iterativeSingleSigmaRCR(std::vector<double> &w, std::vector<double> &y)
    {
        bool stop = false;
        int maxIndex, trueCount = (int) w.size();
        double mu = -1, stDev = -1, sigma = -1, hold, max; max = -99999;

        std::vector<bool> localFlags = result.flags;
        std::vector<int> indices;
        std::vector<double> diff, trueW, trueY, trueWHold, trueYHold;

        while (!stop)
        {
            //this->parametricModel->setTrueVec(localFlags, w, y);
            //this->parametricModel->buildModelSpace();
            //getMu(parametricModel);
            //trueWHold = trueW;
            //trueYHold = trueY;

            //trueY = handleMuTechSelect();
            //mu = 0;
            

            if (this->muType == PARAMETRIC)
            {
                parametricModel->setTrueVec(localFlags, w, y);
                mu = 0;
                parametricModel->buildModelSpace();
                this->delta = parametricModel->parameterSpace.size();
                trueY = handleMuTechSelect();
                indices = this->parametricModel->indices;
                trueW = parametricModel->trueW;
                trueYHold = trueY;
                trueWHold = trueW;

            }
            else if (this->muType == NONPARAMETRIC)
            {
                this->nonParametricModel->muFunc(localFlags, indices, w, y, trueW, trueY);

                trueYHold = trueY;
                trueWHold = trueW;

                mu = handleMuTechSelect(trueW, trueY);

                trueW = trueWHold;
                trueY = trueYHold;
                //indices = this->nonParametricModel->indices;

            }
            else
            {
                setTrueVec(localFlags, indices, w, y, trueW, trueY);

                trueYHold = trueY;
                trueWHold = trueW;

                mu = handleMuTechSelect(trueW, trueY);

                trueW = trueWHold;
                trueY = trueYHold;
            }

            trueCount = (int) trueW.size();
            max = -99999;
            trueW = trueWHold;
            trueY = trueYHold;
            diff.reserve(trueCount);
            for (size_t i = 0; i < trueY.size(); i++)
            {
                hold = std::abs(trueY[i] - mu);
                diff.push_back(hold);
                if (hold > max)
                {
                    max = hold;
                    maxIndex = indices[i];
                }
            }

            sort(trueW, diff);

            stDev = handleSigmaTechSelect(trueCount, trueW, diff);
            
            sigma = stDev * nCorrect(trueCount, trueW);

            diff.clear();
            indices.clear();

            if (this->muType == PARAMETRIC)
            {
                stop = reject((int) this->parametricModel->parameterSpace.size(), trueCount, maxIndex, max / sigma, localFlags, trueY);
            }
            else
            {
                stop = reject(trueCount, maxIndex, max / sigma, localFlags, y);
            }
        }
        result.flags = localFlags;

        this->result.mu = mu;
        this->result.stDev = stDev;
        this->result.sigma = sigma;
    }
    void RCR::iterativeSingleSigmaRCR(std::vector<double> &y)
    {
        bool stop = false;
        int maxIndex, trueCount = (int) y.size();
        double mu = -1, stDev = -1, sigma = -1, hold, max; max = -99999;

        std::vector<bool> localFlags = result.flags;
        std::vector<int> indices;

        std::vector<double> diff, trueY, trueYHold;

        while (!stop)
        {
            //this->parametricModel->setTrueVec(localFlags, y);
            //this->parametricModel->buildModelSpace();
            //getMu(parametricModel);
            //trueWHold = trueW;
            //trueYHold = trueY;
                    
            if (this->muType == PARAMETRIC)
            {
                parametricModel->setTrueVec(localFlags, y);
                mu = 0;
                parametricModel->buildModelSpace();
                this->delta = parametricModel->parameterSpace.size();
                trueY = handleMuTechSelect();
                indices = this->parametricModel->indices;

            }
            else if (this->muType == NONPARAMETRIC)
            {
                this->nonParametricModel->muFunc(localFlags, indices, y, trueY);

                trueYHold = trueY;

                mu = handleMuTechSelect((int) trueY.size(), trueY);

                trueY = trueYHold;
                //indices = this->nonParametricModel->indices;

            }
            else
            {
                setTrueVec(localFlags, indices, y, trueY);
                trueYHold = trueY;

                mu = handleMuTechSelect((int) trueY.size(), trueY);

                trueY = trueYHold;
            }
            //trueY = handleMuTechSelect();
            //mu = 0;

            trueCount = (int) trueY.size();
            //std::cout << trueCount << "\n";

            max = -99999;

            //indices = this->parametricModel->indices;
            diff.reserve(trueCount);
            for (size_t i = 0; i < trueY.size(); i++)
            {
                hold = std::abs(trueY[i] - mu);
                diff.push_back(hold);
                if (hold > max)
                {
                    max = hold;
                    maxIndex = indices[i];
                }
            }

            sort(diff);

            stDev = handleSigmaTechSelect(trueCount, diff);
            sigma =  stDev * nCorrect(trueCount);

            diff.clear();
            indices.clear();
            if (this->muType == PARAMETRIC)
            {
                stop = reject((int) this->parametricModel->parameterSpace.size(), trueCount, maxIndex, max / sigma, localFlags, trueY);
            }
            else
            {
                stop = reject(trueCount, maxIndex, max / sigma, localFlags, y);
            }
        }
        result.flags = localFlags;

        this->result.mu = mu;
        this->result.stDev = stDev;
        this->result.sigma = sigma;
    }
    void RCR::iterativeLowerSigmaRCR(std::vector<double> &w, std::vector<double> &y)
    {
        bool stop = false, nonzeroAbove, nonzeroBelow;
        int maxIndex, trueCount;
        double mu = -1, stDev = -1, sigma = -1, hold, max, nCorrection, stDevAbove = -1, stDevBelow = -1; // , sigmaBelow = 0, sigmaAbove = 0,

        std::vector<bool> localFlags = result.flags;
        std::vector<int> indices;

        std::vector<double> diffBelow, diffAbove, wBelow, wAbove, trueW, trueY, trueWHold, trueYHold;

        while (!stop)
        {

            if (this->muType == PARAMETRIC)
            {
                parametricModel->setTrueVec(localFlags, w, y);
                mu = 0;
                parametricModel->buildModelSpace();
                this->delta = parametricModel->parameterSpace.size();

                trueY = handleMuTechSelect();

                indices = this->parametricModel->indices;
                trueW = this->parametricModel->trueW;
                trueYHold = trueY;
                trueWHold = trueW;

            }
            else if (this->muType == NONPARAMETRIC)
            {
                this->nonParametricModel->muFunc(localFlags, indices, w, y, trueW, trueY);

                trueWHold = trueW;
                trueYHold = trueY;

                mu = handleMuTechSelect(trueW, trueY);

                trueW = trueWHold;
                trueY = trueYHold;
                indices = this->nonParametricModel->indices;

            }
            else
            {
                setTrueVec(localFlags, indices, w, y, trueW, trueY);
                trueWHold = trueW;
                trueYHold = trueY;

                mu = handleMuTechSelect(trueW, trueY);

                trueW = trueWHold;
                trueY = trueYHold;
            }


            trueCount = (int) trueW.size();
            max = -99999;
            diffBelow.reserve(trueCount);
            diffAbove.reserve(trueCount);
            wBelow.reserve(trueCount);
            wAbove.reserve(trueCount);
            for (size_t i = 0; i < trueY.size(); i++)
            {
                hold = getDiff(mu, trueY[i]);
                if (hold > max)
                {
                    max = hold;
                    maxIndex = indices[i];
                }
                if (isEqual(trueY[i], mu))
                {
                    diffBelow.push_back(hold);
                    diffAbove.push_back(hold);
                    wBelow.push_back(0.5*trueW[i]);
                    wAbove.push_back(0.5*trueW[i]);
                }
                else if (trueY[i] > mu)
                {
                    diffAbove.push_back(hold);
                    wAbove.push_back(trueW[i]);
                }
                else
                {
                    diffBelow.push_back(hold);
                    wBelow.push_back(trueW[i]);
                }
            }
            nonzeroAbove = diffAbove.size() > 0;
            nonzeroBelow = diffBelow.size() > 0;

            if (nonzeroBelow)
            {
                sort(wBelow, diffBelow);
            }
            if (nonzeroAbove)
            {
                sort(wAbove, diffAbove);
            }
            nCorrection = nCorrect(trueCount, trueW);

            if (nonzeroAbove && nonzeroBelow)
            {
                stDevAbove = handleSigmaTechSelect(trueCount, wAbove, diffAbove);
                stDevBelow = handleSigmaTechSelect(trueCount, wBelow, diffBelow);
                stDev = min(stDevAbove, stDevBelow);
                sigma = stDev * nCorrection;
            }
            else if (nonzeroAbove)
            {
                stDev = handleSigmaTechSelect(trueCount, wAbove, diffAbove);
                stDevAbove = stDev;
                sigma = stDev * nCorrection;
            }
            else if (nonzeroBelow)
            {
                stDev = handleSigmaTechSelect(trueCount, wBelow, diffBelow);
                stDevBelow = stDev;
                sigma = stDev * nCorrection;
            }

            //sort(wBelow, diffBelow);
            //sort(wAbove, diffAbove);
            //sigma = min(handleSigmaTechSelect(trueCount, wBelow, diffBelow), handleSigmaTechSelect(trueCount, wAbove, diffAbove)) * nCorrect(trueCount, trueW);

            diffBelow.clear();
            diffAbove.clear();
            wBelow.clear();
            wAbove.clear();

            if (this->muType == PARAMETRIC)
            {
                stop = reject((int) this->parametricModel->parameterSpace.size(), trueCount, maxIndex, max / sigma, localFlags, trueY);
            }
            else
            {
                stop = reject(trueCount, maxIndex, max / sigma, localFlags, y);
            }
        }
        result.flags = localFlags;

        this->result.mu = mu;
        this->result.sigma = sigma;
        this->result.stDev = stDev;
        this->result.stDevAbove = stDevAbove;
        this->result.stDevBelow = stDevBelow;
    }
    void RCR::iterativeLowerSigmaRCR(std::vector<double> &y)
    {
        bool stop = false, split = false, nonzeroAbove, nonzeroBelow;
        int maxIndex, trueCount, belowSplitIndex = -1, aboveSplitIndex = -1;
        double mu = -1, sigma = -1, hold, max, stDev = -1, nCorrection, stDevAbove = -1, stDevBelow = -1; // sigmaBelow = 0, sigmaAbove = 0,

        std::vector<bool> localFlags = result.flags;
        std::vector<int> indices;

        std::vector<double> diffBelow, diffAbove, wBelow, wAbove, trueY, trueYHold;

        while (!stop)
        {
            //std::cout << trueY.size() << "\n";

            if (this->muType == PARAMETRIC)
            {
                parametricModel->setTrueVec(localFlags, y);
                mu = 0;
                parametricModel->buildModelSpace();
                this->delta = parametricModel->parameterSpace.size();
                trueY = handleMuTechSelect();
                indices = this->parametricModel->indices;

            }
            else if (this->muType == NONPARAMETRIC)
            {
                this->nonParametricModel->muFunc(localFlags, indices, y, trueY);

                trueYHold = trueY;

                mu = handleMuTechSelect((int) trueY.size(), trueY);

                trueY = trueYHold;
                indices = this->nonParametricModel->indices;

            }
            else
            {
                setTrueVec(localFlags, indices, y, trueY);
                trueYHold = trueY;

                mu = handleMuTechSelect((int) trueY.size(), trueY);

                trueY = trueYHold;
            }
            trueCount = (int) trueY.size();
            //std::cout << trueCount << "\n";

            max = -99999;
            split = false;
            diffBelow.reserve(trueCount);
            diffAbove.reserve(trueCount);
            for (size_t i = 0; i < trueY.size(); i++)
            {
                hold = getDiff(mu, trueY[i]);
                if (hold > max)
                {
                    max = hold;
                    maxIndex = indices[i];
                }
                if (isEqual(trueY[i], mu))
                {
                    diffBelow.push_back(hold);
                    diffAbove.push_back(hold);
                    split = true;
                    belowSplitIndex = (int) diffBelow.size() - 1;
                    aboveSplitIndex = (int) diffAbove.size() - 1;
                }
                else if (trueY[i] > mu)
                {
                    diffAbove.push_back(hold);
                }
                else
                {
                    diffBelow.push_back(hold);
                }

            }
            wBelow.reserve(diffBelow.size());
            wAbove.reserve(diffAbove.size());
            wBelow.resize(diffBelow.size(), 1.0);
            wAbove.resize(diffAbove.size(), 1.0);
            if (split)
            {
                wBelow[belowSplitIndex] = 0.5;
                wAbove[aboveSplitIndex] = 0.5;
            }

            nonzeroAbove = diffAbove.size() > 0;
            nonzeroBelow = diffBelow.size() > 0;
            
            if (nonzeroBelow)
            {
                sort(wBelow, diffBelow);
            }
            if (nonzeroAbove)
            {
                sort(wAbove, diffAbove);
            }
            
            nCorrection = nCorrect(trueCount);

            if (nonzeroAbove && nonzeroBelow)
            {
                stDevAbove = handleSigmaTechSelect(trueCount, wAbove, diffAbove);
                stDevBelow = handleSigmaTechSelect(trueCount, wBelow, diffBelow);
                stDev = min(stDevAbove, stDevBelow);
                sigma = stDev * nCorrection;
            }
            else if (nonzeroAbove)
            {
                stDev = handleSigmaTechSelect(trueCount, wAbove, diffAbove);
                stDevAbove = stDev;
                sigma = stDev * nCorrection;
            }
            else if (nonzeroBelow)
            {
                stDev = handleSigmaTechSelect(trueCount, wBelow, diffBelow);
                stDevBelow = stDev;
                sigma = stDev * nCorrection;
            }
            //sort(wBelow, diffBelow);
            //sort(wAbove, diffAbove);
            //sigma = min(handleSigmaTechSelect(trueCount, wBelow, diffBelow), handleSigmaTechSelect(trueCount, wAbove, diffAbove)) * nCorrect(trueCount);

            diffBelow.clear();
            diffAbove.clear();
            wBelow.clear();
            wAbove.clear();

            if (this->muType == PARAMETRIC)
            {
                stop = reject((int) this->parametricModel->parameterSpace.size(), trueCount, maxIndex, max / sigma, localFlags, trueY);
            }
            else
            {
                stop = reject(trueCount, maxIndex, max / sigma, localFlags, y);
            }
        }
        result.flags = localFlags;

        this->result.mu = mu;
        this->result.sigma = sigma;
        this->result.stDev = stDev;
        this->result.stDevAbove = stDevAbove;
        this->result.stDevBelow = stDevBelow;
        
    }
    void RCR::iterativeEachSigmaRCR(std::vector<double> &w, std::vector<double> &y)
    {
        bool stop = false, nonzeroAbove, nonzeroBelow;
        int maxIndex, trueCount;
        double mu = -1, sigmaBelow = -1, sigmaAbove = -1, stDevBelow = -1, stDevAbove = -1, nCorrection, hold, max = -99999;

        std::vector<bool> localFlags = result.flags;
        std::vector<int> indices;

        std::vector<double> diffBelow, diffAbove, wBelow, wAbove, trueW, trueY, trueWHold, trueYHold;

        while (!stop)
        {
            if (this->muType == PARAMETRIC)
            {
                parametricModel->setTrueVec(localFlags, w, y);
                mu = 0;
                parametricModel->buildModelSpace();
                this->delta = parametricModel->parameterSpace.size();

                trueY = handleMuTechSelect();

                indices = this->parametricModel->indices;
                trueW = this->parametricModel->trueW;
                trueYHold = trueY;
                trueWHold = trueW;

            }
            else if (this->muType == NONPARAMETRIC)
            {
                this->nonParametricModel->muFunc(localFlags, indices, w, y, trueW, trueY);

                trueWHold = trueW;
                trueYHold = trueY;

                mu = handleMuTechSelect(trueW, trueY);

                trueW = trueWHold;
                trueY = trueYHold;
                indices = this->nonParametricModel->indices;

            }
            else
            {
                setTrueVec(localFlags, indices, w, y, trueW, trueY);
                trueWHold = trueW;
                trueYHold = trueY;

                mu = handleMuTechSelect(trueW, trueY);

                trueW = trueWHold;
                trueY = trueYHold;
            }

            trueCount = (int) trueW.size();
            diffBelow.reserve(trueCount);
            diffAbove.reserve(trueCount);
            wBelow.reserve(trueCount);
            wAbove.reserve(trueCount);
            for (size_t i = 0; i < trueY.size(); i++)
            {
                hold = getDiff(mu, trueY[i]);
                if (isEqual(trueY[i],mu))
                {
                    diffBelow.push_back(hold);
                    diffAbove.push_back(hold);
                    wBelow.push_back(0.5*trueW[i]);
                    wAbove.push_back(0.5*trueW[i]);
                }
                else if (trueY[i] > mu)
                {
                    diffAbove.push_back(hold);
                    wAbove.push_back(trueW[i]);
                }
                else
                {
                    diffBelow.push_back(hold);
                    wBelow.push_back(trueW[i]);
                }
            }

            //sort(wBelow, diffBelow);
            //sort(wAbove, diffAbove);

            //nCorrection = nCorrect(trueCount, trueW);
            //sigmaBelow = handleSigmaTechSelect(trueCount, wBelow, diffBelow) * nCorrection;
            //sigmaAbove = handleSigmaTechSelect(trueCount, wAbove, diffAbove) * nCorrection;

            max = -99999;
            //nCorrection = nCorrect(trueCount);
            nonzeroAbove = diffAbove.size() > 0;
            nonzeroBelow = diffBelow.size() > 0;
            nCorrection = nCorrect(trueCount, trueW);

            if (nonzeroBelow)
            {
                sort(wBelow, diffBelow);
                stDevBelow = handleSigmaTechSelect(trueCount, wBelow, diffBelow);
                sigmaBelow = stDevBelow * nCorrection;
            }
            if (nonzeroAbove)
            {
                sort(wAbove, diffAbove);
                stDevAbove = handleSigmaTechSelect(trueCount, wAbove, diffAbove);
                sigmaAbove = stDevAbove * nCorrection;
            }


            for (size_t i = 0; i < y.size(); i++)
            {
                if (localFlags[i])
                {
                    hold = std::abs(y[i] - mu);
                    if (y[i] < mu && (hold / sigmaBelow) > max)
                    {
                        max = hold / sigmaBelow;
                        maxIndex = (int) i;
                    }
                    if (y[i] > mu && (hold / sigmaAbove) > max)
                    {
                        max = hold / sigmaAbove;
                        maxIndex = (int) i;
                    }
                }
            }

            diffBelow.clear();
            diffAbove.clear();
            wBelow.clear();
            wAbove.clear();

            if (this->muType == PARAMETRIC)
            {
                stop = reject((int) this->parametricModel->parameterSpace.size(), trueCount, maxIndex, max, localFlags, trueY);
            }
            else
            {
                stop = reject(trueCount, maxIndex, max, localFlags, y);
            }


        }
        result.flags = localFlags;

        this->result.mu = mu;
        this->result.stDevBelow = stDevBelow;
        this->result.stDevAbove = stDevAbove;
        this->result.sigmaBelow = sigmaBelow;
        this->result.sigmaAbove = sigmaAbove;
    }
    void RCR::iterativeEachSigmaRCR(std::vector<double> &y)
    {
        bool stop = false, split = false, nonzeroAbove, nonzeroBelow;
        int maxIndex, trueCount, belowSplitIndex = -1, aboveSplitIndex = -1;
        double mu = -1, sigmaBelow = -1, sigmaAbove = -1, stDevAbove = -1, stDevBelow = -1, nCorrection, hold, max = -99999;

        std::vector<bool> localFlags = result.flags;
        std::vector<int> indices;

        std::vector<double> diffBelow, diffAbove, wBelow, wAbove, trueY, trueYHold;

        while (!stop)
        {
            //std::cout << trueY.size() << "\n";

            if (this->muType == PARAMETRIC)
            {
                parametricModel->setTrueVec(localFlags, y);
                mu = 0;
                parametricModel->buildModelSpace();
                this->delta = parametricModel->parameterSpace.size();
                trueY = handleMuTechSelect();
                indices = this->parametricModel->indices;

            }
            else if (this->muType == NONPARAMETRIC)
            {
                this->nonParametricModel->muFunc(localFlags, indices, y, trueY);

                trueYHold = trueY;

                mu = handleMuTechSelect((int) trueY.size(), trueY);

                trueY = trueYHold;
                indices = this->nonParametricModel->indices;

            }
            else
            {
                setTrueVec(localFlags, indices, y, trueY);
                trueYHold = trueY;

                mu = handleMuTechSelect((int) trueY.size(), trueY);

                trueY = trueYHold;
            }

            trueCount = (int) trueY.size();
            diffBelow.reserve(trueCount);
            diffAbove.reserve(trueCount);
            split = false;
            for (size_t i = 0; i < trueY.size(); i++)
            {
                hold = getDiff(mu, trueY[i]);
                if (isEqual(trueY[i], mu))
                {
                    diffBelow.push_back(hold);
                    diffAbove.push_back(hold);
                    split = true;
                    belowSplitIndex = (int) diffBelow.size() - 1;
                    aboveSplitIndex = (int) diffAbove.size() - 1;
                }
                else if (trueY[i] > mu)
                {
                    diffAbove.push_back(hold);
                }
                else
                {
                    diffBelow.push_back(hold);
                }
            }
            wBelow.resize(diffBelow.size(), 1.0);
            wAbove.resize(diffAbove.size(), 1.0);
            if (split)
            {
                wBelow[belowSplitIndex] = 0.5;
                wAbove[aboveSplitIndex] = 0.5;
            }

            //sort(wBelow, diffBelow);
            //sort(wAbove, diffAbove);

            //nCorrection = nCorrect(trueCount);
            //sigmaBelow = handleSigmaTechSelect(trueCount, wBelow, diffBelow) * nCorrection;
            //sigmaAbove = handleSigmaTechSelect(trueCount, wAbove, diffAbove) * nCorrection;
        
            nonzeroAbove = diffAbove.size() > 0;
            nonzeroBelow = diffBelow.size() > 0;
            nCorrection = nCorrect(trueCount);
            if (nonzeroBelow)
            {
                sort(wBelow, diffBelow);
                stDevBelow = handleSigmaTechSelect(trueCount, wBelow, diffBelow);
                sigmaBelow = stDevBelow * nCorrection;
            }
            if (nonzeroAbove)
            {
                sort(wAbove, diffAbove);
                stDevAbove = handleSigmaTechSelect(trueCount, wAbove, diffAbove);
                sigmaAbove = stDevAbove * nCorrection;
            }

            max = -99999;

            for (size_t i = 0; i < y.size(); i++)
            {
                if (localFlags[i])
                {
                    hold = std::abs(y[i] - mu);
                    if (y[i] < mu && (hold / sigmaBelow) > max)
                    {
                        max = hold / sigmaBelow;
                        maxIndex = (int) i;
                    }
                    if (y[i] > mu && (hold / sigmaAbove) > max)
                    {
                        max = hold / sigmaAbove;
                        maxIndex = (int) i;
                    }
                }
            }

            diffBelow.clear();
            diffAbove.clear();
            wBelow.clear();
            wAbove.clear();

            if (this->muType == PARAMETRIC)
            {
                stop = reject((int) this->parametricModel->parameterSpace.size(), trueCount, maxIndex, max, localFlags, trueY);
            }
            else
            {
                stop = reject(trueCount, maxIndex, max, localFlags, y);
            }

        }
        result.flags = localFlags;

        this->result.mu = mu;
        this->result.stDevBelow = stDevBelow;
        this->result.stDevAbove = stDevAbove;
        this->result.sigmaBelow = sigmaBelow;
        this->result.sigmaAbove = sigmaAbove;
    }


    //bulk chauvenet loops:
    void RCR::bulkSingleSigmaRCR(std::vector<double> &w, std::vector<double> &y)
    {
        bool stop = false;
        int trueCount = (int) w.size();
        double mu, sigma, hold;

        std::vector<bool> localFlags = result.flags;
        std::vector<int> indices;
        std::vector<double> diff, diffHold, trueW, trueY, trueWHold, trueYHold;

        while (!stop)
        {
            //this->parametricModel->setTrueVec(localFlags, w, y);

            //trueWHold = trueW;
            //trueYHold = trueY;

            //mu = handleMuTechSelect(trueW, trueY);

            //trueW = trueWHold;
            //trueY = trueYHold;
            if (this->muType == PARAMETRIC)
            {
                parametricModel->setTrueVec(localFlags, w, y);
                trueW = parametricModel->trueW;
                mu = 0;
                parametricModel->buildModelSpace();
                this->delta = parametricModel->parameterSpace.size();
                trueY = handleMuTechSelect();
                indices = this->parametricModel->indices;

            }
            else if (this->muType == NONPARAMETRIC)
            {
                this->nonParametricModel->muFunc(localFlags, indices, w, y, trueW, trueY);

                trueYHold = trueY;
                trueWHold = trueW;

                mu = handleMuTechSelect(trueW, trueY);

                trueW = trueWHold;
                trueY = trueYHold;
                indices = this->nonParametricModel->indices;

            }
            else
            {
                setTrueVec(localFlags, indices, w, y, trueW,  trueY);

                trueYHold = trueY;
                trueWHold = trueW;

                mu = handleMuTechSelect(trueW, trueY);

                trueW = trueWHold;
                trueY = trueYHold;
            }


            trueCount = (int) trueW.size();
            diff.reserve(trueCount);

            for (size_t i = 0; i < trueY.size(); i++)
            {
                hold = std::abs(trueY[i] - mu);
                diff.push_back(hold);
            }

            diffHold = diff;
            sort(trueW, diff);
            sort(indices, diffHold);

            sigma = handleBulkSigmaTechSelect(trueCount, trueW, diff) * nCorrect(trueCount, trueW);

            for (size_t i = 0; i < diffHold.size(); i++)
            {
                diffHold[i] = diffHold[i] / sigma;
            }

            if (this->muType == PARAMETRIC)
            {
                stop = bulkReject((int) this->parametricModel->parameterSpace.size(), localFlags, indices, diffHold, trueY);
            }
            else
            {
                stop = bulkReject(localFlags, indices, diffHold, y);
            }

            diff.clear();
            indices.clear();
        }
        result.flags = localFlags;
    }
    void RCR::bulkSingleSigmaRCR(std::vector<double> &y)
    {
        bool stop = false;
        int trueCount = (int) y.size();
        double mu, sigma, hold;

        std::vector<bool> localFlags = result.flags;
        std::vector<int> indices;
        std::vector<double> diff, diffHold, trueY, trueYHold;;

        while (!stop)
        {
            //this->parametricModel->setTrueVec(localFlags, y);
            //this->parametricModel->buildModelSpace();
            //getMu(parametricModel);
            //trueWHold = trueW;
            //trueYHold = trueY;

            //trueY = handleMuTechSelect();
            //mu = 0;

            if (this->muType == PARAMETRIC)
            {
                parametricModel->setTrueVec(localFlags, y);
                mu = 0;
                parametricModel->buildModelSpace();
                this->delta = parametricModel->parameterSpace.size();
                trueY = handleMuTechSelect();
                indices = this->parametricModel->indices;

            }
            else if (this->muType == NONPARAMETRIC)
            {
                this->nonParametricModel->muFunc(localFlags, indices, y, trueY);

                trueYHold = trueY;

                mu = handleMuTechSelect((int) trueY.size(), trueY);

                trueY = trueYHold;
                //indices = this->nonParametricModel->indices;

            }
            else
            {
                setTrueVec(localFlags, indices, y, trueY);

                trueYHold = trueY;

                mu = handleMuTechSelect((int) trueY.size(), trueY);

                trueY = trueYHold;
            }

            trueCount = (int) trueY.size();
            //std::cout << trueCount << "\n";

            //indices = this->parametricModel->indices;
            diff.reserve(trueCount);

            for (size_t i = 0; i < trueY.size(); i++)
            {
                hold = std::abs(trueY[i] - mu);
                diff.push_back(hold);
            }

            diffHold = diff;
            sort(diff);
            sort(indices, diffHold);

            sigma = handleBulkSigmaTechSelect(trueCount, diff) * nCorrect(trueCount);

            for (size_t i = 0; i < diffHold.size(); i++)
            {
                diffHold[i] = diffHold[i] / sigma;
            }

            if (this->muType == PARAMETRIC)
            {
                stop = bulkReject((int) this->parametricModel->parameterSpace.size(), localFlags, indices, diffHold, trueY);
            }
            else
            {
                stop = bulkReject(localFlags, indices, diffHold, y);
            }

            diff.clear();
        }
        result.flags = localFlags;
    }
    void RCR::bulkLowerSigmaRCR(std::vector<double> &w, std::vector<double> &y)
    {
        bool stop = false, nonzeroAbove, nonzeroBelow;
        int trueCount = (int) w.size();
        double mu, sigma = -1, hold;

        std::vector<bool> localFlags = result.flags;
        std::vector<int> indices;
        std::vector<double> diffBelow, diffAbove, wBelow, wAbove, diffHold, trueW, trueY, trueWHold, trueYHold;

        while (!stop)
        {

            if (this->muType == PARAMETRIC)
            {
                parametricModel->setTrueVec(localFlags, w, y);
                mu = 0;
                trueW = parametricModel->trueW;
                parametricModel->buildModelSpace();
                this->delta = parametricModel->parameterSpace.size();
                trueY = handleMuTechSelect();
                indices = this->parametricModel->indices;

            }
            else if (this->muType == NONPARAMETRIC)
            {
                this->nonParametricModel->muFunc(localFlags, indices, w, y, trueW, trueY);

                trueWHold = trueW;
                trueYHold = trueY;

                mu = handleMuTechSelect(trueW, trueY);

                trueW = trueWHold;
                trueY = trueYHold;
                indices = this->nonParametricModel->indices;

            }
            else
            {
                setTrueVec(localFlags, indices, w, y, trueW, trueY);
                trueWHold = trueW;
                trueYHold = trueY;

                mu = handleMuTechSelect(trueW, trueY);

                trueW = trueWHold;
                trueY = trueYHold;
            }


            trueCount = (int) trueW.size();
            diffHold.reserve(trueCount);
            diffBelow.reserve(trueCount);
            diffAbove.reserve(trueCount);
            wBelow.reserve(trueCount);
            wAbove.reserve(trueCount);
            for (size_t i = 0; i < trueY.size(); i++)
            {
                hold = getDiff(mu, trueY[i]);
                diffHold.push_back(hold);
                if (isEqual(trueY[i], mu))
                {
                    diffBelow.push_back(hold);
                    diffAbove.push_back(hold);
                    wBelow.push_back(0.5*trueW[i]);
                    wAbove.push_back(0.5*trueW[i]);
                }
                else if (trueY[i] > mu)
                {
                    diffAbove.push_back(hold);
                    wAbove.push_back(trueW[i]);
                }
                else
                {
                    diffBelow.push_back(hold);
                    wBelow.push_back(trueW[i]);
                }
            }

            sort(indices, diffHold);
            //sort(wBelow, diffBelow);
            //sort(wAbove, diffAbove);

            //double nCorrection = nCorrect(trueCount, trueW);
            //std::cout << trueCount << "\n";
            //sigma = min(handleBulkSigmaTechSelect(trueCount, wBelow, diffBelow), handleBulkSigmaTechSelect(trueCount, wAbove, diffAbove)) * nCorrect(trueCount, trueW);
            nonzeroAbove = diffAbove.size() > 0;
            nonzeroBelow = diffBelow.size() > 0;

            if (nonzeroBelow)
            {
                sort(wBelow, diffBelow);
            }
            if (nonzeroAbove)
            {
                sort(wAbove, diffAbove);
            }
            if (nonzeroAbove && nonzeroBelow)
            {
                sigma = min(handleBulkSigmaTechSelect(trueCount, wBelow, diffBelow), handleBulkSigmaTechSelect(trueCount, wAbove, diffAbove)) * nCorrect(trueCount, trueW);
            }
            else if (nonzeroAbove)
            {
                sigma = handleBulkSigmaTechSelect(trueCount, wAbove, diffAbove) * nCorrect(trueCount, trueW);
            }
            else if (nonzeroBelow)
            {
                sigma = handleBulkSigmaTechSelect(trueCount, wBelow, diffBelow) * nCorrect(trueCount, trueW);
            }

            for (size_t i = 0; i < diffHold.size(); i++)
            {
                diffHold[i] = diffHold[i] / sigma;
            }

            if (this->muType == PARAMETRIC)
            {
                stop = bulkReject((int) this->parametricModel->parameterSpace.size(), localFlags, indices, diffHold, trueY);
            }
            else
            {
                stop = bulkReject(localFlags, indices, diffHold, y);
            }

            diffBelow.clear();
            diffAbove.clear();
            wBelow.clear();
            wAbove.clear();
            diffHold.clear();
        }
        result.flags = localFlags;

    }
    void RCR::bulkLowerSigmaRCR(std::vector<double> &y)
    {
        bool stop = false, split = false, nonzeroAbove, nonzeroBelow;
        int trueCount = (int) y.size(), belowSplitIndex = -1, aboveSplitIndex = -1;
        double mu, sigma = -1, hold;

        std::vector<bool> localFlags = result.flags;
        std::vector<int> indices;
        std::vector<double> diffBelow, diffAbove, wBelow, wAbove, diffHold, trueY, trueYHold;

        while (!stop)
        {

            if (this->muType == PARAMETRIC)
            {
                parametricModel->setTrueVec(localFlags, y);
                mu = 0;
                parametricModel->buildModelSpace();
                this->delta = parametricModel->parameterSpace.size();
                trueY = handleMuTechSelect();
                indices = this->parametricModel->indices;

            }
            else if (this->muType == NONPARAMETRIC)
            {
                this->nonParametricModel->muFunc(localFlags, indices, y, trueY);
                trueYHold = trueY;

                mu = handleMuTechSelect((int) trueY.size(), trueY);
                indices = this->nonParametricModel->indices;

                trueY = trueYHold;
            }
            else
            {
                setTrueVec(localFlags, indices, y, trueY);
                trueYHold = trueY;

                mu = handleMuTechSelect((int) trueY.size(), trueY);

                trueY = trueYHold;
            }
            //mu = 0;
            split = false;
            diffBelow.reserve(trueCount);
            diffAbove.reserve(trueCount);
            for (size_t i = 0; i < trueY.size(); i++)
            {
                hold = getDiff(mu, trueY[i]);
                diffHold.push_back(hold);
                if (isEqual(trueY[i], mu))
                {
                    diffBelow.push_back(hold);
                    diffAbove.push_back(hold);
                    belowSplitIndex = (int) diffBelow.size() - 1;
                    aboveSplitIndex = (int) diffAbove.size() - 1;
                    split = true;
                }
                else if (trueY[i] > mu)
                {
                    diffAbove.push_back(hold);
                }
                else
                {
                    diffBelow.push_back(hold);
                }
            }
            wBelow.reserve(diffBelow.size());
            wAbove.reserve(diffAbove.size());
            wBelow.resize(diffBelow.size(), 1.0);
            wAbove.resize(diffAbove.size(), 1.0);

            if (split)
            {
                wBelow[belowSplitIndex] = 0.5;
                wAbove[aboveSplitIndex] = 0.5;
            }
            sort(indices, diffHold);
            //sort(wBelow, diffBelow);
            //sort(wAbove, diffAbove);

            //sigma = min(handleBulkSigmaTechSelect(trueCount, wBelow, diffBelow), handleBulkSigmaTechSelect(trueCount, wAbove, diffAbove)) * nCorrect(trueCount);
            nonzeroAbove = diffAbove.size() > 0;
            nonzeroBelow = diffBelow.size() > 0;

            if (nonzeroBelow)
            {
                sort(wBelow, diffBelow);
            }
            if (nonzeroAbove)
            {
                sort(wAbove, diffAbove);
            }
            if (nonzeroAbove && nonzeroBelow)
            {
                sigma = min(handleBulkSigmaTechSelect(trueCount, wBelow, diffBelow), handleBulkSigmaTechSelect(trueCount, wAbove, diffAbove)) * nCorrect(trueCount);
            }
            else if (nonzeroAbove)
            {
                sigma = handleBulkSigmaTechSelect(trueCount, wAbove, diffAbove) * nCorrect(trueCount);
            }
            else if (nonzeroBelow)
            {
                sigma = handleBulkSigmaTechSelect(trueCount, wBelow, diffBelow) * nCorrect(trueCount);
            }

            for (size_t i = 0; i < diffHold.size(); i++)
            {
                diffHold[i] = diffHold[i] / sigma;
            }
            if (this->muType == PARAMETRIC)
            {
                stop = bulkReject((int) this->parametricModel->parameterSpace.size(), localFlags, indices, diffHold, trueY);
            }
            else
            {
                stop = bulkReject(localFlags, indices, diffHold, y);
            }

            diffBelow.clear();
            diffAbove.clear();
            wBelow.clear();
            wAbove.clear();
            diffHold.clear();
        }
        result.flags = localFlags;

    }
    void RCR::bulkEachSigmaRCR(std::vector<double> &w, std::vector<double> &y)
    {
        bool stop = false, nonzeroBelow, nonzeroAbove;
        int trueCount = (int) w.size();
        double mu, sigmaBelow = -1, sigmaAbove = -1, hold; //, max = -99999;

        std::vector<bool> localFlags = result.flags;
        std::vector<int> indices;

        std::vector<double> diffHold, diffBelow, diffAbove, wBelow, wAbove, trueW, trueY, trueWHold, trueYHold;

        while (!stop)
        {
            /*
            this->parametricModel->setTrueVec(localFlags, w, y);

            trueWHold = trueW;
            trueYHold = trueY;

            mu = handleMuTechSelect(trueW, trueY);

            trueW = trueWHold;
            trueY = trueYHold;
            */
            if (this->muType == PARAMETRIC)
            {
                parametricModel->setTrueVec(localFlags, w, y);
                mu = 0;
                trueW = parametricModel->trueW;
                parametricModel->buildModelSpace();
                this->delta = parametricModel->parameterSpace.size();
                trueY = handleMuTechSelect();
                indices = this->parametricModel->indices;

            }
            else if (this->muType == NONPARAMETRIC)
            {
                this->nonParametricModel->muFunc(localFlags, indices, w, y, trueW, trueY);

                trueYHold = trueY;
                trueWHold = trueW;

                mu = handleMuTechSelect(trueW, trueY);

                trueW = trueWHold;
                trueY = trueYHold;
                indices = this->nonParametricModel->indices;

            }
            else
            {
                setTrueVec(localFlags, indices, w, y, trueW, trueY);

                trueYHold = trueY;
                trueWHold = trueW;

                mu = handleMuTechSelect(trueW, trueY);

                trueW = trueWHold;
                trueY = trueYHold;
            }


            trueCount = (int) trueW.size();
            diffHold.reserve(trueCount);
            diffBelow.reserve(trueCount);
            diffAbove.reserve(trueCount);
            wBelow.reserve(trueCount);
            wAbove.reserve(trueCount);
            for (size_t i = 0; i < trueY.size(); i++)
            {
                hold = getDiff(mu, trueY[i]);
                diffHold.push_back(hold);

                if (isEqual(trueY[i],mu))
                {
                    diffBelow.push_back(hold);
                    diffAbove.push_back(hold);
                    wBelow.push_back(0.5*trueW[i]);
                    wAbove.push_back(0.5*trueW[i]);
                }
                else if (trueY[i] < mu)
                {
                    diffBelow.push_back(hold);
                    wBelow.push_back(trueW[i]);
                }
                else
                {
                    diffAbove.push_back(hold);
                    wAbove.push_back(trueW[i]);
                }
            }

            //sort(wBelow, diffBelow);
            //sort(wAbove, diffAbove);

            //nCorrection = nCorrect(trueCount, trueW);
            //sigmaBelow = handleBulkSigmaTechSelect(trueCount, wBelow, diffBelow) * nCorrection;
            //sigmaAbove = handleBulkSigmaTechSelect(trueCount, wAbove, diffAbove) * nCorrection;
            nonzeroAbove = diffAbove.size() > 0;
            nonzeroBelow = diffBelow.size() > 0;

            if (nonzeroBelow)
            {
                sort(wBelow, diffBelow);
                sigmaBelow = handleBulkSigmaTechSelect(trueCount, wBelow, diffBelow) * nCorrect(trueCount, trueW);
            }
            if (nonzeroAbove)
            {
                sort(wAbove, diffAbove);
                sigmaAbove = handleBulkSigmaTechSelect(trueCount, wAbove, diffAbove) * nCorrect(trueCount, trueW);
            }

            for (size_t i = 0; i < diffHold.size(); i++)
            {
                if (y[indices[i]] < mu)
                {
                    diffHold[i] = diffHold[i] / sigmaBelow;
                }
                if (y[indices[i]] > mu)
                {
                    diffHold[i] = diffHold[i] / sigmaAbove;
                }
            }
            sort(indices, diffHold);

            if (this->muType == PARAMETRIC)
            {
                stop = bulkReject((int) this->parametricModel->parameterSpace.size(), localFlags, indices, diffHold, trueY);
            }
            else
            {
                stop = bulkReject(localFlags, indices, diffHold, y);
            }


            diffBelow.clear();
            diffAbove.clear();
            wBelow.clear();
            wAbove.clear();
            diffHold.clear();
        }
        result.flags = localFlags;
    }
    void RCR::bulkEachSigmaRCR(std::vector<double> &y)
    {
        bool stop = false, split = false, nonzeroAbove = -1, nonzeroBelow = -1;
        int trueCount = (int) y.size(), belowSplitIndex = -1, aboveSplitIndex = -1;
        double mu, sigmaBelow = -1, sigmaAbove = -1, hold; //, max = -99999;

        std::vector<bool> localFlags = result.flags;

        std::vector<int> indices;
        std::vector<double> diffHold, diffBelow, diffAbove, wBelow, wAbove, trueY, trueYHold;

        while (!stop)
        {
            //this->parametricModel->setTrueVec(localFlags, y);
            //this->parametricModel->buildModelSpace();
            //getMu(parametricModel);
            //trueWHold = trueW;
            //trueYHold = trueY;

            //trueY = handleMuTechSelect();
            //mu = 0;

            if (this->muType == PARAMETRIC)
            {
                parametricModel->setTrueVec(localFlags, y);
                mu = 0;
                parametricModel->buildModelSpace();
                this->delta = parametricModel->parameterSpace.size();
                trueY = handleMuTechSelect();
                indices = this->parametricModel->indices;

            }
            else if (this->muType == NONPARAMETRIC)
            {
                this->nonParametricModel->muFunc(localFlags, indices, y, trueY);

                trueYHold = trueY;

                mu = handleMuTechSelect((int) trueY.size(), trueY);

                trueY = trueYHold;
                //indices = this->nonParametricModel->indices;

            }
            else
            {
                setTrueVec(localFlags, indices, y, trueY);

                trueYHold = trueY;

                mu = handleMuTechSelect((int) trueY.size(), trueY);

                trueY = trueYHold;
            }

            trueCount = (int) trueY.size();

            split = false;
            diffHold.reserve(trueCount);
            diffBelow.reserve(trueCount);
            diffAbove.reserve(trueCount);
            for (size_t i = 0; i < trueY.size(); i++)
            {
                hold = getDiff(mu, trueY[i]);
                diffHold.push_back(hold);
                if (isEqual(trueY[i],mu))
                {
                    diffBelow.push_back(hold);
                    diffAbove.push_back(hold);
                    belowSplitIndex = (int) diffBelow.size() - 1;
                    aboveSplitIndex = (int) diffAbove.size() - 1;
                    split = true;
                }
                else if (trueY[i] > mu)
                {
                    diffAbove.push_back(hold);
                }
                else
                {
                    diffBelow.push_back(hold);
                }
            }

            wBelow.resize(diffBelow.size(), 1.0);
            wAbove.resize(diffAbove.size(), 1.0);

            if (split)
            {
                wBelow[belowSplitIndex] = 0.5;
                wAbove[aboveSplitIndex] = 0.5;
            }
            //sort(wBelow, diffBelow);
            //sort(wAbove, diffAbove);

            //nCorrection = nCorrect(trueCount);

            //sigmaBelow = handleBulkSigmaTechSelect(trueCount, wBelow, diffBelow) * nCorrection;
            //sigmaAbove = handleBulkSigmaTechSelect(trueCount, wAbove, diffAbove) * nCorrection;
            
            nonzeroAbove = diffAbove.size() > 0;
            nonzeroBelow = diffBelow.size() > 0;

            if (nonzeroBelow)
            {
                sort(wBelow, diffBelow);
                sigmaBelow = handleBulkSigmaTechSelect(trueCount, wBelow, diffBelow) * nCorrect(trueCount);
            }
            if (nonzeroAbove)
            {
                sort(wAbove, diffAbove);
                sigmaAbove = handleBulkSigmaTechSelect(trueCount, wAbove, diffAbove) * nCorrect(trueCount);
            }

            for (size_t i = 0; i < diffHold.size(); i++)
            {
                if (y[indices[i]] < mu)
                {
                    diffHold[i] = diffHold[i] / sigmaBelow;
                }
                if (y[indices[i]] > mu)
                {
                    diffHold[i] = diffHold[i] / sigmaAbove;
                }
            }
            sort(indices, diffHold);

            if (this->muType == PARAMETRIC)
            {
                stop = bulkReject((int) this->parametricModel->parameterSpace.size(), localFlags, indices, diffHold, trueY);
            }
            else
            {
                stop = bulkReject(localFlags, indices, diffHold, y);
            }
            diffBelow.clear();
            diffAbove.clear();
            wBelow.clear();
            wAbove.clear();
            diffHold.clear();
        }
        result.flags = localFlags;
    }


    //chauvenet loop handlers
    void RCR::handleRCRLoopSelect(bool bulk, std::vector<double> &w, std::vector<double> &y)
    {
        if (bulk)
        {
            switch (sigmaChoice)
            {
            case SINGLE:
                return bulkSingleSigmaRCR(w, y);
            case LOWER:
                return bulkLowerSigmaRCR(w, y);
            case EACH:
                return bulkEachSigmaRCR(w, y);
            }
        }
        else
        {
            switch (sigmaChoice)
            {
            case SINGLE:
                return iterativeSingleSigmaRCR(w, y);
            case LOWER:
                return iterativeLowerSigmaRCR(w, y);
            case EACH:
                return iterativeEachSigmaRCR(w, y);
            }
        }

    }
    void RCR::handleRCRLoopSelect(bool bulk, std::vector<double> &y)
    {
        if (bulk)
        {
            switch (sigmaChoice)
            {
            case SINGLE:
                return bulkSingleSigmaRCR(y);
            case LOWER:
                return bulkLowerSigmaRCR(y);
            case EACH:
                return bulkEachSigmaRCR(y);
            }
        }
        else
        {
            switch (sigmaChoice)
            {
            case SINGLE:
                return iterativeSingleSigmaRCR(y);
            case LOWER:
                return iterativeLowerSigmaRCR(y);
            case EACH:
                return iterativeEachSigmaRCR(y);
            }
        }

    }


    //destructors
    RCR::~RCR()
    {
    }
}
