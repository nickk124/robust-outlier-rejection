//#pragma once
#include <vector>
#include <cmath>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <float.h>
#include <iterator>

#include "FunctionalForm.h"
#include "NonParametric.h"
//#include "Demo.h"
//#include "RCRwebutils.h"

namespace RCRLib {

    const double PI = 3.1415926535897932384626434;
    const double inverfMult = (8.0 * (PI - 3.0)) / (3 * PI*(4.0 - PI));
    const double squareRootOf2 = sqrt(2.0);


    enum MuTechs { MEAN, MEDIAN, MODE };
    enum SigmaTechs { STANDARD_DEVIATION, SIXTY_EIGHTH_PERCENTILE, SINGLE_LINE, DOUBLE_LINE };
    enum SigmaChoices { SINGLE, LOWER, EACH };
    enum RejectionTechs { SS_MEDIAN_DL, LS_MODE_68, LS_MODE_DL, ES_MODE_DL };
    enum MuTypes {VALUE, PARAMETRIC, NONPARAMETRIC};

    struct RCRResults
    {
        double mu;
        double stDev;
        double stDevBelow;
        double stDevAbove;
        double stDevTotal;
        double sigma;
        double sigmaBelow;
        double sigmaAbove;
        std::vector<bool> flags;
        std::vector<int> indices;
        std::vector<double> cleanW;
        std::vector<double> cleanY;
        std::vector<double> rejectedW;
        std::vector<double> rejectedY;
        std::vector<double> originalW;
        std::vector<double> originalY;
    };

    typedef void(*MuCalcPrepWeighted)(std::vector<bool>&, std::vector<int>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
    typedef void(*MuCalcPrep)(std::vector<bool>&, std::vector<int>&, std::vector<double>&, std::vector<double>&);

    class RCR
    {

    public:

        RCR();
        RCR(RejectionTechs);
        RejectionTechs rejectionTech;
        RCRResults result;
        ~RCR();

        //standard procedures
        void setRejectionTech(RejectionTechs);
        void performRejection(std::vector<double> &);
        void performBulkRejection(std::vector<double> &);
        void performRejection(std::vector<double> &, std::vector<double> &);
        void performBulkRejection(std::vector<double> &, std::vector<double> &);
        void setParametricModel(FunctionalForm&);
        void setNonParametricModel(NonParametric&);

        //void revertMuType();
        void setInitialModel(std::vector<double> &);
        void setMuType(MuTypes);
        //void figureFunc(std::vector<double> &);


    private:

        static bool unityTablesLoaded;
        static std::vector<double> ESUnity;
        static std::vector<double> SSUnity;
        static std::vector<double> LSUnity;
        static std::vector<double> ESDLUnityCF;
        static std::vector<double> LSDLUnityCF;
        static std::vector<double> LS68UnityCF;
        static std::vector<double> SSDLUnityCF;
        static std::vector<std::vector<double> > SSConstants;
        std::vector<int> rejectedBy;

        MuTechs muTech;
        SigmaTechs sigmaTech;
        SigmaChoices sigmaChoice;
        FunctionalForm mF;
        FunctionalForm *parametricModel;
        NonParametric *nonParametricModel;
        MuTypes muType;

        bool useModifiedMuForm;
        double delta;


        //standard procedures
        void alignTechniques();
        void loadUnityTables();
        void setFinalVectors(std::vector<double> &);
        void setFinalVectors(std::vector<double> &, std::vector<double> &);
        void testFunc();


        //cf and fn models
        static double getLower68CF(int, std::vector<double>&);
        static double getLower68CF(int);
        static double getLowerDLCF(int, std::vector<double>&);
        static double getLowerDLCF(int);
        static double getSingleDLCF(int, std::vector<double>&);
        static double getSingleDLCF(int);
        static double getEachDLCF(int, std::vector<double>&);
        static double getEachDLCF(int);
        static double getSingleFN(int, std::vector<double>&, std::vector<double>&);
        static double getSingleFN(int, std::vector<double>&);
        static double getLowerFN(int, std::vector<double>&, std::vector<double>&);
        static double getLowerFN(int, std::vector<double>&);
        static double getEachFN(int, std::vector<double>&, std::vector<double>&);
        static double getEachFN(int, std::vector<double>&);

        //statistic tools
        double nCorrect(int, std::vector<double>&);
        double nCorrect(int);
        double getFN(int, std::vector<double>&, std::vector<double>&);
        double getFN(int, std::vector<double> &x);


        //rejectors:
        bool reject(int, int, int&, double, std::vector<bool>&, std::vector<double>&);
        bool reject(int, int&, double, std::vector<bool>&, std::vector<double>&);
        bool bulkReject(int, std::vector<bool>&, std::vector<int>&, std::vector<double>&, std::vector<double>&);
        bool bulkReject(std::vector<bool>&, std::vector<int>&, std::vector<double>&, std::vector<double>&);


        //sigma calculations:
        double fitDL(double, std::vector<double>&, std::vector<double>&, std::vector<double>&);
        double fitDL(double, std::vector<double>&, std::vector<double>&);


        //technique handlers:
        double handleMuTechSelect(std::vector<double>&, std::vector<double>&);
        double handleMuTechSelect(int, std::vector<double>&);
        std::vector<double> handleMuTechSelect();

        double handleSigmaTechSelect(int, std::vector<double>&, std::vector<double>&);
        double handleSigmaTechSelect(int, std::vector<double>&);
        double handleBulkSigmaTechSelect(int, std::vector<double>&, std::vector<double>&);
        double handleBulkSigmaTechSelect(int, std::vector<double>&);


        //iterative chauvenet loops:
        void iterativeSingleSigmaRCR(std::vector<double>&, std::vector<double>&);
        void iterativeSingleSigmaRCR(std::vector<double>&);
        void iterativeLowerSigmaRCR(std::vector<double>&, std::vector<double>&);
        void iterativeLowerSigmaRCR(std::vector<double>&);
        void iterativeEachSigmaRCR(std::vector<double>&, std::vector<double>&);
        void iterativeEachSigmaRCR(std::vector<double>&);


        //bulk chauvenet loops:
        void bulkSingleSigmaRCR(std::vector<double>&, std::vector<double>&);
        void bulkSingleSigmaRCR(std::vector<double>&);
        void bulkLowerSigmaRCR(std::vector<double>&, std::vector<double>&);
        void bulkLowerSigmaRCR(std::vector<double>&);
        void bulkEachSigmaRCR(std::vector<double>&, std::vector<double>&);
        void bulkEachSigmaRCR(std::vector<double>&);


        //chauvenet loop handlers:
        void handleRCRLoopSelect(bool, std::vector<double>&, std::vector<double>&);
        void handleRCRLoopSelect(bool, std::vector<double>&);



    };
}
