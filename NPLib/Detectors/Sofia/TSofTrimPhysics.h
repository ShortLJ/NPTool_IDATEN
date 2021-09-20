#ifndef TSofTrimPHYSICS_H
#define TSofTrimPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : November 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SofTrim Treated data                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// C++ headers 
#include <vector>
#include <map>
#include <string>
using namespace std;

// ROOT headers
#include "TObject.h"
#include "TH1.h"
#include "TVector3.h"
#include "TSpline.h"
// NPTool headers
#include "TSofTrimData.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"



class TSofTrimPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TSofTrimPhysics();
    ~TSofTrimPhysics() {};


  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};


  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    vector<int>      SectionNbr;
    vector<double>   EnergyPair1;
    vector<double>   EnergyPair2;
    vector<double>   EnergyPair3;
    vector<double>   DriftTimePair1;
    vector<double>   DriftTimePair2;
    vector<double>   DriftTimePair3;
    vector<double>   EnergySection;
    vector<double>   Theta;

  /// A usefull method to bundle all operation to add a detector
  void AddDetector(TVector3 POS); 
  void AddDetector(double R, double Theta, double Phi); 
  
  //////////////////////////////////////////////////////////////
  // methods inherited from the VDetector ABC class
  public:
    // read stream from ConfigFile to pick-up detector parameters
    void ReadConfiguration(NPL::InputParser);

    // add parameters to the CalibrationManger
    void AddParameterToCalibrationManager();

    // method called event by event, aiming at extracting the 
    // physical information from detector
    void BuildPhysicalEvent();

    // same as BuildPhysicalEvent() method but with a simpler
    // treatment
    void BuildSimplePhysicalEvent();

    // same as above but for online analysis
    void BuildOnlinePhysicalEvent()  {BuildPhysicalEvent();};

    // activate raw data object and branches from input TChain
    // in this method mother branches (Detector) AND daughter leaves 
    // (fDetector_parameter) have to be activated
    void InitializeRootInputRaw();

    // activate physics data object and branches from input TChain
    // in this method mother branches (Detector) AND daughter leaves 
    // (fDetector_parameter) have to be activated
    void InitializeRootInputPhysics();

    // create branches of output ROOT file
    void InitializeRootOutput();

    // clear the raw and physical data objects event by event
    void ClearEventPhysics() {Clear();}      
    void ClearEventData()    {m_EventData->Clear();}   


  //////////////////////////////////////////////////////////////
  // specific methods to SofTrim array
  public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();

    // clear the pre-treated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TSofTrimData object to TSofTrimPhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TSofTrimData* rawDataPointer) {m_EventData = rawDataPointer;}
   
    void LoadSplinePairAngle();
    void LoadSplinePairDriftTime();
    void LoadSplineSectionDriftTime();
    void LoadSplineSectionAngle();

    void SetBeta(double beta) {m_Beta = beta;}
    double GetBeta() {return m_Beta;}

    double GetMaxEnergySection();

  // objects are not written in the TTree
  private:
    TSofTrimData*         m_EventData;        //!
    TSofTrimData*         m_PreTreatedData;   //!
    TSofTrimPhysics*      m_EventPhysics;     //!

  // getters for raw and pre-treated data object
  public:
    TSofTrimData* GetRawData()        const {return m_EventData;}
    TSofTrimData* GetPreTreatedData() const {return m_PreTreatedData;}

  // parameters used in the analysis
  private:
    double m_E_Threshold;     //!
    double m_Beta;     //!
    double m_BetaNorm;     //!
    string m_SPLINE_PAIR_ANGLE_PATH;     //!
    string m_SPLINE_PAIR_DT_PATH;     //!
    string m_SPLINE_SECTION_DT_PATH;     //!
    string m_SPLINE_SECTION_ANGLE_PATH;     //!
    bool m_IsSplineSectionDriftTime; //!
    bool m_IsSplineSectionAngle; //!

  private:
    TSpline3* fcorr_EvsA[3][3]; //!
    TSpline3* fcorr_EvsDT[3][3]; //!
    TSpline3* fcorr_sec_angle[3]; //!
    TSpline3* fcorr_sec_dt[3]; //!

  // number of detectors
  private:
    int m_NumberOfDetectors;  //!
    int m_NumberOfSections;  //!
    int m_NumberOfAnodesPaired;  //!
    int m_NumberOfAnodesPerSection;  //!

  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TSofTrimPhysics,1)  // SofTrimPhysics structure
};
#endif
