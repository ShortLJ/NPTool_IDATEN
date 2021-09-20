#ifndef TSofTofWPHYSICS_H
#define TSofTofWPHYSICS_H
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
 *  This class hold TofTofW Treated data                                *
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
#include "TRandom3.h"
// NPTool headers
#include "TSofTofWData.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"



class TSofTofWPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TSofTofWPhysics();
    ~TSofTofWPhysics() {};


  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};


  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    vector<int>      PlasticNbr;
    vector<double>   TimeNs;
    vector<double>   RawPosY;
    vector<double>   CalPosY;
    vector<double>   RawTof;
    vector<double>   CalTof;

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

    double CalculateTimeNs(int, int, int, int);
    double GetStartTime() {return m_StartTime;}
    double GetTofAlignedValue() {return m_TofAlignedValue;}
    void SetStartTime(double val) {m_StartTime = val;}
    void SetTofAlignedValue(double val) {m_TofAlignedValue = val;}

  //////////////////////////////////////////////////////////////
  // specific methods to SofTofW array
  public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();

    // clear the pre-treated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TSofTofWData object to TSofTofWPhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TSofTofWData* rawDataPointer) {m_EventData = rawDataPointer;}
    
  // objects are not written in the TTree
  private:
    TSofTofWData*         m_EventData;        //!
    TSofTofWData*         m_PreTreatedData;   //!
    TSofTofWPhysics*      m_EventPhysics;     //!

  // getters for raw and pre-treated data object
  public:
    TSofTofWData* GetRawData()        const {return m_EventData;}
    TSofTofWData* GetPreTreatedData() const {return m_PreTreatedData;}

  // parameters used in the analysis
  private:
    double m_StartTime; //!
    double m_TofAlignedValue; //!
    int m_NumberOfPlastics; //!
    double m_E_RAW_Threshold; //!
    double m_E_Threshold;     //!

    TRandom3 rand; //!
  // number of detectors
  private:
    int m_NumberOfDetectors;  //!

  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TSofTofWPhysics,1)  // SofTofWPhysics structure
};
#endif
