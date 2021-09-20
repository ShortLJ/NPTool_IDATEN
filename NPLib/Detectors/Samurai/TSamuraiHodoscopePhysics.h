#ifndef TSamuraiHodoscopePHYSICS_H
#define TSamuraiHodoscopePHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : April 2021                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SamuraiHodoscope Treated data                            *
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
#include "TSamuraiHodoscopeData.h"
//#include "TSamuraiHodoscopeSpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
#include "NPXmlParser.h"

// forward declaration
//class TSamuraiHodoscopeSpectra;

class TSamuraiHodoscopePhysics : public TObject, public NPL::VDetector {
 //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TSamuraiHodoscopePhysics();
    ~TSamuraiHodoscopePhysics() {};

  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};

  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    vector<int>      ID;
    vector<double>   Charge;
    vector<double>   Time;
 
  private: // XML reading of calibration and position
    void ReadXML(NPL::XmlParser&);
    // calibration parameter original ID vs coeff
    map<unsigned int , double > m_qup_gain;//!
    map<unsigned int , double > m_qup_offset;//!
    map<unsigned int , double > m_tup_gain;//!
    map<unsigned int , double > m_tup_offset;//!
    map<unsigned int , double > m_qdw_gain;//!
    map<unsigned int , double > m_qdw_offset;//!
    map<unsigned int , double > m_tdw_gain;//!
    map<unsigned int , double > m_tdw_offset;//!

    // old ID to new ID map
    map<unsigned int , unsigned int > m_ID;//!

  private: // analysis parameter
    double rawQ_low_threshold;//!
    double rawQ_high_threshold;//!
    double Q_low_threshold;//!
    double Q_high_threshold;//! 
    double rawT_low_threshold;//!
    double rawT_high_threshold;//!
    double T_low_threshold;//!
    double T_high_threshold;//! 
    TRandom3 rand;


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

    // methods related to the TSamuraiHodoscopeSpectra class
    // instantiate the TSamuraiHodoscopeSpectra class and 
    // declare list of histograms
    void InitSpectra();

    // fill the spectra
    void FillSpectra();

    // used for Online mainly, sanity check for histograms and 
    // change their color if issues are found, for example
    void CheckSpectra();

    // used for Online only, clear all the spectra
    void ClearSpectra();

    // write spectra to ROOT output file
    void WriteSpectra();


  //////////////////////////////////////////////////////////////
  // specific methods to SamuraiHodoscope array
  public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();
    
    // clear the pre-treated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}
    
  
    vector <double> GetCharge() {return Charge;}
    vector <double> GetTime() {return Time;}
    vector<int> GetID() {return ID;}

    // give and external TSamuraiHodoscopeData object to TSamuraiHodoscopePhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TSamuraiHodoscopeData* rawDataPointer) {m_EventData = rawDataPointer;}
    
  // objects are not written in the TTree
  private:
    TSamuraiHodoscopeData*         m_EventData;        //!
    TSamuraiHodoscopeData*         m_PreTreatedData;   //!
    TSamuraiHodoscopePhysics*      m_EventPhysics;     //!

  // getters for raw and pre-treated data object
  public:
    TSamuraiHodoscopeData* GetRawData()        const {return m_EventData;}
    TSamuraiHodoscopeData* GetPreTreatedData() const {return m_PreTreatedData;}

  // spectra class
  private:
   // TSamuraiHodoscopeSpectra* m_Spectra; // !

  // spectra getter
  public:
    map<string, TH1*>   GetSpectra(); 

  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TSamuraiHodoscopePhysics,1)  // SamuraiHodoscopePhysics structure
};
#endif
