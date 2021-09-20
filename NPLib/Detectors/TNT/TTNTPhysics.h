#ifndef TTNTPHYSICS_H
#define TTNTPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Dustin Scriven  contact address: scrivend@tamu.edu                        *
 *                                                                           *
 * Creation Date  : July 2018                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TNT Treated data                                *
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
// NPTool headers
#include "TTNTData.h"
#include "TTNTSpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
// forward declaration
class TTNTSpectra;



class TTNTPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TTNTPhysics();
    ~TTNTPhysics() {};


  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};

  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    vector<int>      DetectorNumber;
    vector<double>   Energy;
    vector<double>   Time;    
    vector<double>   Xpos;
    vector<double>   Ypos;
    vector<double>   Zpos;
    vector<double>   Theta;
    vector<double>   Phi;
    vector<double>   Theta_Lab;
    // Template map<DetectorNumber, <Pos[0], Pos[1], Pos[2]>>
    map<int, vector<double>> DetectorInfo; //!
  
  /// A usefull method to bundle all operation to add a detector
  void AddDetector(TVector3 POS, string shape); 
  void AddDetector(double R, double Theta, double Phi, string shape); 
  
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

    // methods related to the TTNTSpectra class
    // instantiate the TTNTSpectra class and 
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
  // specific methods to TNT array
  public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();

    // clear the pre-treated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TTNTData object to TTNTPhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TTNTData* rawDataPointer) {m_EventData = rawDataPointer;}
    
  // objects are not written in the TTree
  private:
    TTNTData*         m_EventData;        //!
    TTNTData*         m_PreTreatedData;   //!
    TTNTPhysics*      m_EventPhysics;     //!

  // getters for raw and pre-treated data object
  public:
    TTNTData* GetRawData()        const {return m_EventData;}
    TTNTData* GetPreTreatedData() const {return m_PreTreatedData;}

  // parameters used in the analysis
  private:
    // thresholds
    double m_E_RAW_Threshold; //!
    double m_E_Threshold;     //!

  // number of detectors
  private:
    double m_NumberOfDetectors;  //!

  // spectra class
  private:
    TTNTSpectra* m_Spectra; // !

  // spectra getter
  public:
    map<string, TH1*>   GetSpectra(); 

  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TTNTPhysics,1)  // TNTPhysics structure
};
#endif
