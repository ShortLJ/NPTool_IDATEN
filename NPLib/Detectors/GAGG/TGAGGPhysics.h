#ifndef TGAGGPHYSICS_H
#define TGAGGPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Valerian Alcindor  contact address: *
 *                                                                           *
 * Creation Date  : October 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold GAGG Treated data                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <map>
#include <string>
#include <vector>
using namespace std;

// ROOT headers
#include "TH1.h"
#include "TObject.h"
#include "TVector3.h"
// NPTool headers
#include "NPCalibrationManager.h"
#include "NPInputParser.h"
#include "NPVDetector.h"
#include "TGAGGData.h"
#include "TGAGGSpectra.h"
// forward declaration
class TGAGGSpectra;

class TGAGGPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
public:
  TGAGGPhysics();
  ~TGAGGPhysics(){};

  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
public:
  void Clear();
  void Clear(const Option_t*){};

  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
public:
  vector<int>    DetectorNumber;
  vector<double> Energy;
  vector<double> Time;
  vector<TVector3> m_Position;


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
  void BuildOnlinePhysicalEvent() { BuildPhysicalEvent(); };

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
  void ClearEventPhysics() { Clear(); }
  void ClearEventData() { m_EventData->Clear(); }

  // methods related to the TGAGGSpectra class
  // instantiate the TGAGGSpectra class and
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
  // specific methods to GAGG array
public:
  // remove bad channels, calibrate the data and apply thresholds
  void PreTreat();

  // clear the pre-treated object
  void ClearPreTreatedData() { m_PreTreatedData->Clear(); }

  // read the user configuration file. If no file is found, load standard one
  void ReadAnalysisConfig();

  // give and external TGAGGData object to TGAGGPhysics.
  // needed for online analysis for example
  void SetRawDataPointer(TGAGGData* rawDataPointer) {
    m_EventData = rawDataPointer;
  }

  // objects are not written in the TTree
private:
  TGAGGData*    m_EventData; //!
  TGAGGData*    m_PreTreatedData; //!
  TGAGGPhysics* m_EventPhysics; //!

  // getters for raw and pre-treated data object
public:
  TGAGGData* GetRawData() const { return m_EventData; }
  TGAGGData* GetPreTreatedData() const { return m_PreTreatedData; }
  TVector3   GetPositionOfInteraction(int& i);

  // parameters used in the analysis
private:
  // thresholds
  double m_E_RAW_Threshold; //!
  double m_E_Threshold; //!

  // number of detectors
private:
  int m_NumberOfDetectors; //!

  // spectra class
private:
  TGAGGSpectra* m_Spectra; // !

  // spectra getter
public:
  map<string, TH1*> GetSpectra();

  // Static constructor to be passed to the Detector Factory
public:
  static NPL::VDetector* Construct();

  ClassDef(TGAGGPhysics, 1) // GAGGPhysics structure
};
#endif
