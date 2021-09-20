#ifndef TStrassePHYSICS_H
#define TStrassePHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: F. Flavigny    contact : flavigny@lpccaen.in2p3.fr       *
 *                                                                           *
 * Creation Date  : July 2020                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Strasse Treated data                                     *
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
#include "TStrasseData.h"
#include "TStrasseSpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
// forward declaration
class TStrasseSpectra;



class TStrassePhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TStrassePhysics();
    ~TStrassePhysics() {};


  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};

  public:
    vector<TVector2> MatchInner();
    vector<TVector2> MatchOuter();
    int CheckEvent();

  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    Int_t EventMultiplicity;
    vector<int>     DetectorNumber;
    vector<double>  E;
    vector<double>  DE;
    vector<int>     InnerStripT;
    vector<int>     InnerStripL;
    vector<int>     OuterStripT;
    vector<int>     OuterStripL;


    vector<double> InnerPosX;
    vector<double> InnerPosY;
    vector<double> InnerPosZ;
    vector<double> OuterPosX;
    vector<double> OuterPosY;
    vector<double> OuterPosZ;


  
  //////////////////////////////////////////////////////////////
  // methods inherited from the VDetector ABC class
  public:
    // read stream from ConfigFile to pick-up detector parameters
    void ReadConfiguration(NPL::InputParser);

    /// A usefull method to bundle all operation to add a detector
    void AddInnerDetector(double R, double Z, double Phi,double Shift,TVector3 Ref); 
    void AddOuterDetector(double R, double Z, double Phi,double Shift,TVector3 Ref); 
 
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

    // methods related to the TStrasseSpectra class
    // instantiate the TStrasseSpectra class and 
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
  // specific methods to Strasse array
  public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();

    // clear the pre-treated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TStrasseData object to TStrassePhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TStrasseData* rawDataPointer) {m_EventData = rawDataPointer;}
 
    double GetNumberOfInnerDetector() const {return m_NumberOfInnerDetectors;}
    double GetNumberOfOuterDetector() const {return m_NumberOfOuterDetectors;}
    int GetEventMultiplicity() const {return EventMultiplicity;}

    double GetInnerStripPositionX(const int N, const int X, const int Y){
      return m_InnerStripPositionX[N-1][X-1][Y-1];
    };
    double GetInnerStripPositionY(const int N, const int X, const int Y){
      return m_InnerStripPositionY[N-1][X-1][Y-1];
    };
    double GetInnerStripPositionZ(const int N, const int X, const int Y){
      return m_InnerStripPositionZ[N-1][X-1][Y-1];
    };
    
    double GetOuterStripPositionX(const int N, const int X, const int Y){
      return m_OuterStripPositionX[N-1][X-1][Y-1];
    };
    double GetOuterStripPositionY(const int N, const int X, const int Y){
      return m_OuterStripPositionY[N-1][X-1][Y-1];
    };
    double GetOuterStripPositionZ(const int N, const int X, const int Y){
      return m_OuterStripPositionZ[N-1][X-1][Y-1];
    };

    TVector3 GetInnerPositionOfInteraction(const int i);
    TVector3 GetOuterPositionOfInteraction(const int i);

    TVector3 GetDetectorNormal(const int i);
    
  // objects are not written in the TTree
  private:
    TStrasseData*         m_EventData;        //!
    TStrasseData*         m_PreTreatedData;   //!
    TStrassePhysics*      m_EventPhysics;     //!

  // getters for raw and pre-treated data object
  public:
    TStrasseData* GetRawData()        const {return m_EventData;}
    TStrasseData* GetPreTreatedData() const {return m_PreTreatedData;}

  // parameters used in the analysis
  private:
    int m_NumberOfInnerDetectors; //!
    int m_NumberOfOuterDetectors; //!

    vector<vector<vector<double>>> m_InnerStripPositionX; //!
    vector<vector<vector<double>>> m_InnerStripPositionY; //!
    vector<vector<vector<double>>> m_InnerStripPositionZ; //!

    vector<vector<vector<double>>> m_OuterStripPositionX; //!
    vector<vector<vector<double>>> m_OuterStripPositionY; //!
    vector<vector<vector<double>>> m_OuterStripPositionZ; //!


    // thresholds
    double m_E_RAW_Threshold; //!
    double m_E_Threshold;     //!

 private:
    unsigned int m_MaximumStripMultiplicityAllowed;//
    double m_StripEnergyMatching;//


  // spectra class
  private:
    TStrasseSpectra* m_Spectra; // !

  // spectra getter
  public:
    map<string, TH1*>   GetSpectra(); 

  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TStrassePhysics,1)  // StrassePhysics structure
  
  private: // geometry
  ////////////////////
  // Inner Detector //
  ////////////////////
  // Wafer parameter
  double Inner_Wafer_Length;
  double Inner_Wafer_Width;
  double Inner_Wafer_Thickness;
  double Inner_Wafer_AlThickness;
  double Inner_Wafer_PADExternal;
  double Inner_Wafer_PADInternal;
  double Inner_Wafer_GuardRing;

  // PCB parameter
  double Inner_PCB_PortWidth;
  double Inner_PCB_StarboardWidth;
  double Inner_PCB_BevelAngle;
  double Inner_PCB_UpstreamWidth;
  double Inner_PCB_DownstreamWidth;
  double Inner_PCB_MidWidth;
  double Inner_PCB_Thickness;
  double Inner_PCB_Ledge;
  double Inner_PCB_Step;
  double Inner_Wafer_TransverseStrips;
  double Inner_Wafer_LongitudinalStrips;

  ////////////////////
  // Outer Detector //
  ////////////////////
  // Wafer parameter
  double Outer_Wafer_Length;
  double Outer_Wafer_Width;
  double Outer_Wafer_Thickness;
  double Outer_Wafer_AlThickness;
  double Outer_Wafer_PADExternal;
  double Outer_Wafer_PADInternal;
  double Outer_Wafer_GuardRing;

  // PCB parameter
  double Outer_PCB_PortWidth;
  double Outer_PCB_StarboardWidth;
  double Outer_PCB_BevelAngle;
  double Outer_PCB_UpstreamWidth;
  double Outer_PCB_DownstreamWidth;
  double Outer_PCB_MidWidth;
  double Outer_PCB_Thickness;
  double Outer_PCB_Ledge;
  double Outer_PCB_Step;
  double Outer_Wafer_TransverseStrips;
  double Outer_Wafer_LongitudinalStrips;


};
#endif
