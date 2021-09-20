#ifndef TNebulaPHYSICS_H
#define TNebulaPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : December 2019                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Nebula Treated data                                      *
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
#include "TNebulaData.h"
#include "TNebulaSpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
#include "NPXmlParser.h"
// forward declaration
class TNebulaSpectra;



class TNebulaPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TNebulaPhysics();
    ~TNebulaPhysics() {};


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
    vector<double>   Charge;
    vector<double>   TOF;
    vector<double>   PosY;
    vector<double>   PosX;
    vector<double>   PosZ;
    vector<bool>     IsVeto;

  public:
    TVector3 GetPos(const unsigned int& i) const{
      return TVector3(PosX[i],PosY[i],PosZ[i]);
    }

    // Return true if one veto fired
    bool HasVeto(){
      unsigned int size = IsVeto.size();
      for(unsigned int i = 0 ; i < size ; i++){
        if(IsVeto[i])
          return true;
      }
      return false;
    };

    /////////// Get index of fastest neutron
    int GetFirstHit(){
      unsigned int size = TOF.size();
      unsigned int index=0;

      if(!size)
        return -1;

      double tof = TOF[0];
      for(unsigned int i = 1 ; i < size ; i++){
        if(tof<TOF[i]){
          tof=TOF[i];
          index=i;
        }
      }
      return index;
    };

  public:
    /// A usefull method to bundle all operation to add a detector
    void ReadXML(NPL::XmlParser); 

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

    // methods related to the TNebulaSpectra class
    // instantiate the TNebulaSpectra class and 
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
    // specific methods to Nebula array
  public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TNebulaData object to TNebulaPhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TNebulaData* rawDataPointer) {m_EventData = rawDataPointer;}

    // objects are not written in the TTree
  private:
    TNebulaData*         m_EventData;        //!
    TNebulaPhysics*      m_EventPhysics;     //!

    // getters for raw and pre-treated data object
  public:
    TNebulaData* GetRawData()        const {return m_EventData;}

    // parameters used in the analysis
  private:
    // thresholds
    double m_Q_RAW_Threshold; //!
    double m_Q_Threshold;     //!
    double m_V_Threshold;     //!

  public: 
    void SetQThreshold(double t) {m_Q_Threshold=t;};
    void SetVThreshold(double t) {m_V_Threshold=t;};
    // number of detectors
  private:
    int m_NumberOfBars;  //!

  private: // offset and inversion 
    std::map<unsigned int, TVector3> m_offset;//!
    std::map<unsigned int, bool> m_invertX;//!
    std::map<unsigned int, bool> m_invertY;//!

  private: // xml calibration
    // position
    std::map<unsigned int , double > PositionX;//!
    std::map<unsigned int , double > PositionY;//!
    std::map<unsigned int , double > PositionZ;//!

    // linear cal
    std::map<unsigned int , double > aQu;//!
    std::map<unsigned int , double > bQu;//!
    std::map<unsigned int , double > aQd;//!
    std::map<unsigned int , double > bQd;//!
    std::map<unsigned int , double > aTu;//!
    std::map<unsigned int , double > bTu;//!
    std::map<unsigned int , double > aTd;//!
    std::map<unsigned int , double > bTd;//!

    // T average offset
    std::map<unsigned int , double > avgT0;//!

    // slew correction T= tcal +slwT/sqrt(Qcal)
    std::map<unsigned int , double > slwTu;//!
    std::map<unsigned int , double > slwTd;//!

    // DT position cal
    std::map<unsigned int , double > DTa;//!
    std::map<unsigned int , double > DTb;//!


    // spectra class
  private:
    TNebulaSpectra* m_Spectra; // !

    // spectra getter
  public:
    map<string, TH1*>   GetSpectra(); 
  
    // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TNebulaPhysics,1)  // NebulaPhysics structure
};
#endif
