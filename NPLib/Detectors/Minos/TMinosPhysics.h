#ifndef TMinosPHYSICS_H
#define TMinosPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Cyril Lenain  contact address: lenain@lpccaen.in2p3.fr   *
 *                                                                           *
 * Creation Date  : 2019                                                     *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Minos Treated data                                       *
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
#include "TVector3.h"


// NPTool headers
#include "TMinosData.h"
#include "MinosUtility.h"
#include "NPLinearRansac3D.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
#include "NPXmlParser.h"

using namespace NPL;
// forward declaration

class TMinosPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TMinosPhysics();
    ~TMinosPhysics() {};
  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   

    void Clear(const Option_t*) {};
    
    //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    // PAD information
    vector<unsigned short>  Ring_Pad;  
    vector<double>  X_Pad;  
    vector<double>  Y_Pad;  
    vector<double>  Z_Pad;  
    vector<double>  Q_Pad;
    vector<double>  T_Pad;
    
    // Tracks information
    vector<TVector3> Tracks_P0;
    vector<TVector3> Tracks_Dir;

    // Vertex information
    double X_Vertex,Y_Vertex,Z_Vertex,Theta_12,Delta_Vertex;

  private: // used for analysis
    double m_TimeBin;//!
    double m_ShapingTime;//!
    double m_Baseline;//!
    unsigned int m_Sampling;//!
    TVector3 m_Position;//!
    double   m_ZRotation;//!


    NPL::MinosUtility m_utility;//! // an utility to fit the pad signal
    NPL::LinearRansac3D m_ransac;//! // a linear ransac to build the 3D tracks

    double m_E_RAW_Threshold; //! adc channels
    double m_E_Threshold;     //! MeV
    std::map<unsigned short,double> m_X;//!
    std::map<unsigned short,double> m_Y;//! 
    std::map<unsigned short,double> m_Period;//!
    std::map<unsigned short,double> m_Offset;//!
    std::map<unsigned short,double> m_Gain;//!


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
    // declare list of histograms

  //////////////////////////////////////////////////////////////
  // specific methods to Minos array
  public:
    
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();


    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TMinosData object to TMinosPhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TMinosData* rawDataPointer) {m_EventData = rawDataPointer;}
     
    // Read and load the XML info 
    void ReadXML(NPL::XmlParser&);
  // objects are not written in the TTree
  private:
    TMinosData*         m_EventData;        //!
    TMinosPhysics*      m_EventPhysics;     //!
  
    // getters for raw and pre-treated data object
  public:
    TMinosData* GetRawData()        const {return m_EventData;}//!
    vector<double> GetPad_X()  {return X_Pad;} //!
    vector<double> GetPad_Y()  {return Y_Pad;} //!
    vector<double> GetPad_Z()  {return Z_Pad;} //!
    vector<double> GetPad_T()  {return T_Pad;} //!
    double GetVertexX()  {return X_Vertex;} //!
    double GetVertexY()  {return Y_Vertex;} //!
    double GetVertexZ()  {return Z_Vertex;} //!
    
    double GetDeltaVertex()  {return Delta_Vertex;} //!
    double GetTheta12()  {return Theta_12;} //!
  
    int GetNbrOfTracks(){return Tracks_P0.size();}
      
    TVector3 GetTracksP0(unsigned int i){return Tracks_P0[i];}//!
    TVector3 GetTracksDir(unsigned int i){return Tracks_Dir[i];}//!
    double Angle(unsigned int i, unsigned j){return Tracks_Dir[i].Angle(Tracks_Dir[j]);};
  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TMinosPhysics,2)  // MinosPhysics structure
};

#endif
