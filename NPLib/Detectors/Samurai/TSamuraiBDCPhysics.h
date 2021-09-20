#ifndef TSamuraiBDCPhysics_H
#define TSamuraiBDCPhysics_H
/*****************************************************************************
 * Copyright (C) 2009-2020    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : May 2021                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SamuraiBDC treated data                                  *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *  
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// STL
#include <vector>
#include <map>
#include <string>

// NPL
#include "TSamuraiBDCData.h"
#include "SamuraiDCIndex.h"
//#include "TSamuraiBDCSpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
#include "NPXmlParser.h"

#if __cplusplus > 199711L && NPMULTITHREADING 
#include "NPDCReconstructionMT.h"
#else
#include "NPDCReconstruction.h"
#endif

// ROOT 
#include "TVector3.h" 

// Forward declaration
//class TSamuraiBDCSpectra;

class DCHit{
  public:
    DCHit(unsigned int detector, unsigned int layer , unsigned int wire, double time, double tot, double drift){
      Detector=detector;
      Layer=layer;
      Wire=wire;
      Time=time;
      ToT=tot; 
      DriftLength=drift;
    }
    DCHit(){};
    ~DCHit(){};

  public:
    unsigned int Detector;
    unsigned int Layer;
    unsigned int Wire;
    double Time;
    double ToT; 
    double DriftLength;
};

class TSamuraiBDCPhysics : public TObject, public NPL::VDetector{
  public:
    TSamuraiBDCPhysics();
    ~TSamuraiBDCPhysics() {};

    friend class DCHit;
  public: 
    void Clear();   
    void Clear(const Option_t*) {};

  public:
    //  map of [bdc number]=vector of hit
    std::map<unsigned int, std::vector<DCHit>> m_DCHit;//!


    // Computed variable
    std::vector<double> PosX;
    std::vector<double> PosY;
    std::vector<double> PosZ;
    std::vector<double> ThetaX;
    std::vector<double> PhiY;
    std::vector<double> devX;
    std::vector<double> devY;
    std::vector<unsigned int> Detector;
    std::vector<TVector3> Dir;
    std::vector<int> PileUp;

  private: // offset and inversion 
    std::map<unsigned int, TVector3> m_offset;//!
    std::map<unsigned int, bool> m_invertX;//!
    std::map<unsigned int, bool> m_invertY;//!
    std::map<unsigned int, bool> m_invertD;//!


  public:
    // Projected position at given Z plan
    TVector3 ProjectedPosition(int Detector, double Z);

  private: // Charateristic of the DC 
    void AddDC(int det, NPL::XmlParser&);//! take the XML file and fill in Wire_X and Layer_Angle
    std::map<SamuraiDCIndex,double> Wire_X;//! X position of the wires
    std::map<SamuraiDCIndex,double> Wire_Z;//! Z position of the wires
    std::map<SamuraiDCIndex,double> Wire_Angle;//! Wire Angle (0 for X, 90 for Y, U and V are typically at +/-30)

  private: // Analysis
    double ToTThreshold_H;//! a ToT Low threshold to remove noise
    double ToTThreshold_L;//! a ToT High threshold to remove noise
    // since the calibration is a sigmoid there quite a few event at the edge 
    double DriftLowThreshold;//! Minimum Drift length to keep the hit 
    double DriftUpThreshold;//! Maximum Drift length to keep the hit
    double PowerThreshold;//! Maximum P2 minimisation value to keep the track   
    // Construct the 2D track and ref position at Z=0 and Z=100 based on X,Z and Radius provided

    // Object use to perform the DC reconstruction
#if __cplusplus > 199711L && NPMULTITHREADING 
    NPL::DCReconstructionMT m_reconstruction;//!
#else
    NPL::DCReconstruction m_reconstruction;//!
#endif

  public: //   Innherited from VDetector Class

    // Read stream at ConfigFile to pick-up parameters of detector (Position,...) using Token
    void ReadConfiguration(NPL::InputParser) ;


    // Add Parameter to the CalibrationManger
    void AddParameterToCalibrationManager() ;      

    // Activated associated Branches and link it to the private member DetectorData address
    // In this method mother Branches (Detector) AND daughter leaf (fDetector_parameter) have to be activated
    void InitializeRootInputRaw() ;

    // Activated associated Branches and link it to the private member DetectorPhysics address
    // In this method mother Branches (Detector) AND daughter leaf (parameter) have to be activated
    void InitializeRootInputPhysics() ;

    // Create associated branches and associated private member DetectorPhysics address
    void InitializeRootOutput() ;

    // This method is called at each event read from the Input Tree. Aime is to build treat Raw dat in order to extract physical parameter. 
    void BuildPhysicalEvent() ;

    // Same as above, but only the simplest event and/or simple method are used (low multiplicity, faster algorythm but less efficient ...).
    // This method aimed to be used for analysis performed during experiment, when speed is requiered.
    // NB: This method can eventually be the same as BuildPhysicalEvent.
    void BuildSimplePhysicalEvent() ;

    // Same as above but for online analysis
    void BuildOnlinePhysicalEvent()  {BuildPhysicalEvent();};

    // Those two method all to clear the Event Physics or Data
    void ClearEventPhysics() {Clear();}      
    void ClearEventData()    {m_EventData->Clear();}   

    // Method related to the TSpectra classes, aimed at providing a framework for online applications
    // Instantiate the Spectra class and the histogramm throught it
    void InitSpectra();
    // Fill the spectra hold by the spectra class
    void FillSpectra();
    // Used for Online mainly, perform check on the histo and for example change their color if issues are found
    void CheckSpectra();
    // Used for Online only, clear all the spectra hold by the Spectra class
    void ClearSpectra();
    // Write Spectra to file
    void WriteSpectra();

  public:      //   Specific to SamuraiBDC Array

    //   Clear The PreTeated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    //   Remove bad channel, calibrate the data and apply threshold
    void PreTreat();

    // Retrieve raw and pre-treated data
    TSamuraiBDCData* GetRawData()        const {return m_EventData;}
    TSamuraiBDCData* GetPreTreatedData() const {return m_PreTreatedData;}

    TVector3 GetPos(unsigned int det);
    double GetPosX(unsigned int i)  {return PosX[i];}
    double GetPosY(unsigned int i)  {return PosY[i];}
    double GetThetaX(unsigned int i){return ThetaX[i];}
    double GetPhiY(unsigned int i)  {return PhiY[i];}
    double GetDevX(unsigned int i)  {return devX[i];}
    double GetDevY(unsigned int i)  {return devY[i];}
    int    GetPileUp(unsigned int i){return PileUp[i];}

  private:   //   Root Input and Output tree classes
    TSamuraiBDCData*         m_EventData;//!
    TSamuraiBDCData*         m_PreTreatedData;//!
    TSamuraiBDCPhysics*      m_EventPhysics;//!

  private: // Spectra Class
    // TSamuraiBDCSpectra* m_Spectra; // !

  public: // Spectra Getter
    std::map< std::string , TH1*> GetSpectra(); 

  public: // Static constructor to be passed to the Detector Factory
    static NPL::VDetector* Construct();
    ClassDef(TSamuraiBDCPhysics,1)  // SamuraiBDCPhysics structure
};

#endif
