#ifndef TPyramidPHYSICS_H
#define TPyramidPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Mohamad Moukaddam                                        *
 * contact address: mohamad.moukaddam@iphc.cnrs.fr                           *
 *                                                                           *
 * Creation Date  : November 2018                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Pyramid Treated data                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
// STL
#include <vector>
#include <map>
#include <string>
using namespace std;
// NPL
#include "TPyramidData.h"
#include "TPyramidSpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
// ROOT 
#include "TVector2.h" 
#include "TVector3.h" 
#include "TObject.h"
#include "TH1.h"

class TPyramidSpectra;

using namespace std ;

class TPyramidPhysics : public TObject, public NPL::VDetector{
  public:
    TPyramidPhysics();
    ~TPyramidPhysics() {};

  public: 
    void Clear();   
    void Clear(const Option_t*) {};

  public:
    //   Provide Physical Multiplicity
    Int_t EventMultiplicity;

    // Detector
    vector<int> Detector_N ;
    vector<int> Outer_Detector_N ;

    // Inner 
    vector<double> Strip_E;
    vector<double> Strip_T;
    vector<int>    Strip_N;
    vector<double> Strip_Pos;
   
    // Control stuff 
    vector<double> DownStream_E;
    vector<double> DownStream_T;
    vector<double> UpStream_E;
    vector<double> UpStream_T;
    vector<double> Back_E;
    vector<double> Back_T;

    // Outter 
    vector<double> Outer_Strip_E;
    vector<double> Outer_Strip_T;
    vector<double> Outer_Strip_N;
    vector<double> Outer_Back_E;
    vector<double> Outer_Back_T;

    double XCoord;
    double YCoord;
    double ZCoord;


  public:      //   Innherited from VDetector Class
    //   Read stream at ConfigFile to pick-up parameters of detector (Position,...) using Token
    void ReadConfiguration(NPL::InputParser) ;

    //   Add Parameter to the CalibrationManger
    void AddParameterToCalibrationManager() ;      

    //   Activated associated Branches and link it to the private member DetectorData address
    //   In this method mother Branches (Detector) AND daughter leaf (fDetector_parameter) have to be activated
    void InitializeRootInputRaw() ;

    //   Activated associated Branches and link it to the private member DetectorPhysics address
    //   In this method mother Branches (Detector) AND daughter leaf (parameter) have to be activated
    void InitializeRootInputPhysics() ;

    //   Create associated branches and associated private member DetectorPhysics address
    void InitializeRootOutput() ;

    //   This method is called at each event read from the Input Tree. Aime is to build treat Raw dat in order to extract physical parameter. 
    void BuildPhysicalEvent() ;

    //   Same as above, but only the simplest event and/or simple method are used (low multiplicity, faster algorythm but less efficient ...).
    //   This method aimed to be used for analysis performed during experiment, when speed is requiered.
    //   NB: This method can eventually be the same as BuildPhysicalEvent.
    void BuildSimplePhysicalEvent() ;

    // Same as above but for online analysis
    void BuildOnlinePhysicalEvent()  {BuildPhysicalEvent();};

    //   Those two method all to clear the Event Physics or Data
    void ClearEventPhysics() {Clear();}      
    void ClearEventData()    {m_EventData->Clear();}   

    // Method related to the TSpectra classes, aimed at providing a framework 
    // for online applications
    // Instantiate the Spectra class and the histogramm throught it
    void InitSpectra();
    // Fill the spectra hold by the spectra class
    void FillSpectra();
    // Used for Online mainly, perform check on the histo and for example change 
    // their color if issues are found
    void CheckSpectra();
    // Used for Online only, clear all the spectra hold by the Spectra class
    void ClearSpectra();
    // Write Spectra to file
    void WriteSpectra();
  public://   Specific to Pyramid Array
    //   Clear The PreTeated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear(); m_PreTreatedMSData->Clear();}

    //   Remove bad channel, calibrate the data and apply threshold
    void PreTreat();

    //   Return false if the channel is disabled by user
    //   Frist argument is either "X","Y","SiLi","CsI"
    bool IsValidChannel(const string DetectorType, const int detector , const int channel);

    //   Initialize the standard parameter for analysis
    //   ie: all channel enable, maximum multiplicity for strip = number of detector
    void InitializeStandardParameter();

    //   Read the user configuration file; if no file found, load standard one
    void ReadAnalysisConfig();

    //   Add a Detector
    void AddDetector( double X,double Y,double Z);
 
    // Give and external TMustData object to TPyramidPhysics. Needed for online analysis for example.
    void SetRawDataPointer(TPyramidData* rawDataPointer) {m_EventData = rawDataPointer;}
    // Retrieve raw and pre-treated data
    TPyramidData* GetRawData()        const {return m_EventData;}
    TPyramidData* GetPreTreatedData() const {return m_PreTreatedData;}

    double GetNumberOfDetector() const { return m_NumberOfDetector; };

    // To be called after a build Physical Event 
    int GetEventMultiplicity() const { return EventMultiplicity; };

    double GetXCoord() const { return XCoord; };
    double GetYCoord() const { return YCoord; };
    double GetZCoord() const { return ZCoord; };

    TVector3 GetPositionOfInteraction(const int i) const;   
    TVector3 GetRandomisedPositionOfInteraction(const int i) const;  
    TVector3 GetDetectorNormal(const int i) const;

  private:   //   Parameter used in the analysis
    // By default take EX and TY.
    bool m_Take_E_Strip;//!
    bool m_Take_T_Back;//!

    //  Threshold
    double m_Strip_E_Threshold ;//!
    double m_Back_E_Threshold ;//!
    double m_OuterBack_E_Threshold ;//!
    double m_Maximum_FrontBack_Difference ;//!
  private:   //   Root Input and Output tree classes
    TPyramidData*         m_EventData;//!
    TPyramidData*         m_PreTreatedData;//!
    TPyramidData*         m_PreTreatedMSData;//! stores the intermediate Matchsticks calibrated Data
    TPyramidPhysics*      m_EventPhysics;//!
    map<int, vector <double> > m_mapU;//! the maps sorts out the data before storing in m_PreTreatedData
    map<int, vector <double> > m_mapD;//! 
    map<int, vector <double> > m_mapB;//! 
    map<int, vector <double> > m_mapO;//! 
    map<int, vector <double> > m_mapMSU;//! 
    map<int, vector <double> > m_mapMSD;//! 
    
  private:   //   Map of activated channel
    map< int, vector<bool> > m_InnerStripUpstreamChannelStatus;//!
    map< int, vector<bool> > m_InnerStripDownstreamChannelStatus;//!
    map< int, vector<bool> > m_OuterStripChannelStatus;//!
    map< int, vector<bool> > m_InnerBackChannelStatus;//!
    map< int, vector<bool> > m_OuterBackChannelStatus;//!

  private:   //   Spatial Position of Strip Calculated on bases of detector position
    int m_NumberOfDetector;//!
    vector< vector<double> > m_StripPositionX;//!
    vector< vector<double> > m_StripPositionY;//!
    vector< vector<double> > m_StripPositionZ;//!

  private:
    bool m_boolChamber;
    vector<bool> m_boolInner;
    vector<bool> m_boolOuter;
    vector<double> m_Z; // shift on the z-axis
    vector<double> m_ANGLE; // angle of rotation around the downstream width 

  private: // Spectra
    TPyramidSpectra*      m_Spectra;//!

  public:
    map< string,TH1* > GetSpectra(); 

  private: // Usefull method
   // Calibrate data
  double Cal_Strip_Upstream_E(const int i);
  double Cal_Strip_Downstream_E(const int i);
  double Cal_Back_E(const int i);
  double Match_Strip_Upstream_E(const int i);
  double Match_Strip_Downstream_E(const int i);

  public: // Static constructor to be passed to the Detector Factory
     static NPL::VDetector* Construct();
     ClassDef(TPyramidPhysics,1)  // PyramidPhysics structure

};
namespace Pyramid_LOCAL{
 string itoa(unsigned int value);
}
#endif
