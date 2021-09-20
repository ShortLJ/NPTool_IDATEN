#ifndef TCOMPTONTELESCOPEPHYSICS_H
#define TCOMPTONTELESCOPEPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : February 2015                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ComptonTelescope treated data                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *  
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// STL
#include <vector>
using namespace std;

// NPL
#include "TComptonTelescopeData.h"
#include "TComptonTelescopeSpectra.h"

#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"

// ROOT 
#include "TVector2.h" 
#include "TVector3.h" 
#include "TObject.h"

// forward declaration
class TComptonTelescopeSpectra;

class TComptonTelescopePhysics : public TObject, public NPL::VDetector
{
   public:
      TComptonTelescopePhysics();
      ~TComptonTelescopePhysics() {};

   public: 
      void Clear();   
      void Clear(const Option_t*) {};


   public:  // data obtained after BuildPhysicalEvent() and stored in ROOT output file
      //   DSSD
      int EventMultiplicity;
      vector<int> EventType;
      vector<int> TowerNumber;//
      vector<int> DetectorNumber;
      vector<int>    Strip_Front;
      vector<int>    Strip_Back;
      vector<double> Strip_E;//
      vector<double> Strip_T;
      vector<double> Front_Energy;
      vector<double> Back_Energy;
      vector<double> Half_Energy;//
      vector<double> Front_Time;
      vector<double> Back_Time;
      vector<bool> Same_FBTime;//
      // Calorimeter
      vector<double> Calor_E;
      vector<double> Calor_T;
      vector<int>    CalorPosX;
      vector<int>    CalorPosY;
      multimap<int, vector<int>> CalorData;
      vector<double> Calor_E_calib;
      // Both
      vector<int> deltaT;
      long int resetCount;
   
 
   public:  // inherited from VDetector
      // Read stream at ConfigFile to pick-up parameters of detector (Position,...) using Token
      void ReadConfiguration(NPL::InputParser);

      // Add parameters to the CalibrationManger
      void AddParameterToCalibrationManager();

      // Activate associated branches and link them to the private member object m_EventData
      void InitializeRootInputRaw();

      // Activate associated branches and link them to the private member m_EventPhysics
      void InitializeRootInputPhysics();

      // Create associated branches and associated private member m_EventPhysics 
      void InitializeRootOutput();

      // This method is called at each event read from the Input Tree. Aime is to build treat Raw dat in order to extract physical parameter. 
      void BuildPhysicalEvent();

      // Same as above, but only the simplest event and/or simple method are used (low multiplicity, faster algorythm but less efficient ...).
      // This method aimed to be used for analysis performed during experiment, when speed is requiered.
      // NB: This method can eventually be the same as BuildPhysicalEvent.
      void BuildSimplePhysicalEvent();

      // Same as above but for online analysis
      void BuildOnlinePhysicalEvent()  {BuildPhysicalEvent();};

      // Clear raw and physics data 
      void ClearEventPhysics()   {Clear();}
      void ClearEventData()      {m_EventData->Clear();}

      // Methods related to the TW1Spectra classes
      // Instantiate the TW1Spectra class and the histograms
      void InitSpectra();
      // Fill the spectra defined in TW1Spectra
      void FillSpectra();
      // Used for Online mainly, perform check on the histo and for example change their color if issues are found
      void CheckSpectra();
      // Used for Online only, clear all the spectra hold by the Spectra class
      void ClearSpectra();
      // Write Spectra to file
      void WriteSpectra();
   

   public:      //   Specific to ComptonTelescope Array
      // Remove bad channel, calibrate the data and apply threshold
      void PreTreat();
   
      // Clear The PreTeated object
      void ClearPreTreatedData()   {m_PreTreatedData->Clear();}
   
      // Return false if the channel is disabled by user
      // Frist argument is either "FRONT" or "BACK"
      bool IsValidChannel(const string DetectorType, const int telescope, const int channel);
   
      // Initialize the standard parameter for analysis
      // ie: all channel enable, maximum multiplicity for strip = number of telescope
      void InitializeStandardParameter();
      
      // Read the user configuration file; if no file found, load standard one
      void ReadAnalysisConfig();
         
      // Add a Detector
//      void AddDetectorDSSSD(TVector3 C_X0_Y0, TVector3 C_X31_Y0, TVector3 C_X0_Y31, TVector3 C_X31_Y31, double size_dsssd, int nb_strip);
      void AddComptonTelescope(TVector3 C_X0_Y0, TVector3 C_X31_Y0, TVector3 C_X0_Y31, TVector3 C_X31_Y31, double size_dsssd, int nb_strip);
      void AddComptonTelescope(double Z);
      
      // Give and external TComptonTelescopeData object to TComptonTelescopePhysics
      // Needed for online analysis for example
      void SetRawDataPointer(TComptonTelescopeData* rawDataPointer) {m_EventData = rawDataPointer;}

      // Retrieve raw and pre-treated data
      TComptonTelescopeData* GetRawData()        const {return m_EventData;}
      TComptonTelescopeData* GetPreTreatedData() const {return m_PreTreatedData;}

      // Use to access the strip position
      double GetStripPositionX(const int N, const int Front, const int Back)  const {return m_StripPositionX[N-1][Front][Back];};
      double GetStripPositionY(const int N, const int Front, const int Back)  const {return m_StripPositionY[N-1][Front][Back];};
      double GetStripPositionZ(const int N, const int Front, const int Back)  const {return m_StripPositionZ[N-1][Front][Back];};

      double GetNumberOfDetectors() const {return m_NumberOfDetectors;};

      // To be called after a build Physical Event 
      int GetEventMultiplicity() const {return EventMultiplicity;};
      int GetTowerNumber(const int i) {return TowerNumber[i];};
      int GetDetectorNumber(const int i) {return DetectorNumber[i];};
      int GetFrontStrip(const int i) {return Strip_Front[i];};
      int GetBackStrip(const int i) {return Strip_Back[i];};
      double GetFrontEnergy(const int i) {return Front_Energy[i];};
      double GetBackEnergy(const int i) {return Back_Energy[i];};
      double GetHalfEnergy(const int i) {return Half_Energy[i];};
      double GetFrontTime(const int i) {return Front_Time[i];};
      double GetBackTime(const int i) {return Back_Time[i];};

      double GetCalorEnergy(const int i) {return Calor_E_calib[i];};
      
      TVector3 GetPositionOfInteractionDSSSD(const int i) const;   
      TVector3 GetDetectorNormal(const int i) const;

   private:   // Parameter used in the analysis
      // By default take EX and TY.
      bool m_Take_E_Front;      //!
      const int m_NPixels;            //!

      // If multiplicity is greater than m_MaximumStripMultiplicityAllowed 
      // after PreTreat(), event is not treated
      unsigned int m_MaximumStripMultiplicityAllowed; //!
      // Give the allowance in percent of the difference in energy between X and Y
      double m_StripEnergyMatchingSigma;              //!
      double m_StripEnergyMatchingNumberOfSigma;      //!
      bool m_MultOneOnly; //!

      //  Threshold
      double m_StripFront_E_RAW_Threshold;   //!
      double m_StripFront_E_Threshold;       //!
      double m_StripBack_E_RAW_Threshold;    //!
      double m_StripBack_E_Threshold;        //!

   public:  // methods used in event treatment 
      vector<TVector2> Match_Front_Back();
      Int_t CheckEvent();
   

   private: // objects used in event treatment 
      TComptonTelescopeData*         m_EventData;        //!
      TComptonTelescopeData*         m_PreTreatedData;   //!
      TComptonTelescopePhysics*      m_EventPhysics;     //!


   private:   // Map of activated channel
      map< int, vector<bool> > m_FrontChannelStatus;     //!
      map< int, vector<bool> > m_BackChannelStatus;      //! 

   private:   // Spatial Position of Strip Calculated on bases of detector position
      int m_NumberOfDetectors;   //!
      int m_NumberOfStrips;      //!
      vector< vector < vector < double > > >   m_StripPositionX;  //!
      vector< vector < vector < double > > >   m_StripPositionY;  //!
      vector< vector < vector < double > > >   m_StripPositionZ;  //!

   private: // Spectra Class   
      TComptonTelescopeSpectra*      m_Spectra; //! 

   public: // Spectra Getter
      map<string, TH1*> GetSpectra();


   public: // Static constructor to be passed to the Detector Factory
     static NPL::VDetector* Construct();
     ClassDef(TComptonTelescopePhysics,1)  // ComptonTelescopePhysics structure


    public: // Add data for analysis
      // counters
      Int_t m_nCounterEvt; //!
      Int_t m_nCounterHit; //!
      Int_t m_CounterEvt[50]; //!
      Int_t m_CounterHit[50]; //!
      // physical events
      //vector<int> TowerNumber;
      //vector<double> Half_Energy;
      //vector<bool> Same_FBTime; 

};


namespace ComptonTelescope_LOCAL
{
   //   DSSD
   //   Front
   double fStrip_Front_E(const TComptonTelescopeData* Data, const int i);
   double fStrip_Front_T(const TComptonTelescopeData* Data, const int i);

   //   Back   
   double fStrip_Back_E(const TComptonTelescopeData* Data, const int i);
   double fStrip_Back_T(const TComptonTelescopeData* Data, const int i);

   // Calorimeter
   double fCalorimeter_E(double charge, int detectorNumber);
   double fCalorimeter_Q(const TComptonTelescopeData* Data, const int i);
   double fCalorimeter_ped(const TComptonTelescopeData* Data, const int i);
   double fCalorimeter_coeffAlign(const TComptonTelescopeData* Data, const int i);
   double fCalorimeter_calibE(double sumADC, int detectorNumber);
}


#endif
