#ifndef TBIGRIPSPLASTICPHYSICS_H
#define TBIGRIPSPLASTICPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: F. Flavigny contact address: flavigny@lpccaen.in2p3.fr   *
 *                                                                           *
 * Creation Date  : April 2020                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold RIBF Plastic treated data                                *
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
#include <iostream>

// NPL
#include "TBigRIPSPlasticData.h"
//#include "BigRIPSPlasticVariables.h"
//#include "TBigRIPSPlasticSpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
#include "NPXmlParser.h"
#include "NPDCReconstruction.h"
// ROOT 
#include "TVector3.h" 
// Forward declaration
//class TBigRIPSPlasticSpectra;


using namespace std ;

class TBigRIPSPlasticPhysics : public TObject, public NPL::VDetector{
  public:
    TBigRIPSPlasticPhysics();
    ~TBigRIPSPlasticPhysics() {};

  public: 
    void Clear();   
    void Clear(const Option_t*) {};
    void Print();   

  public:
    std::vector<int> ID;
    std::vector<int> FP;
    std::vector<double> T;
    std::vector<double> TL;
    std::vector<double> TR;
    std::vector<double> TSlew;
    std::vector<double> TLSlew;
    std::vector<double> TRSlew;
    std::vector<double> QL;
    std::vector<double> QR;
    std::vector<int> multiHit;
    std::vector<bool> fired;

  public:

    // Projected position at given Z plan
    TVector3 ProjectedPosition(double Z);

  private: // Xml file read to add Plastics and their parameters 
    void AddPlastics(string name, NPL::XmlParser&);//! take the XML file and fill in parameters of each Plastic
    map<int,double> RawUpperLimit;//! Upper Value of TDC range considered for a Plastic
    map<int,double> RawLowerLimit;//! Lower Value of TDC range considered for a Plastic 
    map<int,int>  IDtoFP;//! Focal plane where the Plastic is located
    map<int,double> tcal_left; 
    map<int,double> tcal_right; 
    map<int,double> tslew_left_a; 
    map<int,double> tslew_left_b; 
    map<int,double> tslew_right_a; 
    map<int,double> tslew_right_b; 
  
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

  public:     

    //   Clear The PreTeated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    //   Remove bad channel, calibrate the data and apply threshold
    void PreTreat();

    // Retrieve raw and pre-treated data
    TBigRIPSPlasticData* GetRawData()        const {return m_EventData;}
    TBigRIPSPlasticData* GetPreTreatedData() const {return m_PreTreatedData;}

  private:   //   Root Input and Output tree classes
    TBigRIPSPlasticData*         m_EventData;//!
    TBigRIPSPlasticData*         m_PreTreatedData;//!
    TBigRIPSPlasticPhysics*      m_EventPhysics;//!


  private: // Spectra Class
   // TBigRIPSPlasticSpectra* m_Spectra; // !

  public: // Spectra Getter
    map< string , TH1*> GetSpectra(); 

  public: // Static constructor to be passed to the Detector Factory
    static NPL::VDetector* Construct();
    ClassDef(TBigRIPSPlasticPhysics,1)  // BigRIPSPlasticPhysics structure
};


/*---------------------------------------------------------------------------*
* Comment:                                                                  *
*                                                                           *  
*  Intermediate class necessary to hold all variables per detector per event*
*  Different from TPlasticData whose variable (vectors) are independent     *
*                                                                           *
*****************************************************************************/

class BigRIPSPlasticVariables{
  public:
   BigRIPSPlasticVariables(){Clear();};  
   ~BigRIPSPlasticVariables(){};  

  public:
    std::vector<double> FTL;
    std::vector<double> FTR;
    std::vector<double> FQL;
    std::vector<double> FQR;
    int FmultiHit[4];

    void Clear(){
        FTL.clear();
        FTR.clear();
        FQL.clear();
        FQR.clear();
        for(int i=0; i<4; i++) FmultiHit[i]=0;
    };

    void Print(){
        //cout << "XXXXXXXXXXXXXXXXXXXXXXXX Plastic Event XXXXXXXXXXXXXXXXX" << endl;
        cout << "FTL_Mult = " << FTL.size();
        for (UShort_t i = 0; i < FTL.size(); i++){cout << "\tFTL: " << FTL[i] << endl;}
        cout << "FTR_Mult = " << FTR.size();
        for (UShort_t i = 0; i < FTR.size(); i++){cout << "\tFTR: " << FTR[i] << endl;}
        cout << "FQL_Mult = " << FQL.size();
        for (UShort_t i = 0; i < FQL.size(); i++){cout << "\tFQL: " << FQL[i] << endl;}
        cout << "FQR_Mult = " << FQR.size();
        for (UShort_t i = 0; i < FQR.size(); i++){cout << "\tFQR: " << FQR[i] << endl;}
        cout << "MultHit = " <<endl;
        for (UShort_t i = 0; i <4; i++){cout << FmultiHit[i] << endl;}
    }

    bool HasTLandQL(){
        if(FTL.size()==1 && FQL.size()==1){return true;}
        else{return false;}
    }
    bool HasTRandQR(){
        if(FTR.size()==1 && FQR.size()==1){return true;}
        else{return false;}
    }
    bool HasTLandTR(){
        if(FTL.size()==1 && FTR.size()==1){return true;}
        else{return false;}
    }
    bool HasQLandQR(){
        if(FQL.size()==1 && FQR.size()==1){return true;}
        else{return false;}
    }
    bool HasEverything(){
        if(FTL.size()==1 && FTR.size()==1 &&
           FQL.size()==1 && FQR.size()==1){
            return true;
        }else{return false;}
    }
  };

#endif
