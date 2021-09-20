#ifndef TBIGRIPSPPACPHYSICS_H
#define TBIGRIPSPPACPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : October 2020                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SamuraiFDC2 treated data                                 *
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
#include "TBigRIPSPPACData.h"
//#include "BigRIPSPPACVariables.h"
//#include "TBigRIPSPPACSpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
#include "NPXmlParser.h"
#include "NPDCReconstruction.h"
// ROOT 
#include "TVector3.h" 
// Forward declaration
//class TBigRIPSPPACSpectra;

using namespace std ;

class TBigRIPSPPACPhysics : public TObject, public NPL::VDetector{
  public:
    TBigRIPSPPACPhysics();
    ~TBigRIPSPPACPhysics() {};

  public: 
    void Clear();   
    void Clear(const Option_t*) {};
    void Print();   

  public:
    std::vector<int> ID;
    std::vector<int> FP;
    std::vector<double> TX1;
    std::vector<double> TX2;
    std::vector<double> TY1;
    std::vector<double> TY2;
    std::vector<double> TA;
    std::vector<double> TSumX;
    std::vector<double> TDiffX;
    std::vector<double> X;
    std::vector<double> TSumY;
    std::vector<double> TDiffY;
    std::vector<double> Y;
    std::vector<int> multiHit;
    //map<int,vector<double>> Data ;

    int PileUp;

  public:
    // Projected position at given Z plan
    TVector3 ProjectedPosition(double Z);

  private: // Xml file read to add PPACs and their parameters 
    void AddPPACs(string name, NPL::XmlParser&);//! take the XML file and fill in parameters of each PPAC
    map<int,double> RawUpperLimit;//! Upper Value of TDC range considered for a PPAC
    map<int,double> RawLowerLimit;//! Lower Value of TDC range considered for a PPAC 
    map<int,int>  IDtoFP;//! Focal plane where the PPAC is located
    map<int,double> ch2ns_TX1;//! 
    map<int,double> ch2ns_TX2;//!
    map<int,double> ch2ns_TY1;//!
    map<int,double> ch2ns_TY2;//!
    map<int,double> ch2ns_TA;//!
    map<int,double> xns_off;//!
    map<int,double> yns_off;//!
    map<int,double> x_offset;//!
    map<int,double> xpos_offset;//!
    map<int,double> x_ns2mm;//!
    map<int,bool> ignore_txsum_cut;//!
    map<int,double> txsum_min;//!
    map<int,double> txsum_max;//!
    map<int,double> y_offset;//!
    map<int,double> ypos_offset;//!
    map<int,double> y_ns2mm;//!
    map<int,bool> ignore_tysum_cut;//!
    map<int,double> tysum_min;//!
    map<int,double> tysum_max;//!
  
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

  public:      //   Specific

    //   Clear The PreTeated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    //   Remove bad channel, calibrate the data and apply threshold
    void PreTreat();

    // Retrieve raw and pre-treated data
    TBigRIPSPPACData* GetRawData()        const {return m_EventData;}
    TBigRIPSPPACData* GetPreTreatedData() const {return m_PreTreatedData;}

  private:   //   Root Input and Output tree classes
    TBigRIPSPPACData*         m_EventData;//!
    TBigRIPSPPACData*         m_PreTreatedData;//!
    TBigRIPSPPACPhysics*      m_EventPhysics;//!


  private: // Spectra Class
   // TBigRIPSPPACSpectra* m_Spectra; // !

  public: // Spectra Getter
    map< string , TH1*> GetSpectra(); 

  public: // Static constructor to be passed to the Detector Factory
    static NPL::VDetector* Construct();
    ClassDef(TBigRIPSPPACPhysics,1)  // BigRIPSPPACPhysics structure
};

/*---------------------------------------------------------------------------*
* Comment:                                                                  *
*                                                                           *  
*  Intermediate class necessary to hold all variables per detector per event*
*  Different from TPPACData whose variable (vectors) are independent     *
*                                                                           *
*****************************************************************************/

class BigRIPSPPACVariables{
  public:
   BigRIPSPPACVariables(){Clear();};  
   ~BigRIPSPPACVariables(){};  

  public:
    std::vector<double> FTX1;
    std::vector<double> FTX2;
    std::vector<double> FTY1;
    std::vector<double> FTY2;
    std::vector<double> FTA;
    int FmultiHit[5];

    void Clear(){
        FTX1.clear();
        FTX2.clear();
        FTY1.clear();
        FTY2.clear();
        FTA.clear();
        for(int i=0; i<5; i++) FmultiHit[i]=0;
    };

    void Print(){
        //cout << "XXXXXXXXXXXXXXXXXXXXXXXX PPAC Event XXXXXXXXXXXXXXXXX" << endl;
        cout << "FTX1_Mult = " << FTX1.size();
        for (UShort_t i = 0; i < FTX1.size(); i++){cout << "\tFTX1: " << FTX1[i] << endl;}
        cout << "FTX2_Mult = " << FTX2.size();
        for (UShort_t i = 0; i < FTX2.size(); i++){cout << "\tFTX2: " << FTX2[i] << endl;}
        cout << "FTY1_Mult = " << FTY1.size();
        for (UShort_t i = 0; i < FTY1.size(); i++){cout << "\tFTY1: " << FTY1[i] << endl;}
        cout << "FTY2_Mult = " << FTY2.size();
        for (UShort_t i = 0; i < FTY2.size(); i++){cout << "\tFTY2: " << FTY2[i] << endl;}
        cout << "FTA_Mult = " << FTA.size();
        for (UShort_t i = 0; i < FTA.size(); i++){cout << "\tFTA: " << FTA[i] << endl;}
        cout << "MultHit = " <<endl;
        for (UShort_t i = 0; i <5; i++){cout << FmultiHit[i] << endl;}
    }

    bool HasTXs(){
        if(FTX1.size()==1 && FTX2.size()==1){return true;}
        else{return false;}
    }
    bool HasTYs(){
        if(FTY1.size()==1 && FTY2.size()==1){return true;}
        else{return false;}
    }
    bool HasTA(){
        if(FTA.size()==1){return true;}
        else{return false;}
    }
    bool HasEverything(){
        if(FTX1.size()==1 && FTX2.size()==1 &&
           FTY1.size()==1 && FTY2.size()==1 &&
           FTA.size()==1){
            return true;
        }else{return false;}
    }
};

#endif
