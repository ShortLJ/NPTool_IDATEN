#ifndef TBIGRIPSICPHYSICS_H
#define TBIGRIPSICPHYSICS_H
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
 *  This class hold RIBF IC treated data                                *
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
#include "TBigRIPSICData.h"
//#include "BigRIPSICVariables.h"
//#include "TBigRIPSICSpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
#include "NPXmlParser.h"
#include "NPDCReconstruction.h"
// ROOT 
#include "TVector3.h" 
// Forward declaration
//class TBigRIPSICSpectra;



using namespace std ;

class TBigRIPSICPhysics : public TObject, public NPL::VDetector{
  public:
    TBigRIPSICPhysics();
    ~TBigRIPSICPhysics() {};

  public: 
    void Clear();   
    void Clear(const Option_t*) {};
    void Print();   

  public:
    //vectors with a size equals to the number of IC in dataset
    std::vector<int> ID;
    std::vector<int> FP;
    std::vector<double> RawAvSum;
    std::vector<double> RawSqSum;
    std::vector<double> CalAvSum;
    std::vector<double> CalSqSum;
    std::vector<int> NLayerFired;

    //vectors size = number of IC in dataset * number of layers
    std::vector<double> E;
    std::vector<int> E_Layer;
    std::vector<int> E_ID;


  public:

    // Projected position at given Z plan
    TVector3 ProjectedPosition(double Z);

  private: // Xml file read to add ICs and their parameters 
    void AddICs(string name, NPL::XmlParser&);//! take the XML file and fill in parameters of each IC
    map<int,int>  IDtoFP;//! Focal plane where the IC is located
    map<int,double> ch2mev_0; 
    map<int,double> ch2mev_1; 
    map<int,map<int,int>> pedestal; 
  
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
    TBigRIPSICData* GetRawData()        const {return m_EventData;}
    TBigRIPSICData* GetPreTreatedData() const {return m_PreTreatedData;}

  private:   //   Root Input and Output tree classes
    TBigRIPSICData*         m_EventData;//!
    TBigRIPSICData*         m_PreTreatedData;//!
    TBigRIPSICPhysics*      m_EventPhysics;//!


  private: // Spectra Class
   // TBigRIPSICSpectra* m_Spectra; // !

  public: // Spectra Getter
    map< string , TH1*> GetSpectra(); 

  public: // Static constructor to be passed to the Detector Factory
    static NPL::VDetector* Construct();
    ClassDef(TBigRIPSICPhysics,1)  // BigRIPSICPhysics structure
};


/*****************************************************************************
* Comment:                                                                  *
*                                                                           *  
*  Intermediate class useful to hold all variables per detector per event  *
*  Different from TICData whose variable (vectors) are independent          *
*                                                                           *
*****************************************************************************/

class BigRIPSICVariables{
  public:
   BigRIPSICVariables(){Clear();};  
   ~BigRIPSICVariables(){};  

  public:
    std::vector<double> FT;
    std::vector<double> FE;
    std::vector<double> FE_Layer;
    int FmultiHit[2];

    void Clear(){
        FE.clear();
        FT.clear();
        for(int i=0; i<2; i++) FmultiHit[i]=0;
    };

    void Print(){
        //cout << "XXXXXXXXXXXXXXXXXXXXXXXX IC Event XXXXXXXXXXXXXXXXX" << endl;
        cout << "FE_Mult = " << FE.size();
        for (UShort_t i = 0; i < FE.size(); i++){cout << "\tFE: " << FE[i] << endl;}
        cout << "FT_Mult = " << FT.size();
        for (UShort_t i = 0; i < FT.size(); i++){cout << "\tFT: " << FT[i] << endl;}
        cout << "MultHit = " <<endl;
        for (UShort_t i = 0; i <2; i++){cout << FmultiHit[i] << endl;}
    }

    bool HasEverything(){
        if(FE.size()==1 && FT.size()==1){
            return true;
        }else{return false;}
    }
  };

#endif
