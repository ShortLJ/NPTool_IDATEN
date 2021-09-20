#ifndef __MUST2DATA__
#define __MUST2DATA__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : febuary 2009                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold must2 Raw data                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include <vector>
#include "TObject.h"

using namespace std ;
class TMust2Data : public TObject {
   private:
      // DSSD
      // X strips
      // Energy
      vector<UShort_t>   fMM_StripXE_DetectorNbr;
      vector<UShort_t>   fMM_StripXE_StripNbr;
      vector<Double_t>   fMM_StripXE_Energy;
      // Time
      vector<UShort_t>   fMM_StripXT_DetectorNbr;
      vector<UShort_t>   fMM_StripXT_StripNbr;
      vector<Double_t>   fMM_StripXT_Time;
      // Y strips
      // Energy
      vector<UShort_t>   fMM_StripYE_DetectorNbr;
      vector<UShort_t>   fMM_StripYE_StripNbr;
      vector<Double_t>   fMM_StripYE_Energy;
      // Time
      vector<UShort_t>   fMM_StripYT_DetectorNbr;
      vector<UShort_t>   fMM_StripYT_StripNbr;
      vector<Double_t>   fMM_StripYT_Time;

      // SiLi
      // Energy
      vector<UShort_t>   fMM_SiLiE_DetectorNbr;
      vector<UShort_t>   fMM_SiLiE_PadNbr;
      vector<Double_t>   fMM_SiLiE_Energy;
      // Time
      vector<UShort_t>   fMM_SiLiT_DetectorNbr;
      vector<UShort_t>   fMM_SiLiT_PadNbr;
      vector<Double_t>   fMM_SiLiT_Time;

      // CsI
      // Energy
      vector<UShort_t>   fMM_CsIE_DetectorNbr;
      vector<UShort_t>   fMM_CsIE_CristalNbr;
      vector<Double_t>   fMM_CsIE_Energy;
      // Time
      vector<UShort_t>   fMM_CsIT_DetectorNbr;
      vector<UShort_t>   fMM_CsIT_CristalNbr;
      vector<Double_t>   fMM_CsIT_Time;

   public:
      TMust2Data();
      virtual ~TMust2Data();

      void   Clear();
      void  Clear(const Option_t*) {};
      void   Dump() const;

      /////////////////////           SETTERS           ////////////////////////
      // DSSD
      // (X,E)
      inline void   SetStripXE(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Energy){
        fMM_StripXE_DetectorNbr.push_back(DetNbr);
        fMM_StripXE_StripNbr.push_back(StripNbr);
        fMM_StripXE_Energy.push_back(Energy);
      }
      
      // (X,T)
     inline void   SetStripXT(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Time){   
      fMM_StripXT_DetectorNbr.push_back(DetNbr);  
      fMM_StripXT_StripNbr.push_back(StripNbr);       
      fMM_StripXT_Time.push_back(Time);  
     } 
       // (Y,E)
      inline void   SetStripYE(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Energy){
        fMM_StripYE_DetectorNbr.push_back(DetNbr);
        fMM_StripYE_StripNbr.push_back(StripNbr);
        fMM_StripYE_Energy.push_back(Energy);
      }
      
      // (Y,T)
     inline void   SetStripYT(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Time){   
      fMM_StripYT_DetectorNbr.push_back(DetNbr);  
      fMM_StripYT_StripNbr.push_back(StripNbr);       
      fMM_StripYT_Time.push_back(Time);  
     } 
      
      // SiLi
      // E
      inline void   SetSiLiE(const UShort_t& DetNbr, const UShort_t& PadNbr, const Double_t& Energy){
        fMM_SiLiE_DetectorNbr.push_back(DetNbr);
        fMM_SiLiE_PadNbr.push_back(PadNbr);
        fMM_SiLiE_Energy.push_back(Energy);
      }
      
      // T
     inline void   SetSiLiT(const UShort_t& DetNbr, const UShort_t& PadNbr, const Double_t& Time){   
      fMM_SiLiT_DetectorNbr.push_back(DetNbr);  
      fMM_SiLiT_PadNbr.push_back(PadNbr);       
      fMM_SiLiT_Time.push_back(Time);  
     } 
 
      // CsI
      // E
      inline void   SetCsIE(const UShort_t& DetNbr, const UShort_t& CristalNbr, const Double_t& Energy){
        fMM_CsIE_DetectorNbr.push_back(DetNbr);
        fMM_CsIE_CristalNbr.push_back(CristalNbr);
        fMM_CsIE_Energy.push_back(Energy);
      }
      
      // T
     inline void   SetCsIT(const UShort_t& DetNbr, const UShort_t& CristalNbr, const Double_t& Time){   
      fMM_CsIT_DetectorNbr.push_back(DetNbr);  
      fMM_CsIT_CristalNbr.push_back(CristalNbr);       
      fMM_CsIT_Time.push_back(Time);  
     } 
 
      /////////////////////           GETTERS           ////////////////////////
      // DSSD
      // (X,E)
      UShort_t   GetMMStripXEMult()                      const {return fMM_StripXE_DetectorNbr.size();}
      UShort_t   GetMMStripXEDetectorNbr(const Int_t& i) const {return fMM_StripXE_DetectorNbr[i];}
      UShort_t   GetMMStripXEStripNbr(const Int_t& i)    const {return fMM_StripXE_StripNbr[i];}
      Double_t   GetMMStripXEEnergy(const Int_t& i)      const {return fMM_StripXE_Energy[i];}
      // (X,T)
      UShort_t   GetMMStripXTMult()                      const {return fMM_StripXT_DetectorNbr.size();}
      UShort_t   GetMMStripXTDetectorNbr(const Int_t& i) const {return fMM_StripXT_DetectorNbr[i];}
      UShort_t   GetMMStripXTStripNbr(const Int_t& i)    const {return fMM_StripXT_StripNbr[i];}
      Double_t   GetMMStripXTTime(const Int_t& i)        const {return fMM_StripXT_Time[i];}
      // (Y,E)
      UShort_t   GetMMStripYEMult()                      const {return fMM_StripYE_DetectorNbr.size();}
      UShort_t   GetMMStripYEDetectorNbr(const Int_t& i) const {return fMM_StripYE_DetectorNbr[i];}
      UShort_t   GetMMStripYEStripNbr(const Int_t& i)    const {return fMM_StripYE_StripNbr[i];}
      Double_t   GetMMStripYEEnergy(const Int_t& i)      const {return fMM_StripYE_Energy[i];}
      // (Y,T)
      UShort_t   GetMMStripYTMult()                      const {return fMM_StripYT_DetectorNbr.size();}
      UShort_t   GetMMStripYTDetectorNbr(const Int_t& i) const {return fMM_StripYT_DetectorNbr[i];}
      UShort_t   GetMMStripYTStripNbr(const Int_t& i)    const {return fMM_StripYT_StripNbr[i];}
      Double_t   GetMMStripYTTime(const Int_t& i)        const {return fMM_StripYT_Time[i];}

      // SiLi
      //(E)
      UShort_t   GetMMSiLiEMult()                       const {return fMM_SiLiE_DetectorNbr.size();}
      UShort_t   GetMMSiLiEDetectorNbr(const Int_t& i)  const {return fMM_SiLiE_DetectorNbr[i];}
      UShort_t   GetMMSiLiEPadNbr(const Int_t& i)       const {return fMM_SiLiE_PadNbr[i];}
      Double_t   GetMMSiLiEEnergy(const Int_t& i)       const {return fMM_SiLiE_Energy[i];}
      //(T)
      UShort_t   GetMMSiLiTMult()                       const {return fMM_SiLiT_DetectorNbr.size();}
      UShort_t   GetMMSiLiTDetectorNbr(const Int_t& i)  const {return fMM_SiLiT_DetectorNbr[i];}
      UShort_t   GetMMSiLiTPadNbr(const Int_t& i)       const {return fMM_SiLiT_PadNbr[i];}
      Double_t   GetMMSiLiTTime(const Int_t& i)         const {return fMM_SiLiT_Time[i];}

      // CsI 
      //(E)
      UShort_t   GetMMCsIEMult()                        const {return fMM_CsIE_DetectorNbr.size();}
      UShort_t   GetMMCsIEDetectorNbr(const Int_t& i)   const {return fMM_CsIE_DetectorNbr[i];}
      UShort_t   GetMMCsIECristalNbr(const Int_t& i)    const {return fMM_CsIE_CristalNbr[i];}
      Double_t   GetMMCsIEEnergy(const Int_t& i)        const {return fMM_CsIE_Energy[i];}
      //(T)
      UShort_t   GetMMCsITMult()                        const {return fMM_CsIT_DetectorNbr.size();}
      UShort_t   GetMMCsITDetectorNbr(const Int_t& i)   const {return fMM_CsIT_DetectorNbr[i];}
      UShort_t   GetMMCsITCristalNbr(const Int_t& i)    const {return fMM_CsIT_CristalNbr[i];}
      Double_t   GetMMCsITTime(const Int_t& i)          const {return fMM_CsIT_Time[i];}

      ClassDef(TMust2Data,2)  // Must2Data structure
};

#endif
