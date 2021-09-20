#ifndef SAMURAIDCIndex_H
#define SAMURAIDCIndex_H

/*****************************************************************************
 * Copyright (C) 2009-2020    this file is part of the NPTool Project        *
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
 *  This class hold SamuraiFDC0 treated data                                 *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *  
 *  little class to index each of the DC wire                                *
 *                                                                           *
 *****************************************************************************/


class SamuraiDCIndex{
  public:
   SamuraiDCIndex(){};  
   ~SamuraiDCIndex(){};  
   SamuraiDCIndex(unsigned int det, unsigned int layer, unsigned int wire){
     m_det=det;
     m_layer=layer;
     m_wire=wire;
     m_norme=Norme();
   };  
  
  private:
    unsigned int m_det;
    unsigned int m_layer;
    unsigned int m_wire;
    unsigned int m_norme;
    
  inline int Norme() const {return (m_det*1000000000+m_layer*1000000+m_wire);} ;

  bool operator<(const SamuraiDCIndex i2){
    return this->Norme()<i2.Norme();
    }
  
  friend bool operator<(const SamuraiDCIndex i1,const SamuraiDCIndex i2){
    return i1.Norme()<i2.Norme();
    }
  
  friend bool operator==(const SamuraiDCIndex i1,const SamuraiDCIndex i2){
   return i1.Norme()==i2.Norme();
   }
  };

#endif
