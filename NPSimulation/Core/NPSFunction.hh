#ifndef NPSFunction_h
#define NPSFunction_h 1
#include"TVector3.h"
#include"G4ThreeVector.hh"

namespace NPS{
////////////////////////////////////////////////////////////////////////////////
inline G4ThreeVector ConvertVector(TVector3 vec){
  return G4ThreeVector(vec.X(),vec.Y(),vec.Z());
  }
////////////////////////////////////////////////////////////////////////////////
inline TVector3 ConvertVector(G4ThreeVector vec){
  return TVector3(vec.x(),vec.y(),vec.z());
  }
}
#endif 
