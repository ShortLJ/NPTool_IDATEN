//----------! Written by C. Santamaria / CEA Saclay !----------
#ifndef TMINOSRESULT_H
#define TMINOSRESULT_H

//#include "TClonesArray.h"
#include "TObject.h"

class TMinosResult : public TObject {
 public:
  TMinosResult();
  TMinosResult(Double_t xmm, Double_t ymm, Double_t zmm , Double_t chargemax, Double_t ncluster, Double_t npads, Double_t zmax);
  void Set(Double_t xmm, Double_t ymm, Double_t zmm , Double_t chargemax, Double_t ncluster, Double_t npads, Double_t zmax);
  virtual ~TMinosResult();
 public:
  Double_t x_mm;
  Double_t y_mm;
  Double_t z_mm;
  Double_t Chargemax;
  Double_t n_Cluster;
  Double_t n_Pads;
  Double_t z_max;

  ClassDef(TMinosResult,1);
};
#endif // end of #ifndef TMINOSRESULT_H
