//----------! Written by C. Santamaria / CEA Saclay !----------
#ifndef TMINOSCLUST_H
#define TMINOSCLUST_H

#include "TObject.h"

class TMinosClust : public TObject {
 public:
  TMinosClust();
  TMinosClust(Double_t xmm, Double_t ymm, Double_t tns, Double_t zmm , Double_t chargemax, Double_t ncluster, Double_t npads, Double_t chi2);
  void Set(Double_t xmm, Double_t ymm, Double_t tns, Double_t zmm , Double_t chargemax, Double_t ncluster, Double_t npads, Double_t chi2);
  void SetZ(Double_t zmm);
  virtual ~TMinosClust();
 public:
  Double_t x_mm;
  Double_t y_mm;
  Double_t t_ns;
  Double_t z_mm;
  Double_t Chargemax;
  Double_t n_Cluster;
  Double_t n_Pads;
  Double_t Chi2;

  ClassDef(TMinosClust,1);
};
#endif // end of #ifndef TMINOSCLUST_H
