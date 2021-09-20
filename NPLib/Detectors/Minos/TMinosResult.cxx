//----------! Written by C. Santamaria / CEA Saclay !----------
#include "TMinosResult.h"

ClassImp(TMinosResult)

TMinosResult::TMinosResult()
: TObject(), x_mm(0), y_mm(0), z_mm(0), Chargemax(0), n_Cluster(0), n_Pads(0), z_max(0){
// Constructor
}
TMinosResult::TMinosResult(Double_t xmm, Double_t ymm, Double_t zmm , Double_t chargemax, Double_t ncluster, Double_t npads, Double_t zmax)
  : TObject()
{
  // Constructor
  x_mm = xmm;
  y_mm= ymm;
  z_mm = zmm;
  Chargemax = chargemax;
  n_Cluster=ncluster;
  n_Pads = npads;
  z_max = zmax;
}
TMinosResult::~TMinosResult()
{
// Destructor
}

void TMinosResult::Set(Double_t xmm, Double_t ymm, Double_t zmm, Double_t chargemax, Double_t ncluster, Double_t npads, Double_t zmax)
{
  x_mm = xmm;
  y_mm= ymm;
  z_mm = zmm;
  Chargemax = chargemax;
  n_Cluster=ncluster;
  n_Pads = npads;
  z_max = zmax;
}

