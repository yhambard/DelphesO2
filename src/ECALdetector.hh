/// @author: Roberto Preghenella
/// @email: preghenella@bo.infn.it

#ifndef _DelphesO2_ECALdetector_h_
#define _DelphesO2_ECALdetector_h_

#include "classes/DelphesClasses.h"

namespace o2
{
namespace delphes
{

class ECALdetector
{

 public:
  ECALdetector() = default;
  ~ECALdetector() = default;

  enum { kBarrelHighResolution, kBarrelCoarse, kEndcap}; // Sectors of ECAL detector

  void setup(float resoEA, float resoEB, float resoEC, float resoPosA, float resoPosB);
  bool hasECAL(const Track& track) const;
  bool makeSignal(const GenParticle& particle, TLorentzVector& pECAL, float& posZ, float& posPhi);
  bool makeChargedSignal(const Track& track, TLorentzVector& pECAL, float& posZ, float& posPhi);

 protected:
  Double_t smearPhotonE(const Double_t& eTrue);
  Double_t sigmaX(const Double_t& eTrue);
  TLorentzVector smearPhotonP4(const TLorentzVector& pTrue, float& Z, float& phi);

  TLorentzVector smearMIPP4(const TLorentzVector& pTrue, float& Z, float& phi); //// MIPS

  int kType = kBarrelCoarse;

  float mRadius = 115.; // ECAL barrel inner radius [cm]
  float mLength = 270.; // ECAL half-length along beam axis [cm]
  float mLengthHighRes = 64.; // ECAL high resolution segment half-length along beam axis [cm]

  float zEndcap = 435.; // Endcap in cm
  float rMaxEndcap = 180.; //cm
  float rMinEndcap = 16.; //cm

  float mEnergyResolutionA = 0.02;  // parameter A of energy resolution in GeV
  float mEnergyResolutionB = 0.095;   // parameter B of energy resolution in GeV^{1/2}
  float mEnergyResolutionC = 0.01;   // parameter C of energy resolution


  float mEnergyLowResolutionA = 0.02;  // parameter A of energy resolution in GeV
  float mEnergyLowResolutionB = 0.095;   // parameter B of energy resolution in GeV^{1/2}

  float mEnergyHighResolutionA = 0.002;  // parameter A of energy resolution in GeV
  float mEnergyHighResolutionB = 0.02;   // parameter B of energy resolution in GeV^{1/2}

  float mPositionResolutionA = 0.15; // parameter A of coordinate resolution in cm
  float mPositionResolutionB = 0.30; // parameter B of coordinate resolution in cm*GeV^{1/2}
};

} // namespace delphes
} // namespace o2

#endif /** _DelphesO2_ECALdetector_h_ **/

// test comment