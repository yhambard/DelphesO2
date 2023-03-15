/// @author: Yuri Kharlov
/// @author: Nicolo' Jacazio <nicolo.jacazio@cern.ch>
/// @since 04/09/2021

#include "ECALdetector.hh"
#include "TDatabasePDG.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include <iostream>

namespace o2
{
namespace delphes
{

/*****************************************************************/

void ECALdetector::setup(float resoEA, float resoEB, float resoEC, float resoXA, float resoXB)
{
  mEnergyResolutionA = resoEA;
  mEnergyResolutionB = resoEB;
  mEnergyResolutionC = resoEC;
  mPositionResolutionA = resoXA;
  mPositionResolutionB = resoXB;
}

/*****************************************************************/

bool ECALdetector::hasECAL(const Track& track) const
{
  auto x = track.XOuter * 0.1; // [cm]
  auto y = track.YOuter * 0.1; // [cm]
  auto z = track.ZOuter * 0.1; // [cm]
  /** check if hit **/
  bool ishit = false;
  auto r = hypot(x, y);
  ishit = (fabs(r - mRadius) < 0.001 && fabs(z) < mLength);
  if (!ishit)
    return false;
  auto particle = (GenParticle*)track.Particle.GetObject();
  return true;
}

/*****************************************************************/
bool ECALdetector::makeSignal(const GenParticle& particle,
                              TLorentzVector& p4ECAL,
                              float& posZ,
                              float& posPhi)
{
  // Simulate fast response of ECAL to photons:
  // take generated particle as input and calculate its smeared 4-momentum p4ECAL
  // and hit coordinates posZ, posPhi

  TLorentzVector p4True = particle.P4(); // true 4-momentum
  if (TMath::Abs(p4True.Eta()) > 1.59) {    // ECAL acceptance is rougly limited by |eta|<4
    return false;
  }
  posPhi = p4True.Phi(); // azimuth angle of a hit
  posZ = -1e6;
  Double_t tanTheta = TMath::Tan(p4True.Theta());
  if (tanTheta != 0.) {
    posZ = mRadius / tanTheta; // z-coodrinate of a hit
  }
  if(abs(posZ) <= 64) setup( mEnergyHighResolutionA, mEnergyHighResolutionB, mEnergyResolutionC, mPositionResolutionA, mPositionResolutionB);
  else setup( mEnergyResolutionA, mEnergyResolutionB, mEnergyResolutionC, mPositionResolutionA, mPositionResolutionB);

  p4ECAL = smearPhotonP4(p4True, posZ, posPhi);
  return true;
}

bool ECALdetector::makeChargedSignal(const Track& track,
                                    TLorentzVector& p4ECAL,
                                    float& posZ,
                                    float& posPhi)
{
  const int pid = track.PID;
  if (abs(pid) != 11 && abs(pid) != 13 && abs(pid) != 211 && abs(pid) != 321 && abs(pid) != 2212 ) { // e+-. MIPs will be added later.   /// added mips pions +-211 kaons 321 protons 2212 
    return false;
  }
  if (TMath::Abs(track.EtaOuter) > 1.59) {    // ECAL acceptance is rougly limited by |eta|<4
    return false;
  }
  
  if(abs(track.ZOuter * 0.1) <= 64) setup( mEnergyHighResolutionA, mEnergyHighResolutionB, mEnergyResolutionC, mPositionResolutionA, mPositionResolutionB);
  else setup( mEnergyResolutionA, mEnergyResolutionB, mEnergyResolutionC, mPositionResolutionA, mPositionResolutionB);
  
  //TLorentzVector p4Tracker = track.P4();
  TLorentzVector p4Tracker;
  p4Tracker.SetPtEtaPhiM(track.PT, track.EtaOuter, track.PhiOuter, 0.0);   // OUTERS USED
  
  posZ = track.ZOuter * 0.1;
  posPhi = track.PhiOuter;
  if(abs(pid) == 11) 
    p4ECAL = smearPhotonP4(p4Tracker, posZ, posPhi);
  else 
    p4ECAL = smearMIPP4(p4Tracker, posZ, posPhi);
  return true;
}

/*****************************************************************/
TLorentzVector ECALdetector::smearPhotonP4(const TLorentzVector& pTrue, 
                                          float& Z,
                                          float& phi)
{
  // This function smears the photon 4-momentum from the true one via applying
  // parametrized energy and coordinate resolution

  // Get true energy from true 4-momentum and smear this energy
  Double_t eTrue = pTrue.E();
  Double_t eSmeared = smearPhotonE(eTrue);
  // Smear direction of 3-vector
  phi += gRandom->Gaus(0., sigmaX(eTrue) / mRadius);
  Z += gRandom->Gaus(0., sigmaX(eTrue));
  Double_t theta = TMath::Pi()/2;
  if (Z!=0.)
  {
    theta = TMath::ATan(mRadius/Z); 
    if(theta<0) theta+= TMath::Pi();
  }
  // Calculate smeared components of 3-vector
  Double_t pxSmeared = eSmeared * TMath::Cos(phi) * TMath::Sin(theta);
  Double_t pySmeared = eSmeared * TMath::Sin(phi) * TMath::Sin(theta);
  Double_t pzSmeared = eSmeared * TMath::Cos(theta);
  std::cout<<"Pz = "<<pzSmeared<<"\t theta = "<<theta<<std::endl;
  // Construct new 4-momentum from smeared energy and 3-momentum
  TLorentzVector pSmeared(pxSmeared, pySmeared, pzSmeared, eSmeared);
  return pSmeared;
}
/*****************************************************************/
Double_t ECALdetector::sigmaX(const Double_t& eTrue)
{
  // Calculate sigma of photon coordinate smearing [cm]
  // E is the photon energy
  Double_t dX = sqrt(mPositionResolutionA * mPositionResolutionA + mPositionResolutionB * mPositionResolutionB / eTrue);
  return dX;
}
/*****************************************************************/
Double_t ECALdetector::smearPhotonE(const Double_t& eTrue)
{
  // Smear a photon energy eTrue according to a Gaussian distribution with energy resolution parameters
  // sigma of Gaussian smearing is calculated from parameters A,B,C and true energy

  const Double_t sigmaE = eTrue * sqrt(mEnergyResolutionA * mEnergyResolutionA / eTrue / eTrue +
                                       mEnergyResolutionB * mEnergyResolutionB / eTrue +
                                       mEnergyResolutionC * mEnergyResolutionC);
  Double_t eSmeared = gRandom->Gaus(eTrue, sigmaE);
  if (eSmeared < 0)
    eSmeared = 0;
  return eSmeared;
}
/*****************************************************************/

TLorentzVector ECALdetector::smearMIPP4(const TLorentzVector& pTrue,  
                                        float& Z,
                                        float& phi)
{ 
  Double_t E_mpv=.28; //[GeV]
  Double_t landausigma=0.05;
  Double_t eSmeared = gRandom->Landau(E_mpv,landausigma);  
  while (eSmeared > pTrue.P())
  {
    eSmeared = gRandom->Landau(E_mpv,landausigma);
  }

  phi = phi + gRandom->Gaus(0., sigmaX(eSmeared) / mRadius);
  Z = Z + gRandom->Gaus(0., sigmaX(eSmeared));
  Double_t theta = TMath::Pi()/2;                               ///better spacial resolution if mip?
  if (Z!=0.)
  {
    theta = TMath::ATan(mRadius/Z); 
    if(theta<0) theta+= TMath::Pi();
  }
  
  //considering mass is 0;;
  Double_t pxSmeared = eSmeared * TMath::Cos(phi) * TMath::Sin(theta);
  Double_t pySmeared = eSmeared * TMath::Sin(phi) * TMath::Sin(theta);
  Double_t pzSmeared = eSmeared * TMath::Cos(theta);
  TLorentzVector pSmeared(pxSmeared, pySmeared, pzSmeared, eSmeared);  //reconstructed MIP p
  return pSmeared;
}

/*****************************************************************/

/*

*/


} // namespace delphes
} // namespace o2

