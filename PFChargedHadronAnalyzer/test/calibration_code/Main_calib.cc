// #include <vector>
// #include <TROOT.h>
// #include <TChain.h>
// #include <TFile.h>
// #include <TF1.h>
// #include <TF2.h>
// #include "TH2F.h"
// #include "TLegend.h"
// #include "TProfile.h"
// #include "TProfile2D.h"
// #include "TGraph.h"
// #include "TMath.h"
// #include "TGraphErrors.h"
// #include "Math/SMatrix.h"
// #include "Math/SVector.h"
// #include "TCanvas.h"
// #include "TStyle.h"
// #include "TLine.h"
// #include <string>
// #include <iostream>
// #include <math.h>

#include "Main_calib.h"
//#include "CrystalBall.C"

using namespace std;


bool freezeparameters = true;
bool useMean = false;
bool useMedian = false;
bool changeRange =false;
bool old_logic = false;
bool drawpT = false;
bool drawResoFit = false;
bool saveCanvas = true;
bool payload=false;
bool useP_reco=false;//Instead of using etrue, the calibration code uses p_reco.
bool drawRespPlots = false;
int hadrons_eta_symbol = 0;//This int variable is 0 when we take all the eta values (positives and negatives); +1 when we only take hadrons with positive eta; -1 when we only take hadrons with negatives eta.
bool PFEnergyCalibrationFunction=false; //This bool variable is true for using PFEnergyCalibration function from CMSSW and flase for using the default PFHC calibration function getCalibratedEnergy
bool WriteNTupleFile = false;
bool Parameters24Above500GeV = false; //This bool variable is true for using the parameters for the calibration function for the energy above 500 GeV and false for using the parameters for the calibration function for the energy below 500 GeV (first one)
//char* _region_ = (char*)"EC_outside_tracker";
//char* _region_ = (char*)"EC_within_tracker";
//char* _region_ = (char*)"barrel";
char* _region_ = (char*)"Full";

float _etaMin_ = 0.0;
float _etaMax_ = 0.0;
/*************************************NUEVO BY MIKKO (Ading EnergyCalibration, cmssw) *****************/
PFEnergyCalibration::PFEnergyCalibration() {
  //calibChrisClean.C calibration parameters bhumika Nov, 2018
  faBarrel = std::make_unique<TF1>("faBarrel","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1., sampleRangeHigh);
  faBarrel->SetParameter(0,13.033);
  faBarrel->SetParameter(1,87.2668);
  faBarrel->SetParameter(2,-699.24);
  faBarrel->SetParameter(3,0.304668);
  faBarrel->SetParameter(4,13.8154);
  faBarrel->SetParameter(5,0.266523);
  faBarrel->SetParameter(6,0.0171292);
  faBarrel->SetParameter(7,-0.725741);
  fbBarrel = std::make_unique<TF1>("fbBarrel","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1., sampleRangeHigh);
  fbBarrel->SetParameter(0,1.75412);
  fbBarrel->SetParameter(1,-0.413335);
  fbBarrel->SetParameter(2,-2.08127);
  fbBarrel->SetParameter(3,126.351);
  fbBarrel->SetParameter(4,0.770695);
  fbBarrel->SetParameter(5,0.00404635);
  fbBarrel->SetParameter(6,1.12044);
  fbBarrel->SetParameter(7,-1.38901);
  fcBarrel = std::make_unique<TF1>("fcBarrel","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1., sampleRangeHigh);
  fcBarrel->SetParameter(0,11.384);
  fcBarrel->SetParameter(1,23.1406);
  fcBarrel->SetParameter(2,-28.9497);
  fcBarrel->SetParameter(3,0.97494);
  fcBarrel->SetParameter(4,19.647);
  fcBarrel->SetParameter(5,1.6272);
  fcBarrel->SetParameter(6,-0.0108519);
  fcBarrel->SetParameter(7,-0.423267);
  faEtaBarrelEH = std::make_unique<TF1>("faEtaBarrelEH","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1., sampleRangeHigh);
  faEtaBarrelEH->SetParameter(0,-0.0249299);
  faEtaBarrelEH->SetParameter(1,41.0626);
  faEtaBarrelEH->SetParameter(2,-0.756772);
  faEtaBarrelEH->SetParameter(3,9.82297e-07);
  faEtaBarrelEH->SetParameter(4,41.0056);
  faEtaBarrelEH->SetParameter(5,1.02733e-06);
  faEtaBarrelEH->SetParameter(6,-3.70842);
  faEtaBarrelEH->SetParameter(7,-3.69656);
  fbEtaBarrelEH = std::make_unique<TF1>("fbEtaBarrelEH","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1., sampleRangeHigh);
  fbEtaBarrelEH->SetParameter(0,-0.176791);
  fbEtaBarrelEH->SetParameter(1,0.60533);
  fbEtaBarrelEH->SetParameter(2,-0.891364);
  fbEtaBarrelEH->SetParameter(3,0.033817);
  fbEtaBarrelEH->SetParameter(4,0.4556);
  fbEtaBarrelEH->SetParameter(5,0.134414);
  fbEtaBarrelEH->SetParameter(6,-61.9498);
  fbEtaBarrelEH->SetParameter(7,-0.520039);
  faEtaBarrelH = std::make_unique<TF1>("faEtaBarrelH","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1., sampleRangeHigh);
  faEtaBarrelH->SetParameter(0,-1.0774);
  faEtaBarrelH->SetParameter(1,40.5678);
  faEtaBarrelH->SetParameter(2,1.13601);
  faEtaBarrelH->SetParameter(3,0.0878601);
  faEtaBarrelH->SetParameter(4,39.5163);
  faEtaBarrelH->SetParameter(5,0.0785619);
  faEtaBarrelH->SetParameter(6,-1.39265);
  faEtaBarrelH->SetParameter(7,-1.43368);
  fbEtaBarrelH = std::make_unique<TF1>("fbEtaBarrelH","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1., sampleRangeHigh);
  fbEtaBarrelH->SetParameter(0,-26.2747);
  fbEtaBarrelH->SetParameter(1,26.5685);
  fbEtaBarrelH->SetParameter(2,14.4212);
  fbEtaBarrelH->SetParameter(3,2.15401);
  fbEtaBarrelH->SetParameter(4,0.592654);
  fbEtaBarrelH->SetParameter(5,0.622331);
  fbEtaBarrelH->SetParameter(6,-0.45401);
  fbEtaBarrelH->SetParameter(7,-0.0735657);

  faEndcap = std::make_unique<TF1>("faEndcap","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1., sampleRangeHigh);
  faEndcap->SetParameter(0,29.5353);
  faEndcap->SetParameter(1,-327.967);
  faEndcap->SetParameter(2,-830.992);
  faEndcap->SetParameter(3,0.298631);
  faEndcap->SetParameter(4,15.6682);
  faEndcap->SetParameter(5,0.242381);
  faEndcap->SetParameter(6,-0.00308457);
  faEndcap->SetParameter(7,-0.654633);
  fbEndcap = std::make_unique<TF1>("fbEndcap","[0]+([4]*(x-[5])*exp(-(x*[7])))+(([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))",1., sampleRangeHigh);
  fbEndcap->SetParameter(0,-0.180939);
  fbEndcap->SetParameter(1,32.2565);
  fbEndcap->SetParameter(2,4245.12);
  fbEndcap->SetParameter(3,0.121009);
  fbEndcap->SetParameter(4,0.0736027);
  fbEndcap->SetParameter(5,18.954);
  fbEndcap->SetParameter(6,-0.0734254);
  fbEndcap->SetParameter(7,0.0771871);
  fcEndcap = std::make_unique<TF1>("fcEndcap","([0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5])))))",1., sampleRangeHigh);
  fcEndcap->SetParameter(0,2.09235);
  fcEndcap->SetParameter(1,0.328094);
  fcEndcap->SetParameter(2,-3.98841);
  fcEndcap->SetParameter(3,34.8701);
  fcEndcap->SetParameter(4,1.17097);
  fcEndcap->SetParameter(5,0.030813);
  fcEndcap->SetParameter(6,0.66908);
  fcEndcap->SetParameter(7,-1.26556);
  faEtaEndcapEH = std::make_unique<TF1>("faEtaEndcapEH","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1., sampleRangeHigh);
  faEtaEndcapEH->SetParameter(0,-74.1561);
  faEtaEndcapEH->SetParameter(1,165.092);
  faEtaEndcapEH->SetParameter(2,-35.3997);
  faEtaEndcapEH->SetParameter(3,42.5886);
  faEtaEndcapEH->SetParameter(4,86.787);
  faEtaEndcapEH->SetParameter(5,2.24242);
  faEtaEndcapEH->SetParameter(6,0.00751736);
  faEtaEndcapEH->SetParameter(7,-0.530692);
  fbEtaEndcapEH = std::make_unique<TF1>("fbEtaEndcapEH","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1., sampleRangeHigh);
  fbEtaEndcapEH->SetParameter(0,-244.796);
  fbEtaEndcapEH->SetParameter(1,244.906);
  fbEtaEndcapEH->SetParameter(2,12.0756);
  fbEtaEndcapEH->SetParameter(3,17.4169);
  fbEtaEndcapEH->SetParameter(4,0.115151);
  fbEtaEndcapEH->SetParameter(5,0.00235061);
  fbEtaEndcapEH->SetParameter(6,-0.535085);
  fbEtaEndcapEH->SetParameter(7,-1.31019);
  faEtaEndcapH = std::make_unique<TF1>("faEtaEndcapH","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1., sampleRangeHigh);
  faEtaEndcapH->SetParameter(0,-20.7829);
  faEtaEndcapH->SetParameter(1,1.40346);
  faEtaEndcapH->SetParameter(2,-0.75444);
  faEtaEndcapH->SetParameter(3,-0.394878);
  faEtaEndcapH->SetParameter(4,-19.0334);
  faEtaEndcapH->SetParameter(5,1.72721);
  faEtaEndcapH->SetParameter(6,-0.162634);
  faEtaEndcapH->SetParameter(7,-0.294783);
  fbEtaEndcapH = std::make_unique<TF1>("fbEtaEndcapH","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1., sampleRangeHigh);
  fbEtaEndcapH->SetParameter(0,0.0310929);
  fbEtaEndcapH->SetParameter(1,62.8732);
  fbEtaEndcapH->SetParameter(2,136.556);
  fbEtaEndcapH->SetParameter(3,0.0583275);
  fbEtaEndcapH->SetParameter(4,63.047);
  fbEtaEndcapH->SetParameter(5,0.0603511);
  fbEtaEndcapH->SetParameter(6,-0.612309);
  fbEtaEndcapH->SetParameter(7,-0.650172);

  //added by Bhumika on 2 august 2018

  fcEtaBarrelH = std::make_unique<TF1>("fcEtaBarrelH", "[3]*((x-[0])^[1])+[2]", 0., sampleRangeHigh);
  fcEtaBarrelH->SetParameter(0, 0);
  fcEtaBarrelH->SetParameter(1, 2);
  fcEtaBarrelH->SetParameter(2, 0);
  fcEtaBarrelH->SetParameter(3, 1);

  fcEtaEndcapH = std::make_unique<TF1>("fcEtaEndcapH", "[3]*((x-[0])^[1])+[2]", 0., sampleRangeHigh);
  fcEtaEndcapH->SetParameter(0, 0);
  fcEtaEndcapH->SetParameter(1, 0);
  fcEtaEndcapH->SetParameter(2, 0.05);
  fcEtaEndcapH->SetParameter(3, 0);

  fdEtaEndcapH = std::make_unique<TF1>("fdEtaEndcapH", "[3]*((x-[0])^[1])+[2]", 0., sampleRangeHigh);
  fdEtaEndcapH->SetParameter(0, 1.5);
  fdEtaEndcapH->SetParameter(1, 4);
  fdEtaEndcapH->SetParameter(2, -1.1);
  fdEtaEndcapH->SetParameter(3, 1.0);

  fcEtaBarrelEH = std::make_unique<TF1>("fcEtaBarrelEH", "[3]*((x-[0])^[1])+[2]", 0., sampleRangeHigh);
  fcEtaBarrelEH->SetParameter(0, 0);
  fcEtaBarrelEH->SetParameter(1, 2);
  fcEtaBarrelEH->SetParameter(2, 0);
  fcEtaBarrelEH->SetParameter(3, 1);

  fcEtaEndcapEH = std::make_unique<TF1>("fcEtaEndcapEH", "[3]*((x-[0])^[1])+[2]", 0., sampleRangeHigh);
  fcEtaEndcapEH->SetParameter(0, 0);
  fcEtaEndcapEH->SetParameter(1, 0);
  fcEtaEndcapEH->SetParameter(2, 0);
  fcEtaEndcapEH->SetParameter(3, 0);

  fdEtaEndcapEH = std::make_unique<TF1>("fdEtaEndcapEH", "[3]*((x-[0])^[1])+[2]", 0., sampleRangeHigh);
  fdEtaEndcapEH->SetParameter(0, 1.5);
  fdEtaEndcapEH->SetParameter(1, 2.0);
  fdEtaEndcapEH->SetParameter(2, 0.6);
  fdEtaEndcapEH->SetParameter(3, 1.0);
}

// PFEnergyCalibration::CalibratedEndcapPFClusterEnergies PFEnergyCalibration::calibrateEndcapClusterEnergies(
//     reco::PFCluster const& eeCluster,
//     std::vector<reco::PFCluster const*> const& psClusterPointers,
//     ESChannelStatus const& channelStatus,
//     bool applyCrackCorrections) const {
//   double ps1_energy_sum = 0.;
//   double ps2_energy_sum = 0.;
//   bool condP1 = true;
//   bool condP2 = true;

//   for (auto const& psclus : psClusterPointers) {
//     bool cond = true;
//     for (auto const& recH : psclus->recHitFractions()) {
//       auto strip = recH.recHitRef()->detId();
//       if (strip != ESDetId(0)) {
//         //getStatusCode() == 0 => active channel
//         // apply correction if all recHits are dead
//         if (channelStatus.getMap().find(strip)->getStatusCode() == 0) {
//           cond = false;
//           break;
//         }
//       }
//     }

//     if (psclus->layer() == PFLayer::PS1) {
//       ps1_energy_sum += psclus->energy();
//       condP1 &= cond;
//     } else if (psclus->layer() == PFLayer::PS2) {
//       ps2_energy_sum += psclus->energy();
//       condP2 &= cond;
//     }
//   }

//   double ePS1 = condP1 ? -1. : 0.;
//   double ePS2 = condP2 ? -1. : 0.;

//   double cluscalibe = energyEm(eeCluster, ps1_energy_sum, ps2_energy_sum, ePS1, ePS2, applyCrackCorrections);

//   return {cluscalibe, ePS1, ePS2};
// }

void PFEnergyCalibration::energyEmHad(double t, double& e, double& h, double eta, double phi) const {
  // Use calorimetric energy as true energy for neutral particles
  const double tt = t;
  const double ee = e;
  const double hh = h;
  double etaCorrE = 1.;
  double etaCorrH = 1.;
  auto absEta = std::abs(eta);
  t = min(999.9, max(tt, e + h));
  if (t < 1.)
    return;

  // Barrel calibration
  if (absEta < 1.48) {
    // The energy correction
    double a = e > 0. ? aBarrel(t) : 1.;
    double b = e > 0. ? bBarrel(t) : cBarrel(t);
    double thresh = e > 0. ? threshE : threshH;

    // Protection against negative calibration
    if (a < -0.25 || b < -0.25) {
      a = 1.;
      b = 1.;
      thresh = 0.;
    }

    // The new estimate of the true energy
    t = min(999.9, max(tt, thresh + a * e + b * h));

    // The angular correction
    if (e > 0. && thresh > 0.) {
      etaCorrE = 1.0 + aEtaBarrelEH(t) + 1.3 * bEtaBarrelEH(t) * cEtaBarrelEH(absEta);
      etaCorrH = 1.0;
    } else {
      etaCorrE = 1.0 + aEtaBarrelH(t) + 1.3 * bEtaBarrelH(t) * cEtaBarrelH(absEta);
      etaCorrH = 1.0 + aEtaBarrelH(t) + bEtaBarrelH(t) * cEtaBarrelH(absEta);
    }
    if (e > 0. && thresh > 0.)
      e = h > 0. ? threshE - threshH + etaCorrE * a * e : threshE + etaCorrE * a * e;
    if (h > 0. && thresh > 0.) {
      h = threshH + etaCorrH * b * h;
    }

    // Endcap calibration
  } else {
    // The energy correction
    double a = e > 0. ? aEndcap(t) : 1.;
    double b = e > 0. ? bEndcap(t) : cEndcap(t);
    double thresh = e > 0. ? threshE : threshH;

    // Protection against negative calibration
    if (a < -0.25 || b < -0.25) {
      a = 1.;
      b = 1.;
      thresh = 0.;
    }

    // The new estimate of the true energy
    t = min(999.9, max(tt, thresh + a * e + b * h));

    // The angular correction
    const double dEta = std::abs(absEta - 1.5);
    const double etaPow = dEta * dEta * dEta * dEta;

    if (e > 0. && thresh > 0.) {
      if (absEta < 2.5) {
        etaCorrE = 1. + aEtaEndcapEH(t) + bEtaEndcapEH(t) * cEtaEndcapEH(absEta);
      } else {
        etaCorrE = 1. + aEtaEndcapEH(t) + 1.3 * bEtaEndcapEH(t) * dEtaEndcapEH(absEta);
      }

      etaCorrH = 1. + aEtaEndcapEH(t) + bEtaEndcapEH(t) * (0.04 + etaPow);
    } else {
      etaCorrE = 1.;
      if (absEta < 2.5) {
        etaCorrH = 1. + aEtaEndcapH(t) + bEtaEndcapH(t) * cEtaEndcapH(absEta);
      } else {
        etaCorrH = 1. + aEtaEndcapH(t) + bEtaEndcapH(t) * dEtaEndcapH(absEta);
      }
    }

    //t = min(999.9,max(tt, thresh + etaCorrE*a*e + etaCorrH*b*h));

    if (e > 0. && thresh > 0.)
      e = h > 0. ? threshE - threshH + etaCorrE * a * e : threshE + etaCorrE * a * e;
    if (h > 0. && thresh > 0.) {
      h = threshH + etaCorrH * b * h;
    }
  }

  // Protection
  if (e < 0. || h < 0.) {
    // Some protection against crazy calibration
    if (e < 0.)
      e = ee;
    if (h < 0.)
      h = hh;
  }

  // And that's it !
}

// The calibration functions
double PFEnergyCalibration::aBarrel(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfa_BARREL, point);

  // } else {
    return faBarrel->Eval(x);
  // }
}

double PFEnergyCalibration::bBarrel(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfb_BARREL, point);

  // } else {
    return fbBarrel->Eval(x);
  // }
}

double PFEnergyCalibration::cBarrel(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfc_BARREL, point);

  // } else {
    return fcBarrel->Eval(x);
  // }
}

double PFEnergyCalibration::aEtaBarrelEH(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfaEta_BARRELEH, point);

  // } else {
    return faEtaBarrelEH->Eval(x);
  // }
}

double PFEnergyCalibration::bEtaBarrelEH(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfbEta_BARRELEH, point);

  // } else {
    return fbEtaBarrelEH->Eval(x);
  // }
}

double PFEnergyCalibration::aEtaBarrelH(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfaEta_BARRELH, point);

  // } else {
    return faEtaBarrelH->Eval(x);
  // }
}

double PFEnergyCalibration::bEtaBarrelH(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfbEta_BARRELH, point);

  // } else {
    return fbEtaBarrelH->Eval(x);
  // }
}

double PFEnergyCalibration::aEndcap(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfa_ENDCAP, point);

  // } else {
    return faEndcap->Eval(x);
  // }
}

double PFEnergyCalibration::bEndcap(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfb_ENDCAP, point);

  // } else {
    return fbEndcap->Eval(x);
  // }
}

double PFEnergyCalibration::cEndcap(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfc_ENDCAP, point);

  // } else {
  return fcEndcap->Eval(x);
  // }
}

double PFEnergyCalibration::aEtaEndcapEH(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfaEta_ENDCAPEH, point);

  // } else {
  return faEtaEndcapEH->Eval(x);
  // }
}

double PFEnergyCalibration::bEtaEndcapEH(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfbEta_ENDCAPEH, point);

  // } else {
  return fbEtaEndcapEH->Eval(x);
  // }
}

double PFEnergyCalibration::aEtaEndcapH(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfaEta_ENDCAPH, point);

  // } else {
  return faEtaEndcapH->Eval(x);
  // }
}

double PFEnergyCalibration::bEtaEndcapH(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfbEta_ENDCAPH, point);
  // } else {
  return fbEtaEndcapH->Eval(x);
  // }
}

//added by Bhumika Kansal on 3 august 2018

double PFEnergyCalibration::cEtaBarrelH(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfcEta_BARRELH, point);

  // } else {
  return fcEtaBarrelH->Eval(x);
  // }
}
double PFEnergyCalibration::cEtaEndcapH(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfcEta_ENDCAPH, point);

  // } else {
  return fcEtaEndcapH->Eval(x);
  // }
}

double PFEnergyCalibration::dEtaEndcapH(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfdEta_ENDCAPH, point);

  // } else {
  return fdEtaEndcapH->Eval(x);
  // }
}

double PFEnergyCalibration::cEtaBarrelEH(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfcEta_BARRELEH, point);

  // } else {
  return fcEtaBarrelEH->Eval(x);
  // }
}

double PFEnergyCalibration::cEtaEndcapEH(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfcEta_ENDCAPEH, point);

  // } else {
  return fcEtaEndcapEH->Eval(x);
  // }
}

double PFEnergyCalibration::dEtaEndcapEH(double x) const {
  // if (pfCalibrations) {
  //   BinningPointByMap point;
  //   point.insert(BinningVariables::JetEt, x);
  //   return pfCalibrations->getResult(PerformanceResult::PFfdEta_ENDCAPEH, point);

  // } else {
  return fdEtaEndcapEH->Eval(x);
  // }
}

// double PFEnergyCalibration::energyEm(const reco::PFCluster& clusterEcal,
//                                      double ePS1,
//                                      double ePS2,
//                                      bool crackCorrection) const {
//   return Ecorr(clusterEcal.energy(), ePS1, ePS2, clusterEcal.eta(), clusterEcal.phi(), crackCorrection);
// }

// double PFEnergyCalibration::energyEm(const reco::PFCluster& clusterEcal,
//                                      double ePS1,
//                                      double ePS2,
//                                      double& ps1,
//                                      double& ps2,
//                                      bool crackCorrection) const {
//   return Ecorr(clusterEcal.energy(), ePS1, ePS2, clusterEcal.eta(), clusterEcal.phi(), ps1, ps2, crackCorrection);
// }

std::ostream& operator<<(std::ostream& out, const PFEnergyCalibration& calib) {
  if (!out)
    return out;

  out << "PFEnergyCalibration -- " << endl;

  // if (calib.pfCalibrations) {
  //   static const std::map<std::string, PerformanceResult::ResultType> functType = {
  //       {"PFfa_BARREL", PerformanceResult::PFfa_BARREL},
  //       {"PFfa_ENDCAP", PerformanceResult::PFfa_ENDCAP},
  //       {"PFfb_BARREL", PerformanceResult::PFfb_BARREL},
  //       {"PFfb_ENDCAP", PerformanceResult::PFfb_ENDCAP},
  //       {"PFfc_BARREL", PerformanceResult::PFfc_BARREL},
  //       {"PFfc_ENDCAP", PerformanceResult::PFfc_ENDCAP},
  //       {"PFfaEta_BARRELH", PerformanceResult::PFfaEta_BARRELH},
  //       {"PFfaEta_ENDCAPH", PerformanceResult::PFfaEta_ENDCAPH},
  //       {"PFfbEta_BARRELH", PerformanceResult::PFfbEta_BARRELH},
  //       {"PFfbEta_ENDCAPH", PerformanceResult::PFfbEta_ENDCAPH},
  //       {"PFfaEta_BARRELEH", PerformanceResult::PFfaEta_BARRELEH},
  //       {"PFfaEta_ENDCAPEH", PerformanceResult::PFfaEta_ENDCAPEH},
  //       {"PFfbEta_BARRELEH", PerformanceResult::PFfbEta_BARRELEH},
  //       {"PFfbEta_ENDCAPEH", PerformanceResult::PFfbEta_ENDCAPEH},
  //       {"PFfcEta_BARRELH", PerformanceResult::PFfcEta_BARRELH},
  //       {"PFfcEta_ENDCAPH", PerformanceResult::PFfcEta_ENDCAPH},
  //       {"PFfdEta_ENDCAPH", PerformanceResult::PFfdEta_ENDCAPH},
  //       {"PFfcEta_BARRELEH", PerformanceResult::PFfcEta_BARRELEH},
  //       {"PFfcEta_ENDCAPEH", PerformanceResult::PFfcEta_ENDCAPEH},
  //       {"PFfdEta_ENDCAPEH", PerformanceResult::PFfdEta_ENDCAPEH}

  //   };

  //   for (std::map<std::string, PerformanceResult::ResultType>::const_iterator func = functType.begin();
  //        func != functType.end();
  //        ++func) {
  //     cout << "Function: " << func->first << endl;
  //     PerformanceResult::ResultType fType = func->second;
  //     calib.pfCalibrations->printFormula(fType);
  //   }

  // } else {
    std::cout << "Default calibration functions : " << std::endl;

    calib.faBarrel->Print();
    calib.fbBarrel->Print();
    calib.fcBarrel->Print();
    calib.faEtaBarrelEH->Print();
    calib.fbEtaBarrelEH->Print();
    calib.faEtaBarrelH->Print();
    calib.fbEtaBarrelH->Print();
    calib.faEndcap->Print();
    calib.fbEndcap->Print();
    calib.fcEndcap->Print();
    calib.faEtaEndcapEH->Print();
    calib.fbEtaEndcapEH->Print();
    calib.faEtaEndcapH->Print();
    calib.fbEtaEndcapH->Print();
    //
  // }

  return out;
}


/*********************************************FIN NUEVO BY MIKKO *******************************/

double Calibration::getCalibratedEnergy(double ETrue, double ecalEnergy, 
                                        double hcalEnergy, int DifferentTerms)//This DifferentTerms int variable aims to take into account only the EcalEnergy term (value 1), the HcalEnergy term (value 2), the sum of EcalEnergy and HcalEnergy terms (value 3) and the total correction EcalEnergy + HcalEnergy + independent terms (value 0). To do this, it ignores unnecessary parameters.
{
  double a = functionA_->Eval(ETrue);
  double b = functionB_->Eval(ETrue);
  double c = functionC_->Eval(ETrue);

  if(DifferentTerms==1){
    a = 0.0;
    c = 1.0;
  }else if(DifferentTerms==2){
    a = 0.0;
    b = 1.0;
  }else if(DifferentTerms==3){
    a = 0.0;
  }

  return a+ b*ecalEnergy + c*hcalEnergy;
}

//eta-formula
double Calibration::getCalibratedEnergy(double ETrue, double ecalEnergy, 
                                        double hcalEnergy, double eta, int DifferentTerms)//This DifferentTerms int variable aims to take into account only the contribution of alpha and beta parameters. So, for taking only alpha (value 1), for taking only beta (value 2), and for taking both, alpha and beta (value 0). To do this, it ignores unnecessary parameters.
{
   double etaPow;
   double factor_;
   double a = functionA_->Eval(ETrue);
   double b = functionB_->Eval(ETrue);
   double c = functionC_->Eval(ETrue);
   double alpha = functionAlpha_->Eval(ETrue);
   double beta = functionBeta_->Eval(ETrue);
   double counterAlpha = 0;
   double counterBeta = 0;

   if(isBarrel_) 
   {
      etaPow = eta*eta;
      factor_ = factorB;
      if(ecalEnergy>0) {  //shubham
	counterAlpha = alpha;
	counterBeta = beta;
      }
   
   }
   else 
     {//endcap
       
       if(ecalEnergy > 0) {//EH-hadrons
      //  counterAlpha = alpha;
      //  counterBeta = beta;
	 if( fabs(eta)>2.5) {//EH-hadrons + ENDCAP II region
           
	   //           etaPow = -0.3 + 1.3*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5);

	   etaPow = -0.5 + 1.3*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5);
     //etaPow = 0.04+(fabs(eta)-1.5)*(fabs(eta)-1.5)*(fabs(eta)-1.5)*(fabs(eta)-1.5);//Using 2018 functions of eta AN2022_015
	   //	   etaPow = -0.3 + 1.3*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5);//*(fabs(eta) - 1.5) ; //change for UL2017

	 }
	 else {//EH-hadrons + ENDCAP I region
	   etaPow=0; 
     //etaPow = -0.8 + 2.2*(fabs(eta)-1.5)*(fabs(eta)-1.5);//Using 2018 functions of eta AN2022_015
	   //etaPow = 0.8 - 2*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5);  
	   //	   etaPow = -0.1 + 0.5*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta)-1.5);//*(fabs(eta)-1.5);
	   //	   etaPow = -0.8 + 0.5*(fabs(eta) - 1.5)*(fabs(eta) - 1.5);
	   //	   etaPow = 0.8 - 0.6*(fabs(eta) - 1.5)*(fabs(eta) - 1.5) + 1.1*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5);  
	 }
       }
       else  // H hadrons here
	 {
	   
	   if( fabs(eta)<2.5) {//H-hadrons + ENDCAP I region
	     etaPow=0; //for UL2016 H
       //etaPow = 0.05;//Using 2018 functions of eta AN2022_015
	     //etaPow = 0.08 - 0.01*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5); //current
	     
	   }
	   else  {//H-hadrons + ENDCAP II region
	     //  etaPow =1.2*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5);
	     etaPow = 1.1*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5) ; //current
	     //etaPow = -2.0+(fabs(eta)-1.5)*(fabs(eta)-1.5);//Using 2018 funtions of eta AN2022_015
       //Giving better result
	     //etaPow = -0.6*(fabs(eta) - 1.5)*(fabs(eta) - 1.5) + 1.1*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5)*(fabs(eta) - 1.5);
	   }
	   
	 }
       factor_ = factorE;
       
     }
   

   /*
   if(fabs(eta)>2.5 && fabs(eta)<2.51 && ecalEnergy>0)
     cout<<alpha<<"  "<<beta<<"   "<<ETrue<<endl;
   */
  
  if(DifferentTerms==1){
    beta = 0.0;
  }else if(DifferentTerms==2){
    alpha = 0.0;
  }

   return a + (1.0 + alpha + factor_*beta*etaPow)*b*ecalEnergy + 
      (1.0 + alpha + beta*etaPow - counterAlpha - counterBeta*etaPow)*
      c*hcalEnergy;
}


///////////////////////////////////////////////////////////////////////////////
//Global functions used in the main of the code.
///////////////////////////////////////////////////////////////////////////////


//Takes apart a TH2 and creates response and resolution plots from it. Note: 
//the TGraphs will not draw correctly without passing the TGraphs as references
//to the function.
double CalculateMedian(TH1F* h)
{
  // TFile f("histos.root");
  // TH1F *h = (TH1F*)f.Get("hgaus");

  int numBins = h->GetXaxis()->GetNbins();
  Double_t* x = new Double_t[numBins];
  Double_t* y = new Double_t[numBins];
  for (int i = 0; i < numBins; i++) {
    x[i] = h->GetBinCenter(i);
    y[i] = h->GetBinContent(i);
  }
  double MedianOfHisto = TMath::Median(numBins, &x[0], &y[0]);

  return MedianOfHisto;

  // cout<<"Median -----------> \t"<<MedianOfHisto;
  // return 0;
}




void drawGausFit(TH2F* inHisto, TGraph& response, TGraph& resolution)
{
  cout << "ENTRO EN DRAWGAUSFIT" << endl;
   if(inHisto->GetEntries() == 0) return;
   
   vector<TH1F*> ETrueBin;
   TF1* gaus; 
   string name;
   int digitsSampleRangeHigh = int(std::log10(sampleRangeHigh) + 1);
   char num[digitsSampleRangeHigh];//This value must be at least equal to the number of digits in the sampleRangeHigh value.
   float rebin = .5;//.5
   TCanvas* canvas;
   TCanvas* temp = new TCanvas();
   //TLine *line = new TLine(0.0,0.0,sampleRangeHigh,0.0);
   int rangelow_ = 0.2, rangehigh_ = sampleRangeHigh, bins_ = sampleRangeHigh;

   if (drawpT) {
     rangehigh_ = 20;
     bins_ = 10;
     rebin = rebin/2;
   }
   TLine *line = new TLine(0.0,0.0,rangehigh_,0.0);
   
   // TH2F* respHisto = new TH2F("respHisto", "", sampleRangeHigh, 0, sampleRangeHigh, 100, -0.5, 0.5);
   // //TH2F* resoHisto = new TH2F("resoHisto", "", sampleRangeHigh, 0, sampleRangeHigh, 100, 0.0, 0.5);
   // TH2F* resoHisto = new TH2F("resoHisto", "", sampleRangeHigh, 0, sampleRangeHigh, 200, 0.0, 1.0);
  double Ymax=0.50;
  double Ymin=-0.50;
  if(drawpT){
    Ymax=0.3;
    Ymin=-0.3;
  } 
   
  TH2F* respHisto = new TH2F("respHisto", "", bins_, rangelow_, rangehigh_, 100, Ymin, Ymax);
  //  if(inHisto->GetName()=="corrBarrelEcalHcal_ErawEcal"){
  //   TH2F* respHisto = new TH2F("respHisto", "", bins_, rangelow_, rangehigh_, 100, -1, 15);
  //  }
   //TH2F* resoHisto = new TH2F("resoHisto", "", sampleRangeHigh, 0, sampleRangeHigh, 100, 0.0, 0.5);
   TH2F* resoHisto = new TH2F("resoHisto", "", bins_, rangelow_, rangehigh_, 200, 0.0, 1.0);



  //  TGraph averages;
  //  TGraph rmss;

   vector<double> ETrue;
   vector<double> gausMean; 
   vector<double> gausSigma;
   vector<double> average;
   vector<double> rms;

   char* fileName = new char[1000];
   sprintf(fileName,"projections_%s_%s.root",_region_,inHisto->GetName());
   TFile* file1=new TFile(fileName,"recreate");

   
   temp->cd();//This TCanvas is only used since when we do the fit down below 
              //it creates an unwanted TCanvas. We will get rid of it later on 
              //in the function.

   // TCanvas* cccc=new TCanvas("balda","bacla");
   //cout<<"ETrue.back(), gausMean[0].back()"<<endl;
   //cout<<"**********Draw Gaus**********"<<endl;
   for(unsigned bin = 0.2; bin <= rangehigh_; )
   {
     double x_min = -1.0, x_max = 1.0;
     if(strcmp(inHisto->GetName(),"corrEtaEndcapEcalHcal") == 0 && strcmp(_region_,"EC_outside_tracker") == 0 && false)
       x_min = -0.70;
      name = "histcorhybrid";
      sprintf(num, "%i", bin);
      name += num;
      //Split up the TH2 into many TH1's for each ETrue bin.
     
      ETrueBin.push_back((TH1F*)inHisto->ProjectionY(name.c_str(),bin, 
                                                     bin + 4*rebin));
      cout <<"bin to  bin + 4*rebin: "<<bin<<" to "<<(bin + 4*rebin)<<"   "<<ETrueBin.back()->GetEntries()<<endl;
      if(ETrueBin.back()->GetEntries() > 5)
	{
	  //Fit each ETrue bin to a gaus (iteratively done to get better fit)
	  //cout<<"ETrueBin.back()->GetEntries():"<<ETrueBin.back()->GetEntries()<<endl;
	  if(bin > 2) {

	    gaus =new TF1("gaus","gaus",-3,3);
	    gaus->SetParameters(1000.,0.,0.2);
	    //ETrueBin.back()->Fit("gaus", "Q", "", -1.0, 1.0);
	    ETrueBin.back()->Fit("gaus", "Q", "", x_min, x_max);
	    //ETrueBin.back()->Fit("gaus", "Q", "", -0.7, 0.7);
	    
	    
	    gaus = ETrueBin.back()->GetFunction("gaus");
	    //	    cout<<bin<<" "<<gaus->GetParameter(1)<<"   "<<gaus->GetParameter(2)<<endl;
	    	    
	    if(gaus->GetParameter(2) < 0)
	      goto here1;
	    x_min = gaus->GetParameter(1)- gaus->GetParameter(2);
      
	    if((strcmp(inHisto->GetName(),"corrEtaEndcapEcalHcal") == 0 || strcmp(inHisto->GetName(),"corrEndcapEcalHcal") == 0) && strcmp(_region_,"EC_outside_tracker") == 0 && false)
	      x_min = (gaus->GetParameter(1) - gaus->GetParameter(2)) < -0.7 ? -0.7 : (gaus->GetParameter(1) - gaus->GetParameter(2));

		 
	    //x_max = (gaus->GetParameter(1)+ 2* gaus->GetParameter(2))>1.0 ? 1.0 : (gaus->GetParameter(1)+ 2* gaus->GetParameter(2));
	    x_max = (gaus->GetParameter(1)+ 2* gaus->GetParameter(2));
	    ETrueBin.back()->Fit("gaus", "Q", "", x_min, x_max);

	    // ETrueBin.back()->Fit("gaus", "Q", "",
	    // 			 gaus->GetParameter(1) - 2*gaus->
	    // 			 GetParameter(2), 1.0);  

	    gaus = ETrueBin.back()->GetFunction("gaus");


            x_min = gaus->GetParameter(1)- gaus->GetParameter(2);
            x_max = (gaus->GetParameter(1)+ 2* gaus->GetParameter(2));

            if((strcmp(inHisto->GetName(),"corrEtaEndcapEcalHcal") == 0 || strcmp(inHisto->GetName(),"corrEndcapEcalHcal") == 0) && strcmp(_region_,"EC_outside_tracker") == 0)
	      {
		x_min=-0.5;
		x_max=0.2;
		if(ETrue.back()<114) x_min=-0.85;
		if(ETrue.back()>154) x_max=0.1;
		if(strcmp(inHisto->GetName(),"corrEndcapEcalHcal") == 0) x_max=0.6;
	      }
            //x_max = (gaus->GetParameter(1)+ 2* gaus->GetParameter(2))>1.0 ? 1.0 : (gaus->GetParameter(1)+ 2* gaus->GetParameter(2));
	    //            x_max = (gaus->GetParameter(1)+ 2* gaus->GetParameter(2));



      // if(strcmp(inHisto->GetName(),"corrEtaBarrelEcalHcal") == 0 && strcmp(_region_,"barrel") == 0 && bin==32)
      // 	{
      // 	  x_min = (gaus->GetParameter(1) - gaus->GetParameter(2)) < -0.4 ? -0.4 : (gaus->GetParameter(1) - gaus->GetParameter(2));
      // 	  x_max = (gaus->GetParameter(1) + gaus->GetParameter(2)) > 0.38 ? 0.38 : (gaus->GetParameter(1) + gaus->GetParameter(2));
      // 	  cout<<"x_min :  "<<x_min<<"      x_max:  "<<x_max<<endl;	 
      // 	}
	    
	    ETrueBin.back()->Fit("gaus", "Q", "", x_min, x_max);


            // ETrueBin.back()->Fit("gaus", "Q", "",
            //                      gaus->GetParameter(1) - 2*gaus->
            //                      GetParameter(2), 1.0);
            gaus = ETrueBin.back()->GetFunction("gaus");

	  here1:
	    // if(bin<=16)
	    // gausMean.push_back(ETrueBin.back()->GetMean());
	    // else

	    
            //gausSigma.push_back(gaus->GetParameter(2)/(1.0 + min(0.0, gaus->GetParameter(1))));

	    if (useMedian) {
	      gausMean.push_back(CalculateMedian(ETrueBin.back()));
	      gausSigma.push_back(CalculateMedian(ETrueBin.back())/(1.0 + min(0.0, ETrueBin.back()->GetMean())));

	    }
	    else if (useMean) {
	      gausMean.push_back(ETrueBin.back()->GetMean());
	      gausSigma.push_back(ETrueBin.back()->GetRMS()/(1.0 + min(0.0, ETrueBin.back()->GetMean())));
	    }

	    else {
	      gausMean.push_back(gaus->GetParameter(1));
	      gausSigma.push_back(gaus->GetParameter(2)/(1.0 + min(0.0, gaus->GetParameter(1))));
	    }
	    //gausMean.push_back(gaus->GetParameter(1));
       	  }
	   else {
	   
	     gaus =new TF1("gaus","gaus",-3,3);
	     gaus->SetParameters( 500, 10, 5, 0, 0.20 );
	     gaus->FixParameter(2,5);
	     ETrueBin.back()->Fit("gaus", "QN0", "", -1.0, 1.0);
	     ETrueBin.back()->Fit("gaus", "QN0", "", -1.0, 1.0);
	     ETrueBin.back()->Fit("gaus", "Q", "", -1.0, 1.0);

	     gausMean.push_back(gaus->GetParameter(1));
	     gausSigma.push_back(gaus->GetParameter(2)/
				 (1.0 + min(0.0, gaus->GetParameter(1))));
         //cout << endl << "PARAMETRO0: " << gaus->GetParameter(0) << "\tPARAMETRO1: " << gaus->GetParameter(1) << "\tPARAMETRO2: " << gaus->GetParameter(2) << "\tPARAMETRO3: " << gaus->GetParameter(3) << "\tPARAMETRO4: " << gaus->GetParameter(4) << endl;
	   }


	  // cout<<bin<<"   "<<median1(ETrueBin.back())<<endl;

	  // TFile oFile( ("tmp/"+name+".root").c_str() ,"RECREATE");
	  // ETrueBin.back()->Write();
	  // gaus->Write();
	  // oFile.Close();

	  // cccc->cd();
	  // ETrueBin.back()->Draw();
	  // // //   gaus->Draw("same");
	  // cccc->SaveAs( ("tmp/"+name+".png").c_str() );
	  // cccc->SaveAs( ("tmp/"+name+".C").c_str() );

            ETrue.push_back(bin + 2.0*rebin);
	    //cout<<"bin:"<<bin<<", rebin:"<<rebin<<", bin + 2*rebin:"<<(bin + 2*rebin)<<endl;
	    //cout<<ETrue.back()<<", "<<gausMean.back()<<endl;
	    //shubham
	    //cout<<ETrue.back()<<" ";
	    //  if(bin<=16)
	    //  gausMean.push_back(ETrueBin.back()->GetMean());
	    // else
	    //   gausMean.push_back(gaus->GetParameter(1));

	    
	   
            average.push_back(ETrueBin.back()->GetMean());
            rms.push_back(ETrueBin.back()->GetRMS());
	    // if((!ETrueBin.back()->GetMean())) {
      //   bin += 1*rebin;
      
      // //Increase bin size with increasing ETrue since there are fewer high 
      // //energy events than low energy ones.
      // if(bin > 20) rebin = 2.0;
      // if(bin > 100) rebin = 5.0; //20
      // if(bin > 1000) rebin = 20.0; //50
      //   continue;
      // }
      cout<<bin<<"   "<<ETrue.back()<<"   "<<ETrueBin.back()->GetMean()<<" <> "<<gausMean.back()<<"   "<<ETrueBin.back()->GetRMS()<<"   "<<gausSigma.back()<<endl;

	    //cout<<bin<<"   "<<ETrue.back()<<"   "<<ETrueBin.back()->GetMean()<<" <> "<<gausMean.back()<<"   "<<ETrueBin.back()->GetMeanError()<<"   "<<gaus->GetParError(1)<<endl;
	    //cout<<bin<<"   "<<ETrue.back()<<"   "<<ETrueBin.back()->GetMean()<<" <> "<<gausMean.back()<<"   "<<ETrueBin.back()->GetMeanError()<<endl;

	    if (false)
	      gaus->Delete();

	    (ETrueBin.back())->Write();


	}

      
      bin += 2*rebin;
      
      //Increase bin size with increasing ETrue since there are fewer high 
      //energy events than low energy ones.
      if(bin > 10) rebin = 2.;//2.0;
      if(bin > 100) rebin = 5.;//5.0; //20
      if(bin > 400) rebin = 20.0; //50
      if(bin > 1000) rebin = 50.;
      if(bin > 2000) rebin = 100.;
      //delete gaus;
      
   }

   file1->Close();
   // delete cccc;
   //Added by bhumika 1 april 2019
   sprintf(fileName,"resp_%s_%s.root",_region_,inHisto->GetName());
   TFile* file2=new TFile(fileName,"recreate");
   file2->cd();
   response = TGraph(ETrue.size(), &ETrue[1], &gausMean[1]); //Fill the graphs
   response.SetName("response");
   resolution = TGraph(ETrue.size(),&ETrue[1], &gausSigma[1]);
   resolution.SetName("resolution");
   //averages =  TGraph(ETrue.size(), &ETrue[1], &average[1]);
   //rmss = TGraph(ETrue.size(), &ETrue[1], &rms[1]);
   //
   // response = TGraph(ETrue.size(), &ETrue[0], &gausMean[0]); //Fill the graphs
   // //response = TGraph(ETrue.size(), &ETrue[0], &average[0]); //Fill the graphs
   // resolution = TGraph(ETrue.size(),&ETrue[0], &gausSigma[0]);
   // averages =  TGraph(ETrue.size(), &ETrue[0], &average[0]);
   // rmss = TGraph(ETrue.size(), &ETrue[0], &rms[0]);

   //Set up the graphs to look how you want them to.
   response.SetMarkerStyle(22);
   response.SetMarkerSize(0.8);
   response.SetMarkerColor(4);

   resolution.SetMarkerStyle(22);
   resolution.SetMarkerSize(0.8);
   resolution.SetMarkerColor(4);

  //  averages.SetMarkerStyle(22);
  //  averages.SetMarkerSize(0.8);
  //  averages.SetMarkerColor(4);

  //  rmss.SetMarkerStyle(22);
  //  rmss.SetMarkerSize(0.8);
  //  rmss.SetMarkerColor(4);

   line->SetLineStyle(1);
   line->SetLineWidth(2);
   line->SetLineColor(2);


   //  gStyle->SetOptStat(0); 
   //gStyle->SetOptFit(0);
   //canvas = new TCanvas(("canvas "+ (string)(inHisto->GetName()) ).c_str(), ("Response and Resolution "+ (string)(inHisto->GetName())).c_str(), 1000, 500);
   //spandey
   //canvas = new TCanvas(("canvas "+ (string)(inHisto->GetName()) ).c_str(), ("Response and Resolution "+ (string)(inHisto->GetName())).c_str(), 800, 400);
   canvas = new TCanvas(("canvas "+ (string)(inHisto->GetName()) ).c_str(), ("Response"+ (string)(inHisto->GetName())).c_str(), 1600, 900);


 
   //canvas->Divide(2, 1);
   temp->~TCanvas();  //destroy the TCanvas 

   //canvas->cd(1);
   gPad->SetGridx();
   gPad->SetGridy();
   canvas->SetLogx();
   respHisto->SetStats(0);
   respHisto->SetTitle(("Response "+ (string)(inHisto->GetName()) ).c_str());
   respHisto->Draw();
   response.Draw("P");
   line->Draw();

  //  canvas->cd(2);
  //  gPad->SetGridx();
  //  gPad->SetGridy();
  //  resoHisto->SetStats(0);
  //  resoHisto->SetTitle("Resolution");
  //  resoHisto->Draw();
  //  resolution.Draw("P");


   respHisto->GetYaxis()->SetTitle("(E_{cor}-E_{true})/E_{true}");
   respHisto->GetXaxis()->SetTitle("E_{true} [GeV]");

  //  resoHisto->GetYaxis()->SetTitle("#sigma(E)/E_{true}");
  //  resoHisto->GetXaxis()->SetTitle("E_{true} [GeV]");

   if(drawResoFit) {
     //MM Fit Resolution
     TF1* f=new TF1( ("ResoFit"+ (string)(inHisto->GetName())).c_str(),"sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))",20,1000);// 3.5*4
     f->SetParameters(0.06,1.20,0.);
     f->SetParLimits(0,0,10);
     f->SetParLimits(1,0,10);
     f->SetParLimits(2,0,10);
     resolution.Fit(("ResoFit"+ (string)(inHisto->GetName())).c_str(),"QR");
     resolution.Fit(("ResoFit"+ (string)(inHisto->GetName())).c_str(),"QR");
     resolution.Fit(("ResoFit"+ (string)(inHisto->GetName())).c_str(),"R");
     
     
     
     string legend;
     int fres0 = (int)(f->GetParameter(0)*100.);
     int fres1 = (int)(10.*(f->GetParameter(0)*100.-fres0));
     int fres2 = (int)(f->GetParameter(1)*100.);
     // char text[100];
     // sprintf(text,"#sigma/E = %i%/#sqrt{E} + %i.%i%",fres2,fres0,fres1);
     TString text = "#sigma/E = ";
     text+=(int)fres2;
     text+="%/#sqrt{E} + ";
     text+=(int)fres0;
     text+=".";
     text+=(int)fres1;
     
     legend += text;
     TLegend *leg=new TLegend(0.30,0.75,0.85,0.85);
     leg->AddEntry((&resolution),legend.c_str(),"lp");
     leg->SetTextSize(0.04);
     leg->Draw();
   }
   if (saveCanvas) {
     //string cname = ((string)(inHisto->GetName()) ) + ".gif";
     char  cname[200];
     //sprintf(cname,  "%s_%s_updatedCode.gif",inHisto->GetName(),_region_);
     //canvas->Print(cname);
     //canvas->SaveAs(cname);

     sprintf(cname,  "%s_%s_updatedCode.png",inHisto->GetName(),_region_);
     canvas->Print(cname);
     canvas->SaveAs(cname);
   }
   //Added by Bhumika 1 April 2019
 respHisto->Write();
   //  line->Write();
 resoHisto->Write();
   response.Write();
   resolution.Write();
   file2->cd();
   file2->Write();
   file2->Close();
   cout << "SALGO DE DRAWGAUSFIT" << endl;

}

void drawEtaDependence(TH2F* inHisto, TGraph& responseEta)
{
   if(inHisto->GetEntries() == 0) return;

   vector<TH1F*> etaBin;
   TF1* gaus; 
   TString name;
   //char num[4];

   TCanvas* canvas;
   TCanvas* temp = new TCanvas();
   TLine* line = new TLine(0, 0, 3, 0);

   TH2F* respHisto = new TH2F("respHisto", "", 30, 0.0, 3.00, 10000, -100.0, 100.0);

   TGraph averages;
   TGraph rmss;

   vector<double> etaAverage;
   vector<double> gausMean; 
   vector<double> gausSigma;
   vector<double> average;
   vector<double> etaRms;

   char* fileName2 = new char[1000];
   sprintf(fileName2,"projections_%s_%s_eta.root",_region_,inHisto->GetName());
   //   TFile* file1=new TFile(fileName,"recreate");

   TFile* file1=new TFile(fileName2,"recreate");
   temp->cd();//This TCanvas is only used since when we do the fit down below 
              //it creates an unwanted TCanvas. We will get rid of it later on 
              //in the function.

   float mR=0;
   float mR2=0;
   int N=0;

   // //  TCanvas* cccc=new TCanvas("bala","bala");

   // 
for(unsigned bin = 1; bin < (unsigned)inHisto->GetNbinsX(); bin = bin + 1)
   {
      name = "histEta";
      //   sprintf(num, "%i", bin);
      TString name2 = name;
      name2 += bin;
      //Split up the TH2 into many TH1's for each eta bin.
   
      etaBin.push_back((TH1F*)inHisto->ProjectionY(name,bin, bin + 1));
      
      name += inHisto->GetXaxis()->GetBinCenter(bin);
      

      if(etaBin.back()->GetEntries() > 0)
      {
         //Fit each eta bin to a gaus (iteratively done to get better fit)
	double x_min=0, x_max=0;
       
	if(strcmp(inHisto->GetName(),"corrEtaDependenceEH") == 0  && (strcmp(_region_,"EC_outside_tracker")|| strcmp(_region_,"Full") == 0))
	  {
	    x_min=-0.2;
	    x_max=0.15;
	  }
	
	gaus =new TF1("gaus","gaus(0)",x_min,x_max);
	//gaus->SetParameters(500.,0.,0.2);
	gaus->SetParameters(500.,etaBin.back()->GetMean(),etaBin.back()->GetRMS());

	float rrms  = etaBin.back()->GetRMS();
	float mmean = etaBin.back()->GetMean();
	etaBin.back()->Fit("gaus", "Q", "", -1.0, 1.0);
	    //	etaBin.back()->Fit("gaus", "Q", "", mmean-rrms, mmean+rrms);
	//etaBin.back()->Fit("gaus", "Q", "",x_min,x_max);
	    // gaus = etaBin.back()->GetFunction("gaus");
	
	int nsig = 2.0;
	if(inHisto->GetXaxis()->GetBinLowEdge(bin) > 1.5) nsig = 1;
	etaBin.back()->Fit("gaus", "Q", "",
			   gaus->GetParameter(1) - nsig*gaus->
			   GetParameter(2), 1.0);  
	// gaus = etaBin.back()->GetFunction("gaus");
	
	etaBin.back()->Fit("gaus", "Q", "",
			   gaus->GetParameter(1) - nsig*gaus->
			   GetParameter(2), 1.0);

	if(inHisto->GetXaxis()->GetBinLowEdge(bin) > 1.5 && inHisto->GetXaxis()->GetBinLowEdge(bin) < 2.7)
	  etaBin.back()->Fit("gaus", "Q", "",x_min,x_max);

	
         // etaBin.back()->Fit("gaus", "Q", "",
         //                    gaus->GetParameter(1) - gaus->
         //                    GetParameter(2), 1.0);  
	 // // gaus = etaBin.back()->GetFunction("gaus");
         
         // etaBin.back()->Fit("gaus", "Q", "",
         //                    gaus->GetParameter(1) - gaus->
         //                    GetParameter(2), 1.0);
	 // gaus = etaBin.back()->GetFunction("gaus");
         
         etaAverage.push_back(inHisto->GetXaxis()->GetBinCenter(bin));
         etaRms.push_back(0.1);
	 
	 if (useMean) 
	   gausMean.push_back(etaBin.back()->GetMean());
	 else 
	   gausMean.push_back( gaus->GetParameter(1) );
         
         gausSigma.push_back(etaBin.back()->GetRMS());

	 if(etaAverage.back()>1.6) {
	   mR += gausMean.back();
	   mR2 += gausMean.back()*gausMean.back();
	   N++;
	 }

	 // cccc->cd();
	 // etaBin.back()->Draw();
	 // gaus->Draw("same");
	 // cccc->SaveAs( ("tmp/Eta"+name+".png") );
	 // cccc->SaveAs( ("tmp/Eta"+name+".root") );
	 //cout<<gaus->GetParameter(1)<<endl;


	 (etaBin.back())->Write();
      }
            
      //Increase bin size with increasing eta since there are fewer high 
      //energy events than low energy ones.
   }

   //  delete cccc;
   file1->Close();
   char* fileName = new char[1000];
   sprintf(fileName,"resp_%s_wrtEta.root",inHisto->GetName());
   TFile* file2=new TFile(fileName,"recreate");

   //   TFile* file2=new TFile("output2.root","recreate");
   file2->cd();

   responseEta = TGraph(etaAverage.size(), &etaAverage[0], &gausMean[0]); 
//&etaRms[0], &gausSigma[0]); 

   responseEta.SetMarkerStyle(22);
   responseEta.SetMarkerSize(1);
   responseEta.SetMarkerColor(4);   

   if(changeRange) {
     responseEta.SetMinimum(-1.0);
     responseEta.SetMaximum(1.0);
   }

   line->SetLineStyle(1);
   line->SetLineWidth(2);
   line->SetLineColor(2);

   canvas = new TCanvas( ("canvas"+ (string)(inHisto->GetName()) ).c_str(), ("Response"+ (string)(inHisto->GetName()) ).c_str(),1600, 900);

   
   temp->~TCanvas();  //destroy the TCanvas 
   
   canvas->cd();
   canvas->SetFillColor(0);

   gPad->SetGridx();
   gPad->SetGridy();
   respHisto->SetStats(0);
   respHisto->SetTitle(("Response "+ (string)(inHisto->GetName()) ).c_str());
   respHisto->Draw();
   responseEta.Draw("P");
   line->Draw();   
   respHisto->GetXaxis()->SetRangeUser(0,3);
   if(changeRange)  respHisto->GetYaxis()->SetRangeUser(-0.8,0.4);
   else respHisto->GetYaxis()->SetRangeUser(-0.4,0.1);
   respHisto->GetXaxis()->SetTitle("|#eta|");
   respHisto->GetYaxis()->SetTitle("(E_{cor}-E_{true})/E_{true}");


   // TF1*  f_eta = new TF1("etaFit", "[0]*(x - 1.5)*(x - 1.5) + [1]*(x - 1.5)*(x - 1.5)*(x - 1.5)*(x - 1.5) + [2]" , 1.5, 3.0);
   // f_eta->SetParameter(0,0.18);//,-0.16,0.0);
   // f_eta->SetParameter(1,-0.16);//,-0.16,0.0);
   // f_eta->SetParameter(2,0.0);//,-0.16,0.0);
   // responseEta.Fit("etaFit","","", 1.5, 2.85);
   
   //Spread
   //cout<<" Endcap spread and mean for "<<inHisto->GetName()<<endl;

   mR/=N;
   mR2/=N;
   // cout<<" mean = "<<mR<<endl;
   // cout<<" rms = "<<sqrt(mR2- pow(mR,2) )<<endl;


  char  cname[200];
  //sprintf(cname,  "%s_updatedCode.gif",inHisto->GetName());
  //canvas->Print(cname);
  sprintf(cname,  "%s_updatedCode.png",inHisto->GetName());
  canvas->Print(cname);

  respHisto->Write();
   //  line->Write();
  // resoHisto->Write();
   responseEta.Write();
   // resolution.Write();
   file2->cd();
   file2->Write();
   file2->Close();
   
}

void drawCompare(TGraph& response1, TGraph& response2, TGraph& resolution1, TGraph& resolution2)
{
   

   TCanvas* Compare = new TCanvas("Compare" ,"", 1000, 500);
   TH2F * respHisto = new TH2F("respHisto", "", 100, 0, 1000, 100, -1, 1);
   TH2F * resoHisto = new TH2F("resoHisto", "", 100, 0, 1000, 100, 0, 1);
   TLegend * legend1 = new TLegend(0.75, 0.75, 0.95, 0.9);
   TLegend * legend2 = new TLegend(0.75, 0.75, 0.95, 0.9);

   response1.SetMarkerColor(4);
   response1.SetMarkerStyle(22);
   response1.SetMarkerSize(0.8);

   resolution1.SetMarkerColor(4);
   resolution1.SetMarkerStyle(22);
   resolution1.SetMarkerSize(0.8);

   response2.SetMarkerColor(2);
   response2.SetMarkerStyle(22);
   response2.SetMarkerSize(0.8);

   resolution2.SetMarkerColor(2);
   resolution2.SetMarkerStyle(22);
   resolution2.SetMarkerSize(0.8);
   
   legend1->AddEntry(&response1, "Raw");
   legend1->AddEntry(&response2, "Corrected");
   legend2->AddEntry(&resolution1, "Raw");
   legend2->AddEntry(&resolution2, "Corrected");

   Compare->Divide(2,1);
  
   Compare->cd(1);
   gPad->SetGridx();
   gPad->SetGridy();
   respHisto->SetStats(0);
   respHisto->SetTitle("Response");
   respHisto->Draw();
   response1.Draw("P");
   response2.Draw("P");
   legend1->Draw();

   Compare->cd(2);
   gPad->SetGridx();
   gPad->SetGridy();
   resoHisto->SetStats(0);
   resoHisto->SetTitle("Resolution");
   resoHisto->Draw();
   resolution1.Draw("P");
   resolution2.Draw("P");
   legend2->Draw();

}

vector<float> assignvalues(vector<float> *pfcID_, vector<float> *Ecalenergy_, 
			   vector<float> *Hcalenergy_, vector<float> *dr) {

  vector<float> energies;
  float e = 0.0 , h = 0.0;
  for(unsigned ii = 0; ii < pfcID_->size(); ii++) {
    //cout<<" pfcID_:" << pfcID_->at(ii) << endl;
    if (old_logic) {
      if (pfcID_->at(ii) == 4 && dr->at(ii) < 0.2) e += Ecalenergy_->at(ii);
      if (pfcID_->at(ii) == 5 && dr->at(ii) < 0.4) h += Hcalenergy_->at(ii);
      
    }
    else {
      if (pfcID_->at(ii) == 5 && dr->at(ii) < 0.4) {
	e += Ecalenergy_->at(ii);
	h += Hcalenergy_->at(ii);
      }
    }
  }
  energies.push_back(e);
  energies.push_back(h);
  return energies;

}

//Recursive function to obtain all the Trees coming from different .root files
//located in a specific directory and in all its subdirectories.
//The function first searches the path provided as the second argument.
//If the element is a folder, it calls the function again but with the path of said folder,
//until obtaining all the Trees of all the subfolders of the original path.
void add_root_files_to_a_chain(TChain *chain, const char *path) {
     //int contador = 0;
	TSystemDirectory dir(path, path);
	TList *files = dir.GetListOfFiles();
	if (files) {
		TSystemFile *file;
		TString filename;
		TIter next(files);
		while ((file = (TSystemFile *)next())) {
			filename = file->GetName();
			if (!file->IsDirectory() && filename.EndsWith(".root")) {
				TString filepath = TString::Format("%s/%s", path, filename.Data());
				chain->Add(filepath);
				//  contador++;
				//  if(contador==1) break;
			} else if (file->IsDirectory() && TString(filename) != "." && TString(filename) != "..") {
				TString subpath = TString::Format("%s/%s", path, filename.Data());
				add_root_files_to_a_chain(chain, subpath);
			}
		}
	}
}


//Takes apart a TTree from a root file and puts the wanted information into 
//vectors. 
void getValuesFromTree(TTree* tree, vector<double>& ETrueEnergies, 
                       vector<double>& ecalEnergies, 
                       vector<double>& hcalEnergies, vector<double>& etas, 
                       vector<double>& phis, vector<double>& genE, vector<double>& genP, 
                       vector<double>& genEta, vector<double>& genPhi, vector<double>& trkP,
                       vector<double>& trkEta, vector<double>& trkPhi, vector<double>& momentum, 
                       vector<double>& ho, vector<double>& charge)//offline
{
   Float_t         true_;
   Float_t         p_;
   Float_t         ecal_;
   Float_t         hcal_;
   Float_t         ho_;
   Float_t         eta_;
   Float_t         phi_;
   Float_t         charge_;
   Float_t         genE_;
   Float_t         genP_;
   Float_t         genEta_;
   Float_t         genPhi_;
   Float_t         trkP_;
   Float_t         trkEta_;
   Float_t         trkPhi_;
   TBranch        *b_true; 
   TBranch        *b_p;   
   TBranch        *b_ecal;   
   TBranch        *b_hcal;
   TBranch        *b_ho;   
   TBranch        *b_eta;    
   TBranch        *b_phi;
   TBranch        *b_charge;
   TBranch        *b_genE;
   TBranch        *b_genP;
   TBranch        *b_genEta;
   TBranch        *b_genPhi;
   TBranch        *b_trkP;
   TBranch        *b_trkEta;
   TBranch        *b_trkPhi;   

   vector<float>        *pfcID_;
   vector<float>        *E_ecal_;
   vector<float>        *E_hcal_;
   vector<float>        *dr_;
   TBranch        *b_pfcID;   
   TBranch        *b_E_ecal;
   TBranch        *b_E_hcal;
   TBranch        *b_dr;

   pfcID_ = 0;
   E_ecal_ = 0;
   E_hcal_ = 0;
   dr_ = 0;

   //TChain* sTree = new TChain("s;1");
   //add_root_files_to_a_chain(sTree, "/eos/home-c/cmunozdi/step3_ana/");


   tree->SetMakeClass(1);
   
   tree->SetBranchStatus("*", 0);
   tree->SetBranchStatus("true", 1);
   tree->SetBranchStatus("p", 1);
   tree->SetBranchStatus("ecal", 1);
   tree->SetBranchStatus("hcal", 1);
   tree->SetBranchStatus("ho", 1);
   tree->SetBranchStatus("eta", 1);
   tree->SetBranchStatus("phi", 1);
   tree->SetBranchStatus("charge", 1);
   tree->SetBranchStatus("pfcID", 1);
   tree->SetBranchStatus("Eecal", 1);
   tree->SetBranchStatus("Ehcal", 1);
   tree->SetBranchStatus("dr", 1);
   tree->SetBranchStatus("genE", 1);
   tree->SetBranchStatus("genP", 1);
   tree->SetBranchStatus("genEta", 1);
   tree->SetBranchStatus("genPhi", 1);
   tree->SetBranchStatus("trkP", 1);
   tree->SetBranchStatus("trkEta", 1);
   tree->SetBranchStatus("trkPhi", 1);
   
   
   
   if(tree->GetBranchStatus("true"))
      tree->SetBranchAddress("true", &true_, &b_true);
   tree->SetBranchAddress("p", &p_, &b_p);
   tree->SetBranchAddress("ecal", &ecal_, &b_ecal);
   tree->SetBranchAddress("hcal", &hcal_, &b_hcal);
   tree->SetBranchAddress("eta", &eta_, &b_eta);
   tree->SetBranchAddress("phi", &phi_, &b_phi);
   tree->SetBranchAddress("ho", &ho_, &b_ho);
   tree->SetBranchAddress("charge", &charge_, &b_charge);
   tree->SetBranchAddress("pfcID", &pfcID_, &b_pfcID);
   tree->SetBranchAddress("Eecal", &E_ecal_, &b_E_ecal);
   tree->SetBranchAddress("Ehcal", &E_hcal_, &b_E_hcal);
   tree->SetBranchAddress("dr", &dr_, &b_dr);
   tree->SetBranchAddress("genE", &genE_, &b_genE);
   tree->SetBranchAddress("genP", &genP_, &b_genP);
   tree->SetBranchAddress("genEta", &genEta_, &b_genEta);
   tree->SetBranchAddress("genPhi", &genPhi_, &b_genPhi);
   tree->SetBranchAddress("trkP", &trkP_, &b_trkP);
   tree->SetBranchAddress("trkEta", &trkEta_, &b_trkEta);
   tree->SetBranchAddress("trkPhi", &trkPhi_, &b_trkPhi);

   double sigmaEcalHcal=1;
   long veto = 0 ;
   //int count = 0;
   bool flag[10] = {0,0,0,0,0,0,0,0,0,0};
	  for(int entry = 0; entry < tree->GetEntries(); entry++) {
       tree->GetEntry(entry);
       
       // if(ecal_<0.4) continue; //FIXME MM
       //if (fabs(eta_) < 2.4 && p_ == 0) continue;
       if((hadrons_eta_symbol==+1)&&(eta_<0)) continue;
       else if((hadrons_eta_symbol==-1)&&(eta_>0)) continue;

       //if (fabs(eta_) > 2.5 && (true_/cosh(eta_) < 5)) { continue;}
       //if (true_ < 48 || true_ > 52 ) continue;
       //if (true_ < 178 || true_ > 182 ) continue;
       //if (true_/cosh(eta_) < 20 ) continue;
       //////  HEP17 Veto
       //if (phi_ < -0.4 && phi_ > -1.0 && eta_ < 3.0 && eta_ > 1.5) { veto++; continue; } 
       //if (fabs(eta_) > 1.0)  continue;
       //if (true_>50 ) continue;
       /*if(tree->GetBranchStatus("true"))
	 ETrueEnergies.push_back(true_);
       else
	 ETrueEnergies.push_back(p_);
       if(pfcID_->size() != 0) {
	 vector<float> tmp = assignvalues(pfcID_, E_ecal_, E_hcal_, dr_);
	 ecalEnergies.push_back(tmp.at(0));
	 hcalEnergies.push_back(tmp.at(1));
       }
       else {
	 ecalEnergies.push_back(ecal_);
	 hcalEnergies.push_back(hcal_);
       }*/

       momentum.push_back(p_);
       ho.push_back(ho_);
       charge.push_back(charge_);

       //Gen and trk branches
        genE.push_back(genE_);
        genP.push_back(genP_);
        genEta.push_back(genEta_);
        genPhi.push_back(genPhi_);
        trkP.push_back(trkP_);
        trkEta.push_back(trkEta_);
        trkPhi.push_back(trkPhi_);
      
      double etrue, ecal, hcal, p_reco;
      bool saveData = false;

      if(tree->GetBranchStatus("true")){
          if(pfcID_->size() != 0) {
            vector<float> tmp = assignvalues(pfcID_, E_ecal_, E_hcal_, dr_);
            etrue = true_;
            p_reco = p_;
            ecal = tmp.at(0);
            hcal = tmp.at(1);
            //if((ecal+hcal)/etrue<=2){
              if(useP_reco) ETrueEnergies.push_back(p_reco);
              else ETrueEnergies.push_back(etrue);
              ecalEnergies.push_back(tmp.at(0));
	            hcalEnergies.push_back(tmp.at(1));
              saveData = true;
            //}
          }else{
            //if((ecal+hcal)/etrue<=2){
              if(useP_reco) ETrueEnergies.push_back(p_);
              else ETrueEnergies.push_back(true_);
              ecalEnergies.push_back(ecal_);
	            hcalEnergies.push_back(hcal_);
              saveData = true;
            //}
          }
      }else{
        if(pfcID_->size() != 0) {
            vector<float> tmp = assignvalues(pfcID_, E_ecal_, E_hcal_, dr_);
            etrue = true_;
            ecal = tmp.at(0);
            hcal = tmp.at(1);
            //if((ecal+hcal)/etrue<=2){
              ETrueEnergies.push_back(p_);
              ecalEnergies.push_back(tmp.at(0));
	            hcalEnergies.push_back(tmp.at(1));
              saveData = true;
            //}
          }else{
            //if((ecal+hcal)/etrue<=2){
              ETrueEnergies.push_back(p_);
              ecalEnergies.push_back(ecal_);
	            hcalEnergies.push_back(hcal_);
              saveData = true;
            //}
          }
      }





      if(saveData){
        etas.push_back(eta_);
        phis.push_back(phi_);


        if(fabs(eta_)<1.5) 
    sigmaEcalHcal = sqrt(0.08*0.08 + 1.04*1.04*(std::max((double)(ecal_ + hcal_), 1.0)));
        else
    sigmaEcalHcal = sqrt(0.04*0.04 + 1.80*1.80*(std::max((double)(ecal_ + hcal_), 1.0)));

        sigmas.push_back(sigmaEcalHcal);


        if(fabs(eta_) > 2.5 && ecalEnergies.back() != 0 && false) {
    cout<<"***************"<<endl;
    cout<<fabs(eta_)<<" "<<ecalEnergies.back()<<endl;
      
        }
      }
       //cout<< "**************" << endl;
       
       // cout<<" pfcID_.size(): " << pfcID_->size() << " Eecal->size(): " << E_ecal_->size()
       // 	   << " Ehcal->size(): " << E_hcal_->size() << " dr: " << dr_->size() << endl;
       // for(int ii = 0; ii < pfcID_->size() && entry < 20; ii++) {
       // 	 cout<<" pfcID_:" << pfcID_->at(ii) << endl;
       // }


       
       unsigned N = tree->GetEntriesFast();
       int frac = ((double)entry/N)*100;
       switch(frac) {
       case 10 : if (!flag[0]) { cout<<"10%"<<endl; flag[0] = 1; } break;
       case 20 : if (!flag[1]) { cout<<"20%"<<endl; flag[1] = 1; } break;
       case 30 : if (!flag[2]) { cout<<"30%"<<endl; flag[2] = 1; } break;
       case 40 : if (!flag[3]) { cout<<"40%"<<endl; flag[3] = 1; } break;
       case 50 : if (!flag[4]) { cout<<"50%"<<endl; flag[4] = 1; } break;
       case 60 : if (!flag[5]) { cout<<"60%"<<endl; flag[5] = 1; } break;
       case 70 : if (!flag[6]) { cout<<"70%"<<endl; flag[6] = 1; } break;
       case 80 : if (!flag[7]) { cout<<"80%"<<endl; flag[7] = 1; } break;
       case 90 : if (!flag[8]) { cout<<"90%"<<endl; flag[8] = 1; } break;
       case 99 : if (!flag[9]) { cout<<"100%"<<endl; flag[9] = 1; } break;
       default : break;
       
       }
       
     }

   cout<<" Entries "<<ecalEnergies.size()<<endl;
   cout<<" Vetoed events "<<veto<<endl;
   //exit(0);

}

// void getValuesFromTree(vector<double>& ETrueEnergies, //ONLINE
//     vector<double>& ecalEnergies, 
//     vector<double>& hcalEnergies, vector<double>& etas, 
//     vector<double>& phis) //TTree* sTree
// {
//   vector<float>*         true_=0;
//   vector<float>*         p_=0;
//   vector<float>*         eta_=0;
//   vector<float>*         phi_=0;
//   TBranch        *b_true; 
//   TBranch        *b_p;   
//   TBranch        *b_eta;    
//   TBranch        *b_phi;    

//   vector<float>        *pfcs_=0;
//   vector<float>        *E_ecal_=0;
//   vector<float>        *E_hcal_=0;
//   vector<float>        *dr_=0;
//   TBranch        *b_pfcs;   
//   TBranch        *b_E_ecal;
//   TBranch        *b_E_hcal;
//   TBranch        *b_dr;

//   vector<float>  *pfc_eta_=0;
//   TBranch        *b_pfc_eta_=0;   
//   vector<float>  *pfc_phi_=0;
//   TBranch        *b_pfc_phi_=0;   

//   // TFile *ftemp = new TFile("SinglePion/PFHadCalibration.root");
//   // TTree* sTree=(TTree*)ftemp->Get("pfHadCalibNTuple/Candidates");
//   TChain* sTree = new TChain ("pfHadCalibNTuple/Candidates");
//   //sTree->Add("/eos/user/d/dosite/online_04_07/E0to500/*.root");
//   //sTree->Add("/eos/user/d/dosite/online_06_07_no_whitelist/SinglePionE_2_200_run3_13_0_0/SinglePionGun_E0p2to200/crab_SinglePion_E_2to200_PFHadCalib_run3_2023/230706_092151/0000/*.root");
//       //sTree = (TTree*)chain;

//   //sTree->Add("/eos/user/d/dosite/online_06_07_no_whitelist/Combined_0_500/*.root");
//   sTree->Add("/eos/user/d/dosite/PFHadCalib_online_Conrado/*.root");



//   //sTree->SetMakeClass1);

//   sTree->SetBranchStatus("*", 0);
//   sTree->SetBranchStatus("true_energy", 1);
//   sTree->SetBranchStatus("true_eta", 1);
//   sTree->SetBranchStatus("true_phi", 1);
//   sTree->SetBranchStatus("pfc_trackRef_p", 1);
//   sTree->SetBranchStatus("pfc_id", 1);
//   sTree->SetBranchStatus("pfc_trackRef_eta", 1);
//   sTree->SetBranchStatus("pfc_trackRef_phi", 1);
//   sTree->SetBranchStatus("pfc_ecal", 1);
//   sTree->SetBranchStatus("pfc_hcal", 1);

//   sTree->SetBranchAddress("true_energy", &true_);
//   sTree->SetBranchAddress("true_eta", &eta_);
//   sTree->SetBranchAddress("true_phi", &phi_);
//   sTree->SetBranchAddress("pfc_trackRef_p", &p_);
//   sTree->SetBranchAddress("pfc_id", &pfcs_);
//   sTree->SetBranchAddress("pfc_trackRef_eta", &pfc_eta_);
//   sTree->SetBranchAddress("pfc_trackRef_phi", &pfc_phi_);
//   sTree->SetBranchAddress("pfc_ecal", &E_ecal_);
//   sTree->SetBranchAddress("pfc_hcal", &E_hcal_);

//   double sigmaEcalHcal=1;
//   long veto = 0 ;
//   //int count = 0;
//   bool flag[10] = {0,0,0,0,0,0,0,0,0,0};
//   float e_true, ecal, hcal, minDr, Dr, true_index;

//   cout << sTree->GetEntries() << endl;

//   //for( unsigned entry = 0; entry < std::min((unsigned)500000000,(unsigned)(sTree->GetEntriesFast()) ); entry++) {
//   for(int entry = 0; entry < sTree->GetEntries(); entry++) { 
//     sTree->GetEntry(entry);
//     //sTree->GetEntry(1266);

//     if(pfcs_->size() == 0 || true_->size() == 0) continue;

//     for(int i = 0; i < (int)pfcs_->size(); ++i){
//       minDr = 99.;  Dr = 999.; true_index = -1;
//       for(int j = 0; j < (int)true_->size(); ++j){
//         Dr = (pfc_eta_->at(i)-eta_->at(j))*(pfc_eta_->at(i)-eta_->at(j)) + (pfc_phi_->at(i)-phi_->at(j))*(pfc_phi_->at(i)-phi_->at(j));
//         if((minDr > Dr) && ((fabs(eta_->at(j)) <= 2.4 && p_->at(i) > 0) || fabs(eta_->at(j)) > 2.4)) { minDr = Dr; true_index = j; }
//       }

//       if(true_index >= 0) {
//         if(Dr < 1.){
//           e_true = true_->at(true_index); //Izm
//           //e_true = p_->at(true_index); //Izm
//           if(e_true > 500) continue;
//           ecal = E_ecal_->at(i);
//           hcal = E_hcal_->at(i);
//           ETrueEnergies.push_back(e_true); 
//           etas.push_back(eta_->at(true_index));
//           phis.push_back(phi_->at(true_index));
//           ecalEnergies.push_back(ecal);
//           hcalEnergies.push_back(hcal);
//           if(fabs(eta_->at(true_index))<1.5) sigmaEcalHcal = sqrt(0.08*0.08 + 1.04*1.04*(std::max((double)(ecal + hcal), 1.0)));
//           else																  	sigmaEcalHcal = sqrt(0.04*0.04 + 1.80*1.80*(std::max((double)(ecal + hcal), 1.0)));
//           sigmas.push_back(sigmaEcalHcal);

//         }
//       }
//     }

//     /*
//         if(sTree->GetBranchStatus("true"))
//         ETrueEnergies.push_back(true_);
//         else
//         ETrueEnergies.push_back(p_);
//     //if(pfcs_->size() != 0) {
//     if(pfcs_->size() == 0) {
//     //vector<float> tmp = assignvalues(pfcs_, E_ecal_, E_hcal_, dr_);
//     vector<float> tmp = assignvalues(pfcs_, E_ecal_, E_hcal_, 0);
//     ecalEnergies.push_back(tmp.at(0));
//     hcalEnergies.push_back(tmp.at(1));
//     }
//     else {
//     ecalEnergies.push_back(ecal_);
//     hcalEnergies.push_back(hcal_);
//     }
//     etas.push_back(eta_);
//     phis.push_back(phi_);


//     if(fabs(eta_)<1.5) 
//     sigmaEcalHcal = sqrt(0.08*0.08 + 1.04*1.04*(std::max((double)(ecal_ + hcal_), 1.0)));
//     else
//     sigmaEcalHcal = sqrt(0.04*0.04 + 1.80*1.80*(std::max((double)(ecal_ + hcal_), 1.0)));

//     sigmas.push_back(sigmaEcalHcal);


//     if(fabs(eta_) > 2.5 && ecalEnergies.back() != 0 && false) {
//     cout<<"***************"<<endl;
//     cout<<fabs(eta_)<<" "<<ecalEnergies.back()<<endl;

//     }
//     //cout<< "**************" << endl;

//     // cout<<" pfcID_.size(): " << pfcID_->size() << " Eecal->size(): " << E_ecal_->size()
//     // 	   << " Ehcal->size(): " << E_hcal_->size() << " dr: " << dr_->size() << endl;
//     // for(int ii = 0; ii < pfcID_->size() && entry < 20; ii++) {
//     // 	 cout<<" pfcID_:" << pfcID_->at(ii) << endl;
//     // }
//     */


//     unsigned N = sTree->GetEntriesFast();
//     int frac = ((double)entry/N)*100;
//     switch(frac) {
//       case 10 : if (!flag[0]) { cout<<"10%"<<endl; flag[0] = 1; } break;
//       case 20 : if (!flag[1]) { cout<<"20%"<<endl; flag[1] = 1; } break;
//       case 30 : if (!flag[2]) { cout<<"30%"<<endl; flag[2] = 1; } break;
//       case 40 : if (!flag[3]) { cout<<"40%"<<endl; flag[3] = 1; } break;
//       case 50 : if (!flag[4]) { cout<<"50%"<<endl; flag[4] = 1; } break;
//       case 60 : if (!flag[5]) { cout<<"60%"<<endl; flag[5] = 1; } break;
//       case 70 : if (!flag[6]) { cout<<"70%"<<endl; flag[6] = 1; } break;
//       case 80 : if (!flag[7]) { cout<<"80%"<<endl; flag[7] = 1; } break;
//       case 90 : if (!flag[8]) { cout<<"90%"<<endl; flag[8] = 1; } break;
//       case 99 : if (!flag[9]) { cout<<"100%"<<endl; flag[9] = 1; } break;
//       default : break;

//     }

//   }

// 		cout<<" Entries "<<ecalEnergies.size()<<endl;
// 		cout<<" Vetoed events "<<veto<<endl;
// 		//exit(0);

// 		}

///////////////////////////////////////////////////////////////////////////////
//This is the main of the macro. Everything that you want output must be added 
//in here. I have it so that all the variables that I use were defined above 
//since it looks neater.
///////////////////////////////////////////////////////////////////////////////
//void calibChris()
int main() 

{
   gROOT->Reset();
   gStyle->SetCanvasColor(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);


   InitBarrelAlpha();
   LoadOldThresholds();
   //LoadNewThresholds();

   gStyle->SetOptFit(0);

   //Open the file, get the tree and fill of the vectors of values you need.
   //inputFile = TFile::Open("/eos/home-c/cmunozdi/step3_ana/PGun_step3_RECO_1264_2_200_usingGTRun3v2_noPU/SinglePionGun_E0p2to200/crab_PGun_step3_RECO_1264_2_200_usingGTRun3v2_noPU_v4-v2/230522_081801/0000/step3_999.root");
   //inputFile = TFile::Open("IsolatedChargedHadronsFromQCD.root");
   //inputFile = TFile::Open("pfcalibTestTag_all.root");
   // inputFile = TFile::Open("IsolatedChargedHadronsFromMinBias.root");
    //sTree = (TTree*)inputFile->Get("s;1");
   

   
   TChain* chain= new TChain("s");
   
   /// HCAL issue SAMPLE - MC-v2
   //chain->Add("./input_sample/PGun_930pre1_JetMET_HCALScaleStudies_2_500.root");


   ///JME_checks sample
   //chain->Add("/home/shubham/work/PFCalibration/samples/10_0_2_NO_CUT/PGun_2_200_10_0_2_upgrade2018_NO_CUT.root");

   // chain->Add("/home/shubham/work/PFCalibration/samples/10_0_2_NO_CUT/PGun_2_500_10_0_2_upgrade2018_NO_CUT_new.root");

   //chain->Add("./root_files/PGun_2_500_10_6_0_pre2_UL2018.root");
   // chain->Add("./root_files/PGun_2_500_10_0_3_upgrade2018_ECAL_pfB.root");
   //  chain->Add("./home/bhumika/work/PFCalibration/for_10_0_2/calib_codes/PGun_2_500_10_0_2_upgrade2018_NO_CUT_new.root");
   //chain->Add("./root_files/PGun_Singlepion_10_6_0_UL2016.root");
   //   chain->Add("/Volumes/SSD/bhumi/work/Run3/rootfile/PGun_step3_RECO_1100_2021_Run3.root");
   //   chain->Add("/Volumes/SSD/bhumi/work/Run3/rootfile/PGun_step3_RECO_1248_2_500_usingGTEEleak.root");
   //   chain->Add("./rootfile/PGun_step3_RECO_1264_2_500_withPU.root");
//   chain->Add("/eos/home-c/cmunozdi/step3_ana/PGun_step3_RECO_1264_2_200_usingGTRun3v2_noPU/SinglePionGun_E0p2to200/crab_PGun_step3_RECO_1264_2_200_usingGTRun3v2_noPU_v4-v2/230522_081801/0000/*.root");
 
   //add_root_files_to_a_chain(chain, "/eos/user/c/cmunozdi/OFFLINE_NTUPLES/2024_Merged_3Attempt");//_merged/");//_proof/");//NewNTuplizerVersion/");/2024_Merged/OfflineNTuples_2024GT0_merged

   chain->Add("/eos/home-c/cmunozdi/OFFLINE_NTUPLES/OfflineNTuples_2024GT0_merged/0_500.root");//2024_Merged_3Attempt/rawFromNTuplizer/*.root");
   sTree = (TTree*)chain;
   cout<<"Reading input tree..."<<endl;
   getValuesFromTree(sTree, ETrueEnergies, ecalEnergies, 
                     hcalEnergies, etas, phis, genE, genP, 
                     genEta, genPhi, trkP, trkEta, trkPhi, 
                     momentums, ho_energies, charges);




   if (strcmp(_region_, "barrel") == 0) {
     _etaMin_ = 0.0;
     _etaMax_ = 1.5;
   }
   else if (strcmp(_region_, "EC_within_tracker") == 0 ) {
      _etaMin_ = 1.55;
     _etaMax_ = 2.5;
     //     _etaMin_ = 1.8;
     // _etaMax_ = 2.0;

   }
   
   else if (strcmp(_region_, "EC_outside_tracker") == 0 ) {
     _etaMin_ = 2.5;
     //_etaMax_ = 3.0; //update on 29 Aug 2019
     _etaMax_ = 3.0;//2.75;
   }
   
   else if (strcmp(_region_,  "Full") == 0 ) {
     _etaMin_ = 1.55;
     _etaMax_ = 3.0;
   }
 
   // cout<< " _region_: " << _region_<< " , (_region_ == EC_outside_tracker): " 
   //     << (strcmp(_region_ , "EC_outside_tracker") == 0) << " _etaMax_: " 
   //     << _etaMax_ << " ,_etaMin:_ " << _etaMin_ << endl;

   
   //Create all the ABC objects you need with increasing bin size
   //since there are fewer events at higher energies. 
   cout<<"Creating abc and alphabeta objects..."<<endl;
  
   BinsETrue.clear();
   BinsETrueEta.clear();

   for(double bin = 0.0; bin < 10.0; bin = bin + lBs)
     {
       barrelABCEcalHcal.push_back(new ABC(bin, bin + lBs, true));
       barrelABCEcal.push_back(new ABC(bin, bin + lBs, true));
       barrelABCHcal.push_back(new ABC(bin, bin + lBs, true));
       endcapABCEcalHcal.push_back(new ABC(bin, bin + lBs, false));
       endcapABCEcal.push_back(new ABC(bin, bin + lBs,false));
       endcapABCHcal.push_back(new ABC(bin, bin + lBs, false));
       BinsETrue.push_back(bin);
     }
   
   
   
   for(double bin = 10.0; bin < 100.0 ; bin = bin + mBs) //2
     {
       barrelABCEcalHcal.push_back(new ABC(bin, bin + mBs, true));
       barrelABCEcal.push_back(new ABC(bin, bin + mBs, true));
      barrelABCHcal.push_back(new ABC(bin, bin + mBs, true));
      endcapABCEcalHcal.push_back(new ABC(bin, bin + mBs, false));
      endcapABCEcal.push_back(new ABC(bin, bin + mBs,false));
      endcapABCHcal.push_back(new ABC(bin, bin + mBs, false));
      BinsETrue.push_back(bin);
   }
   
   
   for(double bin = 100.0; bin < sampleRangeHigh ; bin = bin + hBs) //10
   {
     barrelABCEcalHcal.push_back(new ABC(bin, bin + hBs, true));
     barrelABCEcal.push_back(new ABC(bin, bin + hBs, true));
     barrelABCHcal.push_back(new ABC(bin, bin + hBs, true));
     endcapABCEcalHcal.push_back(new ABC(bin, bin + hBs, false));
     endcapABCEcal.push_back(new ABC(bin, bin + hBs,false));
     endcapABCHcal.push_back(new ABC(bin, bin + hBs, false));  
     BinsETrue.push_back(bin);
   }
   BinsETrue.push_back( BinsETrue.back() + hBs );



   // cout<<"barrelABCEcalHcal size: "<<barrelABCEcalHcal.size()<<endl;
   // cout<<"barrelABCEcal size: "<<barrelABCEcal.size()<<endl;
   // cout<<"barrelABCHcal size: "<<barrelABCHcal.size()<<endl;
   /*for(int i = 0; i < barrelABCEcalHcal.size(); i++) {
     if((barrelABCEcalHcal.at(i)->getSize()) != 0 )
       cout<<"EcalHcal: found one!! at "<<i<<endl;
     if((barrelABCEcal.at(i)->getSize()) != 0 )
       cout<<"Ecal: found one!! at "<<i<<endl;
     if((barrelABCHcal.at(i)->getSize()) != 0 )
       cout<<"Hcal: found one!! at "<<i<<endl;

       }*/
   //cout<<"(barrelABCEcalHcal.at(300)->getBinLowEdge()) : "<<(barrelABCEcalHcal.at(300)->getBinLowEdge())<<endl;
   //cout<<"(barrelABCEcalHcal.at(300)->getBinHighEdge()) : "<<(barrelABCEcalHcal.at(300)->getBinHighEdge())<<endl;
   //cout<<"(barrelABCEcalHcal.at(0)->getETrue(0)) : "<<(barrelABCEcalHcal.at(0)->getETrue(0))<<endl;




   
   for(double bin = 0.0; bin < 10.0; bin = bin + lBs*RBE)
     {
       barrelAlphaBetaEcalHcal.push_back(new AlphaBeta(bin, bin + lBs*RBE, true));
       barrelAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + lBs*RBE, true));
       endcapAlphaBetaEcalHcal.push_back(new AlphaBeta(bin, bin + lBs*RBE, false));
       endcapAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + lBs*RBE, false));
       BinsETrueEta.push_back(bin);
       // barrelAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + lBs, true));
       // endcapAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + lBs, false));

     }
   for(double bin = 10.0; bin < 100.0 ; bin = bin + mBs*RBE)
     {
       barrelAlphaBetaEcalHcal.push_back(new AlphaBeta(bin, bin + mBs*RBE, true));
       barrelAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + mBs*RBE, true));
       endcapAlphaBetaEcalHcal.push_back(new AlphaBeta(bin, bin + mBs*RBE, false));
       endcapAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + mBs*RBE, false));
       BinsETrueEta.push_back(bin);
     //   barrelAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + mBs, true));
     //   endcapAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + mBs, false));
     }
   
   for(double bin = 100.0; bin < sampleRangeHigh ; bin = bin + hBs*RBE)
     {
       barrelAlphaBetaEcalHcal.push_back(new AlphaBeta(bin, bin + hBs*RBE, true));
       barrelAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + hBs*RBE, true));
       endcapAlphaBetaEcalHcal.push_back(new AlphaBeta(bin, bin + hBs*RBE, false));
       endcapAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + mBs*RBE, false));
       BinsETrueEta.push_back(bin);
       // barrelAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + hBs, true));
       // endcapAlphaBetaHcal.push_back(new AlphaBeta(bin, bin + hBs, false));
     }
   BinsETrueEta.push_back( BinsETrueEta.back() + hBs*RBE );
   
   //Fill all the ABC Objects with their respective events. They are all 
   //divided up into the six possible case ( (endcap or barrel)x(ecalhcal or 
   //ecal or hcal))
   

   TH1F* EcalSpectrum=new TH1F("EcalSpectrum","EcalSpectrum",1000,0,100);

   cout<<"Filling abc objects..."<<endl;
   for( unsigned bin = 0; bin < barrelABCEcal.size(); ++bin)
     {
       barrelABCEcalHcal[bin]->computeA(aEH);
       barrelABCEcal[bin]->computeA(aE);
       barrelABCHcal[bin]->computeA(aH);

       endcapABCEcalHcal[bin]->computeA(aEHe);
       endcapABCEcal[bin]->computeA(aEe);
       endcapABCHcal[bin]->computeA(aHe);
     }
   

   // cout<<"ETrueEnergies size: "<<ETrueEnergies.size()<<endl;
   // //cout<<"GetETrueBinEta(100): "<<GetETrueBinEta(99.9858)<<endl;
   // cout<<"GetETrueBinEta(100): "<<GetETrueBinEta(100)<<endl;

   //Filling ==============================
   {

     unsigned bin = 0;
     for(unsigned entry = 0; entry < ETrueEnergies.size(); entry++)
       {
	 etrue = ETrueEnergies[entry];
	 ecal = ecalEnergies[entry];
	 hcal = hcalEnergies[entry];
	 eta = etas[entry];

	 if(hcal == 0.0) continue;
	 if( etrue <1 ) continue;
	 if( etrue >sampleRangeHigh ) continue;
	 // if( etrue <10 ) continue;
	 // if( etrue >12 ) continue;

	 bin = GetETrueBin( etrue );

	 if( ecal > 0.0 && hcal > 0.0)
	   {
	     barrelABCEcalHcal[bin]->addEntry(etrue, ecal, hcal, eta);
	     endcapABCEcalHcal[bin]->addEntry(etrue, ecal, hcal, eta);
            
	     if(eta<1.3)
	       EcalSpectrum->Fill(ecal);

	   }
	 else if(ecal > 0.0)
	   {
	     barrelABCEcal[bin]->addEntry(etrue, ecal, hcal ,eta);
	     endcapABCEcal[bin]->addEntry(etrue, ecal, hcal ,eta);
	   }
	 else if(hcal > 0.0)
	   {
	     barrelABCHcal[bin]->addEntry(etrue, ecal, hcal, eta);
	     endcapABCHcal[bin]->addEntry(etrue, ecal, hcal, eta);
	   }
         
	 bin = GetETrueBinEta( etrue );

	 if(bin < barrelAlphaBetaEcalHcal.size())
	   {

	    
	     if( ecal > 0.0 && hcal >= 0.0 )
	       {
		 endcapAlphaBetaEcalHcal[bin]->addEntry(etrue, ecal, hcal, eta);
		 barrelAlphaBetaEcalHcal[bin]->addEntry(etrue, ecal, hcal, eta);
	       }
	     else {
	       endcapAlphaBetaHcal[bin]->addEntry(etrue, ecal, hcal, eta);
	       barrelAlphaBetaHcal[bin]->addEntry(etrue, ecal, hcal, eta);
	     }
	   }
       }
   }


   //   cout<<"#####################################################################"<<endl;
   //for(unsigned bin = 2; bin < barrelABCEcalHcal.size() - 1; ++bin) {
   // cout<<"barrelABCEcalHcal["<<bin<<"]->isEmptyInFitRange(): "<<barrelABCEcalHcal[bin]->isEmptyInFitRange()<<endl;
   //}
   //cout<<"#####################################################################"<<endl;

   //Filling ==============================  
   /*
   TFile* file=new TFile("output.root","recreate");
       EcalSpectrum->Write();
   file->Close();
   */
   
   //Compute the calibration constants along with their uncertainties for each
   //ETrue bin, then add their values to a Calibration object.
   cout<<"Computing a, b, c coefficients..."<<endl;

   for(unsigned bin = 2; bin < barrelABCEcalHcal.size() - 1; ++bin)
   {
      
      if(!barrelABCEcalHcal[bin]->isEmptyInFitRange())
      { 
         barrelABCEcalHcal[bin]->computeETrueAverage();
         barrelABCEcalHcal[bin]->computeETrueRMS();
         barrelABCEcalHcal[bin]->computeA(aEH);
         barrelABCEcalHcal[bin]->computeBC();
	 //exit(0);
      }
      
      if(!barrelABCEcal[bin]->isEmptyInFitRange())
      { 
         barrelABCEcal[bin]->computeETrueAverage();
         barrelABCEcal[bin]->computeETrueRMS();
         barrelABCEcal[bin]->computeA(aEH);
         barrelABCEcal[bin]->computeB();
      }
      if(!barrelABCHcal[bin]->isEmptyInFitRange())
      { 
         barrelABCHcal[bin]->computeETrueAverage();
         barrelABCHcal[bin]->computeETrueRMS();
         barrelABCHcal[bin]->computeA(aH);
         barrelABCHcal[bin]->computeC();
      }
      if(!endcapABCEcalHcal[bin]->isEmptyInFitRange())
      {
         endcapABCEcalHcal[bin]->computeETrueAverage();
         endcapABCEcalHcal[bin]->computeETrueRMS();
         endcapABCEcalHcal[bin]->computeA(aEHe);
         endcapABCEcalHcal[bin]->computeBC();
      }
      if(!endcapABCEcal[bin]->isEmptyInFitRange())
      {
         endcapABCEcal[bin]->computeETrueAverage();
         endcapABCEcal[bin]->computeETrueRMS();
         endcapABCEcal[bin]->computeA(aEe);
         endcapABCEcal[bin]->computeB();
      }
      if(!endcapABCHcal[bin]->isEmptyInFitRange())
      {
         endcapABCHcal[bin]->computeETrueAverage();
         endcapABCHcal[bin]->computeETrueRMS();
         endcapABCHcal[bin]->computeA(aHe);
         endcapABCHcal[bin]->computeC();
      }
      

      if(!barrelABCEcalHcal[bin]->isEmpty() && 
         barrelABCEcalHcal[bin]->getBinHighEdge() >
         barrelWithEcalHcalCalib->getETrueMax())
      {
         barrelWithEcalHcalCalib->setETrueMax(
            barrelABCEcalHcal[bin]->getBinHighEdge());
      }
      if(!barrelABCEcal[bin]->isEmpty() && 
         barrelABCEcal[bin]->getBinHighEdge() >
         barrelWithEcalCalib->getETrueMax())
      {
         barrelWithEcalCalib->setETrueMax(
            barrelABCEcal[bin]->getBinHighEdge());
      }
      if(!barrelABCHcal[bin]->isEmpty() && 
         barrelABCHcal[bin]->getBinHighEdge() >
         barrelWithHcalCalib->getETrueMax())
      {
         barrelWithHcalCalib->setETrueMax(
            barrelABCHcal[bin]->getBinHighEdge());
      }
      if(!endcapABCEcalHcal[bin]->isEmpty() && 
         endcapABCEcalHcal[bin]->getBinHighEdge() >
         endcapWithEcalHcalCalib->getETrueMax())
      {
         endcapWithEcalHcalCalib->setETrueMax(
            endcapABCEcalHcal[bin]->getBinHighEdge());
      }
      if(!endcapABCEcal[bin]->isEmpty() && 
         endcapABCEcal[bin]->getBinHighEdge() >
         endcapWithEcalCalib->getETrueMax())
      {
         endcapWithEcalCalib->setETrueMax(
            endcapABCEcal[bin]->getBinHighEdge());
      }
      if(!endcapABCHcal[bin]->isEmpty() && 
         endcapABCHcal[bin]->getBinHighEdge() >
         endcapWithHcalCalib->getETrueMax())
      {
         endcapWithHcalCalib->setETrueMax(
            endcapABCHcal[bin]->getBinHighEdge());
      }
 

      barrelWithEcalHcalCalib->addGraphPoints(barrelABCEcalHcal[bin]); 
      barrelWithEcalCalib->addGraphPoints(barrelABCEcal[bin]); 
      barrelWithHcalCalib->addGraphPoints(barrelABCHcal[bin]); 
      endcapWithEcalHcalCalib->addGraphPoints(endcapABCEcalHcal[bin]); 
      endcapWithEcalCalib->addGraphPoints(endcapABCEcal[bin]); 
      endcapWithHcalCalib->addGraphPoints(endcapABCHcal[bin]);                 


   }
   
   cout<<"Fitting a, b, c coefficients..."<<endl;
   //Initialize all the ABC graphs in the calibration objects.
   barrelWithEcalHcalCalib->initializeGraphs("abc");
   barrelWithEcalCalib->initializeGraphs("abc");
   barrelWithHcalCalib->initializeGraphs("abc");   
  
   endcapWithEcalHcalCalib->initializeGraphs("abc");
   endcapWithEcalCalib->initializeGraphs("abc");
   endcapWithHcalCalib->initializeGraphs("abc");   

   //Define the functions that you will fit your ABC calibration constants to.
   functionBarrelEcalHcalA = new TF1("functionBarrelEcalHcalA","[0]", 0, sampleRangeHigh);
   // functionBarrelEcalHcalB = new TF1("functionBarrelEcalHcalB","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])", 0, 1000);

   functionBarrelEcalHcalB = new TF1("functionBarrelEcalHcalB","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, sampleRangeHigh);//-[8]*exp(-x^[9]/[10])", 0, 1000);//-[8]*exp(-x^[9]/[10]))", 0, 1000);
   //  functionBarrelEcalHcalB = new TF1("functionBarrelEcalHcalB","[0]+((([1]+([2]/(x^[5])))*exp(-(x^[4]/[3]))))", 0, 1000);
   //functionBarrelEcalHcalC = new TF1("functionBarrelEcalHcalC","[0]+(([1]+([2]/sqrt(x)))*exp(-(x^[4]/[3])))",0,1000); //[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])", 0, 1000);
   functionBarrelEcalHcalC = new TF1("functionBarrelEcalHcalC","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",0,sampleRangeHigh);
  
   functionEndcapEcalHcalA = new TF1("functionEndcapEcalHcalA","[0]", 0, sampleRangeHigh);
   // functionEndcapEcalHcalC = new TF1("functionEndcapEcalHcalC","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])", 0, 1000);
   functionEndcapEcalHcalB = new TF1("functionEndcapEcalHcalB","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, sampleRangeHigh);//-[8]*exp(-x^[9]/[10])", 0, 1000); 
   //functionEndcapEcalHcalB = new TF1("functionEndcapEcalHcalB","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[4]/[3]))))", 0, 1000);
   //functionEndcapEcalHcalC = new TF1("functionEndcapEcalHcalC","([0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5])))))", 0, 1000);
   
   //Offline
   functionEndcapEcalHcalC = new TF1("functionEndcapEcalHcalC","[0]+([4]*(x-[5])*exp(-(x*[7])))+(([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))", 0, sampleRangeHigh);
  //Online
   //functionEndcapEcalHcalC = new TF1("functionEndcapEcalHcalC","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, 1000);


   functionBarrelHcalA = new TF1("functionBarrelHcalA","[0]", 0, sampleRangeHigh);
   functionBarrelHcalB = new TF1("functionBarrelHcalB","[0]", 0, sampleRangeHigh);
   // functionBarrelHcalC = new TF1("functionBarrelHcalC","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])", 0, 1000);
   functionBarrelHcalC = new TF1("functionBarrelHcalC","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, sampleRangeHigh);//-[8]*exp(-x^[9]/[10])", 0, 1000);
   //spandey
   //functionBarrelHcalC = new TF1("functionBarrelHcalC","1.03*([0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5])))))", 0, 1000);
  
   functionEndcapHcalA = new TF1("functionEndcapHcalA","[0]", 0, sampleRangeHigh);
   functionEndcapHcalB = new TF1("functionEndcapHcalB","[0]", 0, sampleRangeHigh);
   functionEndcapHcalC = new TF1("functionEndcapHcalC","([0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5])))))", 0, sampleRangeHigh); //[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])


   if(freezeparameters) {
    if(!Parameters24Above500GeV){//Parameters to fix for 2024 pfhc
      functionBarrelEcalHcalA->FixParameter(0, aEH);

        functionEndcapEcalHcalA->FixParameter(0, aEHe);

      functionBarrelHcalA->FixParameter(0, aH);
      functionBarrelHcalB->FixParameter(0, 0.0);

      functionEndcapHcalA->FixParameter(0, aHe);
      functionEndcapHcalB->FixParameter(0, 0.0);
        //2024 current parameters
        functionBarrelEcalHcalB->FixParameter(0, 14.9081);
        functionBarrelEcalHcalB->FixParameter(1, -92.531);
        functionBarrelEcalHcalB->FixParameter(2, -586.723);
        functionBarrelEcalHcalB->FixParameter(3, 0.281367);
        functionBarrelEcalHcalB->FixParameter(4, 13.0608);
        functionBarrelEcalHcalB->FixParameter(5, 0.450973);
        functionBarrelEcalHcalB->FixParameter(6, 0.03687);
        functionBarrelEcalHcalB->FixParameter(7, -0.583429);


        functionBarrelEcalHcalC->FixParameter(0, 2.414);
        functionBarrelEcalHcalC->FixParameter(1, -2.99257);
        functionBarrelEcalHcalC->FixParameter(2, -3.10022);
        functionBarrelEcalHcalC->FixParameter(3, 2.4884);
        functionBarrelEcalHcalC->FixParameter(4, 1.49647);
        functionBarrelEcalHcalC->FixParameter(5, 0.0591164);
        functionBarrelEcalHcalC->FixParameter(6, 0.401639);
        functionBarrelEcalHcalC->FixParameter(7, -0.848485);


        functionBarrelHcalC->FixParameter(0, 10.7719);
        functionBarrelHcalC->FixParameter(1, 6.36096);
        functionBarrelHcalC->FixParameter(2, -23.8131);
        functionBarrelHcalC->FixParameter(3, 1.77669);
        functionBarrelHcalC->FixParameter(4, 12.6614);
        functionBarrelHcalC->FixParameter(5, 0.722518);
        functionBarrelHcalC->FixParameter(6, 0.0447024);
        functionBarrelHcalC->FixParameter(7, -0.636043);


        functionEndcapEcalHcalB->FixParameter(0, 37.7103);
        functionEndcapEcalHcalB->FixParameter(1, -318.409);
        functionEndcapEcalHcalB->FixParameter(2, -1095.14);
        functionEndcapEcalHcalB->FixParameter(3, 0.298041);
        functionEndcapEcalHcalB->FixParameter(4, 23.6435);
        functionEndcapEcalHcalB->FixParameter(5, 0.326262);
        functionEndcapEcalHcalB->FixParameter(6, -0.0049515);
        functionEndcapEcalHcalB->FixParameter(7, -0.631655);


        functionEndcapEcalHcalC->FixParameter(0, -2.24814);
        functionEndcapEcalHcalC->FixParameter(1, 3.15142);
        functionEndcapEcalHcalC->FixParameter(2, 3.76944);
        functionEndcapEcalHcalC->FixParameter(3, 1.06815);
        functionEndcapEcalHcalC->FixParameter(4, 0.0726289);
        functionEndcapEcalHcalC->FixParameter(5, 24.8456);
        functionEndcapEcalHcalC->FixParameter(6, -0.609985);
        functionEndcapEcalHcalC->FixParameter(7, 0.0968119);


        functionEndcapHcalC->FixParameter(0, 1.63161);
        functionEndcapHcalC->FixParameter(1, 6.36717);
        functionEndcapHcalC->FixParameter(2, -33.0955);
        functionEndcapHcalC->FixParameter(3, 0.501949);
        functionEndcapHcalC->FixParameter(4, 0.856495);
        functionEndcapHcalC->FixParameter(5, 0.0255776);
        functionEndcapHcalC->FixParameter(6, 0.0809049);
        functionEndcapHcalC->FixParameter(7, -1.41804);
      }
      else{
        functionBarrelEcalHcalA->FixParameter(0, aEH);

        functionBarrelEcalHcalB->FixParameter(0,13.033);
        functionBarrelEcalHcalB->FixParameter(1,87.2668);
        functionBarrelEcalHcalB->FixParameter(2,-699.24);
        functionBarrelEcalHcalB->FixParameter(3,0.304668);
        functionBarrelEcalHcalB->FixParameter(4,13.8154);
        functionBarrelEcalHcalB->FixParameter(5,0.266523);
        functionBarrelEcalHcalB->FixParameter(6,0.0171292);
        functionBarrelEcalHcalB->FixParameter(7,-0.725741);

        functionBarrelEcalHcalC->FixParameter(0,1.75412);
        functionBarrelEcalHcalC->FixParameter(1,-0.413335);
        functionBarrelEcalHcalC->FixParameter(2,-2.08127);
        functionBarrelEcalHcalC->FixParameter(3,126.351);
        functionBarrelEcalHcalC->FixParameter(4,0.770695);
        functionBarrelEcalHcalC->FixParameter(5,0.00404635);
        functionBarrelEcalHcalC->FixParameter(6,1.12044);
        functionBarrelEcalHcalC->FixParameter(7,-1.38901);

        functionBarrelHcalC->FixParameter(0,11.384);
        functionBarrelHcalC->FixParameter(1,23.1406);
        functionBarrelHcalC->FixParameter(2,-28.9497);
        functionBarrelHcalC->FixParameter(3,0.97494);
        functionBarrelHcalC->FixParameter(4,19.647);
        functionBarrelHcalC->FixParameter(5,1.6272);
        functionBarrelHcalC->FixParameter(6,-0.0108519);
        functionBarrelHcalC->FixParameter(7,-0.423267);

        functionEndcapEcalHcalA->FixParameter(0, aEHe);

        functionEndcapEcalHcalB->FixParameter(0,29.5353);
        functionEndcapEcalHcalB->FixParameter(1,-327.967);
        functionEndcapEcalHcalB->FixParameter(2,-830.992);
        functionEndcapEcalHcalB->FixParameter(3,0.298631);
        functionEndcapEcalHcalB->FixParameter(4,15.6682);
        functionEndcapEcalHcalB->FixParameter(5,0.242381);
        functionEndcapEcalHcalB->FixParameter(6,-0.00308457);
        functionEndcapEcalHcalB->FixParameter(7,-0.654633);

        functionEndcapEcalHcalC->FixParameter(0,-0.180939);
        functionEndcapEcalHcalC->FixParameter(1,32.2565);
        functionEndcapEcalHcalC->FixParameter(2,4245.12);
        functionEndcapEcalHcalC->FixParameter(3,0.121009);
        functionEndcapEcalHcalC->FixParameter(4,0.0736027);
        functionEndcapEcalHcalC->FixParameter(5,18.954);
        functionEndcapEcalHcalC->FixParameter(6,-0.0734254);
        functionEndcapEcalHcalC->FixParameter(7,0.0771871);

        functionBarrelHcalA->FixParameter(0, aH);
        functionBarrelHcalB->FixParameter(0, 0.0);
        functionEndcapHcalA->FixParameter(0, aHe);
        functionEndcapHcalB->FixParameter(0, 0.0);

        functionEndcapHcalC->FixParameter(0,2.09235);
        functionEndcapHcalC->FixParameter(1,0.328094);
        functionEndcapHcalC->FixParameter(2,-3.98841);
        functionEndcapHcalC->FixParameter(3,34.8701);
        functionEndcapHcalC->FixParameter(4,1.17097);
        functionEndcapHcalC->FixParameter(5,0.030813);
        functionEndcapHcalC->FixParameter(6,0.66908);
        functionEndcapHcalC->FixParameter(7,-1.26556);
      }

        
   }
   else {
      functionBarrelEcalHcalA->FixParameter(0, aEH);
      functionBarrelEcalHcalB->SetParameters(13.033, 87.2666, -699.239, 0.304668, 13.8154, 0.266522, 0.0171293, -0.725742);
      functionBarrelEcalHcalC->SetParameters(1.75412, -0.413337, -2.08127, 126.351, 0.770696, 0.00404633, 1.12044, -1.38901);
      functionBarrelHcalC->SetParameters(11.384, 23.1406, -28.9496, 0.97494, 19.647, 1.6272, -0.0108519, -0.423267);
      functionEndcapEcalHcalA->FixParameter(0, aEHe);
      functionEndcapEcalHcalB->SetParameters(29.5353, -327.967, -830.992, 0.298631, 15.6682, 0.242381, -0.00308458, -0.654633);
      functionEndcapEcalHcalC->SetParameters(-0.139322, 19.7095, 2354.6, 0.12935, 0.0736876, 18.8147, -0.0801738, 0.0771942);
      functionBarrelHcalA->FixParameter(0, aH);
      functionBarrelHcalB->FixParameter(0, 0.0);
      functionEndcapHcalA->FixParameter(0, aHe);
      functionEndcapHcalB->FixParameter(0, 0.0);
      functionEndcapHcalC->SetParameters(2.09235, 0.328094, -3.98841, 34.8701, 1.17097, 0.0308129, 0.669081, -1.26556);
   }
   barrelWithEcalHcalCalib->fitAsToFunction(functionBarrelEcalHcalA);
   //Printing parameters:
   barrelWithEcalHcalCalib->fitBsToFunction(functionBarrelEcalHcalB);
   barrelWithEcalHcalCalib->fitBsToFunction();

   barrelWithEcalHcalCalib->fitCsToFunction(functionBarrelEcalHcalC);

   barrelWithEcalHcalCalib->fitCsToFunction();

   barrelWithEcalHcalCalib->fitCsToFunction();

   endcapWithEcalHcalCalib->fitAsToFunction(functionEndcapEcalHcalA);

   // cout<<"********************************************"<<endl;
   // cout<<"Fit Parameters, functionEndcapEcalHcalB, First"<<endl;
   // for ( unsigned i = 0; i < 10; ++i ) {
   //   double barrelEcalHcalC_spandey = functionEndcapEcalHcalB->GetParameter(i);
   //   if ( barrelEcalHcalC_spandey != 0. )
   //     cout<<"  functionEndcapEcalHcalB,Parameter("<<i<<","<<barrelEcalHcalC_spandey<<");"<<endl;
   // }
   endcapWithEcalHcalCalib->fitBsToFunction(functionEndcapEcalHcalB);

   // cout<<"********************************************"<<endl;
   // cout<<"Fit Parameters, functionEndcapEcalHcalB, First"<<endl;
   // for ( unsigned i = 0; i < 10; ++i ) {
   //   double barrelEcalHcalC_spandey = functionEndcapEcalHcalB->GetParameter(i);
   //   if ( barrelEcalHcalC_spandey != 0. )
   //     cout<<"  functionEndcapEcalHcalB,Parameter("<<i<<","<<barrelEcalHcalC_spandey<<");"<<endl;
   // }
      
   endcapWithEcalHcalCalib->fitBsToFunction();

   // cout<<"********************************************"<<endl;
   // cout<<"Fit Parameters, functionEndcapEcalHcalB, First"<<endl;
   // for ( unsigned i = 0; i < 10; ++i ) {
   //   double barrelEcalHcalC_spandey = functionEndcapEcalHcalB->GetParameter(i);
   //   if ( barrelEcalHcalC_spandey != 0. )
   //     cout<<"  functionEndcapEcalHcalB,Parameter("<<i<<","<<barrelEcalHcalC_spandey<<");"<<endl;
   // }

   //   cout<<"********************************************"<<endl;
   endcapWithEcalHcalCalib->fitCsToFunction(functionEndcapEcalHcalC);
   endcapWithEcalHcalCalib->fitCsToFunction();
   endcapWithEcalHcalCalib->fitCsToFunction();
   //   cout<<"Fit check 1\n";
   barrelWithHcalCalib->fitAsToFunction(functionBarrelHcalA);
   barrelWithHcalCalib->fitBsToFunction(functionBarrelHcalB);
   //   cout<<"Fit check 1.a\n";
   barrelWithHcalCalib->fitBsToFunction();
   //   cout<<"Fit check 1.b\n";
   barrelWithHcalCalib->fitCsToFunction(functionBarrelHcalC);
   barrelWithHcalCalib->fitCsToFunction();
   barrelWithHcalCalib->fitCsToFunction();
   //   cout<<"Fit check 2\n";
   endcapWithHcalCalib->fitAsToFunction(functionEndcapHcalA);
   endcapWithHcalCalib->fitBsToFunction(functionEndcapHcalB);
   endcapWithHcalCalib->fitBsToFunction();

      //cout<<"1 FITTING ENDCAP HCAL C#######"<<endl;
   endcapWithHcalCalib->fitCsToFunction(functionEndcapHcalC);
      // cout<<"2 FITTING ENDCAP HCAL C#######"<<endl;
   endcapWithHcalCalib->fitCsToFunction();
      // cout<<"3 FITTING ENDCAP HCAL C#######"<<endl;
   endcapWithHcalCalib->fitCsToFunction();
   //   cout<<"Fit check 3\n";
   
   //exit(0);
   //Here we fill up the AlphaBeta objects, compute alpha and beta, then add 
   //them to the Calibration objects. 
   cout<<"Computing alpha and beta coefficients..."<<endl;
   for(unsigned bin = 2; bin < barrelAlphaBetaEcalHcal.size() - 1; bin++)
   {
     for(unsigned entry = 0; entry < barrelAlphaBetaEcalHcal[bin]->getSize(); entry++)
      {
         
         etrue = barrelAlphaBetaEcalHcal[bin]->getETrue(entry);
         ecal = barrelAlphaBetaEcalHcal[bin]->getEcal(entry);
         hcal = barrelAlphaBetaEcalHcal[bin]->getHcal(entry);
         bpar = 1.0;
         cpar = 1.0;
         

         if(ecal > 0 && hcal > 0)
         {
            bpar = barrelWithEcalHcalCalib->getFunctionB()->Eval(etrue);
            cpar = barrelWithEcalHcalCalib->getFunctionC()->Eval(etrue);
         }
         else if(ecal > 0)
            bpar = barrelWithEcalHcalCalib->getFunctionB()->Eval(etrue);
       
         
         barrelAlphaBetaEcalHcal[bin]->correctEcal(entry, bpar);
         barrelAlphaBetaEcalHcal[bin]->correctHcal(entry, cpar);
      }

   
     for(unsigned entry = 0; entry < barrelAlphaBetaHcal[bin]->getSize(); entry++)
       {
         
	 etrue = barrelAlphaBetaHcal[bin]->getETrue(entry);
	 ecal = barrelAlphaBetaHcal[bin]->getEcal(entry);
	 hcal = barrelAlphaBetaHcal[bin]->getHcal(entry);
	 bpar = 1.0;
	 cpar = 1.0;
	 
	 if(hcal > 0 && ecal==0)
	   cpar = barrelWithHcalCalib->getFunctionC()->Eval(etrue);
      
	 // if(etrue<10)
	 //   cout<<etrue<<"   "<<cpar<<endl;

	 barrelAlphaBetaHcal[bin]->correctEcal(entry, bpar);
         barrelAlphaBetaHcal[bin]->correctHcal(entry, cpar);

       }     

     for(unsigned entry = 0; entry < endcapAlphaBetaEcalHcal[bin]->getSize(); entry++)
	{

	  etrue = endcapAlphaBetaEcalHcal[bin]->getETrue(entry);
	  ecal = endcapAlphaBetaEcalHcal[bin]->getEcal(entry);
	  hcal = endcapAlphaBetaEcalHcal[bin]->getHcal(entry);
	  bpar = 1.0;
	  cpar = 1.0;
         
	  if(ecal > 0 && hcal > 0)
	    {
	      bpar = endcapWithEcalHcalCalib->getFunctionB()->Eval(etrue);
	      cpar = endcapWithEcalHcalCalib->getFunctionC()->Eval(etrue);
	    }
	  else if(ecal > 0)
            bpar = endcapWithEcalHcalCalib->getFunctionB()->Eval(etrue);
	  
	  endcapAlphaBetaEcalHcal[bin]->correctEcal(entry, bpar);
	  endcapAlphaBetaEcalHcal[bin]->correctHcal(entry, cpar);
	}

      for(unsigned entry = 0; entry < endcapAlphaBetaHcal[bin]->getSize(); entry++)
	{

	  etrue = endcapAlphaBetaHcal[bin]->getETrue(entry);
	  ecal = endcapAlphaBetaHcal[bin]->getEcal(entry);
	  hcal = endcapAlphaBetaHcal[bin]->getHcal(entry);
	  bpar = 1.0;
	  cpar = 1.0;
         
	  if(ecal == 0 && hcal > 0)
	    {
	      cpar = endcapWithHcalCalib->getFunctionC()->Eval(etrue);
	    }
	  endcapAlphaBetaHcal[bin]->correctEcal(entry, bpar);
	  endcapAlphaBetaHcal[bin]->correctHcal(entry, cpar);
	}
      
      
      barrelAlphaBetaEcalHcal[bin]->computeSigmaEcalHcal();
      barrelAlphaBetaEcalHcal[bin]->computeETrueAverage();
      barrelAlphaBetaEcalHcal[bin]->computeETrueRMS();

      barrelAlphaBetaHcal[bin]->computeSigmaEcalHcal();
      barrelAlphaBetaHcal[bin]->computeETrueAverage();
      barrelAlphaBetaHcal[bin]->computeETrueRMS();
      
      endcapAlphaBetaEcalHcal[bin]->computeSigmaEcalHcal();
      endcapAlphaBetaEcalHcal[bin]->computeETrueAverage();
      endcapAlphaBetaEcalHcal[bin]->computeETrueRMS();

      endcapAlphaBetaHcal[bin]->computeSigmaEcalHcal();
      endcapAlphaBetaHcal[bin]->computeETrueAverage();
      endcapAlphaBetaHcal[bin]->computeETrueRMS();

      if(barrelAlphaBetaEcalHcal[bin]->computeAlphaBeta())
      {
	barrelWithEcalHcalCalib->addGraphPoints(barrelAlphaBetaEcalHcal[bin]);
	barrelWithEcalCalib->addGraphPoints(barrelAlphaBetaEcalHcal[bin]);
      }
      if(barrelAlphaBetaHcal[bin]->computeAlphaBeta()) {
	barrelWithHcalCalib->addGraphPoints(barrelAlphaBetaHcal[bin]);
      }

      if(endcapAlphaBetaEcalHcal[bin]->computeAlphaBeta())
	{
	  endcapWithEcalHcalCalib->addGraphPoints(endcapAlphaBetaEcalHcal[bin]);
	  endcapWithEcalCalib->addGraphPoints(endcapAlphaBetaEcalHcal[bin]);
	}

     
      if(endcapAlphaBetaHcal[bin]->computeAlphaBeta()) //FIXME
	{
	  endcapWithHcalCalib->addGraphPoints(endcapAlphaBetaHcal[bin]);
	}
   }
   
   cout<<"Fitting alpha, beta coefficients..."<<endl;
   barrelWithEcalHcalCalib->initializeGraphs("alphabeta");
   barrelWithEcalCalib->initializeGraphs("alphabeta");
   barrelWithHcalCalib->initializeGraphs("alphabeta");   
   endcapWithEcalHcalCalib->initializeGraphs("alphabeta");
   endcapWithEcalCalib->initializeGraphs("alphabeta");
   endcapWithHcalCalib->initializeGraphs("alphabeta");   

    //Offline
    /*
   functionBarrelAlphaEcalHcal = new TF1("functionBarrelAlphaEcalHcal","[0]+[1]*exp(-x*[3]/[2])", 0, 1000);//exp(-x*[2])/(x^[2]
   //trial1
	 functionBarrelBetaEcalHcal = new TF1("functionBarrelBetaEcalHcal","[0]+((([1]+([2]/(x^[5])))*exp(-(x^[4]/[3]))))", 0, 1000);          
   functionBarrelAlphaHcal = new TF1("functionBarrelAlphaHcal","[0]+[1]*x", 0, 1000);
      //UL2018                                                                                                                                                                
   // functionBarrelAlphaHcal = new TF1("functionBarrelAlphaHcal","[0]+[1]*exp(-x/[2])", 0, 1000);// +[1]*exp(-x/[2])
   functionBarrelBetaHcal = new TF1("functionBarrelBetaHcal","[0]+[1]*exp(-x/[2])", 0, 1000);


   //faEtaEndCapEH
   //UL2018
   //   functionEndcapAlphaEcalHcal = new TF1("functionEndcapAlphaEcalHcal","[0]+([1]*x^[2]*exp(-x))", 0, 1000);
   //   functionEndcapBetaEcalHcal = new TF1("functionEndcapBetaEcalHcal","[0]+[1]*x^[3]*exp(-x/[2])",0,1000); //+[3]*[3]*exp(-x*x/([4]*[4]))
   //trial1
   //   functionEndcapAlphaEcalHcal = new TF1("functionEndcapAlphaEcalHcal","[0]+[1]*exp(-x/[2])+[3]*[3]*exp(-x*x/([4]*[4]))", 0, 1000);
   //   functionEndcapAlphaEcalHcal = new TF1("functionEndcapAlphaEcalHcal","[0]+([1]*x^[2]*exp(-x/[3]))", 0, 1000);     
   functionEndcapAlphaEcalHcal = new TF1("functionEndcapAlphaEcalHcal","[0]+([1]*x)", 0, 1000);  

   //   functionEndcapAlphaEcalHcal = new TF1("functionEndcapAlphaEcalHcal","
   //functionEndcapAlphaEcalHcal = new TF1("functionEndcapAlphaEcalHcal","[0]+((([1]+([2]/(x^[5])))*exp(-(x^[4]/[3]))))",0,1000);
 	 functionEndcapBetaEcalHcal = new TF1("functionEndcapBetaEcalHcal","[0]+((([1]+([2]/(x^[5])))*exp(-(x^[4]/[3]))))",0,1000);
   
   functionEndcapAlphaHcal = new TF1("functionEndcapAlphaHcal","[0]+[1]*x", 0, 1000);// +[1]*exp(-x/[2])
   //   functionEndcapAlphaHcal = new TF1("functionEndcapAlphaHcal","[0]+[1]*x+[3]*exp(-x/[2])", 0, 1000);// +[1]*exp(-x/[2])

   // functionEndcapAlphaHcal = new TF1("functionEndcapAlphaHcal","[0]+[1]*exp(-x/[2])+[3]*[3]*exp(-x*x/([4]*[4]))", 0, 1000);// +[1]*exp(-x/[2])
   // functionEndcapBetaHcal = new TF1("functionEndcapBetaHcal","[0]+[1]*exp(-x/[2])+[3]*[3]*exp(-x*x/([4]*[4]))",0,1000);
   //UL2018
   //   functionEndcapAlphaHcal = new TF1("functionEndcapAlphaHcal","[0]+[1]*(x^[3])*exp(-x/[2])", 0, 1000);// +[1]*exp(-x/[2])
	functionEndcapAlphaHcal = new TF1("functionEndcapAlphaHcal","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[4]/[3]))))", 0, 1000);
	functionEndcapBetaHcal = new TF1("functionEndcapBetaHcal","[0]+[1]*x*exp(-x/[2])",0,1000);
   //functionEndcapBetaHcal = new TF1("functionEndcapBetaHcal","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[4]/[3]))))",0,1000);
   functionEndcapGammaHcal = new TF1("functionEndcapGammaHcal","[0]+[1]*exp(-x/[2])+[3]*[3]*exp(-x*x/([4]*[4]))",0,1000);
  */
  //Online			
      functionEndcapGammaHcal = new TF1("functionEndcapGammaHcal","[0]+[1]*exp(-x/[2])+[3]*[3]*exp(-x*x/([4]*[4]))",0,sampleRangeHigh);
			functionBarrelAlphaEcalHcal = new TF1("functionBarrelAlphaEcalHcal","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, sampleRangeHigh);
			functionBarrelBetaEcalHcal = new TF1("functionBarrelBetaEcalHcal","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, sampleRangeHigh);
      //functionBarrelBetaEcalHcal = new TF1("functionBarrelBetaEcalHcal","[0] + [1]*( [2]*x+[3]/( [4]*sqrt([5]*x) ) )*exp([6]*x) + [7]*exp([8]*x*x)", 0, 1000); 
			functionBarrelAlphaHcal = new TF1("functionBarrelAlphaHcal","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, sampleRangeHigh);
			functionBarrelBetaHcal = new TF1("functionBarrelBetaHcal","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, sampleRangeHigh);//+[8]*exp(-x^[9]/[10])", 0, 1000);
			functionEndcapAlphaEcalHcal = new TF1("functionEndcapAlphaEcalHcal","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, sampleRangeHigh);
			functionEndcapBetaEcalHcal = new TF1("functionEndcapBetaEcalHcal","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, sampleRangeHigh);//+[8]*exp(-x^[9]/[10])", 0, 1000);
			functionEndcapAlphaHcal = new TF1("functionEndcapAlphaHcal","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, sampleRangeHigh);
			functionEndcapBetaHcal = new TF1("functionEndcapBetaHcal","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0, sampleRangeHigh);


   if(freezeparameters) {
    if(!Parameters24Above500GeV){//Parameters to fix for 2024 pfhc
       //2024 current parameters

      functionBarrelAlphaEcalHcal->FixParameter(0, -0.0337473);
      functionBarrelAlphaEcalHcal->FixParameter(1, 40.0744);
      functionBarrelAlphaEcalHcal->FixParameter(2, -1.07833);
      functionBarrelAlphaEcalHcal->FixParameter(3, 1.7682e-05);
      functionBarrelAlphaEcalHcal->FixParameter(4, 39.9921);
      functionBarrelAlphaEcalHcal->FixParameter(5, 1.81186e-05);
      functionBarrelAlphaEcalHcal->FixParameter(6, -2.83026);
      functionBarrelAlphaEcalHcal->FixParameter(7, -2.82365);


      functionBarrelBetaEcalHcal->FixParameter(0, 0.114298);
      functionBarrelBetaEcalHcal->FixParameter(1, 0.260535);
      functionBarrelBetaEcalHcal->FixParameter(2, -0.76894);
      functionBarrelBetaEcalHcal->FixParameter(3, 0.690161);
      functionBarrelBetaEcalHcal->FixParameter(4, 0.454889);
      functionBarrelBetaEcalHcal->FixParameter(5, 0.149286);
      functionBarrelBetaEcalHcal->FixParameter(6, -19.212);
      functionBarrelBetaEcalHcal->FixParameter(7, -0.44035);


      functionBarrelAlphaHcal->FixParameter(0, -5.8631);
      functionBarrelAlphaHcal->FixParameter(1, 42.9711);
      functionBarrelAlphaHcal->FixParameter(2, 0.647296);
      functionBarrelAlphaHcal->FixParameter(3, 0.380114);
      functionBarrelAlphaHcal->FixParameter(4, 37.1128);
      functionBarrelAlphaHcal->FixParameter(5, 0.300629);
      functionBarrelAlphaHcal->FixParameter(6, -1.22073);
      functionBarrelAlphaHcal->FixParameter(7, -1.26553);


      functionBarrelBetaHcal->FixParameter(0, -26.3294);
      functionBarrelBetaHcal->FixParameter(1, 26.5531);
      functionBarrelBetaHcal->FixParameter(2, 19.5048);
      functionBarrelBetaHcal->FixParameter(3, 1.6943);
      functionBarrelBetaHcal->FixParameter(4, 0.116759);
      functionBarrelBetaHcal->FixParameter(5, 0.0158723);
      functionBarrelBetaHcal->FixParameter(6, -0.435307);
      functionBarrelBetaHcal->FixParameter(7, -0.7014);


      functionEndcapAlphaEcalHcal->FixParameter(0, -74.1841);
      functionEndcapAlphaEcalHcal->FixParameter(1, 165.064);
      functionEndcapAlphaEcalHcal->FixParameter(2, -35.0902);
      functionEndcapAlphaEcalHcal->FixParameter(3, 42.2696);
      functionEndcapAlphaEcalHcal->FixParameter(4, 86.8144);
      functionEndcapAlphaEcalHcal->FixParameter(5, 2.28983);
      functionEndcapAlphaEcalHcal->FixParameter(6, 0.00469653);
      functionEndcapAlphaEcalHcal->FixParameter(7, -0.523986);


      functionEndcapBetaEcalHcal->FixParameter(0, -244.815);
      functionEndcapBetaEcalHcal->FixParameter(1, 244.887);
      functionEndcapBetaEcalHcal->FixParameter(2, 13.679);
      functionEndcapBetaEcalHcal->FixParameter(3, 13.778);
      functionEndcapBetaEcalHcal->FixParameter(4, 0.143508);
      functionEndcapBetaEcalHcal->FixParameter(5, 0.00139485);
      functionEndcapBetaEcalHcal->FixParameter(6, -0.573462);
      functionEndcapBetaEcalHcal->FixParameter(7, -1.62602);


      functionEndcapAlphaHcal->FixParameter(0, -21.2045);
      functionEndcapAlphaHcal->FixParameter(1, 1.61199);
      functionEndcapAlphaHcal->FixParameter(2, -0.924482);
      functionEndcapAlphaHcal->FixParameter(3, -0.356028);
      functionEndcapAlphaHcal->FixParameter(4, -18.8272);
      functionEndcapAlphaHcal->FixParameter(5, 0.957782);
      functionEndcapAlphaHcal->FixParameter(6, -0.144766);
      functionEndcapAlphaHcal->FixParameter(7, -0.318);


      functionEndcapBetaHcal->FixParameter(0, 0.00907483);
      functionEndcapBetaHcal->FixParameter(1, 62.9639);
      functionEndcapBetaHcal->FixParameter(2, -12.057);
      functionEndcapBetaHcal->FixParameter(3, 0.0878058);
      functionEndcapBetaHcal->FixParameter(4, 62.8559);
      functionEndcapBetaHcal->FixParameter(5, 0.0854849);
      functionEndcapBetaHcal->FixParameter(6, -0.678145);
      functionEndcapBetaHcal->FixParameter(7, -0.675762);
    }
    else{
      functionBarrelAlphaEcalHcal->FixParameter(0,-0.0249299);
      functionBarrelAlphaEcalHcal->FixParameter(1,41.0626);
      functionBarrelAlphaEcalHcal->FixParameter(2,-0.756772);
      functionBarrelAlphaEcalHcal->FixParameter(3,9.82297e-07);
      functionBarrelAlphaEcalHcal->FixParameter(4,41.0056);
      functionBarrelAlphaEcalHcal->FixParameter(5,1.02733e-06);
      functionBarrelAlphaEcalHcal->FixParameter(6,-3.70842);
      functionBarrelAlphaEcalHcal->FixParameter(7,-3.69656);

      functionBarrelBetaEcalHcal->FixParameter(0,-0.176791);
      functionBarrelBetaEcalHcal->FixParameter(1,0.60533);
      functionBarrelBetaEcalHcal->FixParameter(2,-0.891364);
      functionBarrelBetaEcalHcal->FixParameter(3,0.033817);
      functionBarrelBetaEcalHcal->FixParameter(4,0.4556);
      functionBarrelBetaEcalHcal->FixParameter(5,0.134414);
      functionBarrelBetaEcalHcal->FixParameter(6,-61.9498);
      functionBarrelBetaEcalHcal->FixParameter(7,-0.520039);

      functionBarrelAlphaHcal->FixParameter(0,-1.0774);
      functionBarrelAlphaHcal->FixParameter(1,40.5678);
      functionBarrelAlphaHcal->FixParameter(2,1.13601);
      functionBarrelAlphaHcal->FixParameter(3,0.0878601);
      functionBarrelAlphaHcal->FixParameter(4,39.5163);
      functionBarrelAlphaHcal->FixParameter(5,0.0785619);
      functionBarrelAlphaHcal->FixParameter(6,-1.39265);
      functionBarrelAlphaHcal->FixParameter(7,-1.43368);

      functionBarrelBetaHcal->FixParameter(0,-26.2747);
      functionBarrelBetaHcal->FixParameter(1,26.5685);
      functionBarrelBetaHcal->FixParameter(2,14.4212);
      functionBarrelBetaHcal->FixParameter(3,2.15401);
      functionBarrelBetaHcal->FixParameter(4,0.592654);
      functionBarrelBetaHcal->FixParameter(5,0.622331);
      functionBarrelBetaHcal->FixParameter(6,-0.45401);
      functionBarrelBetaHcal->FixParameter(7,-0.0735657);

      functionEndcapAlphaEcalHcal->FixParameter(0,-74.1561);
      functionEndcapAlphaEcalHcal->FixParameter(1,165.092);
      functionEndcapAlphaEcalHcal->FixParameter(2,-35.3997);
      functionEndcapAlphaEcalHcal->FixParameter(3,42.5886);
      functionEndcapAlphaEcalHcal->FixParameter(4,86.787);
      functionEndcapAlphaEcalHcal->FixParameter(5,2.24242);
      functionEndcapAlphaEcalHcal->FixParameter(6,0.00751736);
      functionEndcapAlphaEcalHcal->FixParameter(7,-0.530692);

      functionEndcapBetaEcalHcal->FixParameter(0,-244.796);
      functionEndcapBetaEcalHcal->FixParameter(1,244.906);
      functionEndcapBetaEcalHcal->FixParameter(2,12.0756);
      functionEndcapBetaEcalHcal->FixParameter(3,17.4169);
      functionEndcapBetaEcalHcal->FixParameter(4,0.115151);
      functionEndcapBetaEcalHcal->FixParameter(5,0.00235061);
      functionEndcapBetaEcalHcal->FixParameter(6,-0.535085);
      functionEndcapBetaEcalHcal->FixParameter(7,-1.31019);

      functionEndcapAlphaHcal->FixParameter(0,-20.7829);
      functionEndcapAlphaHcal->FixParameter(1,1.40346);
      functionEndcapAlphaHcal->FixParameter(2,-0.75444);
      functionEndcapAlphaHcal->FixParameter(3,-0.394878);
      functionEndcapAlphaHcal->FixParameter(4,-19.0334);
      functionEndcapAlphaHcal->FixParameter(5,1.72721);
      functionEndcapAlphaHcal->FixParameter(6,-0.162634);
      functionEndcapAlphaHcal->FixParameter(7,-0.294783);

      functionEndcapBetaHcal->FixParameter(0,0.0310929);
      functionEndcapBetaHcal->FixParameter(1,62.8732);
      functionEndcapBetaHcal->FixParameter(2,136.556);
      functionEndcapBetaHcal->FixParameter(3,0.0583275);
      functionEndcapBetaHcal->FixParameter(4,63.047);
      functionEndcapBetaHcal->FixParameter(5,0.0603511);
      functionEndcapBetaHcal->FixParameter(6,-0.612309);
      functionEndcapBetaHcal->FixParameter(7,-0.650172);

    }

    
   }

   else {
      functionBarrelAlphaEcalHcal->SetParameters(-0.024931, 41.0626, -0.75677, 9.82294e-07, 41.0056, 1.02733e-06, -3.70842, -3.69656);
      functionBarrelBetaEcalHcal->SetParameters(-0.176792, 0.605329, -0.891352, 0.033817, 0.455599, 0.134413, -61.9498, -0.520039);
      functionBarrelAlphaHcal->SetParameters(-1.07738, 40.5678, 1.13599, 0.08786, 39.5163, 0.0785619, -1.39265, -1.43368);
      functionBarrelBetaHcal->SetParameters(-26.2747, 26.5685, 14.4212, 2.15401, 0.592654, 0.622343, -0.45401, -0.073567);
      functionEndcapAlphaEcalHcal->SetParameters(-74.1563, 165.092, -35.3301, 42.5849, 86.7871, 2.2473, 0.00750041, -0.530646);
      functionEndcapBetaEcalHcal->SetParameters(-244.796, 244.906, 12.0756, 17.4149, 0.1152, 0.00234601, -0.535109, -1.31063);
      functionEndcapAlphaHcal->SetParameters(-20.7829, 1.40346, -0.754448, -0.394878, -19.0334, 1.72721, -0.162634, -0.294782);
      functionEndcapBetaHcal->SetParameters(0.0310899, 62.8732, 136.555, 0.0583275, 63.047, 0.0603511, -0.612309, -0.650172);
   }


   barrelWithEcalHcalCalib->fitAlphasToFunction(functionBarrelAlphaEcalHcal);
   barrelWithEcalHcalCalib->fitAlphasToFunction();
   barrelWithEcalHcalCalib->fitBetasToFunction(functionBarrelBetaEcalHcal);
   barrelWithEcalHcalCalib->fitBetasToFunction();
   endcapWithEcalHcalCalib->fitAlphasToFunction(functionEndcapAlphaEcalHcal);
   endcapWithEcalHcalCalib->fitAlphasToFunction();
   endcapWithEcalHcalCalib->fitBetasToFunction(functionEndcapBetaEcalHcal);
   endcapWithEcalHcalCalib->fitBetasToFunction();

   barrelWithHcalCalib->fitAlphasToFunction(functionBarrelAlphaHcal);
   barrelWithHcalCalib->fitAlphasToFunction();
   barrelWithHcalCalib->fitBetasToFunction(functionBarrelBetaHcal);
   barrelWithHcalCalib->fitBetasToFunction();
   endcapWithHcalCalib->fitAlphasToFunction(functionEndcapAlphaHcal);
   endcapWithHcalCalib->fitAlphasToFunction();
   endcapWithHcalCalib->fitBetasToFunction(functionEndcapBetaHcal);
   endcapWithHcalCalib->fitBetasToFunction();
   endcapWithHcalCalib->fitGammasToFunction(functionEndcapGammaHcal);
   endcapWithHcalCalib->fitGammasToFunction();

   
   
   //Fill all the TH2's that can be put into drawGausFit in order to produce 
   //response and resolution plots.
  //if(drawRespPlots){
    cout<<"Making response and resolution plots..."<<endl;

    
    int contador = -1;
    unsigned N = ETrueEnergies.size();
    PFEnergyCalibration* pec = nullptr;  // Inicializar fuera del bucle si es necesario

    bool usePFEnergyCalibration = PFEnergyCalibrationFunction;
    if (usePFEnergyCalibration||WriteNTupleFile) {
        pec = new PFEnergyCalibration();  // Crear solo una vez si es necesario
    }
    double etrue, momentum,ecal,hcal,ho,eta, abseta,phi, charge,gE,gP,gEta,gPhi,tP,tEta,tPhi,correctedEta, correctedE, PFHCclosure, PFECclosure;
    
    TFile *outFileN = new TFile("/eos/user/c/cmunozdi/OFFLINE_NTUPLES/2024_Merged_3Attempt/2024_0p2to5000GeV_withCorrections.root","recreate");
    TTree *Ntree = new TTree("s","NTuple for energy btw 0.2 and 5000 GeV and PFHC and PFEC energies");
    if (WriteNTupleFile){

      Ntree->Branch("true",&etrue);
      Ntree->Branch("p", &momentum);
      Ntree->Branch("ecal",&ecal);
      Ntree->Branch("hcal",&hcal);
      Ntree->Branch("ho", &ho);
      Ntree->Branch("eta",&eta);
      Ntree->Branch("phi",&phi);
      Ntree->Branch("charge", &charge);
      Ntree->Branch("genE", &gE);
      Ntree->Branch("genP", &gP);
      Ntree->Branch("genEta", &gEta);
      Ntree->Branch("genPhi", &gPhi);
      Ntree->Branch("trkP", &tP);
      Ntree->Branch("trkEta", &tEta);
      Ntree->Branch("trkPhi", &tPhi);
      Ntree->Branch("PFHC_energy",&correctedEta);
      Ntree->Branch("PFEC_energy",&correctedE);
      // Ntree->Branch("PFHC_closure",&PFHCclosure);
      // Ntree->Branch("PFEC_closure",&PFECclosure);
    }

    for (unsigned entry = 0; entry < N; ++entry) {
        contador++;
        etrue = ETrueEnergies[entry];
        momentum = momentums[entry];
        ecal = ecalEnergies[entry];
        hcal = hcalEnergies[entry];
        ho = ho_energies[entry];
        eta = etas[entry];
        abseta = std::abs(eta);
        phi = phis[entry];
        charge = charges[entry];
        gE = genE[entry];
        gP = genP[entry];
        gEta = genEta[entry];
        gPhi = genPhi[entry];
        tP = trkP[entry];
        tEta = trkEta[entry];
        tPhi = trkPhi[entry];

        // Condiciones de filtro
        //if ((ecal + hcal) < 0.5 || etrue < 1.0 || hcal == 0) continue;

        correctedEta = 0;
        double eecalcorr = ecal;
        double ehcalcorr = hcal;
        double etrue_org=etrue; //etrue for charge hadrons (max{etrue, ecal+hcal}) or -1 for neutral hadrons (max{-1, ecal+hcal})


        if (!usePFEnergyCalibration) {
            if (abseta < 1.3) {
                if(ecal > 0) correctedEta = barrelWithEcalHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, abseta, 0);
                else correctedEta = barrelWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, abseta, 0);
                corrEta_range1->Fill(etrue, (correctedEta - etrue) / etrue);
            } else {
                if (ecal > 0 ) correctedEta = endcapWithEcalHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, abseta, 0);
                else correctedEta = endcapWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, abseta, 0);
                if (abseta < 2.4) corrEta_range2->Fill(etrue, (correctedEta - etrue) / etrue);
                else if (abseta < 2.7) corrEta_range3->Fill(etrue, (correctedEta - etrue) / etrue);
                else if (abseta < 3.0) corrEta_range4->Fill(etrue, (correctedEta - etrue) / etrue);
            }
        } else {
            pec->energyEmHad(etrue_org, eecalcorr, ehcalcorr, abseta, phi);
            correctedEta = eecalcorr + ehcalcorr;
            if (abseta < 1.3) {
                corrEta_range1->Fill(etrue, (correctedEta - etrue) / etrue);
            } else if (abseta < 2.4) {
                corrEta_range2->Fill(etrue, (correctedEta - etrue) / etrue);
            } else if (abseta < 2.7) {
                corrEta_range3->Fill(etrue, (correctedEta - etrue) / etrue);
            } else if (abseta < 3.0) {
                corrEta_range4->Fill(etrue, (correctedEta - etrue) / etrue);
            }
        }

        if (WriteNTupleFile){
          pec->energyEmHad(etrue, eecalcorr, ehcalcorr, abseta, phi);
          correctedE = eecalcorr + ehcalcorr;
          PFECclosure = (correctedE - etrue) / etrue;
          if (abseta<1.5) {
            if (ecal > 0) correctedEta = barrelWithEcalHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, abseta, 0);
            else correctedEta = barrelWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, abseta, 0);
          }else{
            if (ecal > 0) correctedEta = endcapWithEcalHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, abseta, 0);
            else correctedEta = endcapWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, abseta, 0);
          }
          PFHCclosure = (correctedEta - etrue) / etrue;
          Ntree->Fill();
        }
        

        // Impresión del progreso
        int frac = static_cast<int>(static_cast<double>(entry) / N * 100);
        if (frac % 10 == 0) {
            static int lastReportedFrac = -1;
            if (frac != lastReportedFrac) {
                std::cout << frac << "%" << std::endl;
                lastReportedFrac = frac;
            }
        }

        //}
        // if(fabs(eta) < 1.5){//alpha beta fit range for barrel
        //     raw->Fill(etrue, (ecal + hcal-etrue)/etrue);
            
        //     if(ecal > 0){//EH-hadrons


        //         correctedEta = barrelWithEcalHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, eta, 0);

        //         correctedE = barrelWithEcalHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, 0);

        //         if(drawpT){
        //             // etrue = etrue/cosh(eta);
        //             // correctedEta = correctedEta/cosh(eta);
        //             // ecal = ecal/cosh(eta);
        //             // hcal = hcal/cosh(eta);
        //             // correctedE = correctedE/cosh(eta);
        //         }else{
        //           correctedE_ErawEcal_EH = barrelWithEcalHcalCalib-> getCalibratedEnergy(etrue, ecal, hcal, 1);
        //           correctedE_ErawHcal_EH = barrelWithEcalHcalCalib-> getCalibratedEnergy(etrue, ecal, hcal, 2);
        //           correctedE_ErawEcalHcal_EH = barrelWithEcalHcalCalib-> getCalibratedEnergy(etrue, ecal, hcal, 3);
        //           correctedEta_org=correctedEta;
        //           correctedE_org=correctedE;
        //           correctedEta_Alpha = barrelWithEcalHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, eta, 1);
        //           correctedEta_Beta = barrelWithEcalHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, eta, 2);
                  
                  
        //           corrEtaBarrelEcalHcal->Fill(etrue, (correctedEta-etrue)/etrue);
        //           corrEtaBarrelEcalHcal_Alpha->Fill(etrue, (correctedEta_Alpha-etrue)/etrue);
        //           corrEtaBarrelEcalHcal_Beta->Fill(etrue, (correctedEta_Beta-etrue)/etrue);

        //           EtaCorrEtaDependenceEH->Fill(eta, (correctedEta-etrue)/etrue);
        //           EtaCorrEtaDependenceEH_Alpha->Fill(eta, (correctedEta_Alpha-etrue)/etrue);
        //           EtaCorrEtaDependenceEH_Beta->Fill(eta, (correctedEta_Beta-etrue)/etrue);
                  
        //           rawEtaDependenceEH->Fill(eta, (ecal + hcal-etrue)/etrue);
                  
        //           corrEtaDependenceEH->Fill(eta, (correctedE-etrue)/etrue);
        //           corrEtaDependenceEH_ErawEcal->Fill(eta, (correctedE_ErawEcal_EH-etrue)/etrue);
        //           corrEtaDependenceEH_ErawHcal->Fill(eta, (correctedE_ErawHcal_EH-etrue)/etrue);
        //           corrEtaDependenceEH_ErawEcalHcal->Fill(eta, (correctedE_ErawEcalHcal_EH-etrue)/etrue);
        //         }
        //         corrEta->Fill(etrue, (correctedEta-etrue)/etrue);
        //         corrEtaBarrel->Fill(etrue, (correctedEta-etrue)/etrue);

        //         EtaCorrEtaDependence->Fill(eta, (correctedEta-etrue)/etrue);

        //         rawEtaDependence->Fill(eta, (ecal + hcal-etrue)/etrue);

        //         corrEtaDependence->Fill(eta, (correctedE-etrue)/etrue);

        //         //if (etrue > 20) {
        //         //corrEtaDependenceEH->Fill(eta, (correctedEta - etrue)/etrue);
        //         //hcorrEtaDependenceEH->Fill(eta, (correctedE - etrue)/etrue);
        //         //}
        //         //corrEtaDependenceProfEH->Fill(etrue, eta, (correctedEta - etrue)/etrue);


        //         h_trueE_vs_mod_eta_response_normalized->Fill(eta,etrue, (correctedEta-etrue)/etrue);
        //         h_trueE_vs_mod_eta_response->Fill(eta,etrue);
        //         if(drawpT) {
        //             etrue = etrue_org;
        //             correctedEta = correctedEta_org;
        //             ecal = ecal_org;
        //             hcal = hcal_org;
        //             correctedE = correctedE_org;
        //         }

        //     }
        //     else{

        //         correctedEta = barrelWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, eta, 0);




        //         correctedE = barrelWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, 0);

        //         correctedEta_org=correctedEta;
        //         correctedE_org=correctedE;
        //         if(drawpT){
        //             // etrue = etrue/cosh(eta);
        //             // correctedEta = correctedEta/cosh(eta);
        //             // ecal = ecal/cosh(eta);
        //             // hcal = hcal/cosh(eta);
        //             // correctedE = correctedE/cosh(eta);
        //         }else{
        //             correctedEta_Alpha = barrelWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, eta, 1);
        //             correctedEta_Beta = barrelWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, eta, 2);
                    
        //             correctedE_ErawHcal_H = barrelWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, 2);

        //             corrEtaDependenceH->Fill(eta, (correctedE-etrue)/etrue);
        //             corrEtaDependenceH_ErawHcal->Fill(eta, (correctedE_ErawHcal_H-etrue)/etrue);

        //             corrEtaBarrelHcal->Fill(etrue, (correctedEta-etrue)/etrue);
        //             corrEtaBarrelHcal_Alpha->Fill(etrue, (correctedEta_Alpha-etrue)/etrue);
        //             corrEtaBarrelHcal_Beta->Fill(etrue, (correctedEta_Beta-etrue)/etrue);

        //             EtaCorrEtaDependenceH->Fill(eta, (correctedEta-etrue)/etrue);
        //             EtaCorrEtaDependenceH_Alpha->Fill(eta, (correctedEta_Alpha-etrue)/etrue);
        //             EtaCorrEtaDependenceH_Beta->Fill(eta, (correctedEta_Beta-etrue)/etrue);    
                    
        //             rawEtaDependenceH->Fill(eta, (ecal + hcal-etrue)/etrue);
        //         }

        //         if(etrue>7 && etrue<9) {
        //             for(int k=0;k<1000;k++) {
        //                 float step=k/500.-0.99995;
        //                 float b= step;
        //                 float a = (etrue - 3.5 )/hcal - 1 -b*eta*eta;
        //                 bcplot->Fill(b,a);
        //             }
        //         }



        //         corrEta->Fill(etrue, (correctedEta-etrue)/etrue);
        //         corrEtaDependence->Fill(eta, (correctedE-etrue)/etrue);

        //         corrEtaBarrel->Fill(etrue, (correctedEta-etrue)/etrue);


        //         EtaCorrEtaDependence->Fill(eta, (correctedEta-etrue)/etrue);

        //         rawEtaDependence->Fill(eta, (ecal + hcal-etrue)/etrue);
        //         //if((fabs(eta) < 1.5) && (correctedEta != correctedE)) cout<<"yolo "<<fabs(eta)<<", correctedEta:"<<correctedEta<<", correctedE:"<<correctedE<<", (correctedEta != correctedE):"
        //         //<<(correctedEta != correctedE)<<endl;
        //         //corrEtaDependenceH->Fill(eta, (correctedEta - etrue)/etrue);
        //         //hcorrEtaDependenceH->Fill(eta, (correctedE - etrue)/etrue);
        //         //corrEtaDependenceProfH->Fill(etrue, eta, (correctedEta - etrue)/etrue);
        //         if(drawpT) {
        //             etrue = etrue_org;
        //             correctedEta = correctedEta_org;
        //             ecal = ecal_org;
        //             hcal = hcal_org;
        //             correctedE = correctedE_org;
        //         }


        //     }

        //     //if(fabs(eta) < 1.0) //b, c fit range

        //     if(fabs(eta) < 1.5){ //b, c fit range //shubham Mar 27
        //         rawBarrel->Fill(etrue, (ecal + hcal-etrue)/etrue);

        //         if(ecal > 0){
        //             correctedE = barrelWithEcalHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, 0);

        //             correctedE_org=correctedE;
        //             if(drawpT) {
        //                 // etrue = etrue/cosh(eta);
        //                 // correctedEta = correctedEta/cosh(eta);
        //                 // ecal = ecal/cosh(eta);
        //                 // hcal = hcal/cosh(eta);
        //                 // correctedE = correctedE/cosh(eta);
        //             }else{
        //                 correctedE_ErawEcal_EH = barrelWithEcalHcalCalib-> getCalibratedEnergy(etrue, ecal, hcal, 1);
        //                 correctedE_ErawHcal_EH = barrelWithEcalHcalCalib-> getCalibratedEnergy(etrue, ecal, hcal, 2);
        //                 correctedE_ErawEcalHcal_EH = barrelWithEcalHcalCalib-> getCalibratedEnergy(etrue, ecal, hcal, 3);

        //                 corrBarrelEcalHcal->Fill(etrue, (correctedE-etrue)/etrue);
        //                 corrBarrelEcalHcal_ErawEcal->Fill(etrue, (correctedE_ErawEcal_EH-etrue)/etrue);
        //                 corrBarrelEcalHcal_ErawHcal->Fill(etrue, (correctedE_ErawHcal_EH-etrue)/etrue);
        //                 corrBarrelEcalHcal_ErawEcalHcal->Fill(etrue, (correctedE_ErawEcalHcal_EH-etrue)/etrue);
        //             }

        //             rawBarrelEcalHcal->Fill(etrue, (ecal + hcal -etrue)/etrue );
        //             corrBarrel->Fill(etrue, (correctedE-etrue)/etrue);


        //             // hcorrEtaDependence->Fill(eta, (correctedE - etrue)/etrue);

        //             //rawEtaDependence->Fill(eta, (ecal + hcal - etrue)/etrue);
        //             // corrEtaDependence->Fill(eta, (correctedEta - etrue)/etrue);

        //             // if(entry<5000) 
        //             // 	 cout<<entry<<"   "<<eta<<"   "<<etrue<<"   "<<ecal+hcal<<"   "<<correctedE<<"   "<<correctedEta<<endl;


        //             //h_response_vs_phi_barrel_EH->Fill(phi, (correctedE - etrue)/etrue); //shuham

        //             if(drawpT) {
        //                 etrue = etrue_org;
        //                 correctedEta = correctedEta_org;
        //                 ecal = ecal_org;
        //                 hcal = hcal_org;
        //                 correctedE = correctedE_org;
        //             }

        //         }
        //         else{
        //             correctedE = barrelWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, 0);
        //             correctedE_org=correctedE;		

        //             if(drawpT) {
        //                 // etrue = etrue/cosh(eta);
        //                 // correctedEta = correctedEta/cosh(eta);
        //                 // ecal = ecal/cosh(eta);
        //                 // hcal = hcal/cosh(eta);
        //                 // correctedE = correctedE/cosh(eta);
        //             }else{
        //                 correctedE_ErawHcal_H = barrelWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, 2);   
                        
        //                 rawBarrelHcal->Fill(etrue, ( ecal + hcal-etrue)/etrue );// (etrue-3.0)/(ecal+hcal) );//, 11936/(3917*sigmas[entry]*sigmas[entry]) );// ( ecal + hcal - etrue)/etrue );

        //                 corrBarrelHcal->Fill(etrue, (correctedE- etrue -etrue)/etrue);
        //                 corrBarrelHcal_ErawHcal->Fill(etrue, (correctedE_ErawHcal_H-etrue)/etrue);
        //             }


        //             corrBarrel->Fill(etrue, (correctedE- etrue -etrue)/etrue);

        //             //h_response_vs_phi_barrel_H->Fill(phi, (correctedE - etrue)/etrue); //shuham
        //             if(drawpT) {
        //                 etrue = etrue_org;
        //                 correctedEta = correctedEta_org;
        //                 ecal = ecal_org;
        //                 hcal = hcal_org;
        //                 correctedE = correctedE_org;
        //             }


        //         }
        //     }
        // }
          
        //   //if(fabs(eta) < 2.5 && fabs(eta) > 1.55) //WITHIN TRACKER alpha beta fit range for endcap 
        //   //if(fabs(eta) < 3.0 && fabs(eta) > 1.55) //FULL EndCap alpha beta fit range for endcap   //shubham
        // //if(fabs(eta) < 3.0 && fabs(eta) > 2.5) //OUTSIDE TRACKER alpha beta fit range for endcap   //shubham
        // if(fabs(eta) < _etaMax_ && fabs(eta) > _etaMin_){
        //     //if (fabs(eta) > 2.7) cout<<"yolo "<<fabs(eta)<<endl;
        //     raw->Fill(etrue, (ecal + hcal-etrue)/etrue);

        //     ////////////////////////
        //     // RAW Proxy
        //     double etrue_proxy;
        //     if (fabs(eta) > 2.5) etrue_proxy = etrue;//ecal + hcal;
        //     else etrue_proxy = etrue;

        //     if(ecal > 0){
        //         correctedEta = endcapWithEcalHcalCalib->getCalibratedEnergy(etrue_proxy, ecal, hcal, eta, 0);




        //         correctedE = endcapWithEcalHcalCalib->getCalibratedEnergy(etrue_proxy, ecal, hcal, 0);
        //         correctedEta_org=correctedEta;
        //         correctedE_org=correctedE;
        //         if(drawpT) {
        //             // etrue = etrue/cosh(eta);
        //             // etrue_proxy = etrue_proxy/cosh(eta);
        //             // correctedEta = correctedEta/cosh(eta);
        //             // ecal = ecal/cosh(eta);
        //             // hcal = hcal/cosh(eta);
        //             // correctedE = correctedE/cosh(eta);
        //         }else{
        //             correctedEta_Alpha = endcapWithEcalHcalCalib->getCalibratedEnergy(etrue_proxy, ecal, hcal, eta, 1);
        //             correctedEta_Beta = endcapWithEcalHcalCalib->getCalibratedEnergy(etrue_proxy, ecal, hcal, eta, 2);

        //             correctedE_ErawEcal_EH = endcapWithEcalHcalCalib->getCalibratedEnergy(etrue_proxy, ecal, hcal, 1);
        //             correctedE_ErawHcal_EH = endcapWithEcalHcalCalib->getCalibratedEnergy(etrue_proxy, ecal, hcal, 2);
        //             correctedE_ErawEcalHcal_EH = endcapWithEcalHcalCalib->getCalibratedEnergy(etrue_proxy, ecal, hcal, 3);

        //             corrEtaDependenceEH->Fill(eta, (correctedE-etrue)/etrue);
        //             corrEtaDependenceEH_ErawEcal->Fill(eta, (correctedE_ErawEcal_EH-etrue)/etrue);
        //             corrEtaDependenceEH_ErawHcal->Fill(eta, (correctedE_ErawHcal_EH-etrue)/etrue);
        //             corrEtaDependenceEH_ErawEcalHcal->Fill(eta, (correctedE_ErawEcalHcal_EH-etrue)/etrue);
                    
        //             corrEtaEndcapEcalHcal->Fill(etrue, (correctedEta-etrue)/etrue);
        //             corrEtaEndcapEcalHcal_Alpha->Fill(etrue, (correctedEta_Alpha-etrue)/etrue);
        //             corrEtaEndcapEcalHcal_Beta->Fill(etrue, (correctedEta_Beta-etrue)/etrue);

        //             EtaCorrEtaDependenceEH->Fill(eta, (correctedEta-etrue)/etrue);
        //             EtaCorrEtaDependenceEH_Alpha->Fill(eta, (correctedEta_Alpha-etrue)/etrue);
        //             EtaCorrEtaDependenceEH_Beta->Fill(eta, (correctedEta_Beta-etrue)/etrue);

        //             rawEtaDependenceEH->Fill(eta, (ecal + hcal-etrue)/etrue);
        //         }

        //         corrEta->Fill(etrue, (correctedEta-etrue)/etrue);
        //         corrEtaDependence->Fill(eta, (correctedE-etrue)/etrue);

        //         corrEtaEndcap->Fill(etrue, (correctedEta-etrue)/etrue);

        //         EtaCorrEtaDependence->Fill(eta, (correctedEta-etrue)/etrue);


        //         //////changed changed changed 30 Apr 
        //         //corrEtaEndcapEcalHcal->Fill((ecal+hcal), (correctedEta - etrue)/etrue);
        //         rawEtaDependence->Fill(eta, (ecal + hcal-etrue)/etrue);
        //         //if (etrue > 20) {
        //         //corrEtaDependenceEH->Fill(eta, (correctedEta - etrue)/etrue);
        //         //hcorrEtaDependenceEH->Fill(eta, (correctedE - etrue)/etrue); //FIXME
        //         //}
        //         //corrEtaDependenceProfEH->Fill(etrue, eta, (correctedEta - etrue)/etrue);

        //         //h_trueE_vs_mod_eta_response_normalized->Fill(eta,etrue, (correctedEta - etrue)/etrue);
        //         //h_trueE_vs_mod_eta_response->Fill(eta,etrue);
        //         if(drawpT) {
        //             etrue = etrue_org;
        //             correctedEta = correctedEta_org;
        //             ecal = ecal_org;
        //             hcal = hcal_org;
        //             correctedE = correctedE_org;
        //         }

        //     }
        //     else{
        //         correctedEta = endcapWithHcalCalib->getCalibratedEnergy(etrue_proxy, ecal, hcal, eta, 0);

        //         correctedE = endcapWithHcalCalib->getCalibratedEnergy(etrue_proxy, ecal, hcal, 0);

        //         correctedEta_org=correctedEta;
        //         correctedE_org=correctedE;

        //         if(drawpT) {
        //             // etrue = etrue/cosh(eta);
        //             // etrue_proxy = etrue_proxy/cosh(eta);
        //             // correctedEta = correctedEta/cosh(eta);
        //             // ecal = ecal/cosh(eta);
        //             // hcal = hcal/cosh(eta);
        //             // correctedE = correctedE/cosh(eta);
        //         }else{
        //             correctedEta_Alpha = endcapWithHcalCalib->getCalibratedEnergy(etrue_proxy, ecal, hcal, eta, 1);
        //             correctedEta_Beta = endcapWithHcalCalib->getCalibratedEnergy(etrue_proxy, ecal, hcal, eta, 2);

        //             correctedE_ErawHcal_H = endcapWithHcalCalib->getCalibratedEnergy(etrue_proxy, ecal, hcal, 2);

        //             corrEtaDependenceH->Fill(eta, (correctedE-etrue)/etrue);
        //             corrEtaDependenceH_ErawHcal->Fill(eta, (correctedE_ErawHcal_H-etrue)/etrue); 

        //             corrEtaEndcapHcal->Fill(etrue, (correctedEta-etrue)/etrue);
        //             corrEtaEndcapHcal_Alpha->Fill(etrue, (correctedEta_Alpha-etrue)/etrue);
        //             corrEtaEndcapHcal_Beta->Fill(etrue, (correctedEta_Beta-etrue)/etrue);

        //             EtaCorrEtaDependenceH->Fill(eta, (correctedEta-etrue)/etrue);
        //             EtaCorrEtaDependenceH_Alpha->Fill(eta, (correctedEta_Alpha-etrue)/etrue);
        //             EtaCorrEtaDependenceH_Beta->Fill(eta, (correctedEta_Beta-etrue)/etrue);

        //             rawEtaDependenceH->Fill(eta, (ecal + hcal-etrue)/etrue);
        //         }

        //         corrEta->Fill(etrue, (correctedEta-etrue)/etrue);  
        //         corrEtaDependence->Fill(eta, (correctedE-etrue)/etrue);

        //         corrEtaEndcap->Fill(etrue, (correctedEta-etrue)/etrue);

        //         EtaCorrEtaDependence->Fill(eta, (correctedEta-etrue)/etrue);


        //         // corrEtaEndcapEcalHcal->Fill(etrue, (correctedEta - etrue)/etrue);
        //         // corrEtaEndcapEcalHcal_Alpha->Fill(etrue, (correctedEta_Alpha - etrue)/etrue);
        //         // corrEtaEndcapEcalHcal_Beta->Fill(etrue, (correctedEta_Beta - etrue)/etrue);
        //         rawEtaDependence->Fill(eta, (ecal + hcal-etrue)/etrue);
        //         //corrEtaDependenceH->Fill(eta, (correctedEta - etrue)/etrue);
        //         //hcorrEtaDependenceH->Fill(eta, (correctedE - etrue)/etrue);
        //         //corrEtaDependenceProfH->Fill(etrue, eta, (correctedEta - etrue)/etrue);
        //         if(drawpT) {
        //             etrue = etrue_org;
        //             correctedEta = correctedEta_org;
        //             ecal = ecal_org;
        //             hcal = hcal_org;
        //             correctedE = correctedE_org;
        //         }

        //     }
        //     //if(fabs(eta) < 2.2) //b, c fi trange
        //     if(fabs(eta) < 3.0){ //b, c fi trange   //shubham
                
        //         rawEndcap->Fill(etrue, (ecal + hcal-etrue)/etrue);

        //         if(ecal > 0){

        //             correctedEta = endcapWithEcalHcalCalib->getCalibratedEnergy(etrue_proxy, ecal, hcal, eta);

        //             correctedE = endcapWithEcalHcalCalib->getCalibratedEnergy(etrue_proxy, ecal, hcal, 0);

        //             correctedEta_org=correctedEta;
        //             correctedE_org=correctedE;

        //             if(drawpT) {
        //                 // etrue = etrue/cosh(eta);
        //                 // etrue_proxy = etrue_proxy/cosh(eta);
        //                 // correctedEta = correctedEta/cosh(eta);
        //                 // ecal = ecal/cosh(eta);
        //                 // hcal = hcal/cosh(eta);
        //                 // correctedE = correctedE/cosh(eta);
        //             }else{
        //                 correctedE_ErawEcal_EH = endcapWithEcalHcalCalib-> getCalibratedEnergy(etrue_proxy, ecal, hcal, 1);
        //                 correctedE_ErawHcal_EH = endcapWithEcalHcalCalib-> getCalibratedEnergy(etrue_proxy, ecal, hcal, 2);
        //                 correctedE_ErawEcalHcal_EH = endcapWithEcalHcalCalib-> getCalibratedEnergy(etrue_proxy, ecal, hcal, 3);

        //                 corrEndcapEcalHcal->Fill(etrue, (correctedE-etrue)/etrue);
        //                 corrEndcapEcalHcal_ErawEcal->Fill(etrue, (correctedE_ErawEcal_EH-etrue)/etrue);
        //                 corrEndcapEcalHcal_ErawHcal->Fill(etrue, (correctedE_ErawHcal_EH-etrue)/etrue);
        //                 corrEndcapEcalHcal_ErawEcalHcal->Fill(etrue, (correctedE_ErawEcalHcal_EH-etrue)/etrue);
        //             }

        //             rawEndcapEcalHcal->Fill(etrue, (ecal + hcal-etrue)/etrue);
        //             corrEndcap->Fill(etrue, (correctedE-etrue)/etrue);



        //             //rawEtaDependence->Fill(eta, (ecal + hcal - etrue)/etrue);
        //             // corrEtaDependence->Fill(eta, (correctedEta - etrue)/etrue);
        //             // hcorrEtaDependence->Fill(eta, (correctedE - etrue)/etrue);

        //             //cout<<"yolo, eta:"<<eta<<endl;
        //             if(etas[entry] > 0) h_response_vs_phi_EndCap_EH_posZ->Fill(phi,(correctedE-etrue)/etrue);
        //             else if(etas[entry] < 0) h_response_vs_phi_EndCap_EH_negZ->Fill(phi,(correctedE-etrue)/etrue);

        //             if(drawpT) {
        //                 etrue = etrue_org;
        //                 correctedEta = correctedEta_org;
        //                 ecal = ecal_org;
        //                 hcal = hcal_org;
        //                 correctedE = correctedE_org;
        //             }


        //         }
        //         else{
        //             correctedE = endcapWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, 0);

        //             //            correctedEta_org=correctedEta;
        //             correctedE_org=correctedE;

        //             if(drawpT) {
        //                 // etrue = etrue/cosh(eta);
        //                 // etrue_proxy = etrue_proxy/cosh(eta);
        //                 // correctedEta = correctedEta/cosh(eta);
        //                 // ecal = ecal/cosh(eta);
        //                 // hcal = hcal/cosh(eta);
        //                 // correctedE = correctedE/cosh(eta);
        //             }else{
        //                 correctedE_ErawHcal_H = endcapWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, 2);

        //                 rawEndcapHcal->Fill(etrue, (ecal + hcal-etrue)/etrue);
        //                 corrEndcapHcal->Fill(etrue, (correctedE-etrue)/etrue);
        //                 corrEndcapHcal_ErawHcal->Fill(etrue, (correctedE_ErawHcal_H-etrue)/etrue);
        //             }

        //             corrEndcap->Fill(etrue, (correctedE-etrue)/etrue);

        //             if(etas[entry] > 0) h_response_vs_phi_EndCap_H_posZ->Fill(phi,(correctedE-etrue)/etrue);
        //             else if(etas[entry] < 0) h_response_vs_phi_EndCap_H_negZ->Fill(phi,(correctedE-etrue)/etrue);

        //             correctedEta_org=correctedEta;
        //             correctedE_org=correctedE;


        //         }
        //     }
        //     else{   //shubham
        //         if(ecal > 0) corrEndcapEcalHcal->Fill(etrue, (correctedE-etrue)/etrue);
        //     }
        // }
    }

    if (usePFEnergyCalibration) {
      delete pec;
    }
    if(WriteNTupleFile){
      Ntree->Write();
      outFileN->Close();
      delete outFileN;
    }
    //if (Ntree != nullptr) {
      //delete Ntree;
    //}
    //if (outFileN != nullptr) {
      //delete outFileN;
    //}
  //}

   ////////////////////////////////////////////////////////////////////////////
   //Add all the draw functions that you would like here, as well as any 
   //additional output you would like.
   ////////////////////////////////////////////////////////////////////////////

   cout<<" Now Summary "<<endl;
   
   //   exit(0);
  drawGausFit(corrEta_range1, response, resolution);
  cout << "Ya he salido de drawGausFit1" << endl;
  drawGausFit(corrEta_range2, response, resolution);
  cout << "Ya he salido de drawGausFit2" << endl;
  drawGausFit(corrEta_range3, response, resolution);
  cout << "Ya he salido de drawGausFit3" << endl;
  drawGausFit(corrEta_range4, response, resolution);
  cout << "Ya he salido de drawGausFit4" << endl;
  drawGausFit(corrEta_range5, response, resolution);
  cout << "Ya he salido de drawGausFit5" << endl;
  //drawEtaDependence(EtaCorrEtaDependence, responseEtaEtaEH_and_H);
  //  rawBarrelEcalHcal->Draw("colz");
  //  rawBarrelHcal->Draw("colz");

   if(drawRespPlots){
       rawBarrel->Draw("colz");
      //// raw barrel response for EH-hdarons
      drawGausFit(rawBarrelEcalHcal,responseRaw,resolutionRaw);
      /// E-corrected barrel response for EH-hdarons
      drawGausFit(corrBarrelEcalHcal,responseCor,resolutionCor);
      drawGausFit(corrBarrelEcalHcal_ErawEcal,responseCor,resolutionCor);
      drawGausFit(corrBarrelEcalHcal_ErawHcal,responseCor,resolutionCor);
      drawGausFit(corrBarrelEcalHcal_ErawEcalHcal,responseCor,resolutionCor);
      
      //// Eta-corrected barrel response for EH-hdarons
      drawGausFit(corrEtaBarrelEcalHcal, responseEta, resolutionEta);
      drawGausFit(corrEtaBarrelEcalHcal_Alpha, responseEta, resolutionEta);
      drawGausFit(corrEtaBarrelEcalHcal_Beta, responseEta, resolutionEta);

      
      //// raw barrel response for H-hdarons
      drawGausFit(rawBarrelHcal,responseRaw,resolutionRaw);
      //// E-corrected barrel response for H-hdarons
      drawGausFit(corrBarrelHcal,responseCor,resolutionCor);
      drawGausFit(corrBarrelHcal_ErawHcal, responseCor, resolutionCor);
      //// Eta-corrected barrel response for H-hdarons
      drawGausFit(corrEtaBarrelHcal,responseCor,resolutionCor);
      drawGausFit(corrEtaBarrelHcal_Alpha,responseCor,resolutionCor);
      drawGausFit(corrEtaBarrelHcal_Beta,responseCor,resolutionCor);
      
      
      //// raw endcap response for EH-hdarons 
      rawEndcapEcalHcal->Draw("colz");
      drawGausFit(rawEndcapEcalHcal,responseRaw,resolutionRaw);
      //// E-corrected endcap response for EH-hdarons
      drawGausFit(corrEndcapEcalHcal,responseCor,resolutionCor);
      drawGausFit(corrEndcapEcalHcal_ErawEcal,responseCor,resolutionCor);
      drawGausFit(corrEndcapEcalHcal_ErawHcal,responseCor,resolutionCor);
      drawGausFit(corrEndcapEcalHcal_ErawEcalHcal,responseCor,resolutionCor);
      //// Eta-corrected endcap response for EH-hdarons
      drawGausFit(corrEtaEndcapEcalHcal,responseCor,resolutionCor);
      drawGausFit(corrEtaEndcapEcalHcal_Alpha,responseCor,resolutionCor);
      drawGausFit(corrEtaEndcapEcalHcal_Beta,responseCor,resolutionCor);
      corrEtaEndcapEcalHcal->Draw("colz");
        
            
      //// raw endcap response for H-hdarons
      rawEndcapHcal->Draw("colz");
      drawGausFit(rawEndcapHcal,responseRaw,resolutionRaw);
      ///// E-corrected endcap response for H-hdarons
      drawGausFit(corrEndcapHcal,responseCor,resolutionCor);
      drawGausFit(corrEndcapHcal_ErawHcal, responseCor, resolutionCor);
      //// Eta-corrected endcap response for H-hdarons
      corrEtaEndcapHcal->Draw("colz");
      drawGausFit(corrEtaEndcapHcal, responseEta, resolutionEta);
      drawGausFit(corrEtaEndcapHcal_Alpha, responseEta, resolutionEta);   
      drawGausFit(corrEtaEndcapHcal_Beta, responseEta, resolutionEta);


      drawEtaDependence(EtaCorrEtaDependenceEH, responseEtaEtaEH);
      drawEtaDependence(EtaCorrEtaDependenceEH_Alpha, responseEtaEtaEH);
      drawEtaDependence(EtaCorrEtaDependenceEH_Beta, responseEtaEtaEH);  
      drawEtaDependence(EtaCorrEtaDependenceH, responseEtaEtaH);
      drawEtaDependence(EtaCorrEtaDependenceH_Alpha, responseEtaEtaH);
      drawEtaDependence(EtaCorrEtaDependenceH_Beta, responseEtaEtaH);   


      drawEtaDependence(rawEtaDependenceEH, responseEtaEtaEH);
      //drawEtaDependence(hcorrEtaDependenceEH, responseEtaHCorrEtaEH);
      //drawEtaDependence(corrEtaDependenceEH, responseEtaEtaEH);
      
      drawEtaDependence(rawEtaDependenceH, responseEtaEtaH);
      //drawEtaDependence(hcorrEtaDependenceH, responseEtaHCorrEtaH);
      //drawEtaDependence(corrEtaDependenceH, responseEtaEtaH);
      drawEtaDependence(corrEtaDependenceH, responseEtaEtaH);
      drawEtaDependence(corrEtaDependenceH_ErawHcal, responseEtaEtaH);

      drawEtaDependence(corrEtaDependenceEH, responseEtaEtaEH);
      drawEtaDependence(corrEtaDependenceEH_ErawEcal, responseEtaEtaEH);
      drawEtaDependence(corrEtaDependenceEH_ErawHcal, responseEtaEtaEH);
      drawEtaDependence(corrEtaDependenceEH_ErawEcalHcal, responseEtaEtaEH);
    

   drawEtaDependence(EtaCorrEtaDependence, responseEtaEtaEH_and_H);   
         
   
   // something for overall
   drawGausFit(rawBarrel,responseRaw,resolutionRaw);
   drawGausFit(rawEndcap, responseRaw, resolutionRaw);
   drawEtaDependence(rawEtaDependence, responseEtaEtaEH_and_H);

   drawGausFit(corrBarrel, responseCor, resolutionCor);
   drawGausFit(corrEndcap, responseCor, resolutionCor);

   
   
   //drawGausFit(corrEta,response, resolution);
   drawEtaDependence(corrEtaDependence, responseEtaEtaEH_and_H);


   drawGausFit(corrEtaBarrel, response, resolution);
   drawGausFit(corrEtaEndcap, response, resolution);
   drawCompare(responseRaw, response, resolutionRaw, resolution);

    TCanvas *canvas_off = new TCanvas("canvas_off", "Histograma 2D OFF", 1200, 900);
    canvas_off->SetLogz(); // Aplicar escala logarítmica en el eje z
    //canvas_off->SetLogy();
    canvas_off->SetLogx();
    canvas_off->SetRightMargin(0.15); // Ajustar margen derecho para etiquetas

    corrEta->SetMinimum(1);
    corrEta->SetMaximum(1e4);
    corrEta->GetXaxis()->SetTitle("E_{true} (GeV)");
    corrEta->GetYaxis()->SetTitle("(E_{PFHC}-E_{true})/E_{true}");
    corrEta->SetTitle("Response (offline)");
    corrEta->Draw("COLZ");

    //TLatex* latex = new TLatex();
    /*latex->SetTextFont(42);
    latex->SetTextSize(0.04);
    latex->SetTextAlign(12);
    latex->SetNDC();*/
    //latex->DrawLatex(0.05, 0.02, "#it{Private work} (#bf{CMS} #it{simulation})");

    canvas_off->Update();
    canvas_off->Draw();
    canvas_off->SaveAs("Offline_Etrue_EcalPlusHcalMinusEtrueDivEtrue_histogram.png");

    TFile outputFile_off("Offline_Etrue_EcalPlusHcalMinusEtrueDivEtrue_histogram.root", "RECREATE");
    corrEta->Write(); // Save the histogram as an object
    outputFile_off.Close();



   
   }



      
   // barrel H calibration coefficient

   barrelWithHcalCalib->drawCoeffGraph("C", "H_barrel");
   barrelWithHcalCalib->drawCoeffGraph("Alpha","H_barrel");
   barrelWithHcalCalib->drawCoeffGraph("Beta", "H_barrel");
   
   // endcap H calibration coefficient
   endcapWithHcalCalib->drawCoeffGraph("C", "H_endcap");
   endcapWithHcalCalib->drawCoeffGraph("Alpha","H_endcap");
   endcapWithHcalCalib->drawCoeffGraph("Beta", "H_endcap");

   // barrel EH calibration coefficient
   
   barrelWithEcalHcalCalib->drawCoeffGraph("A","EH_barrel");
   barrelWithEcalHcalCalib->drawCoeffGraph("B", "EH_barrel");
   barrelWithEcalHcalCalib->drawCoeffGraph("Alpha","EH_barrel");
   barrelWithEcalHcalCalib->drawCoeffGraph("Beta", "EH_barrel");
   
   // endcap EH calibration coefficient
      
   endcapWithEcalHcalCalib->drawCoeffGraph("A","EH_endcap");
   endcapWithEcalHcalCalib->drawCoeffGraph("B", "EH_endcap");
   endcapWithEcalHcalCalib->drawCoeffGraph("Alpha","EH_endcap");
   endcapWithEcalHcalCalib->drawCoeffGraph("Beta", "EH_endcap");
   







   //cout<<"Check pt 1"<<endl;
   
   h_trueE_vs_mod_eta_response_normalized->SetXTitle("|#eta|");
   h_trueE_vs_mod_eta_response_normalized->SetYTitle("True Energy");

   h_trueE_vs_mod_eta_response_normalized->Divide(h_trueE_vs_mod_eta_response);

   h_trueE_vs_mod_eta_response_normalized->SetMinimum(-0.2);
   h_trueE_vs_mod_eta_response_normalized->SetMaximum(0.2);
   h_trueE_vs_mod_eta_response_normalized->Draw("colz");
   

   h_response_vs_phi_barrel_EH->Draw(); //shuham






   h_response_vs_phi_EndCap_EH_posZ->SetXTitle("#phi");
   h_response_vs_phi_EndCap_EH_negZ->SetXTitle("#phi");
   h_response_vs_phi_barrel_EH->SetXTitle("#phi"); //shuham
   h_response_vs_phi_EndCap_H_posZ->SetXTitle("#phi");
   h_response_vs_phi_EndCap_H_negZ->SetXTitle("#phi");
   h_response_vs_phi_barrel_H->SetXTitle("#phi");


   h_response_vs_phi_EndCap_EH_posZ->SetYTitle("(E_{corr} - E_{true} / E_{true})");
   h_response_vs_phi_EndCap_EH_negZ->SetYTitle("(E_{corr} - E_{true} / E_{true})");
   h_response_vs_phi_barrel_EH->SetYTitle("(E_{corr} - E_{true} / E_{true})"); //shuham
   h_response_vs_phi_EndCap_H_posZ->SetYTitle("(E_{corr} - E_{true} / E_{true})");
   h_response_vs_phi_EndCap_H_negZ->SetYTitle("(E_{corr} - E_{true} / E_{true})");
   h_response_vs_phi_barrel_H->SetYTitle("(E_{corr} - E_{true} / E_{true})");


   TFile* file=new TFile("output.root","recreate");
   h_response_vs_phi_EndCap_EH_posZ->Write();
   h_response_vs_phi_EndCap_EH_negZ->Write();
   h_response_vs_phi_barrel_EH->Write(); //shuham
   h_response_vs_phi_EndCap_H_posZ->Write();
   h_response_vs_phi_EndCap_H_negZ->Write();
   h_response_vs_phi_barrel_H->Write(); //shuham

   //h_occupancy_correct_response->Write();                                                                                                                                                                 
   file->Close();

   

   // don't know what these are, may be all inclusive
   //barrelWithEcalHcalCalib->drawCoeffGraph("Alpha","EH");
   //barrelWithEcalHcalCalib->drawCoeffGraph("Beta", "EH");
   //endcapWithHcalCalib->drawCoeffGraph("Gamma", "H");



    TCanvas *cded = new TCanvas("cefzf","ceced");
    bcplot->Draw("colz");
    corrEtaDependenceProfH->Draw("colz");

   
   functionBarrelEcalHcalB_e = functionBarrelEcalHcalB->GetTitle();
   functionBarrelEcalHcalC_e = functionBarrelEcalHcalC->GetTitle();
   functionBarrelHcalC_e = functionBarrelHcalC->GetTitle();
   functionBarrelAlphaEH_e = functionBarrelAlphaEcalHcal->GetTitle();
   functionBarrelBetaEH_e = functionBarrelBetaEcalHcal->GetTitle();

   functionBarrelAlphaH_e = functionBarrelAlphaHcal->GetTitle();
   functionBarrelBetaH_e = functionBarrelBetaHcal->GetTitle();


   functionEndcapEcalHcalB_e = functionEndcapEcalHcalB->GetTitle();
   functionEndcapEcalHcalC_e = functionEndcapEcalHcalC->GetTitle();
   functionEndcapHcalC_e = functionEndcapHcalC->GetTitle();
   functionEndcapAlphaEH_e = functionEndcapAlphaEcalHcal->GetTitle();
   functionEndcapBetaEH_e = functionEndcapBetaEcalHcal->GetTitle();
 
   functionEndcapAlphaH_e = functionEndcapAlphaHcal->GetTitle();
   functionEndcapBetaH_e = functionEndcapBetaHcal->GetTitle();


   //FUnction printing =============================================
   //takes more place, but easier to read
   //Thresholds first
   cout<<"  threshE = "<<barrelABCEcalHcal[2]->getA()<<";"<<endl;
   cout<<"  threshH = "<<barrelABCHcal[2]->getA()<<";"<<endl; 

   //Now functions : Ecal coef barrel
   cout<<"  faBarrel = std::make_unique<TF1>(\"faBarrel\",\""<<
     functionBarrelEcalHcalB_e<<"\",1., sampleRangeHigh);"<<endl;
   for ( int i = 0; i < 100; ++i ) 
     {
       if ( i > (functionBarrelEcalHcalB->GetNpar()-1)) break;
       barrelEcalHcalB = functionBarrelEcalHcalB->GetParameter(i);
       if(payload) cout<<barrelEcalHcalB<<",";
       else cout<<"  faBarrel->SetParameter("<<i<<","<<barrelEcalHcalB<<");"<<endl;
     }
   cout<<endl;
   cout<<"  faBarrel at x = 10 ;  "<<","<<functionBarrelEcalHcalB->Eval(10.)<<endl<<endl;
   // Hcal coef barrel for ecalHcal
   cout<<"  fbBarrel = std::make_unique<TF1>(\"fbBarrel\",\""<<
     functionBarrelEcalHcalC_e<<"\",1., sampleRangeHigh);"<<endl;
   for ( int i = 0; i < 100; ++i ) 
     {
       if ( i > (functionBarrelEcalHcalC->GetNpar()-1)) break;
       barrelEcalHcalC = functionBarrelEcalHcalC->GetParameter(i);
       if(payload) cout<<barrelEcalHcalC<<",";
       else cout<<"  fbBarrel->SetParameter("<<i<<","<<barrelEcalHcalC<<");"<<endl;
     }
     cout<<endl;
   cout<<"  fbBarrel at x = 10 ;  "<<","<<functionBarrelEcalHcalC->Eval(10.)<<endl<<endl;

   // Hcal coef barrel for Hcal
   cout<<"  fcBarrel = std::make_unique<TF1>(\"fcBarrel\",\""<<
     functionBarrelHcalC_e<<"\",1., sampleRangeHigh);"<<endl;

   for ( int i = 0; i < 100; ++i ) 
     {
       if ( i > (functionBarrelHcalC->GetNpar()-1)) break;

       barrelHcalC = functionBarrelHcalC->GetParameter(i);
       if(payload) cout<<barrelHcalC<<",";
       else cout<<"  fcBarrel->SetParameter("<<i<<","<<barrelHcalC<<");"<<endl;
     }
   cout<<endl;
   cout<<"  fcBarrel at x = 10 ;  "<<","<<functionBarrelHcalC->Eval(10.)<<endl<<endl;

   //alpha function EcalHcal
   cout<<"  faEtaBarrelEH = std::make_unique<TF1>(\"faEtaBarrelEH\",\""<<
     functionBarrelAlphaEH_e<<"\",1., sampleRangeHigh);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionBarrelAlphaEcalHcal->GetNpar()-1)) break;
       barrelAlpha = functionBarrelAlphaEcalHcal->GetParameter(i);
       if(payload) cout<<barrelAlpha<<",";
       else cout<<"  faEtaBarrelEH->SetParameter("<<i<<","<<barrelAlpha<<");"<<endl;
     }
   cout<<endl;
   cout<<"  faEtaBarrelEH at x = 10 ;  "<<","<<functionBarrelAlphaEcalHcal->Eval(10.)<<endl<<endl;

   //beta function EcalHcal
   cout<<"  fbEtaBarrelEH = std::make_unique<TF1>(\"fbEtaBarrelEH\",\""<<
     functionBarrelBetaEH_e<<"\",1., sampleRangeHigh);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionBarrelBetaEcalHcal->GetNpar()-1)) break;
       barrelBeta = functionBarrelBetaEcalHcal->GetParameter(i);
       if(payload) cout<<barrelBeta<<",";
       else cout<<"  fbEtaBarrelEH->SetParameter("<<i<<","<<barrelBeta<<");"<<endl;
     }
   cout<<endl;
   cout<<"  fbEtaBarrelEH at x = 10 ;  "<<","<<functionBarrelBetaEcalHcal->Eval(10.)<<endl<<endl;

  //alpha function Hcal
   cout<<"  faEtaBarrelH = std::make_unique<TF1>(\"faEtaBarrelH\",\""<<
     functionBarrelAlphaH_e<<"\",1., sampleRangeHigh);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionBarrelAlphaHcal->GetNpar()-1)) break;
       barrelAlpha = functionBarrelAlphaHcal->GetParameter(i);
       if(payload) cout<<barrelAlpha<<",";
       else cout<<"  faEtaBarrelH->SetParameter("<<i<<","<<barrelAlpha<<");"<<endl;
     }
   cout<<endl;
   cout<<"  faEtaBarrelH at x = 10 ;  "<<","<<functionBarrelAlphaHcal->Eval(10.)<<endl<<endl;

   //beta function Hcal
   cout<<"  fbEtaBarrelH = std::make_unique<TF1>(\"fbEtaBarrelH\",\""<<
     functionBarrelBetaH_e<<"\",1., sampleRangeHigh);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionBarrelBetaHcal->GetNpar()-1)) break;
       barrelBeta = functionBarrelBetaHcal->GetParameter(i);
       if(payload) cout<<barrelBeta<<",";
       else cout<<"  fbEtaBarrelH->SetParameter("<<i<<","<<barrelBeta<<");"<<endl;
     }
   cout<<endl;
   cout<<"  fbEtaBarrelH at x = 10 ;  "<<","<<functionBarrelBetaHcal->Eval(10.)<<endl<<endl;

   //Now endcaps (just a copy)
   //Now functions : Ecal coef Endcap
   cout<<"  faEndcap = std::make_unique<TF1>(\"faEndcap\",\""<<
     functionEndcapEcalHcalB_e<<"\",1., sampleRangeHigh);"<<endl;
   for ( int i = 0; i < 100; ++i ) 
     {
       if ( i > (functionEndcapEcalHcalB->GetNpar()-1)) break;
       endcapEcalHcalB = functionEndcapEcalHcalB->GetParameter(i);
       if(payload) cout<<endcapEcalHcalB<<",";
       else cout<<"  faEndcap->SetParameter("<<i<<","<<endcapEcalHcalB<<");"<<endl;
     }
   cout<<endl;
   cout<<"  faEndcap at x = 10 ;  "<<","<<functionEndcapEcalHcalB->Eval(10.)<<endl<<endl;

   // Hcal coef endcap for ecalHcal
   cout<<"  fbEndcap = std::make_unique<TF1>(\"fbEndcap\",\""<<
     functionEndcapEcalHcalC_e<<"\",1., sampleRangeHigh);"<<endl;
   for ( int i = 0; i < 100; ++i ) 
     {
       if ( i > (functionEndcapEcalHcalC->GetNpar()-1)) break;
       endcapEcalHcalC = functionEndcapEcalHcalC->GetParameter(i);
       if(payload) cout<<endcapEcalHcalC<<",";
       else cout<<"  fbEndcap->SetParameter("<<i<<","<<endcapEcalHcalC<<");"<<endl;
     }
   cout<<endl;
   cout<<"  fbEndcap at x = 10 ;  "<<","<<functionEndcapEcalHcalC->Eval(10.)<<endl<<endl;
   // Hcal coef endcap for Hcal
   cout<<"  fcEndcap = std::make_unique<TF1>(\"fcEndcap\",\""<<
     functionEndcapHcalC_e<<"\",1., sampleRangeHigh);"<<endl;
   for ( int i = 0; i < 100; ++i ) 
     {
       if ( i > (functionEndcapHcalC->GetNpar()-1)) break;
       endcapHcalC = functionEndcapHcalC->GetParameter(i);
       if(payload) cout<<endcapHcalC<<",";
       else cout<<"  fcEndcap->SetParameter("<<i<<","<<endcapHcalC<<");"<<endl;
     }
   cout<<endl;
   cout<<"  fcEndcap at x = 10 ;  "<<","<<functionEndcapHcalC->Eval(10.)<<endl<<endl;

   //alpha function for EcalHcal
   cout<<"  faEtaEndcapEH = std::make_unique<TF1>(\"faEtaEndcapEH\",\""<<
     functionEndcapAlphaEH_e<<"\",1., sampleRangeHigh);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionEndcapAlphaEcalHcal->GetNpar()-1)) break;
       endcapAlpha = functionEndcapAlphaEcalHcal->GetParameter(i);
       if(payload) cout<<endcapAlpha<<",";
       else cout<<"  faEtaEndcapEH->SetParameter("<<i<<","<<endcapAlpha<<");"<<endl;
     }
   cout<<endl;
   cout<<"  faEtaEndcapEH at x = 10 ;  "<<","<<functionEndcapAlphaEcalHcal->Eval(10.)<<endl<<endl;

   //beta function for EcalHcal
   cout<<"  fbEtaEndcapEH = std::make_unique<TF1>(\"fbEtaEndcapEH\",\""<<
     functionEndcapBetaEH_e<<"\",1., sampleRangeHigh);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionEndcapBetaEcalHcal->GetNpar()-1)) break;
       endcapBeta = functionEndcapBetaEcalHcal->GetParameter(i);
       if(payload) cout<<endcapBeta<<",";
       else cout<<" fbEtaEndcapEH->SetParameter("<<i<<","<<endcapBeta<<");"<<endl;
     }
   cout<<endl;
   cout<<"  fbEtaEndcapEH at x = 10 ;  "<<","<<functionEndcapBetaEcalHcal->Eval(10.)<<endl<<endl;

 //alpha function for Hcal
   cout<<"  faEtaEndcapH = std::make_unique<TF1>(\"faEtaEndcapH\",\""<<
     functionEndcapAlphaH_e<<"\",1., sampleRangeHigh);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionEndcapAlphaHcal->GetNpar()-1)) break;
       endcapAlpha = functionEndcapAlphaHcal->GetParameter(i);
       if(payload) cout<<endcapAlpha<<",";
       else cout<<"  faEtaEndcapH->SetParameter("<<i<<","<<endcapAlpha<<");"<<endl;

     }
   cout<<endl;
   cout<<"  faEtaEndcapH at x = 10 ;  "<<","<<functionEndcapAlphaHcal->Eval(10.)<<endl<<endl;

   //beta function for Hcal
   cout<<"  fbEtaEndcapH = std::make_unique<TF1>(\"fbEtaEndcapH\",\""<<
     functionEndcapBetaH_e<<"\",1., sampleRangeHigh);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionEndcapBetaHcal->GetNpar()-1)) break;
       endcapBeta = functionEndcapBetaHcal->GetParameter(i);
       if(payload) cout<<endcapBeta<<",";
       else cout<<"  fbEtaEndcapH->SetParameter("<<i<<","<<endcapBeta<<");"<<endl;
     }
   cout<<endl;
   cout<<"  fbEtaEndcapH at x = 10 ;  "<<","<<functionEndcapBetaHcal->Eval(10.)<<endl;


   cout<<"Ndf for  faBarrel ;  "<<functionBarrelEcalHcalB->GetNDF()<<endl;
   cout<<"Chi Square for  faBarrel ;  "<<functionBarrelEcalHcalB->GetChisquare()<<endl;
   cout<<"Chi Square/Ndf for  faBarrel ;  "<<","<<functionBarrelEcalHcalB->GetChisquare()/functionBarrelEcalHcalB->GetNDF()<<endl<<endl;

   cout<<"Ndf for  fbBarrel ;  "<<functionBarrelEcalHcalC->GetNDF()<<endl;
   cout<<"Chi Square for  fbBarrel ;  "<<functionBarrelEcalHcalC->GetChisquare()<<endl;
   cout<<"Chi Square/Ndf for  fbBarrel ;  "<<","<<functionBarrelEcalHcalC->GetChisquare()/functionBarrelEcalHcalC->GetNDF()<<endl<<endl;

   cout<<"Ndf for  fcBarrel ;  "<<functionBarrelEcalHcalB->GetNDF()<<endl;
   cout<<"Chi Square for  fcBarrel ;  "<<functionBarrelEcalHcalB->GetChisquare()<<endl;
   cout<<"Chi Square/Ndf for  fcBarrel ;  "<<","<<functionBarrelHcalC->GetChisquare()/functionBarrelHcalC->GetNDF()<<endl<<endl;

   cout<<"Ndf for  faEndcap ;  "<<functionEndcapEcalHcalB->GetNDF()<<endl;
   cout<<"Chi Square for  faEndcap ;  "<<functionEndcapEcalHcalB->GetChisquare()<<endl;
   cout<<"Chi Square/Ndf for  faEndcap ;  "<<","<<functionEndcapEcalHcalB->GetChisquare()/functionEndcapEcalHcalB->GetNDF()<<endl<<endl;

   cout<<"Ndf for  fbEndcap ;  "<<functionEndcapEcalHcalC->GetNDF()<<endl;
   cout<<"Chi Square for  fbEndcap ;  "<<functionEndcapEcalHcalC->GetChisquare()<<endl;
   cout<<"Chi Square/Ndf for  fbEndcap ;  "<<","<<functionEndcapEcalHcalC->GetChisquare()/functionEndcapEcalHcalC->GetNDF()<<endl<<endl;

   cout<<"Ndf for  fcEndcap ;  "<<functionEndcapHcalC->GetNDF()<<endl;
   cout<<"Chi Square for  fcEndcap ;  "<<functionEndcapHcalC->GetChisquare()<<endl;
   cout<<"Chi Square/Ndf for  fcEndcap ;  "<<","<<functionEndcapHcalC->GetChisquare()/functionEndcapHcalC->GetNDF()<<endl<<endl;

   cout<<"Ndf for  faEtaEndcapEH ;  "<<functionEndcapAlphaEcalHcal->GetNDF()<<endl;
   cout<<"Chi Square for  faEtaEndcapEH ;  "<<functionEndcapAlphaEcalHcal->GetChisquare()<<endl;
   cout<<"Chi Square/Ndf for   faEtaEndcapEH ;  "<<","<<functionEndcapAlphaEcalHcal->GetChisquare()/functionEndcapAlphaEcalHcal->GetNDF()<<endl<<endl;

   cout<<"Ndf for  fbEtaEndcapEH ;  "<<functionEndcapBetaEcalHcal->GetNDF()<<endl;
   cout<<"Chi Square for  fbEtaEndcapEH ;  "<<functionEndcapBetaEcalHcal->GetChisquare()<<endl;
   cout<<"Chi Square/Ndf for   fbEtaEndcapEH ;  "<<","<<functionEndcapBetaEcalHcal->GetChisquare()/functionEndcapBetaEcalHcal->GetNDF()<<endl<<endl;

   cout<<"Ndf for  faEtaEndcapH ;  "<<functionEndcapAlphaHcal->GetNDF()<<endl;
   cout<<"Chi Square for  faEtaEndcapH ;  "<<functionEndcapAlphaHcal->GetChisquare()<<endl;
   cout<<"Chi Square/Ndf for   faEtaEndcapH ;  "<<","<<functionEndcapAlphaHcal->GetChisquare()/functionEndcapAlphaHcal->GetNDF()<<endl<<endl;

   cout<<"Ndf for  fbEtaEndcapH ;  "<<functionEndcapBetaHcal->GetNDF()<<endl;
   cout<<"Chi Square for  fbEtaEndcapH ;  "<<functionEndcapBetaHcal->GetChisquare()<<endl;
   cout<<"Chi Square/Ndf for   fbEtaEndcapH ;  "<<","<<functionEndcapBetaHcal->GetChisquare()/functionEndcapBetaHcal->GetNDF()<<endl<<endl;

   cout<<"Ndf for  faEtaBarrelEH ;  "<<functionBarrelAlphaEcalHcal->GetNDF()<<endl;
   cout<<"Chi Square for  faEtaBarrelEH ;  "<<functionBarrelAlphaEcalHcal->GetChisquare()<<endl;
   cout<<"Chi Square/Ndf for   faEtaBarrelEH ;  "<<","<<functionBarrelAlphaEcalHcal->GetChisquare()/functionBarrelAlphaEcalHcal->GetNDF()<<endl<<endl;

   cout<<"Ndf for  fbEtaBarrelEH ;  "<<functionBarrelBetaEcalHcal->GetNDF()<<endl;
   cout<<"Chi Square for  fbEtaBarrelEH ;  "<<functionBarrelBetaEcalHcal->GetChisquare()<<endl;
   cout<<"Chi Square/Ndf for   fbEtaBarrelEH ;  "<<","<<functionBarrelBetaEcalHcal->GetChisquare()/functionBarrelBetaEcalHcal->GetNDF()<<endl<<endl;

   cout<<"Ndf for  faEtaBarrelH ;  "<<functionBarrelAlphaHcal->GetNDF()<<endl;
   cout<<"Chi Square for  faEtaBarrelH ;  "<<functionBarrelAlphaHcal->GetChisquare()<<endl;
   cout<<"Chi Square/Ndf for   faEtaBarrelH ;  "<<","<<functionBarrelAlphaHcal->GetChisquare()/functionBarrelAlphaHcal->GetNDF()<<endl<<endl;

   cout<<"Ndf for  fbEtaBarrelH ;  "<<functionBarrelBetaHcal->GetNDF()<<endl;
   cout<<"Chi Square for  fbEtaBarrelH ;  "<<functionBarrelBetaHcal->GetChisquare()<<endl;
   cout<<"Chi Square/Ndf for   fbEtaBarrelH ;  "<<","<<functionBarrelBetaHcal->GetChisquare()/functionBarrelBetaHcal->GetNDF()<<endl<<endl;

      // Printing all the coefficients just to copy-paste the SetParameters or FixParameters
    cout << "***************ENERGY COEFFICIENTS SET PARAMETERS***************" << endl << endl;
    cout << "functionBarrelEcalHcalA->FixParameter(0, aEH);" << endl;
    cout << "functionBarrelEcalHcalB->SetParameters(";
    for (int i = 0; i < functionBarrelEcalHcalB->GetNpar(); ++i) {
        barrelEcalHcalB = functionBarrelEcalHcalB->GetParameter(i);
        cout << barrelEcalHcalB;
        if (i < functionBarrelEcalHcalB->GetNpar() - 1) cout << ", ";
    }
    cout << ");" << endl;

    // Hcal coef barrel para ecalHcal
    cout << "functionBarrelEcalHcalC->SetParameters(";
    for (int i = 0; i < functionBarrelEcalHcalC->GetNpar(); ++i) {
        barrelEcalHcalC = functionBarrelEcalHcalC->GetParameter(i);
        cout << barrelEcalHcalC;
        if (i < functionBarrelEcalHcalC->GetNpar() - 1) cout << ", ";
    }
    cout << ");" << endl;

    // Hcal coef barrel para Hcal
    cout << "functionBarrelHcalC->SetParameters(";
    for (int i = 0; i < functionBarrelHcalC->GetNpar(); ++i) {
        barrelHcalC = functionBarrelHcalC->GetParameter(i);
        cout << barrelHcalC;
        if (i < functionBarrelHcalC->GetNpar() - 1) cout << ", ";
    }
    cout << ");" << endl;
    cout << "functionEndcapEcalHcalA->FixParameter(0, aEHe);" << endl;

    cout << "functionEndcapEcalHcalB->SetParameters(";
    for (int i = 0; i < functionEndcapEcalHcalB->GetNpar(); ++i) {
        endcapEcalHcalB = functionEndcapEcalHcalB->GetParameter(i);
        cout << endcapEcalHcalB;
        if (i < functionEndcapEcalHcalB->GetNpar() - 1) cout << ", ";
    }
    cout << ");" << endl;

    cout << "functionEndcapEcalHcalC->SetParameters(";
    for (int i = 0; i < functionEndcapEcalHcalC->GetNpar(); ++i) {
        endcapEcalHcalC = functionEndcapEcalHcalC->GetParameter(i);
        cout << endcapEcalHcalC;
        if (i < functionEndcapEcalHcalC->GetNpar() - 1) cout << ", ";
    }
    cout << ");" << endl;

    cout << "functionBarrelHcalA->FixParameter(0, aH);" << endl;
    cout << "functionBarrelHcalB->FixParameter(0, 0.0);" << endl;
    cout << "functionEndcapHcalA->FixParameter(0, aHe);" << endl;
    cout << "functionEndcapHcalB->FixParameter(0, 0.0);" << endl;

    cout << "functionEndcapHcalC->SetParameters(";
    for (int i = 0; i < functionEndcapHcalC->GetNpar(); ++i) {
        endcapHcalC = functionEndcapHcalC->GetParameter(i);
        cout << endcapHcalC;
        if (i < functionEndcapHcalC->GetNpar() - 1) cout << ", ";
    }
    cout << ");" << endl << endl;

    cout << "***************ETA COEFFICIENTS SET PARAMETERS***************" << endl << endl;

    // Alpha function EcalHcal
    cout << "functionBarrelAlphaEcalHcal->SetParameters(";
    for (int i = 0; i < functionBarrelAlphaEcalHcal->GetNpar(); ++i) {
        barrelAlpha = functionBarrelAlphaEcalHcal->GetParameter(i);
        cout << barrelAlpha;
        if (i < functionBarrelAlphaEcalHcal->GetNpar() - 1) cout << ", ";
    }
    cout << ");" << endl;

    // Beta function EcalHcal
    cout << "functionBarrelBetaEcalHcal->SetParameters(";
    for (int i = 0; i < functionBarrelBetaEcalHcal->GetNpar(); ++i) {
        barrelBeta = functionBarrelBetaEcalHcal->GetParameter(i);
        cout << barrelBeta;
        if (i < functionBarrelBetaEcalHcal->GetNpar() - 1) cout << ", ";
    }
    cout << ");" << endl;

    cout << "functionBarrelAlphaHcal->SetParameters(";
    for (int i = 0; i < functionBarrelAlphaHcal->GetNpar(); ++i) {
        barrelAlpha = functionBarrelAlphaHcal->GetParameter(i);
        cout << barrelAlpha;
        if (i < functionBarrelAlphaHcal->GetNpar() - 1) cout << ", ";
    }
    cout << ");" << endl;

    cout << "functionBarrelBetaHcal->SetParameters(";
    for (int i = 0; i < functionBarrelBetaHcal->GetNpar(); ++i) {
        barrelBeta = functionBarrelBetaHcal->GetParameter(i);
        cout << barrelBeta;
        if (i < functionBarrelBetaHcal->GetNpar() - 1) cout << ", ";
    }
    cout << ");" << endl;

    cout << "functionEndcapAlphaEcalHcal->SetParameters(";
    for (int i = 0; i < functionEndcapAlphaEcalHcal->GetNpar(); ++i) {
        endcapAlpha = functionEndcapAlphaEcalHcal->GetParameter(i);
        cout << endcapAlpha;
        if (i < functionEndcapAlphaEcalHcal->GetNpar() - 1) cout << ", ";
    }
    cout << ");" << endl;

    cout << "functionEndcapBetaEcalHcal->SetParameters(";
    for (int i = 0; i < functionEndcapBetaEcalHcal->GetNpar(); ++i) {
        endcapBeta = functionEndcapBetaEcalHcal->GetParameter(i);
        cout << endcapBeta;
        if (i < functionEndcapBetaEcalHcal->GetNpar() - 1) cout << ", ";
    }
    cout << ");" << endl;

    cout << "functionEndcapAlphaHcal->SetParameters(";
    for (int i = 0; i < functionEndcapAlphaHcal->GetNpar(); ++i) {
        endcapAlpha = functionEndcapAlphaHcal->GetParameter(i);
        cout << endcapAlpha;
        if (i < functionEndcapAlphaHcal->GetNpar() - 1) cout << ", ";
    }
    cout << ");" << endl;

    cout << "functionEndcapBetaHcal->SetParameters(";
    for (int i = 0; i < functionEndcapBetaHcal->GetNpar(); ++i) {
        endcapBeta = functionEndcapBetaHcal->GetParameter(i);
        cout << endcapBeta;
        if (i < functionEndcapBetaHcal->GetNpar() - 1) cout << ", ";
    }
    cout << ");" << endl;

    cout << endl << "***************ENERGY COEFFICIENTS FIX PARAMETERS***************" << endl << endl;
    cout << "  functionBarrelEcalHcalA->FixParameter(0, aEH);" << endl << endl;
   //Now functions : Ecal coef barrel
  //  cout<<"  faBarrel = new TF1(\"faBarrel\",\""<<
  //    functionBarrelEcalHcalB_e<<"\",1.,1000.);"<<endl;
   for ( int i = 0; i < 100; ++i ) 
     {
       if ( i > (functionBarrelEcalHcalB->GetNpar()-1)) break;
       barrelEcalHcalB = functionBarrelEcalHcalB->GetParameter(i);
       cout<<"  functionBarrelEcalHcalB->FixParameter("<<i<<","<<barrelEcalHcalB<<");"<<endl;
     }
   cout<<endl;
  //  cout<<"  faBarrel at x = 10 ;  "<<","<<functionBarrelEcalHcalB->Eval(10.)<<endl<<endl;
   // Hcal coef barrel for ecalHcal
  //  cout<<"  fbBarrel = new TF1(\"fbBarrel\",\""<<
  //    functionBarrelEcalHcalC_e<<"\",1.,1000.);"<<endl;
   for ( int i = 0; i < 100; ++i ) 
     {
       if ( i > (functionBarrelEcalHcalC->GetNpar()-1)) break;
       barrelEcalHcalC = functionBarrelEcalHcalC->GetParameter(i);
       cout<<"  functionBarrelEcalHcalC->FixParameter("<<i<<","<<barrelEcalHcalC<<");"<<endl;
     }
     cout << endl;
  //  cout<<"  fbBarrel at x = 10 ;  "<<","<<functionBarrelEcalHcalC->Eval(10.)<<endl<<endl;

   // Hcal coef barrel for Hcal
  //  cout<<"  fcBarrel = new TF1(\"fcBarrel\",\""<<
  //    functionBarrelHcalC_e<<"\",1.,1000.);"<<endl;

   for ( int i = 0; i < 100; ++i ) 
     {
       if ( i > (functionBarrelHcalC->GetNpar()-1)) break;

       barrelHcalC = functionBarrelHcalC->GetParameter(i);
       cout<<"  functionBarrelHcalC->FixParameter("<<i<<","<<barrelHcalC<<");"<<endl;
     }
   cout<<endl;
  //  cout<<"  fcBarrel at x = 10 ;  "<<","<<functionBarrelHcalC->Eval(10.)<<endl<<endl;

  cout << "  functionEndcapEcalHcalA->FixParameter(0, aEHe);" << endl << endl;

   //Now functions : Ecal coef Endcap
  //  cout<<"  faEndcap = new TF1(\"faEndcap\",\""<<
  //    functionEndcapEcalHcalB_e<<"\",1.,1000.);"<<endl;
   for ( int i = 0; i < 100; ++i ) 
     {
       if ( i > (functionEndcapEcalHcalB->GetNpar()-1)) break;
       endcapEcalHcalB = functionEndcapEcalHcalB->GetParameter(i);
       cout<<"  functionEndcapEcalHcalB->FixParameter("<<i<<","<<endcapEcalHcalB<<");"<<endl;
     }
   cout<<endl;
  //  cout<<"  faEndcap at x = 10 ;  "<<","<<functionEndcapEcalHcalB->Eval(10.)<<endl<<endl;

   // Hcal coef endcap for ecalHcal
  //  cout<<"  fbEndcap = new TF1(\"fbEndcap\",\""<<
  //    functionEndcapEcalHcalC_e<<"\",1.,1000.);"<<endl;
   for ( int i = 0; i < 100; ++i ) 
     {
       if ( i > (functionEndcapEcalHcalC->GetNpar()-1)) break;
       endcapEcalHcalC = functionEndcapEcalHcalC->GetParameter(i);
       cout<<"  functionEndcapEcalHcalC->FixParameter("<<i<<","<<endcapEcalHcalC<<");"<<endl;
     }
   cout<<endl;
  //  cout<<"  fbEndcap at x = 10 ;  "<<","<<functionEndcapEcalHcalC->Eval(10.)<<endl<<endl;
  cout << "  functionBarrelHcalA->FixParameter(0, aH);" << endl;
  cout << "  functionBarrelHcalB->FixParameter(0, 0.0);" << endl;
  cout << "  functionEndcapHcalA->FixParameter(0, aHe);" << endl;
  cout << "  functionEndcapHcalB->FixParameter(0, 0.0);" << endl << endl;

   // Hcal coef endcap for Hcal
  //  cout<<"  fcEndcap = new TF1(\"fcEndcap\",\""<<
  //    functionEndcapHcalC_e<<"\",1.,1000.);"<<endl;
   for ( int i = 0; i < 100; ++i ) 
     {
       if ( i > (functionEndcapHcalC->GetNpar()-1)) break;
       endcapHcalC = functionEndcapHcalC->GetParameter(i);
       cout<<"  functionEndcapHcalC->FixParameter("<<i<<","<<endcapHcalC<<");"<<endl;
     }
   cout<<endl;


    cout << "***************ETA COEFFICIENTS FIX PARAMETERS***************" << endl << endl;
   //Now eta (just a copy)
   //alpha function EcalHcal
  //  cout<<"  faEtaBarrelEH = new TF1(\"faEtaBarrelEH\",\""<<
  //    functionBarrelAlphaEH_e<<"\",1.,1000.);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionBarrelAlphaEcalHcal->GetNpar()-1)) break;
       barrelAlpha = functionBarrelAlphaEcalHcal->GetParameter(i);
       cout<<"  functionBarrelAlphaEcalHcal->FixParameter("<<i<<","<<barrelAlpha<<");"<<endl;
     }
   cout<<endl;
   //cout<<"  faEtaBarrelEH at x = 10 ;  "<<","<<functionBarrelAlphaEcalHcal->Eval(10.)<<endl<<endl;

   //beta function EcalHcal
  //  cout<<"  fbEtaBarrelEH = new TF1(\"fbEtaBarrelEH\",\""<<
  //    functionBarrelBetaEH_e<<"\",1.,1000.);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionBarrelBetaEcalHcal->GetNpar()-1)) break;
       barrelBeta = functionBarrelBetaEcalHcal->GetParameter(i);
       cout<<"  functionBarrelBetaEcalHcal->FixParameter("<<i<<","<<barrelBeta<<");"<<endl;
     }
   cout<<endl;
   //cout<<"  fbEtaBarrelEH at x = 10 ;  "<<","<<functionBarrelBetaEcalHcal->Eval(10.)<<endl<<endl;

  //alpha function Hcal
  //  cout<<"  faEtaBarrelH = new TF1(\"faEtaBarrelH\",\""<<
  //    functionBarrelAlphaH_e<<"\",1.,1000.);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionBarrelAlphaHcal->GetNpar()-1)) break;
       barrelAlpha = functionBarrelAlphaHcal->GetParameter(i);
       cout<<"  functionBarrelAlphaHcal->FixParameter("<<i<<","<<barrelAlpha<<");"<<endl;
     }
   cout<<endl;
  //  cout<<"  faEtaBarrelH at x = 10 ;  "<<","<<functionBarrelAlphaHcal->Eval(10.)<<endl<<endl;

   //beta function Hcal
  //  cout<<"  fbEtaBarrelH = new TF1(\"fbEtaBarrelH\",\""<<
  //    functionBarrelBetaH_e<<"\",1.,1000.);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionBarrelBetaHcal->GetNpar()-1)) break;
       barrelBeta = functionBarrelBetaHcal->GetParameter(i);
       cout<<"  functionBarrelBetaHcal->FixParameter("<<i<<","<<barrelBeta<<");"<<endl;
     }
   cout<<endl;
  //  cout<<"  fbEtaBarrelH at x = 10 ;  "<<","<<functionBarrelBetaHcal->Eval(10.)<<endl<<endl;

   //alpha function for EcalHcal
  //  cout<<"  faEtaEndcapEH = new TF1(\"faEtaEndcapEH\",\""<<
  //    functionEndcapAlphaEH_e<<"\",1.,1000.);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionEndcapAlphaEcalHcal->GetNpar()-1)) break;
       endcapAlpha = functionEndcapAlphaEcalHcal->GetParameter(i);
       cout<<"  functionEndcapAlphaEcalHcal->FixParameter("<<i<<","<<endcapAlpha<<");"<<endl;
     }
   cout<<endl;
  //  cout<<"  faEtaEndcapEH at x = 10 ;  "<<","<<functionEndcapAlphaEcalHcal->Eval(10.)<<endl<<endl;

   //beta function for EcalHcal
  //  cout<<"  fbEtaEndcapEH = new TF1(\"fbEtaEndcapEH\",\""<<
  //    functionEndcapBetaEH_e<<"\",1.,1000.);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionEndcapBetaEcalHcal->GetNpar()-1)) break;
       endcapBeta = functionEndcapBetaEcalHcal->GetParameter(i);
       cout<<"  functionEndcapBetaEcalHcal->FixParameter("<<i<<","<<endcapBeta<<");"<<endl;
     }
   cout<<endl;
  //  cout<<"  fbEtaEndcapEH at x = 10 ;  "<<","<<functionEndcapBetaEcalHcal->Eval(10.)<<endl<<endl;

 //alpha function for Hcal
  //  cout<<"  faEtaEndcapH = new TF1(\"faEtaEndcapH\",\""<<
  //    functionEndcapAlphaH_e<<"\",1.,1000.);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionEndcapAlphaHcal->GetNpar()-1)) break;
       endcapAlpha = functionEndcapAlphaHcal->GetParameter(i);
       cout<<"  functionEndcapAlphaHcal->FixParameter("<<i<<","<<endcapAlpha<<");"<<endl;

     }
   cout<<endl;
  //  cout<<"  faEtaEndcapH at x = 10 ;  "<<","<<functionEndcapAlphaHcal->Eval(10.)<<endl<<endl;

   //beta function for Hcal
  //  cout<<"  fbEtaEndcapH = new TF1(\"fbEtaEndcapH\",\""<<
  //    functionEndcapBetaH_e<<"\",1.,1000.);"<<endl;
   for ( int i = 0; i < 10; ++i ) 
     {
       if ( i > (functionEndcapBetaHcal->GetNpar()-1)) break;
       endcapBeta = functionEndcapBetaHcal->GetParameter(i);
       cout<<"  functionEndcapBetaHcal->FixParameter("<<i<<","<<endcapBeta<<");"<<endl;
     }
   cout<<endl;



   return 1;
}



