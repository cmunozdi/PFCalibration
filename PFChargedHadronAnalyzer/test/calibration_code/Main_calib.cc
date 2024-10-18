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


bool freezeparameters = false;
bool useMean = false;
bool useMedian = false;
bool changeRange =false;
bool old_logic = false;
bool drawpT = false;
bool drawResoFit = false;
bool saveCanvas = true;
bool payload=true;
bool useP_reco=false;//Instead of using etrue, the calibration code uses p_reco.
bool drawRespPlots = false;
int hadrons_eta_symbol = 0;//This int variable is 0 when we take all the eta values (positives and negatives); +1 when we only take hadrons with positive eta; -1 when we only take hadrons with negatives eta.
bool PFEnergyCalibrationFunction=false; //This bool variable is true for using PFEnergyCalibration function from CMSSW and flase for using the default PFHC calibration function getCalibratedEnergy
bool WriteNTupleFile = false;
//char* _region_ = (char*)"EC_outside_tracker";
//char* _region_ = (char*)"EC_within_tracker";
//char* _region_ = (char*)"barrel";
char* _region_ = (char*)"Full";

float _etaMin_ = 0.0;
float _etaMax_ = 0.0;
/*************************************NUEVO BY MIKKO (Ading EnergyCalibration, cmssw) *****************/
PFEnergyCalibration::PFEnergyCalibration() {
  //calibChrisClean.C calibration parameters bhumika Nov, 2018
  faBarrel = std::make_unique<TF1>(
      "faBarrel", "[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0., 1000.);
  faBarrel->SetParameter(0, 14.9081);
  faBarrel->SetParameter(1, -92.531);
  faBarrel->SetParameter(2, -586.723);
  faBarrel->SetParameter(3, 0.281367);
  faBarrel->SetParameter(4, 13.0608);
  faBarrel->SetParameter(5, 0.450973);
  faBarrel->SetParameter(6, 0.03687);
  faBarrel->SetParameter(7, -0.583429);
  fbBarrel = std::make_unique<TF1>(
      "fbBarrel", "[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0., 1000.);
  fbBarrel->SetParameter(0, 2.414);
  fbBarrel->SetParameter(1, -2.99257);
  fbBarrel->SetParameter(2, -3.10022);
  fbBarrel->SetParameter(3, 2.4884);
  fbBarrel->SetParameter(4, 1.49647);
  fbBarrel->SetParameter(5, 0.0591164);
  fbBarrel->SetParameter(6, 0.401639);
  fbBarrel->SetParameter(7, -0.848485);
  fcBarrel = std::make_unique<TF1>(
      "fcBarrel", "[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0., 1000.);
  fcBarrel->SetParameter(0, 10.7719);
  fcBarrel->SetParameter(1, 6.36096);
  fcBarrel->SetParameter(2, -23.8131);
  fcBarrel->SetParameter(3, 1.77669);
  fcBarrel->SetParameter(4, 12.6614);
  fcBarrel->SetParameter(5, 0.722518);
  fcBarrel->SetParameter(6, 0.0447024);
  fcBarrel->SetParameter(7, -0.636043);
  faEtaBarrelEH = std::make_unique<TF1>("faEtaBarrelEH", "[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0., 1000.);
  faEtaBarrelEH->SetParameter(0, -0.0337473);
  faEtaBarrelEH->SetParameter(1, 40.0744);
  faEtaBarrelEH->SetParameter(2, -1.07833);
  faEtaBarrelEH->SetParameter(3, 1.7682e-05);
  faEtaBarrelEH->SetParameter(4, 39.9921);
  faEtaBarrelEH->SetParameter(5, 1.81186e-05);
  faEtaBarrelEH->SetParameter(6, -2.83026);
  faEtaBarrelEH->SetParameter(7, -2.82365);
  fbEtaBarrelEH = std::make_unique<TF1>("fbEtaBarrelEH", "[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0., 1000.);
  fbEtaBarrelEH->SetParameter(0, 0.114298);
  fbEtaBarrelEH->SetParameter(1, 0.260535);
  fbEtaBarrelEH->SetParameter(2, -0.76894);
  fbEtaBarrelEH->SetParameter(3, 0.690161);
  fbEtaBarrelEH->SetParameter(4, 0.454889);
  fbEtaBarrelEH->SetParameter(5, 0.149286);
  fbEtaBarrelEH->SetParameter(6, -19.212);
  fbEtaBarrelEH->SetParameter(7, -0.44035);
  faEtaBarrelH = std::make_unique<TF1>("faEtaBarrelH", "[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0., 1000.);
  faEtaBarrelH->SetParameter(0, -5.8631);
  faEtaBarrelH->SetParameter(1, 42.9711);
  faEtaBarrelH->SetParameter(2, 0.647296);
  faEtaBarrelH->SetParameter(3, 0.380114);
  faEtaBarrelH->SetParameter(4, 37.1128);
  faEtaBarrelH->SetParameter(5, 0.300629);
  faEtaBarrelH->SetParameter(6, -1.22073);
  faEtaBarrelH->SetParameter(7, -1.26553);
  fbEtaBarrelH = std::make_unique<TF1>("fbEtaBarrelH", "[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0., 1000.);
  fbEtaBarrelH->SetParameter(0, -26.3294);
  fbEtaBarrelH->SetParameter(1, 26.5531);
  fbEtaBarrelH->SetParameter(2, 19.5048);
  fbEtaBarrelH->SetParameter(3, 1.6943);
  fbEtaBarrelH->SetParameter(4, 0.116759);
  fbEtaBarrelH->SetParameter(5, 0.0158723);
  fbEtaBarrelH->SetParameter(6, -0.435307);
  fbEtaBarrelH->SetParameter(7, -0.7014);

  faEndcap = std::make_unique<TF1>(
      "faEndcap", "[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0., 1000.);
  faEndcap->SetParameter(0, 37.7103);
  faEndcap->SetParameter(1, -318.409);
  faEndcap->SetParameter(2, -1095.14);
  faEndcap->SetParameter(3, 0.298041);
  faEndcap->SetParameter(4, 23.6435);
  faEndcap->SetParameter(5, 0.326262);
  faEndcap->SetParameter(6, -0.0049515);
  faEndcap->SetParameter(7, -0.631655);
  fbEndcap = std::make_unique<TF1>(
      "fbEndcap", "[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0., 1000.);
  fbEndcap->SetParameter(0, -2.24814);
  fbEndcap->SetParameter(1, 3.15142);
  fbEndcap->SetParameter(2, 3.76944);
  fbEndcap->SetParameter(3, 1.06815);
  fbEndcap->SetParameter(4, 0.0726289);
  fbEndcap->SetParameter(5, 24.8456);
  fbEndcap->SetParameter(6, -0.609985);
  fbEndcap->SetParameter(7, 0.0968119);
  fcEndcap = std::make_unique<TF1>(
      "fcEndcap", "[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0., 1000.);
  fcEndcap->SetParameter(0, 1.63161);
  fcEndcap->SetParameter(1, 6.36717);
  fcEndcap->SetParameter(2, -33.0955);
  fcEndcap->SetParameter(3, 0.501949);
  fcEndcap->SetParameter(4, 0.856495);
  fcEndcap->SetParameter(5, 0.0255776);
  fcEndcap->SetParameter(6, 0.0809049);
  fcEndcap->SetParameter(7, -1.41804);
  faEtaEndcapEH = std::make_unique<TF1>("faEtaEndcapEH", "[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0., 1000.);
  faEtaEndcapEH->SetParameter(0, -74.1841);
  faEtaEndcapEH->SetParameter(1, 165.064);
  faEtaEndcapEH->SetParameter(2, -35.0902);
  faEtaEndcapEH->SetParameter(3, 42.2696);
  faEtaEndcapEH->SetParameter(4, 86.8144);
  faEtaEndcapEH->SetParameter(5, 2.28983);
  faEtaEndcapEH->SetParameter(6, 0.00469653);
  faEtaEndcapEH->SetParameter(7, -0.523986);
  fbEtaEndcapEH = std::make_unique<TF1>("fbEtaEndcapEH", "[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0., 1000.);
  fbEtaEndcapEH->SetParameter(0, -244.815);
  fbEtaEndcapEH->SetParameter(1, 244.887);
  fbEtaEndcapEH->SetParameter(2, 13.679);
  fbEtaEndcapEH->SetParameter(3, 13.778);
  fbEtaEndcapEH->SetParameter(4, 0.143508);
  fbEtaEndcapEH->SetParameter(5, 0.00139485);
  fbEtaEndcapEH->SetParameter(6, -0.573462);
  fbEtaEndcapEH->SetParameter(7, -1.62602);
  faEtaEndcapH = std::make_unique<TF1>("faEtaEndcapH", "[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0., 1000.);
  faEtaEndcapH->SetParameter(0, -21.2045);
  faEtaEndcapH->SetParameter(1, 1.61199);
  faEtaEndcapH->SetParameter(2, -0.924482);
  faEtaEndcapH->SetParameter(3, -0.356028);
  faEtaEndcapH->SetParameter(4, -18.8272);
  faEtaEndcapH->SetParameter(5, 0.957782);
  faEtaEndcapH->SetParameter(6, -0.144766);
  faEtaEndcapH->SetParameter(7, -0.318);
  fbEtaEndcapH = std::make_unique<TF1>("fbEtaEndcapH", "[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))", 0., 1000.);
  fbEtaEndcapH->SetParameter(0, 0.00907483);
  fbEtaEndcapH->SetParameter(1, 62.9639);
  fbEtaEndcapH->SetParameter(2, -12.057);
  fbEtaEndcapH->SetParameter(3, 0.0878058);
  fbEtaEndcapH->SetParameter(4, 62.8559);
  fbEtaEndcapH->SetParameter(5, 0.0854849);
  fbEtaEndcapH->SetParameter(6, -0.678145);
  fbEtaEndcapH->SetParameter(7, -0.675762);

  //added by Bhumika on 2 august 2018

  fcEtaBarrelH = std::make_unique<TF1>("fcEtaBarrelH", "[3]*((x-[0])^[1])+[2]", 0., 1000.);
  fcEtaBarrelH->SetParameter(0, 0);
  fcEtaBarrelH->SetParameter(1, 2);
  fcEtaBarrelH->SetParameter(2, 0);
  fcEtaBarrelH->SetParameter(3, 1);

  fcEtaEndcapH = std::make_unique<TF1>("fcEtaEndcapH", "[3]*((x-[0])^[1])+[2]", 0., 1000.);
  fcEtaEndcapH->SetParameter(0, 0);
  fcEtaEndcapH->SetParameter(1, 0);
  fcEtaEndcapH->SetParameter(2, 0.05);
  fcEtaEndcapH->SetParameter(3, 0);

  fdEtaEndcapH = std::make_unique<TF1>("fdEtaEndcapH", "[3]*((x-[0])^[1])+[2]", 0., 1000.);
  fdEtaEndcapH->SetParameter(0, 1.5);
  fdEtaEndcapH->SetParameter(1, 4);
  fdEtaEndcapH->SetParameter(2, -1.1);
  fdEtaEndcapH->SetParameter(3, 1.0);

  fcEtaBarrelEH = std::make_unique<TF1>("fcEtaBarrelEH", "[3]*((x-[0])^[1])+[2]", 0., 1000.);
  fcEtaBarrelEH->SetParameter(0, 0);
  fcEtaBarrelEH->SetParameter(1, 2);
  fcEtaBarrelEH->SetParameter(2, 0);
  fcEtaBarrelEH->SetParameter(3, 1);

  fcEtaEndcapEH = std::make_unique<TF1>("fcEtaEndcapEH", "[3]*((x-[0])^[1])+[2]", 0., 1000.);
  fcEtaEndcapEH->SetParameter(0, 0);
  fcEtaEndcapEH->SetParameter(1, 0);
  fcEtaEndcapEH->SetParameter(2, 0);
  fcEtaEndcapEH->SetParameter(3, 0);

  fdEtaEndcapEH = std::make_unique<TF1>("fdEtaEndcapEH", "[3]*((x-[0])^[1])+[2]", 0., 1000.);
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
  double Ymax=0.20;
  double Ymin=-0.5;
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
	    gaus->SetParameters(5000.,0.,0.2);
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
	     gaus->SetParameters( 5000, 10, 5, 0, 0.20 );
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
                       vector<double>& phis)//offline
{
   Float_t         true_;
   Float_t         p_;
   Float_t         ecal_;
   Float_t         hcal_;
   Float_t         eta_;
   Float_t         phi_;
   TBranch        *b_true; 
   TBranch        *b_p;   
   TBranch        *b_ecal;   
   TBranch        *b_hcal;   
   TBranch        *b_eta;    
   TBranch        *b_phi;    

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
   tree->SetBranchStatus("eta", 1);
   tree->SetBranchStatus("phi", 1);
   tree->SetBranchStatus("pfcID", 1);
   tree->SetBranchStatus("Eecal", 1);
   tree->SetBranchStatus("Ehcal", 1);
   tree->SetBranchStatus("dr", 1);
   
   
   
   if(tree->GetBranchStatus("true"))
      tree->SetBranchAddress("true", &true_, &b_true);
   tree->SetBranchAddress("p", &p_, &b_p);
   tree->SetBranchAddress("ecal", &ecal_, &b_ecal);
   tree->SetBranchAddress("hcal", &hcal_, &b_hcal);
   tree->SetBranchAddress("eta", &eta_, &b_eta);
   tree->SetBranchAddress("phi", &phi_, &b_phi);
   tree->SetBranchAddress("pfcID", &pfcID_, &b_pfcID);
   tree->SetBranchAddress("Eecal", &E_ecal_, &b_E_ecal);
   tree->SetBranchAddress("Ehcal", &E_hcal_, &b_E_hcal);
   tree->SetBranchAddress("dr", &dr_, &b_dr);

   double sigmaEcalHcal=1;
   long veto = 0 ;
   //int count = 0;
   bool flag[10] = {0,0,0,0,0,0,0,0,0,0};
	  for(int entry = 0; entry < tree->GetEntries(); entry++) { 
       tree->GetEntry(entry);
       
       // if(ecal_<0.4) continue; //FIXME MM
       if (fabs(eta_) < 2.4 && p_ == 0) continue;
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
 
   add_root_files_to_a_chain(chain, "/eos/home-c/cmunozdi/OFFLINE_NTUPLES/OfflineNTuples_2024GT0_merged");//_merged/");//_proof/");//NewNTuplizerVersion/");/2024_Merged/

//chain->Add("/eos/home-c/cmunozdi/step3_ana/*.root");
   sTree = (TTree*)chain;
   cout<<"Reading input tree..."<<endl;
   getValuesFromTree(sTree, ETrueEnergies, ecalEnergies, 
                     hcalEnergies, etas, phis);




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
   
     //functionEndcapHcalC->FixParameter(0,1.75770);
     //Set the parameters of the function you just defined.
     functionBarrelEcalHcalA->FixParameter(0, aEH);

     //faBarrel

     //for trial1/*
     //Offline Parameters
    //  functionBarrelEcalHcalB->FixParameter(0,14.3654);
    //  functionBarrelEcalHcalB->FixParameter(1,-522.051);
    //  functionBarrelEcalHcalB->FixParameter(2,62.5789);
    //  functionBarrelEcalHcalB->FixParameter(3,0.313092);
    //  functionBarrelEcalHcalB->FixParameter(4,13.3344);
    //  functionBarrelEcalHcalB->FixParameter(5,0.203111);
    //  functionBarrelEcalHcalB->FixParameter(6,0.102733);
    //  functionBarrelEcalHcalB->FixParameter(7,-0.611068);
    
    /*//Online Parameters
    functionBarrelEcalHcalB->FixParameter(0, -2.80407e+01);
    functionBarrelEcalHcalB->FixParameter(1,  2.46818e+01);
    functionBarrelEcalHcalB->FixParameter(2,  3.25377e+01);
    functionBarrelEcalHcalB->FixParameter(3,  8.20102e-01);
    functionBarrelEcalHcalB->FixParameter(4, -4.61369e+00);
    functionBarrelEcalHcalB->FixParameter(5,  1.94869e-01);
    functionBarrelEcalHcalB->FixParameter(6, -7.71434e-01);
    functionBarrelEcalHcalB->FixParameter(7, -4.45407e-01);
     */
     /*
     functionBarrelEcalHcalB->FixParameter(0,19.4031);
     functionBarrelEcalHcalB->FixParameter(1,-288.465);
     functionBarrelEcalHcalB->FixParameter(2,-408.44);
     functionBarrelEcalHcalB->FixParameter(3,0.298629);
     functionBarrelEcalHcalB->FixParameter(4,8.82632);
     functionBarrelEcalHcalB->FixParameter(5,0.273511);
     functionBarrelEcalHcalB->FixParameter(6,0.00277729);
     functionBarrelEcalHcalB->FixParameter(7,-0.590834);
     */
     /*
     functionBarrelEcalHcalB->FixParameter(0,-3.32905);
     functionBarrelEcalHcalB->FixParameter(1,4.47069);
     functionBarrelEcalHcalB->FixParameter(2,8.48508);
     functionBarrelEcalHcalB->FixParameter(3,0.156984);
     functionBarrelEcalHcalB->FixParameter(4,-1.56608);
     functionBarrelEcalHcalB->FixParameter(5,0.6864);
     */
     
     //fbBarrel

     //trial 
     //Offline Parameters
    //  functionBarrelEcalHcalC->FixParameter(0,2.34455);
    //  functionBarrelEcalHcalC->FixParameter(1,-10.7725);
    //  functionBarrelEcalHcalC->FixParameter(2,1.55872);
    //  functionBarrelEcalHcalC->FixParameter(3,0.974501);
    //  functionBarrelEcalHcalC->FixParameter(4,1.4455);
    //  functionBarrelEcalHcalC->FixParameter(5,0.0679553);
    //  functionBarrelEcalHcalC->FixParameter(6,0.277572);
    //  functionBarrelEcalHcalC->FixParameter(7,-0.836124);

    /*//Online Parameters
    functionBarrelEcalHcalC->FixParameter(0, -2.53495e+01);
    functionBarrelEcalHcalC->FixParameter(1,  2.83393e+01);
    functionBarrelEcalHcalC->FixParameter(2,  3.04789e+00);
    functionBarrelEcalHcalC->FixParameter(3,  4.03819e+00);
    functionBarrelEcalHcalC->FixParameter(4,  1.98035e+00);
    functionBarrelEcalHcalC->FixParameter(5,  1.57109e-01);
    functionBarrelEcalHcalC->FixParameter(6, -3.55649e-01);
    functionBarrelEcalHcalC->FixParameter(7, -4.61410e-01);
     */
     /*
     functionBarrelEcalHcalC->FixParameter(0,1.93084);
     functionBarrelEcalHcalC->FixParameter(1,-3.63003);
     functionBarrelEcalHcalC->FixParameter(2,-18.4879);
     functionBarrelEcalHcalC->FixParameter(3,0.573916);
     functionBarrelEcalHcalC->FixParameter(4,0.965864);
     functionBarrelEcalHcalC->FixParameter(5,0.111865);
     functionBarrelEcalHcalC->FixParameter(6,0.10636);
     functionBarrelEcalHcalC->FixParameter(7,-0.626559);
     */

     //fcBarrel

     //for UL2016
     //Offline Parameters
    //  functionBarrelHcalC->FixParameter(0,7.55465);
    //  functionBarrelHcalC->FixParameter(1,8.53437);
    //  functionBarrelHcalC->FixParameter(2,-14.0856);
    //  functionBarrelHcalC->FixParameter(3,5.32547);
    //  functionBarrelHcalC->FixParameter(4,13.1732);
    //  functionBarrelHcalC->FixParameter(5,0.978521);
    //  functionBarrelHcalC->FixParameter(6,0.0427651);
    //  functionBarrelHcalC->FixParameter(7,-0.582561);

    /*//Online Parameters
    functionBarrelHcalC->FixParameter(0, -1.12778e+00);
    functionBarrelHcalC->FixParameter(1,  1.88221e+01);
    functionBarrelHcalC->FixParameter(2,  2.25486e+01);
    functionBarrelHcalC->FixParameter(3,  2.19053e-01);
    functionBarrelHcalC->FixParameter(4,  1.69533e+01);
    functionBarrelHcalC->FixParameter(5,  2.06521e-01);
    functionBarrelHcalC->FixParameter(6, -6.63135e-01);
    functionBarrelHcalC->FixParameter(7, -7.89464e-01);
     */
     /*
       functionBarrelHcalC->FixParameter(0,7.53817);
       functionBarrelHcalC->FixParameter(1,8.48519);
       functionBarrelHcalC->FixParameter(2,-14.9824);
       functionBarrelHcalC->FixParameter(3,5.66342);
       functionBarrelHcalC->FixParameter(4,13.4196);
       functionBarrelHcalC->FixParameter(5,0.836512);
       functionBarrelHcalC->FixParameter(6,0.0244799);
       functionBarrelHcalC->FixParameter(7,-0.575392);
     */

       functionEndcapEcalHcalA->FixParameter(0, aEHe);
  

     //faEndcap
     //trial2
     //Offline Parameters
      //  functionEndcapEcalHcalB->FixParameter(0,19.4940);
      //  functionEndcapEcalHcalB->FixParameter(1,-276.99);
      //  functionEndcapEcalHcalB->FixParameter(2,-538.549);
      //  functionEndcapEcalHcalB->FixParameter(3,0.296178);
      //  functionEndcapEcalHcalB->FixParameter(4,7.15250);
      //  functionEndcapEcalHcalB->FixParameter(5,0.114480);
      //  functionEndcapEcalHcalB->FixParameter(6,-0.00592512);
      //  functionEndcapEcalHcalB->FixParameter(7,-0.763133);

       /*//Online Parameters
      functionEndcapEcalHcalB->FixParameter(0, -7.72477e+00);
      functionEndcapEcalHcalB->FixParameter(1,  8.94414e+00);
      functionEndcapEcalHcalB->FixParameter(2,  9.28640e+01);
      functionEndcapEcalHcalB->FixParameter(3,  3.78652e-01);
      functionEndcapEcalHcalB->FixParameter(4, -2.49208e+00);
      functionEndcapEcalHcalB->FixParameter(5,  8.71209e-03);
      functionEndcapEcalHcalB->FixParameter(6, -2.12572e-01);
      functionEndcapEcalHcalB->FixParameter(7, -2.03657e+00);
       */
       /*
       functionEndcapEcalHcalB->FixParameter(0,1.10850);
       functionEndcapEcalHcalB->FixParameter(1,1666.24);
       functionEndcapEcalHcalB->FixParameter(2,-3984.34);
       functionEndcapEcalHcalB->FixParameter(3,0.206562);
       functionEndcapEcalHcalB->FixParameter(4,0.105031);
       */       
       /*
       functionEndcapEcalHcalB->FixParameter(0,1.12385);
       functionEndcapEcalHcalB->FixParameter(1,340.155);
       functionEndcapEcalHcalB->FixParameter(2,-812.381);
       functionEndcapEcalHcalB->FixParameter(3,0.296042);
       functionEndcapEcalHcalB->FixParameter(4,0.135141);
       */

       /*
       functionEndcapEcalHcalB->FixParameter(0,1.17239);
     functionEndcapEcalHcalB->FixParameter(1,158.189);
     functionEndcapEcalHcalB->FixParameter(2,-378.552);
     functionEndcapEcalHcalB->FixParameter(3,0.36484);
     functionEndcapEcalHcalB->FixParameter(4,0.153998);
       */

     // fbEndcap
      //Offline Parameters
      //  functionEndcapEcalHcalC->FixParameter(0,-2.25169);
      //  functionEndcapEcalHcalC->FixParameter(1,3.15393);
      //  functionEndcapEcalHcalC->FixParameter(2,2.2564);
      //  functionEndcapEcalHcalC->FixParameter(3,0.970055);
      //  functionEndcapEcalHcalC->FixParameter(4,0.0561773);
      //  functionEndcapEcalHcalC->FixParameter(5,21.5052);
      //  functionEndcapEcalHcalC->FixParameter(6,-0.845143);
      //  functionEndcapEcalHcalC->FixParameter(7,0.0898543);

      /*//Online Parameters
      functionEndcapEcalHcalC->FixParameter(0, -2.66848e+01);
      functionEndcapEcalHcalC->FixParameter(1,  2.62733e+01);
      functionEndcapEcalHcalC->FixParameter(2,  4.89762e+00);
      functionEndcapEcalHcalC->FixParameter(3,  6.12116e+00);
      functionEndcapEcalHcalC->FixParameter(4, -1.29643e+00);
      functionEndcapEcalHcalC->FixParameter(5,  5.81421e-03);
      functionEndcapEcalHcalC->FixParameter(6, -5.87004e-01);
      functionEndcapEcalHcalC->FixParameter(7, -2.14092e+00);     
       */
       /*
      functionEndcapEcalHcalC->FixParameter(0,-2.25272);
     functionEndcapEcalHcalC->FixParameter(1,3.15246);
     functionEndcapEcalHcalC->FixParameter(2,2.26428);
     functionEndcapEcalHcalC->FixParameter(3,0.905836);
     functionEndcapEcalHcalC->FixParameter(4,0.107686);
     functionEndcapEcalHcalC->FixParameter(5,18.6568);
     functionEndcapEcalHcalC->FixParameter(6,-0.878386);
     functionEndcapEcalHcalC->FixParameter(7,0.115737);
       */

     functionBarrelHcalA->FixParameter(0, aH);
     functionBarrelHcalB->FixParameter(0, 0.0);

     functionEndcapHcalA->FixParameter(0, aHe);
     functionEndcapHcalB->FixParameter(0, 0.0);

     //fcEndcap 

     //for UL2016
     //Offline Parameters
    //  functionEndcapHcalC->FixParameter(0,1.60863);
    //  functionEndcapHcalC->FixParameter(1,2.56704);
    //  functionEndcapHcalC->FixParameter(2,-15.2573);
    //  functionEndcapHcalC->FixParameter(3,0.772001);
    //  functionEndcapHcalC->FixParameter(4,0.802565);
    //  functionEndcapHcalC->FixParameter(5,0.0215722);
    //  functionEndcapHcalC->FixParameter(6,0.116599);
    //  functionEndcapHcalC->FixParameter(7,-1.48251);

    /*//Online Parameters
    functionEndcapHcalC->FixParameter(0, -4.49522e-02);
    functionEndcapHcalC->FixParameter(1,  3.89030e+01);
    functionEndcapHcalC->FixParameter(2,  2.65705e-01);
    functionEndcapHcalC->FixParameter(3,  1.97489e-02);
    functionEndcapHcalC->FixParameter(4,  3.79233e+01);
    functionEndcapHcalC->FixParameter(5,  1.98860e-02);
    functionEndcapHcalC->FixParameter(6, -2.57695e+00);
    functionEndcapHcalC->FixParameter(7, -2.56589e+00);
   */
   /*
     functionEndcapHcalC->FixParameter(0,1.63389);
     functionEndcapHcalC->FixParameter(1,2.10783);
     functionEndcapHcalC->FixParameter(2,-14.7279);
     functionEndcapHcalC->FixParameter(3,0.792741);
     functionEndcapHcalC->FixParameter(4,0.758436);
     functionEndcapHcalC->FixParameter(5,0.0226573);
     functionEndcapHcalC->FixParameter(6,0.150561);
   functionEndcapHcalC->FixParameter(7,-1.48251);
   */
     //  functionEndcapHcalC->FixParameter(0,1.75770);
     //     functionEndcapHcalC->FixParameter(1,0.866230);
     //     functionEndcapHcalC->FixParameter(2,-5.60717);
     //     functionEndcapHcalC->FixParameter(3,3.08493);
     //     functionEndcapHcalC->FixParameter(4,0.966077);
     //     functionEndcapHcalC->FixParameter(5,0.0315277);
     //     functionEndcapHcalC->FixParameter(6,0.253619);
     //     functionEndcapHcalC->FixParameter(7,-1.33786);


     //functionEndcapHcalC->FixParameter(0,1.32375);
     //functionEndcapHcalC->FixParameter(1,0.191687);
     //functionEndcapHcalC->FixParameter(2,-2.81052);
     //functionEndcapHcalC->FixParameter(3,27.1729);
     //functionEndcapHcalC->FixParameter(4,0.418262);
     //functionEndcapHcalC->FixParameter(5,0.107711);
     //functionEndcapHcalC->FixParameter(6,1.71039);
     //functionEndcapHcalC->FixParameter(7,-1.01447);
     

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
   else {

     functionBarrelEcalHcalA->FixParameter(0, aEH);

     //faBarrel

     //for UL2016
    //  //OFFLINE
    //  functionBarrelEcalHcalB->SetParameter(0,-1.81635);
    //  functionBarrelEcalHcalB->SetParameter(1,2.88284);
    //  functionBarrelEcalHcalB->SetParameter(2,7.12287);
    //  functionBarrelEcalHcalB->SetParameter(3,0.118197);
    //  functionBarrelEcalHcalB->SetParameter(4,-1.49896);
    //  functionBarrelEcalHcalB->SetParameter(5,0.5988);
    //OFFLINE 23
    //  functionBarrelEcalHcalB->SetParameter(0,14.3654);
    //  functionBarrelEcalHcalB->SetParameter(1,-522.051);
    //  functionBarrelEcalHcalB->SetParameter(2,62.5789);
    //  functionBarrelEcalHcalB->SetParameter(3,0.313092);
    //  functionBarrelEcalHcalB->SetParameter(4,13.3344);
    //  functionBarrelEcalHcalB->SetParameter(5,0.203111);
    //  functionBarrelEcalHcalB->SetParameter(6,0.102733);
    //  functionBarrelEcalHcalB->SetParameter(7,-0.611068);
     functionBarrelEcalHcalB->SetParameters(14.917,-125.322,-1065.61,0.276941,13.4689,0.309006,0.0585969,-0.651582);//, 2., 6., 9.);
     

    /*//ONLINE
    functionBarrelEcalHcalB->FixParameter(0, -13.9219);
    functionBarrelEcalHcalB->FixParameter(1, 14.9124);
    functionBarrelEcalHcalB->FixParameter(2, 5.38578);
    functionBarrelEcalHcalB->FixParameter(3, 0.861981);
    functionBarrelEcalHcalB->FixParameter(4, -0.00759275);
    functionBarrelEcalHcalB->FixParameter(5, 3.73563e-23);
    functionBarrelEcalHcalB->FixParameter(6, -1.17946);
    functionBarrelEcalHcalB->FixParameter(7, -13.3644);
    */
     //fbBarrel

     //for UL2016 
    //  //OFFLINE
    //  functionBarrelEcalHcalC->SetParameter(0,1.78346);
    //  functionBarrelEcalHcalC->SetParameter(1,-0.907339);
    //  functionBarrelEcalHcalC->SetParameter(2,-21.1843);
    //  functionBarrelEcalHcalC->SetParameter(3,0.64064);
    //  functionBarrelEcalHcalC->SetParameter(4,0.820368);
    //  functionBarrelEcalHcalC->SetParameter(5,0.0867443);
    //  functionBarrelEcalHcalC->SetParameter(6,0.166651);
    //  functionBarrelEcalHcalC->SetParameter(7,-0.854152);
    //OFFLINE 23
    //  functionBarrelEcalHcalC->SetParameter(0,2.34455);
    //  functionBarrelEcalHcalC->SetParameter(1,-10.7725);
    //  functionBarrelEcalHcalC->SetParameter(2,1.55872);
    //  functionBarrelEcalHcalC->SetParameter(3,0.974501);
    //  functionBarrelEcalHcalC->SetParameter(4,1.4455);
    //  functionBarrelEcalHcalC->SetParameter(5,0.0679553);
    //  functionBarrelEcalHcalC->SetParameter(6,0.277572);
    //  functionBarrelEcalHcalC->SetParameter(7,-0.836124);
     functionBarrelEcalHcalC->SetParameters(2.48124,-5.76933,-3.53787,1.87364,1.55158,0.0386927,0.37346,-0.943844);    
    /*//ONLINE
  		functionBarrelEcalHcalC->SetParameters(1.70114,0.404676,-3.88962,1.2109e+06,
						0.970741,0.0527482,2.60552,-0.8956);
    */

     //fcBarrel

     //for UL2016
    //  //OFFLINE
    //  functionBarrelHcalC->SetParameter(0,7.40006);
    //  functionBarrelHcalC->SetParameter(1,8.34692);
    //  functionBarrelHcalC->SetParameter(2,-16.4046);
    //  functionBarrelHcalC->SetParameter(3,6.11434);
    //  functionBarrelHcalC->SetParameter(4,13.3222);
    //  functionBarrelHcalC->SetParameter(5,0.684788);
    //  functionBarrelHcalC->SetParameter(6,0.0159107);
    //  functionBarrelHcalC->SetParameter(7,-0.584424);
    //OFFLINE 23
    //  functionBarrelHcalC->SetParameter(0,7.55465);
    //  functionBarrelHcalC->SetParameter(1,8.53437);
    //  functionBarrelHcalC->SetParameter(2,-14.0856);
    //  functionBarrelHcalC->SetParameter(3,5.32547);
    //  functionBarrelHcalC->SetParameter(4,13.1732);
    //  functionBarrelHcalC->SetParameter(5,0.978521);
    //  functionBarrelHcalC->SetParameter(6,0.0427651);
    //  functionBarrelHcalC->SetParameter(7,-0.582561);    
    functionBarrelHcalC->SetParameters(8.65119,5.11266,-22.203,3.87654,10.9399,0.492967,0.0698724,-0.613534);//, -2., 6., 9.);
     /*//ONLINE
    functionBarrelHcalC->SetParameters(1.58827e+00,4.06865e-01,-3.69939e+00,1.28926e+03,
        7.13400e-01,2.21925e-02,1.47842e+00,-1.22041e+00);
    */
     functionEndcapEcalHcalA->FixParameter(0, aEHe);
  

     //faEndcap

     //for UL2016
     //OFFLINE
    //  functionEndcapEcalHcalB->SetParameter(0,1.13723);
    //  functionEndcapEcalHcalB->SetParameter(1,10.1688);
    //  functionEndcapEcalHcalB->SetParameter(2,-23.9953);
    //  functionEndcapEcalHcalB->SetParameter(3,1.23784);
    //  functionEndcapEcalHcalB->SetParameter(4,0.278710);

    //OFFLINE 23
    // functionEndcapEcalHcalB->SetParameter(0,19.4940);
    // functionEndcapEcalHcalB->SetParameter(1,-276.99);
    // functionEndcapEcalHcalB->SetParameter(2,-538.549);
    // functionEndcapEcalHcalB->SetParameter(3,0.296178);
    // functionEndcapEcalHcalB->SetParameter(4,7.15250);
    // functionEndcapEcalHcalB->SetParameter(5,0.114480);
    // functionEndcapEcalHcalB->SetParameter(6,-0.00592512);
    // functionEndcapEcalHcalB->SetParameter(7,-0.763133);
    functionEndcapEcalHcalB->SetParameters(29.8453,-324.339,-1039.52,0.29767,15.2922,0.190211,-0.00569525,-0.672166);//, 10., 6., 10.);
    /*//ONLINE
				functionEndcapEcalHcalB->SetParameters(0.930193,11.9536,-30.0337,0.76133,
						0.0776373,7.3809e-10,0.158734,-6.92163);
    */
				//** begin modified by seema  
				//functionEndcapEcalHcalB->FixParameter(0, 0.930193);
				/*functionEndcapEcalHcalB->FixParameter(0, 0.962468);
				functionEndcapEcalHcalB->FixParameter(1, 11.9536);
				//functionEndcapEcalHcalB->FixParameter(2, -30.0337);
				//functionEndcapEcalHcalB->FixParameter(2, -28.6722);
				functionEndcapEcalHcalB->FixParameter(2, -27.7088);
				//functionEndcapEcalHcalB->FixParameter(3, 0.76133);
				//functionEndcapEcalHcalB->FixParameter(3, 0.757575);   
				//functionEndcapEcalHcalB->FixParameter(3, 0.758274);   
				functionEndcapEcalHcalB->FixParameter(3, 0.755474);
				//functionEndcapEcalHcalB->FixParameter(4, 0.0776373);
				functionEndcapEcalHcalB->FixParameter(4, 0.0791012);
				//functionEndcapEcalHcalB->FixParameter(5, 7.3809e-10);
				//functionEndcapEcalHcalB->FixParameter(5, 2.6901e-11);
				functionEndcapEcalHcalB->FixParameter(5, 2.6901e-3);
				functionEndcapEcalHcalB->FixParameter(6, 0.158734);
				//functionEndcapEcalHcalB->FixParameter(7, -6.92163);
				functionEndcapEcalHcalB->FixParameter(7, -0.92163);
*/
     // fbEndcap
     
     //for UL2016
    //  //OFFLINE
    //  functionEndcapEcalHcalC->SetParameter(0,-2.24333e+00);                                      
    //  functionEndcapEcalHcalC->SetParameter(1,3.15985e+00);                                        
    //  functionEndcapEcalHcalC->SetParameter(2,2.55316e+00);                                       
    //  functionEndcapEcalHcalC->SetParameter(3,9.53963e+00);                                       
    //  functionEndcapEcalHcalC->SetParameter(4,1.66999e-01);                                       
    //  functionEndcapEcalHcalC->SetParameter(5,2.17192e+01);                                        
    //  functionEndcapEcalHcalC->SetParameter(6,-2.76862e-01);                                      
    //  functionEndcapEcalHcalC->SetParameter(7,0.117280);     

    //OFFLINE 23
    // functionEndcapEcalHcalC->SetParameter(0,-2.25169);
    // functionEndcapEcalHcalC->SetParameter(1,3.15393);
    // functionEndcapEcalHcalC->SetParameter(2,2.2564);
    // functionEndcapEcalHcalC->SetParameter(3,0.970055);
    // functionEndcapEcalHcalC->SetParameter(4,0.0561773);
    // functionEndcapEcalHcalC->SetParameter(5,21.5052);
    // functionEndcapEcalHcalC->SetParameter(6,-0.845143);
    // functionEndcapEcalHcalC->SetParameter(7,0.0898543);
    functionEndcapEcalHcalC->SetParameters(-2.22929,3.12281,6.21902,0.554492,0.0635686,29.5719,-0.567791,0.0968042);

  /*//ONLINE
      functionEndcapEcalHcalC->SetParameters(-0.436687,2.73698,-3.1509,1.20536,
          -1.39685,0.0180331,0.270058,-2.30372);
      //** begin modified by seema  
      //functionEndcapEcalHcalC->FixParameter(0, -0.436687 );
      functionEndcapEcalHcalC->FixParameter(0, -0.43671);
      //functionEndcapEcalHcalC->FixParameter(1, 2.73698);
      functionEndcapEcalHcalC->FixParameter(1, 2.90096);
      //functionEndcapEcalHcalC->FixParameter(2, -3.1509);
      functionEndcapEcalHcalC->FixParameter(2, -5.10099);
      //functionEndcapEcalHcalC->FixParameter(3, 1.20536);
      functionEndcapEcalHcalC->FixParameter(3, 1.20771);
      //functionEndcapEcalHcalC->FixParameter(4, -1.39685);
      functionEndcapEcalHcalC->FixParameter(4, -1.30656);
      //functionEndcapEcalHcalC->FixParameter(5, 0.0180331);
      functionEndcapEcalHcalC->FixParameter(5, 0.0189607);
      //functionEndcapEcalHcalC->FixParameter(6, 0.270058);
      functionEndcapEcalHcalC->FixParameter(6, 0.270027);
      functionEndcapEcalHcalC->FixParameter(7, -2.30372);			  
*/
     // functionEndcapEcalHcalC->FixParameter(0,-4651.06);
     // functionEndcapEcalHcalC->FixParameter(1,4651.98);
     // functionEndcapEcalHcalC->FixParameter(2,12.1135);
     // functionEndcapEcalHcalC->FixParameter(3,503.253);
     // functionEndcapEcalHcalC->FixParameter(4,-1.39108);
     // functionEndcapEcalHcalC->FixParameter(5,-0.247085);
     // functionEndcapEcalHcalC->FixParameter(6,-0.469361);
     // functionEndcapEcalHcalC->FixParameter(7,0.296653);

     functionBarrelHcalA->FixParameter(0, aH);
     functionBarrelHcalB->FixParameter(0, 0.0);

     functionEndcapHcalA->FixParameter(0, aHe);
     functionEndcapHcalB->FixParameter(0, 0.0);

     //fcEndcap 

     //for UL2016
    //  //OFFLINE
    //  functionEndcapHcalC->FixParameter(0,1.53279);
    //  functionEndcapHcalC->FixParameter(1,3.96794);
    //  functionEndcapHcalC->FixParameter(2,-22.6036);
    //  functionEndcapHcalC->FixParameter(3,0.584018);
    //  functionEndcapHcalC->FixParameter(4,0.719891);
    //  functionEndcapHcalC->FixParameter(5,0.0170555);
    //  functionEndcapHcalC->FixParameter(6,0.0990277);
    //  functionEndcapHcalC->FixParameter(7,-1.55789);
    //OFFLINE 23
    //  functionEndcapHcalC->SetParameter(0,1.60863);
    //  functionEndcapHcalC->SetParameter(1,2.56704);
    //  functionEndcapHcalC->SetParameter(2,-15.2573);
    //  functionEndcapHcalC->SetParameter(3,0.772001);
    //  functionEndcapHcalC->SetParameter(4,0.802565);
    //  functionEndcapHcalC->SetParameter(5,0.0215722);
    //  functionEndcapHcalC->SetParameter(6,0.116599);
    //  functionEndcapHcalC->SetParameter(7,-1.48251);
    functionEndcapHcalC->SetParameters(1.57266,6.33265,-48.722,0.539945,0.732707,0.0152813,0.108066,-1.4727);

/*//ONLINE
    functionEndcapHcalC->SetParameters(1.13795,1.21698,-3.81192,115.409,
        0.673456,0.217077,1.95596,-0.252215);
*/


     //trial
     //functionEndcapHcalC->SetParameter(0,1.32375);
     //functionEndcapHcalC->SetParameter(1,0.191687);
     //functionEndcapHcalC->SetParameter(2,-2.81052);
     //functionEndcapHcalC->SetParameter(3,27.1729);
     //functionEndcapHcalC->SetParameter(4,0.418262);
     //functionEndcapHcalC->SetParameter(5,0.107711);
     //functionEndcapHcalC->SetParameter(6,1.71039);
     //functionEndcapHcalC->SetParameter(7,-1.01447);


     

   }
   barrelWithEcalHcalCalib->fitAsToFunction(functionBarrelEcalHcalA);
   //Printing parameters:
   cout << "1" << endl;
   barrelWithEcalHcalCalib->fitBsToFunction(functionBarrelEcalHcalB);
   cout << "2" << endl;
   barrelWithEcalHcalCalib->fitBsToFunction();

    cout << "3" << endl;
   barrelWithEcalHcalCalib->fitCsToFunction(functionBarrelEcalHcalC);

    cout << "4" << endl;
   barrelWithEcalHcalCalib->fitCsToFunction();

    cout << "5" << endl;
   barrelWithEcalHcalCalib->fitCsToFunction();

   endcapWithEcalHcalCalib->fitAsToFunction(functionEndcapEcalHcalA);

   // cout<<"********************************************"<<endl;
   // cout<<"Fit Parameters, functionEndcapEcalHcalB, First"<<endl;
   // for ( unsigned i = 0; i < 10; ++i ) {
   //   double barrelEcalHcalC_spandey = functionEndcapEcalHcalB->GetParameter(i);
   //   if ( barrelEcalHcalC_spandey != 0. )
   //     cout<<"  functionEndcapEcalHcalB,Parameter("<<i<<","<<barrelEcalHcalC_spandey<<");"<<endl;
   // }
    cout << "6" << endl;
   endcapWithEcalHcalCalib->fitBsToFunction(functionEndcapEcalHcalB);

   // cout<<"********************************************"<<endl;
   // cout<<"Fit Parameters, functionEndcapEcalHcalB, First"<<endl;
   // for ( unsigned i = 0; i < 10; ++i ) {
   //   double barrelEcalHcalC_spandey = functionEndcapEcalHcalB->GetParameter(i);
   //   if ( barrelEcalHcalC_spandey != 0. )
   //     cout<<"  functionEndcapEcalHcalB,Parameter("<<i<<","<<barrelEcalHcalC_spandey<<");"<<endl;
   // }
      
      cout << "7" << endl;
   endcapWithEcalHcalCalib->fitBsToFunction();

   // cout<<"********************************************"<<endl;
   // cout<<"Fit Parameters, functionEndcapEcalHcalB, First"<<endl;
   // for ( unsigned i = 0; i < 10; ++i ) {
   //   double barrelEcalHcalC_spandey = functionEndcapEcalHcalB->GetParameter(i);
   //   if ( barrelEcalHcalC_spandey != 0. )
   //     cout<<"  functionEndcapEcalHcalB,Parameter("<<i<<","<<barrelEcalHcalC_spandey<<");"<<endl;
   // }

   //   cout<<"********************************************"<<endl;
   cout << "8" << endl;
   endcapWithEcalHcalCalib->fitCsToFunction(functionEndcapEcalHcalC);
   cout << "9" << endl;
   endcapWithEcalHcalCalib->fitCsToFunction();
    cout << "10" << endl;
   endcapWithEcalHcalCalib->fitCsToFunction();
   //   cout<<"Fit check 1\n";
   barrelWithHcalCalib->fitAsToFunction(functionBarrelHcalA);
   cout << "11" << endl;
   barrelWithHcalCalib->fitBsToFunction(functionBarrelHcalB);
   //   cout<<"Fit check 1.a\n";
   cout << "12" << endl;
   barrelWithHcalCalib->fitBsToFunction();
   //   cout<<"Fit check 1.b\n";
   cout << "13" << endl;
   barrelWithHcalCalib->fitCsToFunction(functionBarrelHcalC);
   cout << "14" << endl;
   barrelWithHcalCalib->fitCsToFunction();
   cout << "15" << endl;
   barrelWithHcalCalib->fitCsToFunction();
   //   cout<<"Fit check 2\n";
   endcapWithHcalCalib->fitAsToFunction(functionEndcapHcalA);
    cout << "16" << endl;
   endcapWithHcalCalib->fitBsToFunction(functionEndcapHcalB);
    cout << "17" << endl;
   endcapWithHcalCalib->fitBsToFunction();

      //cout<<"1 FITTING ENDCAP HCAL C#######"<<endl;
      cout << "18" << endl;
   endcapWithHcalCalib->fitCsToFunction(functionEndcapHcalC);
      // cout<<"2 FITTING ENDCAP HCAL C#######"<<endl;
      cout << "19" << endl;
   endcapWithHcalCalib->fitCsToFunction();
      // cout<<"3 FITTING ENDCAP HCAL C#######"<<endl;
      cout << "20" << endl;
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


     ////////////////////////////////////////////////////


     //faEtaBarrelEH
     //Offline Parameters
    //  functionBarrelAlphaEcalHcal->FixParameter(0,-0.0153451);
    //  functionBarrelAlphaEcalHcal->FixParameter(1,-0.0582379);
    //  functionBarrelAlphaEcalHcal->FixParameter(2,34.1286);
    //  functionBarrelAlphaEcalHcal->FixParameter(3,0.148752);

     /*
     functionBarrelAlphaEcalHcal->SetParameter(0,0.000597082);
     functionBarrelAlphaEcalHcal->SetParameter(1,-0.0510944);
     functionBarrelAlphaEcalHcal->SetParameter(2,73.4545);
     functionBarrelAlphaEcalHcal->SetParameter(3,0.500292);
     */
     /*//Online Parameters
    functionBarrelAlphaEcalHcal->FixParameter(0, -5.08869e-01);
    functionBarrelAlphaEcalHcal->FixParameter(1,  4.03183e+01);
    functionBarrelAlphaEcalHcal->FixParameter(2, -1.14391e-01);
    functionBarrelAlphaEcalHcal->FixParameter(3,  8.05039e-04);
    functionBarrelAlphaEcalHcal->FixParameter(4,  3.98178e+01);
    functionBarrelAlphaEcalHcal->FixParameter(5,  9.81520e-04);
    functionBarrelAlphaEcalHcal->FixParameter(6, -3.67984e+00);
    functionBarrelAlphaEcalHcal->FixParameter(7, -3.58806e+00);
    */
     //fbEtaBarrelEH
     //trial1
     //Offline Parameters
    //  functionBarrelBetaEcalHcal->FixParameter(0,-0.327691);
    //  functionBarrelBetaEcalHcal->FixParameter(1,0.470381);
    //  functionBarrelBetaEcalHcal->FixParameter(2,35.5676);
    //  functionBarrelBetaEcalHcal->FixParameter(3,0.184173);
    //  functionBarrelBetaEcalHcal->FixParameter(4,-0.278133);
    //  functionBarrelBetaEcalHcal->FixParameter(5,0.678245);

    /*//Online Paramters
    functionBarrelBetaEcalHcal->FixParameter(0,  1.89735e-01);
    functionBarrelBetaEcalHcal->FixParameter(1,  9.78327e-01);
    functionBarrelBetaEcalHcal->FixParameter(2, -2.25719e+00);
    functionBarrelBetaEcalHcal->FixParameter(3,  1.40863e-08);
    functionBarrelBetaEcalHcal->FixParameter(4,  1.20856e+00);
    functionBarrelBetaEcalHcal->FixParameter(5,  2.75584e-01);
    functionBarrelBetaEcalHcal->FixParameter(6, -7.89721e+00);
    functionBarrelBetaEcalHcal->FixParameter(7, -4.89917e-01);
     */
     /*
     functionBarrelBetaEcalHcal->SetParameter(0,-0.426957);
     functionBarrelBetaEcalHcal->SetParameter(1,0.425907);
     functionBarrelBetaEcalHcal->FixParameter(2,7.60554);
     functionBarrelBetaEcalHcal->FixParameter(3,0.226156);
     functionBarrelBetaEcalHcal->FixParameter(4,-0.503740);
     functionBarrelBetaEcalHcal->FixParameter(5,0.638326);
     */
    //faEtaBarrelH
    //Offline Parameters
    //  functionBarrelAlphaHcal->FixParameter(0,0.000004304);
    //  functionBarrelAlphaHcal->FixParameter(1,-1.23718e-05);

    /*//Online Parameters
    functionBarrelAlphaHcal->FixParameter(0, -5.78383e+00);
    functionBarrelAlphaHcal->FixParameter(1,  4.29129e+01);
    functionBarrelAlphaHcal->FixParameter(2,  3.37724e+00);
    functionBarrelAlphaHcal->FixParameter(3,  8.38615e-01);
    functionBarrelAlphaHcal->FixParameter(4,  3.71731e+01);
    functionBarrelAlphaHcal->FixParameter(5,  7.30817e-01);
    functionBarrelAlphaHcal->FixParameter(6, -7.70590e-01);
    functionBarrelAlphaHcal->FixParameter(7, -8.15919e-01);
    */
    //fbEtaBarrelH
   //Offline Parameters
  //  functionBarrelBetaHcal->FixParameter(0,-0.000114709);
  //  functionBarrelBetaHcal->FixParameter(1,0.0001779134);
  //  functionBarrelBetaHcal->FixParameter(2,21.9054);

   /*//Online Parameters
  functionBarrelBetaHcal->FixParameter(0, -2.63360e+01);
  functionBarrelBetaHcal->FixParameter(1,  2.65503e+01);
  functionBarrelBetaHcal->FixParameter(2,  2.28575e+01);
  functionBarrelBetaHcal->FixParameter(3,  1.48954e+00);
  functionBarrelBetaHcal->FixParameter(4,  3.86649e-02);
  functionBarrelBetaHcal->FixParameter(5,  1.78269e-04);
  functionBarrelBetaHcal->FixParameter(6, -4.33317e-01);
  functionBarrelBetaHcal->FixParameter(7, -1.54204e+00);
   */
   /*
     functionBarrelBetaHcal->FixParameter(0,-0.0141891);
     functionBarrelBetaHcal->FixParameter(1,0.0493169);
     functionBarrelBetaHcal->FixParameter(2,17.9909);
   */

    //faEtaEndcapEH
   /*
   functionEndcapAlphaEcalHcal->FixParameter(0,-67.0750);
   functionEndcapAlphaEcalHcal->FixParameter(1,67.0703);
   functionEndcapAlphaEcalHcal->FixParameter(2,8.11031e+08);
   functionEndcapAlphaEcalHcal->FixParameter(3,-1.68302e-01);
   functionEndcapAlphaEcalHcal->FixParameter(4,-1.65204e+02);
   */
   
   //trial1
   //Offline Parameters
  //  functionEndcapAlphaEcalHcal->FixParameter(0,0);
  //  functionEndcapAlphaEcalHcal->FixParameter(1,0);

   /*//Online Parameters
  functionEndcapAlphaEcalHcal->FixParameter(0, -7.35050e+01);
  functionEndcapAlphaEcalHcal->FixParameter(1,  1.65774e+02);
  functionEndcapAlphaEcalHcal->FixParameter(2, -8.65271e+01);
  functionEndcapAlphaEcalHcal->FixParameter(3,  4.65056e+01);
  functionEndcapAlphaEcalHcal->FixParameter(4,  8.60982e+01);
  functionEndcapAlphaEcalHcal->FixParameter(5,  7.60022e-01);
  functionEndcapAlphaEcalHcal->FixParameter(6,  6.13895e-02);
  functionEndcapAlphaEcalHcal->FixParameter(7, -5.97798e-01);
  */
   /*
   functionEndcapAlphaEcalHcal->FixParameter(0,0.0226425);
     functionEndcapAlphaEcalHcal->FixParameter(1,0);
     functionEndcapAlphaEcalHcal->FixParameter(2,7.24241);
     functionEndcapAlphaEcalHcal->FixParameter(3,1.11497);
     */

     //fbEtaEndcapEH
   /*
   fbEtaEndcapEH->FixParameter(0,-0.145124);
   fbEtaEndcapEH->FixParameter(1,8.06545);
   fbEtaEndcapEH->FixParameter(2,1.13674e+06);
   fbEtaEndcapEH->FixParameter(3,0.000418838);
   fbEtaEndcapEH->FixParameter(4,-1.01025);
   fbEtaEndcapEH->FixParameter(5,14.5312);
   */
   //Offline Parameters
  //  functionEndcapBetaEcalHcal->FixParameter(0,-0.233182);
  //  functionEndcapBetaEcalHcal->FixParameter(1,0.170723);
  //  functionEndcapBetaEcalHcal->FixParameter(2,997.592);
  //  functionEndcapBetaEcalHcal->FixParameter(3,0.00882785);
  //  functionEndcapBetaEcalHcal->FixParameter(4,-0.913707);
  //  functionEndcapBetaEcalHcal->FixParameter(5,1.78399);

   /*//Online Parameters
  functionEndcapBetaEcalHcal->FixParameter(0, -2.44502e+02);
  functionEndcapBetaEcalHcal->FixParameter(1,  2.45206e+02);
  functionEndcapBetaEcalHcal->FixParameter(2,  2.00784e+01);
  functionEndcapBetaEcalHcal->FixParameter(3,  1.22865e+01);
  functionEndcapBetaEcalHcal->FixParameter(4,  3.03398e-01);
  functionEndcapBetaEcalHcal->FixParameter(5,  3.52688e-04);
  functionEndcapBetaEcalHcal->FixParameter(6, -4.39338e-01);
  functionEndcapBetaEcalHcal->FixParameter(7, -1.90222e+00);
   */
   /*
    functionEndcapBetaEcalHcal->FixParameter(0,-0.187598);
    functionEndcapBetaEcalHcal->FixParameter(1,0.130653);
    functionEndcapBetaEcalHcal->FixParameter(2,50.0069);
    functionEndcapBetaEcalHcal->FixParameter(3,0.0236769);
    functionEndcapBetaEcalHcal->FixParameter(4,-0.868183);
    functionEndcapBetaEcalHcal->FixParameter(5,1.14030);
   */
    //faEtaEndcapH
    //current
    //    functionEndcapAlphaHcal->FixParameter(0,0.0120625);
    //    functionEndcapAlphaHcal->FixParameter(1,1.14891);
    //    functionEndcapAlphaHcal->FixParameter(2,-2.82297);
    //    functionEndcapAlphaHcal->FixParameter(3,2.22078);
    //    functionEndcapAlphaHcal->FixParameter(4,0.522372);

    //Offline Parameters
    // functionEndcapAlphaHcal->FixParameter(0,0);
    // functionEndcapAlphaHcal->FixParameter(1,0);
    // functionEndcapAlphaHcal->FixParameter(2,0);
    // functionEndcapAlphaHcal->FixParameter(3,2.22078);
    // functionEndcapAlphaHcal->FixParameter(4,0.522372);

    /*//Online Parameters
    functionEndcapAlphaHcal->FixParameter(0, -2.12947e+01);
    functionEndcapAlphaHcal->FixParameter(1,  1.58322e+00);
    functionEndcapAlphaHcal->FixParameter(2, -1.27701e+00);
    functionEndcapAlphaHcal->FixParameter(3, -2.97957e-01);
    functionEndcapAlphaHcal->FixParameter(4, -1.88906e+01);
    functionEndcapAlphaHcal->FixParameter(5,  7.39816e-01);
    functionEndcapAlphaHcal->FixParameter(6, -2.00330e-01);
    functionEndcapAlphaHcal->FixParameter(7, -4.36912e-01);
    */

    //    functionEndcapAlphaHcal->FixParameter(0,0.0220633);
    //functionEndcapAlphaHcal->FixParameter(1,2.53937);
    //functionEndcapAlphaHcal->FixParameter(2,-5.90793);
    //functionEndcapAlphaHcal->FixParameter(3,3.02423);
    //functionEndcapAlphaHcal->FixParameter(4,0.723839);
    
    //fbEtaEndcapH
    //Offline Parameters
    // functionEndcapBetaHcal->FixParameter(0,-0.0524166);
    // functionEndcapBetaHcal->FixParameter(1,-0.223219);
    // functionEndcapBetaHcal->FixParameter(2,1.35631);

    /*//Online Parameters
    functionEndcapBetaHcal->FixParameter(0,  3.67453e-02);
    functionEndcapBetaHcal->FixParameter(1,  6.23461e+01);
    functionEndcapBetaHcal->FixParameter(2, -7.61315e+00);
    functionEndcapBetaHcal->FixParameter(3,  5.64491e-02);
    functionEndcapBetaHcal->FixParameter(4,  6.21069e+01);
    functionEndcapBetaHcal->FixParameter(5,  5.60656e-02);
    functionEndcapBetaHcal->FixParameter(6, -1.02854e+00);
    functionEndcapBetaHcal->FixParameter(7, -1.01863e+00);
    */
    //    functionEndcapBetaHcal->SetParameter(0,-0.0772361);
    //    functionEndcapBetaHcal->SetParameter(1,-0.0646546);
    //    functionEndcapBetaHcal->SetParameter(2,1.12264);
    //    functionEndcapBetaHcal->SetParameter(3,0.0155008);
    //    functionEndcapBetaHcal->SetParameter(4,-1.41743);

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

   else {
    //OFFLINE
    //  functionBarrelAlphaEcalHcal->SetParameters(0.02, -0.1, 500);
    //  functionBarrelBetaEcalHcal->SetParameters(1.17842e-01,1.71167e-01,5.88921e+00);

    //  functionBarrelAlphaHcal->SetParameters(0.02, -0.1, 500);
    //  functionBarrelBetaHcal->SetParameters(1.17842e-01,1.71167e-01,5.88921e+00);
    //OFFLINE 23
    //  functionBarrelAlphaEcalHcal->SetParameter(0,-0.0153451);
    //  functionBarrelAlphaEcalHcal->SetParameter(1,-0.0582379);
    //  functionBarrelAlphaEcalHcal->SetParameter(2,34.1286);
    //  functionBarrelAlphaEcalHcal->SetParameter(3,0.148752);
    //functionBarrelAlphaEcalHcal->SetParameters(-1030.75,1030.71,-802293,0.0609009);
    //ONLINE    
				// functionBarrelAlphaEcalHcal->SetParameter(0, -5.08869e-01);
				// functionBarrelAlphaEcalHcal->SetParameter(1,  4.03183e+01);
				// functionBarrelAlphaEcalHcal->SetParameter(2, -1.14391e-01);
				// functionBarrelAlphaEcalHcal->SetParameter(3,  8.05039e-04);
				// functionBarrelAlphaEcalHcal->SetParameter(4,  3.98178e+01);
				// functionBarrelAlphaEcalHcal->SetParameter(5,  9.81520e-04);
				// functionBarrelAlphaEcalHcal->SetParameter(6, -3.67984e+00);
				// functionBarrelAlphaEcalHcal->SetParameter(7, -3.58806e+00);
        functionBarrelAlphaEcalHcal->SetParameters(-0.0297792,41.0737,-1.54161,0.000420402,40.9805,0.000431157,-2.0854,-2.07823);

     //OFFLINE 23
    //  functionBarrelBetaEcalHcal->SetParameter(0,-0.327691);
    //  functionBarrelBetaEcalHcal->SetParameter(1,0.470381);
    //  functionBarrelBetaEcalHcal->SetParameter(2,35.5676);
    //  functionBarrelBetaEcalHcal->SetParameter(3,0.184173);
    //  functionBarrelBetaEcalHcal->SetParameter(4,-0.278133);
    //  functionBarrelBetaEcalHcal->SetParameter(5,0.678245);
    //functionBarrelBetaEcalHcal->SetParameters(-0.809912,3.35686,250.494,0.166035,-0.13914,0.565185);
    //ONLINE
				// functionBarrelBetaEcalHcal->SetParameter(0,  1.89735e-01);
				// functionBarrelBetaEcalHcal->SetParameter(1,  9.78327e-01);
				// functionBarrelBetaEcalHcal->SetParameter(2, -2.25719e+00);
				// functionBarrelBetaEcalHcal->SetParameter(3,  1.40863e-08);
				// functionBarrelBetaEcalHcal->SetParameter(4,  1.20856e+00);
				// functionBarrelBetaEcalHcal->SetParameter(5,  2.75584e-01);
				// functionBarrelBetaEcalHcal->SetParameter(6, -7.89721e+00);
				// functionBarrelBetaEcalHcal->SetParameter(7, -4.89917e-01);
        functionBarrelBetaEcalHcal->SetParameters(0.158827,0.940327,-1.59628,-1.53431e-09,1.20873,0.489381,-17.5062,-0.367684);




     //OFFLINE 23
    //  functionBarrelAlphaHcal->SetParameter(0,0.000004304);
    //  functionBarrelAlphaHcal->SetParameter(1,-1.23718e-05);
      //functionBarrelAlphaHcal->SetParameters(0.0148632,7.4178e-07);
      //ONLINE
				// functionBarrelAlphaHcal->SetParameter(0, -5.78383e+00);
				// functionBarrelAlphaHcal->SetParameter(1,  4.29129e+01);
				// functionBarrelAlphaHcal->SetParameter(2,  3.37724e+00);
				// functionBarrelAlphaHcal->SetParameter(3,  8.38615e-01);
				// functionBarrelAlphaHcal->SetParameter(4,  3.71731e+01);
				// functionBarrelAlphaHcal->SetParameter(5,  7.30817e-01);
				// functionBarrelAlphaHcal->SetParameter(6, -7.70590e-01);
				// functionBarrelAlphaHcal->SetParameter(7, -8.15919e-01);
          functionBarrelAlphaHcal->SetParameters(-5.79215,42.9144,6.11378,0.763853,37.1711,0.680129,-0.729795,-0.794311);


    //OFFLINE 23
      // functionBarrelBetaHcal->SetParameter(0,-0.000114709);
      // functionBarrelBetaHcal->SetParameter(1,0.0001779134);
      // functionBarrelBetaHcal->SetParameter(2,21.9054);
      //functionBarrelBetaHcal->SetParameters(-0.0170179,-2.13782,0.906182);
      // //ONLINE
			// 	functionBarrelBetaHcal->SetParameter(0, -2.63360e+01);
			// 	functionBarrelBetaHcal->SetParameter(1,  2.65503e+01);
			// 	functionBarrelBetaHcal->SetParameter(2,  2.28575e+01);
			// 	functionBarrelBetaHcal->SetParameter(3,  1.48954e+00);
			// 	functionBarrelBetaHcal->SetParameter(4,  3.86649e-02);
			// 	functionBarrelBetaHcal->SetParameter(5,  1.78269e-04);
			// 	functionBarrelBetaHcal->SetParameter(6, -4.33317e-01);
			// 	functionBarrelBetaHcal->SetParameter(7, -1.54204e+00);
        functionBarrelBetaHcal->SetParameters(-26.2578,26.6346,19.0027,1.739,0.302043,0.0645809,-0.423556,-0.469157);//, 3., 6., 9.);

     //faEtaBarrelEH                                                                                                                                                       
     //functionBarrelAlphaEcalHcal->SetParameter(0,340.572);
     //functionBarrelAlphaEcalHcal->SetParameter(1,-340.606);
     //functionBarrelAlphaEcalHcal->SetParameter(2,6689.73);
     //functionBarrelAlphaEcalHcal->SetParameter(3,0.000886848);

    //fbEtaBarrelEH                                                                                                                                                        
    //functionBarrelBetaEcalHcal->SetParameter(0,0.0532902);
    //functionBarrelBetaEcalHcal->SetParameter(1,0.124218);
    //functionBarrelBetaEcalHcal->SetParameter(2,207.676);

    //faEtaBarrelH                                                                                                                                                         
    //functionBarrelAlphaHcal->SetParameter(0,-0.00157472);
    //functionBarrelAlphaHcal->SetParameter(1,0.00000);

    //fbEtaBarrelH                                                                                                                                                         
    //functionBarrelBetaHcal->SetParameter(0,-0.0154675);
    //functionBarrelBetaHcal->SetParameter(1,0.0728367);
    //functionBarrelBetaHcal->SetParameter(2,56.2748);


   //faEtaEndcapEH
   //OFFLINE
  //  functionEndcapAlphaEcalHcal->SetParameter(0,0.00134831);
  //  functionEndcapAlphaEcalHcal->SetParameter(1,-26.4964);
  //  functionEndcapAlphaEcalHcal->SetParameter(2,0.206669);
  //OFFLINE 23
  //  functionEndcapAlphaEcalHcal->SetParameter(0,0);
  //  functionEndcapAlphaEcalHcal->SetParameter(1,0);
    //functionEndcapAlphaEcalHcal->SetParameters(0.0250812,-3.41435e-05);
  //ONLINE
				// functionEndcapAlphaEcalHcal->SetParameter(0, -7.35050e+01);
				// functionEndcapAlphaEcalHcal->SetParameter(1,  1.65774e+02);
				// functionEndcapAlphaEcalHcal->SetParameter(2, -8.65271e+01);
				// functionEndcapAlphaEcalHcal->SetParameter(3,  4.65056e+01);
				// functionEndcapAlphaEcalHcal->SetParameter(4,  8.60982e+01);
				// functionEndcapAlphaEcalHcal->SetParameter(5,  7.60022e-01);
				// functionEndcapAlphaEcalHcal->SetParameter(6,  6.13895e-02);
				// functionEndcapAlphaEcalHcal->SetParameter(7, -5.97798e-01);
        functionEndcapAlphaEcalHcal->SetParameters(-74.2096,165.037,-35.8544,42.0327,86.8414,2.24809,0.00177795,-0.520057);

   /*//ONLINE
	functionEndcapAlphaEcalHcal->SetParameters(0.02, -0.1, sampleRangeHigh);
*/
   //fbEtaEndcapEH                                                                                                                                
   /*//ONLINE  
	 functionEndcapBetaEcalHcal->SetParameters(0.0399873, -1.51747, 3.22236);   
   */
  //  //OFFLINE
  //  functionEndcapBetaEcalHcal->SetParameter(0,0.0476086);
  //  functionEndcapBetaEcalHcal->SetParameter(1,0.000573041);
  //  functionEndcapBetaEcalHcal->SetParameter(2,133.97);
  //  functionEndcapBetaEcalHcal->SetParameter(3,0.798439);
  //OFFLINE 23
  //  functionEndcapBetaEcalHcal->SetParameter(0,-0.233182);
  //  functionEndcapBetaEcalHcal->SetParameter(1,0.170723);
  //  functionEndcapBetaEcalHcal->SetParameter(2,997.592);
  //  functionEndcapBetaEcalHcal->SetParameter(3,0.00882785);
  //  functionEndcapBetaEcalHcal->SetParameter(4,-0.913707);
  //  functionEndcapBetaEcalHcal->SetParameter(5,1.78399);
   //functionEndcapBetaEcalHcal->SetParameters(-0.0936792,0.0299816,0.319867,1.87256e-07,-4.59142,0.178485);
   //ONLINE
				// functionEndcapBetaEcalHcal->SetParameter(0, -2.44502e+02);
				// functionEndcapBetaEcalHcal->SetParameter(1,  2.45206e+02);
				// functionEndcapBetaEcalHcal->SetParameter(2,  2.00784e+01);
				// functionEndcapBetaEcalHcal->SetParameter(3,  1.22865e+01);
				// functionEndcapBetaEcalHcal->SetParameter(4,  3.03398e-01);
				// functionEndcapBetaEcalHcal->SetParameter(5,  3.52688e-04);
				// functionEndcapBetaEcalHcal->SetParameter(6, -4.39338e-01);
				// functionEndcapBetaEcalHcal->SetParameter(7, -1.90222e+00);
        functionEndcapBetaEcalHcal->SetParameters(-244.79,244.912,14.5169,12.9739,0.189224,0.00310742,-0.569849,-1.45211);//, 2., 6., 9.);
   
   //faEtaEndcapH                                                                                                                                                      //OFFLINE  
  //  functionEndcapAlphaHcal->SetParameter(0,-6.05737e-03);
  //  functionEndcapAlphaHcal->SetParameter(1,-7.20466e-01);
  //  functionEndcapAlphaHcal->SetParameter(2,1.75442e-01);
  //  functionEndcapAlphaHcal->SetParameter(3,1.45679e+01);
  //OFFLINE 23
    // functionEndcapAlphaHcal->SetParameter(0,0);
    // functionEndcapAlphaHcal->SetParameter(1,0);
    // functionEndcapAlphaHcal->SetParameter(2,0);
    // functionEndcapAlphaHcal->SetParameter(3,2.22078);
    // functionEndcapAlphaHcal->SetParameter(4,0.522372);
    //functionEndcapAlphaHcal->SetParameters(0.00844682,0.839335,-2.01159,2.99437,0.606582);
    //ONLINE
				// functionEndcapAlphaHcal->SetParameter(0, -2.12947e+01);
				// functionEndcapAlphaHcal->SetParameter(1,  1.58322e+00);
				// functionEndcapAlphaHcal->SetParameter(2, -1.27701e+00);
				// functionEndcapAlphaHcal->SetParameter(3, -2.97957e-01);
				// functionEndcapAlphaHcal->SetParameter(4, -1.88906e+01);
				// functionEndcapAlphaHcal->SetParameter(5,  7.39816e-01);
				// functionEndcapAlphaHcal->SetParameter(6, -2.00330e-01);
				// functionEndcapAlphaHcal->SetParameter(7, -4.36912e-01);
        functionEndcapAlphaHcal->SetParameters(-21.405,1.64842,-0.832719,-0.358038,-18.5864,0.779407,-0.116406,-0.296534);
  
   /*//ONLINE
  functionEndcapAlphaHcal->SetParameters(0.02, -0.1, sampleRangeHigh, 0.5, 0.6);
*/
   
   //fbEtaEndcapH                              
   //OFFLINE                                                                      
  //  functionEndcapBetaHcal->SetParameter(0,0.0564225);
  //  functionEndcapBetaHcal->SetParameter(1,0.00725947);
  //  functionEndcapBetaHcal->SetParameter(2,23.1278);
  //OFFLINE 23
    // functionEndcapBetaHcal->SetParameter(0,-0.0524166);
    // functionEndcapBetaHcal->SetParameter(1,-0.223219);
    // functionEndcapBetaHcal->SetParameter(2,1.35631);
    //functionEndcapBetaHcal->SetParameters(0.0375576,-6.99123e+06,0.135152);
    //ONLINE
				// functionEndcapBetaHcal->SetParameter(0,  3.67453e-02);
				// functionEndcapBetaHcal->SetParameter(1,  6.23461e+01);
				// functionEndcapBetaHcal->SetParameter(2, -7.61315e+00);
				// functionEndcapBetaHcal->SetParameter(3,  5.64491e-02);
				// functionEndcapBetaHcal->SetParameter(4,  6.21069e+01);
				// functionEndcapBetaHcal->SetParameter(5,  5.60656e-02);
				// functionEndcapBetaHcal->SetParameter(6, -1.02854e+00);
				// functionEndcapBetaHcal->SetParameter(7, -1.01863e+00);
          functionEndcapBetaHcal->SetParameters(0.00191295,62.9369,-0.527866,0.0798405,62.8668,0.0785216,-0.673072,-0.676033);

   /*//ONLINE
 	functionEndcapBetaHcal->SetParameters(0.07, -2.5, 6.0, 0.3, 175.0);
  */
     //////////////
     //   functionEndcapAlphaEcalHcal->SetParameters(0.02, -0.1, 500);
     //   functionEndcapBetaEcalHcal->SetParameters(0.0399873, -1.51747, 3.22236);

     //   functionEndcapAlphaHcal->SetParameters(0.02, -0.1, 500, 0.5, 0.6);
     //   functionEndcapBetaHcal->SetParameters(0.07, -2.5, 6.0, 0.3, 175.0);
      //functionEndcapGammaHcal->SetParameters(-1.71458e+00,9.61337e+00,2.94747e+00);

   }


    cout << "21" << endl;
   barrelWithEcalHcalCalib->fitAlphasToFunction(functionBarrelAlphaEcalHcal);
    cout << "22" << endl;
   barrelWithEcalHcalCalib->fitAlphasToFunction();
    cout << "23" << endl;
   barrelWithEcalHcalCalib->fitBetasToFunction(functionBarrelBetaEcalHcal);
    cout << "24" << endl;
   barrelWithEcalHcalCalib->fitBetasToFunction();
   endcapWithEcalHcalCalib->fitAlphasToFunction(functionEndcapAlphaEcalHcal);
   endcapWithEcalHcalCalib->fitAlphasToFunction();
   endcapWithEcalHcalCalib->fitBetasToFunction(functionEndcapBetaEcalHcal);
   endcapWithEcalHcalCalib->fitBetasToFunction();

    cout << "25" << endl;
   barrelWithHcalCalib->fitAlphasToFunction(functionBarrelAlphaHcal);
    cout << "26" << endl;
   barrelWithHcalCalib->fitAlphasToFunction();
    cout << "27" << endl;
   barrelWithHcalCalib->fitBetasToFunction(functionBarrelBetaHcal);
    cout << "28" << endl;
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
    double etrue,ecal,hcal,eta,phi,correctedEta, correctedE, PFHCclosure, PFECclosure;
    
    TFile *outFileN = new TFile("NTuple_0_500_PFHC_PFEC_24MC24Corrections.root","recreate");
    TTree *Ntree = new TTree("s","NTuple for energy btw 0 and 500 GeV and PFHC and PFEC energies");
    if (WriteNTupleFile){

      Ntree->Branch("true",&etrue);
      Ntree->Branch("ecal",&ecal);
      Ntree->Branch("hcal",&hcal);
      Ntree->Branch("eta",&eta);
      Ntree->Branch("phi",&phi);
      Ntree->Branch("PFHC_energy",&correctedEta);
      Ntree->Branch("PFEC_energy",&correctedE);
      Ntree->Branch("PFHC_closure",&PFHCclosure);
      Ntree->Branch("PFEC_closure",&PFECclosure);
    }

    for (unsigned entry = 0; entry < N; ++entry) {
        contador++;
        etrue = ETrueEnergies[entry];
        ecal = ecalEnergies[entry];
        hcal = hcalEnergies[entry];
        eta = std::abs(etas[entry]);
        phi = phis[entry];

        // Condiciones de filtro
        if ((ecal + hcal) < 0.5 || etrue < 1.0 || hcal == 0) continue;

        correctedEta = 0;
        double eecalcorr = ecal;
        double ehcalcorr = hcal;
        double etrue_org=etrue; //etrue for charge hadrons (max{etrue, ecal+hcal}) or -1 for neutral hadrons (max{-1, ecal+hcal})


        if (!usePFEnergyCalibration) {
            if (eta < 1.3) {
                if(ecal > 0) correctedEta = barrelWithEcalHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, eta, 0);
                else correctedEta = barrelWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, eta, 0);
                corrEta_range1->Fill(etrue, (correctedEta - etrue) / etrue);
            } else {
                if (ecal > 0 ) correctedEta = endcapWithEcalHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, eta, 0);
                else correctedEta = endcapWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, eta, 0);
                if (eta < 2.4) corrEta_range2->Fill(etrue, (correctedEta - etrue) / etrue);
                else if (eta < 2.7) corrEta_range3->Fill(etrue, (correctedEta - etrue) / etrue);
                else if (eta < 3.0) corrEta_range4->Fill(etrue, (correctedEta - etrue) / etrue);
            }
        } else {
            pec->energyEmHad(etrue_org, eecalcorr, ehcalcorr, eta, phi);
            correctedEta = eecalcorr + ehcalcorr;
            if (eta < 1.3) {
                corrEta_range1->Fill(etrue, (correctedEta - etrue) / etrue);
            } else if (eta < 2.4) {
                corrEta_range2->Fill(etrue, (correctedEta - etrue) / etrue);
            } else if (eta < 2.7) {
                corrEta_range3->Fill(etrue, (correctedEta - etrue) / etrue);
            } else if (eta < 3.0) {
                corrEta_range4->Fill(etrue, (correctedEta - etrue) / etrue);
            }
        }

        if (WriteNTupleFile){
          pec->energyEmHad(etrue, eecalcorr, ehcalcorr, eta, phi);
          correctedE = eecalcorr + ehcalcorr;
          PFECclosure = (correctedE - etrue) / etrue;
          if (eta<1.5) {
            if (ecal > 0) correctedEta = barrelWithEcalHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, eta, 0);
            else correctedEta = barrelWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, eta, 0);
          }else{
            if (ecal > 0) correctedEta = endcapWithEcalHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, eta, 0);
            else correctedEta = endcapWithHcalCalib->getCalibratedEnergy(etrue, ecal, hcal, eta, 0);
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
    if (Ntree != nullptr) {
      delete Ntree;
    }
    if (outFileN != nullptr) {
      delete outFileN;
    }
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
   cout<<"  faBarrel = new TF1(\"faBarrel\",\""<<
     functionBarrelEcalHcalB_e<<"\",1.,1000.);"<<endl;
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
   cout<<"  fbBarrel = new TF1(\"fbBarrel\",\""<<
     functionBarrelEcalHcalC_e<<"\",1.,1000.);"<<endl;
   for ( int i = 0; i < 100; ++i ) 
     {
       if ( i > (functionBarrelEcalHcalC->GetNpar()-1)) break;
       barrelEcalHcalC = functionBarrelEcalHcalC->GetParameter(i);
       if(payload) cout<<barrelEcalHcalC<<",";
       else cout<<"  fbBarrel->SetParameter("<<i<<","<<barrelEcalHcalC<<");"<<endl;
     }
   cout<<"  fbBarrel at x = 10 ;  "<<","<<functionBarrelEcalHcalC->Eval(10.)<<endl<<endl;

   // Hcal coef barrel for Hcal
   cout<<"  fcBarrel = new TF1(\"fcBarrel\",\""<<
     functionBarrelHcalC_e<<"\",1.,1000.);"<<endl;

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
   cout<<"  faEtaBarrelEH = new TF1(\"faEtaBarrelEH\",\""<<
     functionBarrelAlphaEH_e<<"\",1.,1000.);"<<endl;
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
   cout<<"  fbEtaBarrelEH = new TF1(\"fbEtaBarrelEH\",\""<<
     functionBarrelBetaEH_e<<"\",1.,1000.);"<<endl;
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
   cout<<"  faEtaBarrelH = new TF1(\"faEtaBarrelH\",\""<<
     functionBarrelAlphaH_e<<"\",1.,1000.);"<<endl;
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
   cout<<"  fbEtaBarrelH = new TF1(\"fbEtaBarrelH\",\""<<
     functionBarrelBetaH_e<<"\",1.,1000.);"<<endl;
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
   cout<<"  faEndcap = new TF1(\"faEndcap\",\""<<
     functionEndcapEcalHcalB_e<<"\",1.,1000.);"<<endl;
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
   cout<<"  fbEndcap = new TF1(\"fbEndcap\",\""<<
     functionEndcapEcalHcalC_e<<"\",1.,1000.);"<<endl;
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
   cout<<"  fcEndcap = new TF1(\"fcEndcap\",\""<<
     functionEndcapHcalC_e<<"\",1.,1000.);"<<endl;
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
   cout<<"  faEtaEndcapEH = new TF1(\"faEtaEndcapEH\",\""<<
     functionEndcapAlphaEH_e<<"\",1.,1000.);"<<endl;
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
   cout<<"  fbEtaEndcapEH = new TF1(\"fbEtaEndcapEH\",\""<<
     functionEndcapBetaEH_e<<"\",1.,1000.);"<<endl;
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
   cout<<"  faEtaEndcapH = new TF1(\"faEtaEndcapH\",\""<<
     functionEndcapAlphaH_e<<"\",1.,1000.);"<<endl;
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
   cout<<"  fbEtaEndcapH = new TF1(\"fbEtaEndcapH\",\""<<
     functionEndcapBetaH_e<<"\",1.,1000.);"<<endl;
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

   return 1;
}



