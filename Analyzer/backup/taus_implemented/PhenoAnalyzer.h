////////////////////////////////////////////////////////////////
//                                                            //
// Author: Andrés Flórez, Universidad de los Andes, Colombia  //
//                                                            //
////////////////////////////////////////////////////////////////

#ifndef _PHENOANALYZER_H_
#define _PHENOANALYZER_H_

#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TEnv.h"
#include "TRandom.h"

#include <iostream>
#include <utility>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"

#include "TString.h"

#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TApplication.h"
#include "DelphesFunctions.h"
#include "TDirectory.h"
#include "TFile.h"
#include <fstream>

using namespace std;

class PhenoAnalysis {
public :
   PhenoAnalysis(TChain&, TFile*, TDirectory* dir[], int nDir);
   ~PhenoAnalysis();
   void crateHistoMasps (int);
   bool overlapingObjects(double, double, double, double, double);
   double calculateE(double, double, double);
   double normalizedDphi(double);
   bool TauID(TLorentzVector TauCand);
   bool TauID2(TLorentzVector TauCand);
   bool TauIDJet(TLorentzVector jet);
   TLorentzVector TauSmearing(TLorentzVector TauCand);

   // For Jets
   std::map<unsigned int, TH1*> _hmap_Nevents;
   std::map<unsigned int, TH1*> _hmap_lead_jet_pT;
   std::map<unsigned int, TH1*> _hmap_lead_jet_eta;
   std::map<unsigned int, TH1*> _hmap_lead_jet_phi;
   std::map<unsigned int, TH1*> _hmap_RM;
   
   std::map<unsigned int, TH1*> _hmap_n_jets;
   // For Taus
   std::map<unsigned int, TH1*> _hmap_tau1_pT;
   std::map<unsigned int, TH1*> _hmap_tau1_eta;
   std::map<unsigned int, TH1*> _hmap_tau1_phi;
   std::map<unsigned int, TH1*> _hmap_tau2_pT;
   std::map<unsigned int, TH1*> _hmap_tau2_eta;
   std::map<unsigned int, TH1*> _hmap_tau2_phi;
   std::map<unsigned int, TH1*> _hmap_n_tau;
   // Topology
   std::map<unsigned int, TH1*> _hmap_jet_met_Dphi;
   std::map<unsigned int, TH1*> _hmap_transverse_mass;
   std::map<unsigned int, TH1*> _hmap_tau1_met_Dphi;
   std::map<unsigned int, TH1*> _hmap_tau2_met_Dphi;
   std::map<unsigned int, TH1*> _hmap_met;
   std::map<unsigned int, TH1*> _hmap_tau1_jet_Dphi;
   std::map<unsigned int, TH1*> _hmap_efective_mass;
   //2D plots
   std::map<unsigned int, TH2*> _hmap_jet_met_metDphi;
private :

};

#endif
