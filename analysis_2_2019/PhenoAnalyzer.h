////////////////////////////////////////////////////////////////
//                                                            //
// Author: Andrés Flórez, Universidad de los Andes, Colombia  //
//                                                            //
////////////////////////////////////////////////////////////////

#ifndef PHENOANALYZER_H
#define PHENOANALYZER_H

#include "ROOTFunctions.h"
#include "DelphesFunctions.h"


using namespace std;

class PhenoAnalysis {
public :
   PhenoAnalysis(TChain&, TFile*, TDirectory* dir[], int nDir);
   ~PhenoAnalysis();
   void crateHistoMasps (int);
   bool overlapingObjects(double, double, double, double, double);
   double calculateE(double, double, double);
   double normalizedDphi(double);
   double calculate_deltaR(TLorentzVector,  Track*);

   // For Jets
   std::map<unsigned int, TH1*> _hmap_Nevents;
   std::map<unsigned int, TH1*> _hmap_lead_jet_pT;
   std::map<unsigned int, TH1*> _hmap_lead_jet_eta;
   std::map<unsigned int, TH1*> _hmap_lead_jet_phi;
   std::map<unsigned int, TH1*> _hmap_n_jets;
   std::map<unsigned int, TH1*> _hmap_slead_jet_pT;
   std::map<unsigned int, TH1*> _hmap_slead_jet_eta;
   std::map<unsigned int, TH1*> _hmap_slead_jet_phi;
   // For Taus
   std::map<unsigned int, TH1*> _hmap_tau1_pT;
   std::map<unsigned int, TH1*> _hmap_tau1_eta;
   std::map<unsigned int, TH1*> _hmap_tau1_phi;
   std::map<unsigned int, TH1*> _hmap_tau2_pT;
   std::map<unsigned int, TH1*> _hmap_tau2_eta;
   std::map<unsigned int, TH1*> _hmap_tau2_phi;
   std::map<unsigned int, TH1*> _hmap_Delta_pT;
   // For Muons
   std::map<unsigned int, TH1*> _hmap_muon1_pT;
   std::map<unsigned int, TH1*> _hmap_muon1_eta;
   std::map<unsigned int, TH1*> _hmap_muon1_phi;
   std::map<unsigned int, TH1*> _hmap_muon2_pT;
   std::map<unsigned int, TH1*> _hmap_muon2_eta;
   std::map<unsigned int, TH1*> _hmap_muon2_phi;
   // For Electrons
   std::map<unsigned int, TH1*> _hmap_elec1_pT;
   std::map<unsigned int, TH1*> _hmap_elec1_eta;
   std::map<unsigned int, TH1*> _hmap_elec1_phi;
   std::map<unsigned int, TH1*> _hmap_elec2_pT;
   std::map<unsigned int, TH1*> _hmap_elec2_eta;
   std::map<unsigned int, TH1*> _hmap_elec2_phi;
   //Gen
   std::map<unsigned int, TH1*> _hmap_Gentau1_pT;
   std::map<unsigned int, TH1*> _hmap_GenDelta_pT;
   std::map<unsigned int, TH1*> _hmap_Gentau1_eta;
   std::map<unsigned int, TH1*> _hmap_Gentau1_phi;
   std::map<unsigned int, TH1*> _hmap_Gentau2_pT;
   std::map<unsigned int, TH1*> _hmap_Gentau2_eta;
   std::map<unsigned int, TH1*> _hmap_Gentau2_phi;
   std::map<unsigned int, TH1*> _hmap_Genst;   

   std::map<unsigned int, TH1*> _hmap_n_tau;
   std::map<unsigned int, TH1*> _hmap_ipTau1;

   // Topology
   std::map<unsigned int, TH1*> _hmap_ht;
   std::map<unsigned int, TH1*> _hmap_st;
   std::map<unsigned int, TH1*> _hmap_dijet_mass;
   std::map<unsigned int, TH1*> _hmap_dijet_deltaEta;
   std::map<unsigned int, TH1*> _hmap_tauMass;
   std::map<unsigned int, TH1*> _hmap_MET;
   std::map<unsigned int, TH1*> _hmap_tau_transmass;
   std::map<unsigned int, TH1*> _hmap_muon_transmass;
   std::map<unsigned int, TH1*> _hmap_elec_transmass;

   //2D plots
   std::map<unsigned int, TH2*> _hmap_jet_met_metDphi;
   std::map<unsigned int, TH2*> _hmap_Dphi_tau_met_with_met;
   std::map<unsigned int, TH2*> _hmap_Dphi_tau_jet_with_jPt;
   std::map<unsigned int, TH2*> _hmap_MT_with_MET;
   std::map<unsigned int, TH2*> _hmap_MT_with_Dphi_tau_met;
   std::map<unsigned int, TH2*> _hmap_MT_with_Dphi_jet_met;
   std::map<unsigned int, TH2*> _hmap_Meff_with_MET;
};

#endif
