////////////////////////////////////////////////////////////////
//                                                            //
// Author: Andrés Flórez, Universidad de los Andes, Colombia  //
//                                                            //
////////////////////////////////////////////////////////////////



#include <iostream>
#include "ROOTFunctions.h"
#include "DelphesFunctions.h"
#include "PhenoAnalyzer.h"
#include <time.h>

int main(int argc, char *argv[]) {
  
  //TApplication app("App",&argc, argv);
  TChain chain("Delphes");
  chain.Add(argv[1]);
  TFile * HistoOutputFile = new TFile(argv[2], "RECREATE");
  int nDir = 23;
  TDirectory *theDirectory[nDir];
  theDirectory[0]  = HistoOutputFile->mkdir("Staus_No_cuts");
  theDirectory[1]  = HistoOutputFile->mkdir("Staus_with_one_Tau");
  theDirectory[2]  = HistoOutputFile->mkdir("Staus_with_two_Taus");
  theDirectory[3]  = HistoOutputFile->mkdir("Staus_After_Tau_pt_min_1Tau");
  theDirectory[4]  = HistoOutputFile->mkdir("Staus_After_Tau_pt_min_2Taus");
  theDirectory[5]  = HistoOutputFile->mkdir("Staus_After_Tau_Eta_min_1Tau");
  theDirectory[6]  = HistoOutputFile->mkdir("Staus_After_Tau_Eta_min_2Taus");
  theDirectory[7]  = HistoOutputFile->mkdir("Staus_After_b_jet_veto_1Tau");
  theDirectory[8]  = HistoOutputFile->mkdir("Staus_After_b_jet_veto_2Taus");
  theDirectory[9]  = HistoOutputFile->mkdir("Staus_After_jet_pt_cut_1Tau");
  theDirectory[10]  = HistoOutputFile->mkdir("Staus_After_jet_pt_cut_2Taus");
  theDirectory[11]  = HistoOutputFile->mkdir("Staus_After_jet_Eta_cut_1Taus");
  theDirectory[12]  = HistoOutputFile->mkdir("Staus_After_jet_Eta_cut_2Taus");
  theDirectory[13]  = HistoOutputFile->mkdir("Staus_Dphi_tau_jet_cut_1Tau");
  theDirectory[14]  = HistoOutputFile->mkdir("Staus_Dphi_tau_jet_cut_2Taus");
  theDirectory[15]  = HistoOutputFile->mkdir("Staus_Transverse_mass_cut_1Tau");
  theDirectory[16]  = HistoOutputFile->mkdir("Staus_Transverse_mass_cut_2Taus");
  theDirectory[17]  = HistoOutputFile->mkdir("Staus_After_MET_1Tau");
  theDirectory[18]  = HistoOutputFile->mkdir("Staus_After_MET_2Taus");
  theDirectory[19]  = HistoOutputFile->mkdir("Staus_Dphi_jet_met_1Tau");
  theDirectory[20]  = HistoOutputFile->mkdir("Staus_Dphi_jet_met_2Taus"); 
  theDirectory[21]  = HistoOutputFile->mkdir("Staus_Efective_mass_cut_1Tau");
  theDirectory[22]  = HistoOutputFile->mkdir("Staus_Efective_mass_cut_2Taus");
  PhenoAnalysis BSM_analysis(chain, HistoOutputFile, theDirectory, nDir);
  
}

using namespace std;
PhenoAnalysis::PhenoAnalysis(TChain& chain, TFile* theFile, TDirectory *cdDir[], int nDir)
{
  ifstream inFile;
  inFile.open ("config.in", ios::in);
  
  if (!inFile)
    {
      cerr << "ERROR: Can't open input file: " << endl;
      exit (1);
    }
  
  string inputType = "";
  
  //This set of lines are used to open and read the "config.in" file. 
  /////////////////////////////////////////////////////////////////////// 
  TEnv *params = new TEnv ("config_file");
  params->ReadFile ("config.in", kEnvChange);
  
  int BGQCD               = params->GetValue ("BGQCD", 0);
  double lead_jet_pt      = params->GetValue ("lead_jet_pt", 100.);
  double jet_eta_max      = params->GetValue ("jet_eta_max", 5.0);
  double tau1_pt_min      = params->GetValue ("tau1_pt_min", 10.);
  double tau2_pt_min      = params->GetValue ("tau2_pt_min", 10.);
  double tau_eta_max      = params->GetValue ("tau_eta_max", 2.1);
  double b_jet_pt_min     = params->GetValue ("b_jet_pt_min", 20.0);
  double DR_jet_elec_max  = params->GetValue ("DR_jet_elec_max", 0.3);
  double DR_taui_tauj_min = params->GetValue ("DR_taui_tauj_min", 0.1);
  double met_min          = params->GetValue ("met_min", 10.);
  double deltajetmetphi   = params->GetValue ("deltajetmetphi",1.5);
  double transversemassmin = params->GetValue ("transversemassmin",0.);
  double efectivemassmin =  params->GetValue ("efectivemassmin",0.);
  double Dphitaujetmax   =  params->GetValue ("Dphitaujetmax",0.); 
  crateHistoMasps(nDir);
  
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");  
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle"); 
  
  MissingET *METpointer; 
  TH1 *Numbertaus_D1     = new TH1F("Numbertaus_D1",  "Numbertaus_D1", 30, -5.5, 24.5); 		 
  TH1 *Numbertaus_D2     = new TH1F("Numbertaus_D2",  "Numbertaus_D2", 10, 0, 10);
  TH1 *Numbertaus_2      = new TH1F("Numbertaus_2",   "Numbertaus_D2", 30, -5.5, 24.5);
  TH1 *Mothers           = new TH1F("Mothers",        "Mothers", 510, -10, 500);
  
  int counter_elec_muon  = 0;
  
  
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    //for(Int_t entry = 0; entry < 10; ++entry)
    {
      //cout <<"-----------Evento: "<<entry<<"-------------"<<endl;
      int pass_cuts[nDir] = {0};
      TLorentzVector Jet_leading_vec(0., 0., 0., 0.);
      TLorentzVector Tau1Had_vec(0., 0., 0., 0.);
      TLorentzVector Tau2Had_vec(0., 0., 0., 0.);
      TLorentzVector Tau3Had_vec(0., 0., 0., 0.);
      TLorentzVector Tau4Had_vec(0., 0., 0., 0.);
      TLorentzVector Tau1HadTLV (0., 0., 0., 0.);
      TLorentzVector Tau2HadTLV (0., 0., 0., 0.);    
      TLorentzVector Tau3HadTLV (0., 0., 0., 0.);
      TLorentzVector Tau4HadTLV (0., 0., 0., 0.);
      TLorentzVector MCnu(0., 0., 0., 0.);
      bool is_b_jet = false;
      treeReader->ReadEntry(entry);
      int tau_conter=0;
      METpointer = (MissingET*) branchMissingET->At(0);
      double MET = METpointer->MET;
      double MET_phi = METpointer->Phi;
      int njets_counter = 0;
      int ntau_counter = 0; 
      int countertau=0;
      //int neutrinos=0;
      //bool passed_jet_tau = false;
      bool fill_tau1 = false;
      bool fill_tau2 = false;
      bool fill_tau3 = false;
      bool fill_tau4 = false;
      double jet_min_pt = 10.;
      
      // We need at least 2 jets, 1 (2) tau jets and the ISR jet.
      
      TLorentzVector jet_i(0., 0., 0., 0.);
      TLorentzVector elec_i(0., 0., 0., 0.);
      double N_muons = 0.;
      for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++){
	Muon *muon = (Muon*) branchMuon->At(muo);
	if ((muon->PT > 10.) && (abs(muon->Eta) < 2.5)){N_muons++;}            
      }
      
      //////Jet-tau fake rate///////////////////
      
      
      for (int j = 0; j < branchJet->GetEntriesFast(); j++) {
	bool passed_jet_tau = false;
	Jet *jet = (Jet*) branchJet->At(j);
	double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
	jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
	passed_jet_tau = TauIDJet(jet_i);
	//cout << "passed_jet_tau "<<passed_jet_tau<<endl;
	if((jet->PT < 20.) || (abs(jet->Eta) > 2.5)){continue;}
	if((passed_jet_tau == true) && (fill_tau1 == false)){
	  Tau1HadTLV.SetPtEtaPhiE(jet_i.Pt(), jet_i.Eta(), jet_i.Phi(), jet_i.E());
	  fill_tau1 = true;
	  continue;
	}
	if((passed_jet_tau == true) && (fill_tau2 == false) && (fill_tau1 == true)){
	  Tau2HadTLV.SetPtEtaPhiE(jet_i.Pt(), jet_i.Eta(), jet_i.Phi(), jet_i.E());
	  fill_tau2 = true;
	  continue;
	}
	if((passed_jet_tau == true) && (fill_tau3 == false) && (fill_tau1 == true) && (fill_tau2 == true)){
	  Tau3HadTLV.SetPtEtaPhiE(jet_i.Pt(), jet_i.Eta(), jet_i.Phi(), jet_i.E());
	  fill_tau3 = true;
	  continue;
	}
	if((passed_jet_tau == true) && (fill_tau4 == false) && (fill_tau1 == true) && (fill_tau2 == true) && 
	   (fill_tau3 == true)){
	  Tau4HadTLV.SetPtEtaPhiE(jet_i.Pt(), jet_i.Eta(), jet_i.Phi(), jet_i.E());
	  fill_tau4 = true;
	}
      }
      
      ////////////////////////////////////
      for (int j = 0; j < branchJet->GetEntriesFast(); j++)
	{
	  //bool passed_jet_tau = false;
          
	  bool is_jet_elec_overlap = false;
	  Jet *jet = (Jet*) branchJet->At(j);
	  double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
	  jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
          
	  // Remove overlaps of all jets with electrons
	  for (int el = 0; el < branchElectron->GetEntriesFast(); el++){
	    Electron *elec = (Electron*) branchElectron->At(el);
	    double elec_energy = calculateE(elec->Eta, elec->PT, 0.000510998902);
	    elec_i.SetPtEtaPhiE(elec->PT, elec->Eta, elec->Phi, elec_energy);
	    double DR_jet_elec = jet_i.DeltaR(elec_i);
	    if ((DR_jet_elec < DR_jet_elec_max) && (elec->PT > 20.) && (abs(elec->Eta) < 2.5)){ 
	      is_jet_elec_overlap = true;
	      break;
	    }   
	  }         
	  if(!is_jet_elec_overlap){
	    if ((jet->PT > b_jet_pt_min) && (jet->BTag == 1)){is_b_jet = true;}
	    if (jet->TauTag == 0){
	      if(jet->PT > jet_min_pt){
		Jet_leading_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		jet_min_pt = jet->PT;
		njets_counter++;
	      }
	    }
	  }
	}
      
      //PRUEBA
      
      //for(int gen_i = 0; gen_i < branchGenParticle->GetEntriesFast(); gen_i++) {
      // GenParticle *gen = (GenParticle*) branchGenParticle->At(gen_i);
      // cout<<"N: "<<gen_i<<", St: "<<gen->Status<<", PID: "<<gen->PID<<", M1: "<<gen->M1<<", M2: "<<gen->M2<<", D1: "<<gen->D1<<", D2: "<<gen->D2<<endl; 
      //}
      //////// 
      
      // Loop over GenTau particles and emulate the identification and 
      // reconstruction of hadronic taus.
      for(int tau_i = 0; tau_i < branchGenParticle->GetEntriesFast(); tau_i++){
	// Look for the particle PDGID
	GenParticle *tau_mother = (GenParticle*) branchGenParticle->At(tau_i);
	if ((abs(tau_mother->PID) == 15) && (tau_mother->Status == 2)){
	  double tau_energy = calculateE(tau_mother->Eta, tau_mother->PT, tau_mother->Mass);
	  //Search for the daughters
	  //loop over the tau daughters 
	  int neutrinos = 0;
	  for(int ii=tau_mother->D1; ii<=tau_mother->D2; ii++) {
	    GenParticle *daughterCand = (GenParticle*) branchGenParticle->At(ii);
	    if((abs(daughterCand->PID) == 12) || (abs(daughterCand->PID) == 14) || (abs(daughterCand->PID) == 16)){
	      neutrinos++;
	      double cand_energy = calculateE(daughterCand->Eta, daughterCand->PT, daughterCand->Mass);
	      MCnu.SetPtEtaPhiE(daughterCand->PT, daughterCand->Eta, daughterCand->Phi, cand_energy);
	    }
	  }
	  if((neutrinos == 1) && (fill_tau1 == false) && (fill_tau2 == false) && (fill_tau3 == false)){
	    Tau1HadTLV.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, tau_energy);
	    Tau1HadTLV = Tau1HadTLV - MCnu;                       
	    fill_tau1 = true;
	    continue;
	  }
	  if((neutrinos == 1) && (fill_tau1 == true) && (fill_tau2 == false) && (fill_tau3 == false)){
	    Tau2HadTLV.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, tau_energy);
	    Tau2HadTLV = Tau2HadTLV - MCnu;
	    fill_tau2 = true;
	    continue;
	  }
	  if((neutrinos == 1) && (fill_tau1 == true) && (fill_tau2 == true) && (fill_tau3 == false)){
	    Tau3HadTLV.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, tau_energy);
	    Tau3HadTLV = Tau3HadTLV - MCnu;
	    fill_tau3 = true;
	    continue;
	  }
	  if((neutrinos == 1) && (fill_tau4 == false) && (fill_tau1 == true) && (fill_tau2 == true) && (fill_tau3 == true)){
	    Tau4HadTLV.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, tau_energy);
	    Tau4HadTLV = Tau4HadTLV - MCnu;
	  } 
	  ///////////////////////////////////////////////////////////////////////////////////////////////////////
	} // Close: if ((abs(tau_mother->PID) == 15) && (tau_mother->Status == 2))
      }// Close: for(int tau_i = 0; tau_i < branchGenParticle->GetEntriesFast(); tau_i++)
      //}// Close: if (branchJet->GetEntriesFast() > 0)
       
	 srand48(entry);
	 if((Tau1HadTLV.Pt() > 1.0) && (Tau2HadTLV.Pt() == 0.0)){
	 bool passed_tau1Cand = false; 
         passed_tau1Cand =  TauID(Tau1HadTLV);
         if (passed_tau1Cand) {
	 Tau1HadTLV = TauSmearing(Tau1HadTLV);
	  }
	 else {
	 Tau1HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
	  }
	 }
         if((Tau1HadTLV.Pt() > 0.0) && (Tau2HadTLV.Pt() > 0.0)){
         srand48(entry+numberOfEntries);
         bool passed_tau1Cand = false;
         passed_tau1Cand =  TauID2(Tau1HadTLV);
         if (passed_tau1Cand) {
         Tau1HadTLV = TauSmearing(Tau1HadTLV);
          }
         else {
         Tau1HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          } 
         srand48(entry+2*numberOfEntries);
    	 bool passed_tau2Cand = false;
	 passed_tau2Cand =  TauID2(Tau2HadTLV);
         if (passed_tau2Cand) { 
	 Tau2HadTLV = TauSmearing(Tau2HadTLV);
	  } 
	 else {
	 Tau2HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
	  }
         }
      /*
      bool passed_tau1Cand = false;
      if((Tau1HadTLV.Pt() > 1.0) && (Tau2HadTLV.Pt() == 0.0)){
	srand48(entry);
	passed_tau1Cand =  TauID(Tau1HadTLV);
	if (passed_tau1Cand) {
	  Tau1HadTLV = TauSmearing(Tau1HadTLV);
	}
            else {
              Tau1HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
            }
      }
      
      bool passed_tau2Cand = false;
      if((Tau1HadTLV.Pt() > 1.0) && (Tau2HadTLV.Pt() > 1.0)){
	srand48(entry);
	passed_tau1Cand =  TauID(Tau1HadTLV);
	if (passed_tau1Cand) {
	  Tau1HadTLV = TauSmearing(Tau1HadTLV);
	  srand48(entry+numberOfEntries);
	  passed_tau2Cand =  TauID(Tau2HadTLV);
	  if (passed_tau2Cand){
	    Tau2HadTLV = TauSmearing(Tau2HadTLV);
	  }else {
	    Tau2HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
	  }
	}
	else {
	  Tau1HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
	  srand48(entry+numberOfEntries);
	  passed_tau2Cand =  TauID(Tau2HadTLV);
	  if (passed_tau2Cand){
	    Tau2HadTLV = TauSmearing(Tau2HadTLV);
	  }else {
	    Tau2HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
	  }                
	}
        
      }
     */ 
      if((Tau1HadTLV.Pt()>2.0) ){ntau_counter++;}
      if((Tau2HadTLV.Pt()>2.0) ){ntau_counter++;}
      
      bool pass_lead_jet_cuts = false;
      // Events with no cuts
      
      //pass_cuts[1] = 1;
      
      if (Jet_leading_vec.Pt() > 2.0){
	pass_cuts[0] = 1;
      } 
      // for events with exactly 1 Tau
      if ((pass_cuts[0] == 1) && (ntau_counter == 1)){
	pass_cuts[1] = 1; 
      }
      ///Exactly 2 Taus
      if ((pass_cuts[0] == 1) && (ntau_counter == 2)){
	pass_cuts[2] = 1;
      }  
      // for events where the tau(s) passes a pT cut
      if ((pass_cuts[1] == 1) && (Tau1HadTLV.Pt() > tau1_pt_min)){
	pass_cuts[3] = 1;
      } 
      if ((pass_cuts[2] == 1) && (Tau1HadTLV.Pt() > tau1_pt_min) && (Tau2HadTLV.Pt() > tau2_pt_min)){
	pass_cuts[4] = 1;
      }
      // for events where the tau(s) passes an eta cut
      if((pass_cuts[3] == 1) && (abs(Tau1HadTLV.Eta()) < tau_eta_max)){
	pass_cuts[5] = 1;
      }
      if((pass_cuts[4] == 1) && (abs(Tau1HadTLV.Eta()) < tau_eta_max) && (abs(Tau2HadTLV.Eta()) < tau_eta_max)){
	pass_cuts[6] = 1;
      }
      // events passing also b-jet veto requirement
      if ( (pass_cuts[5]==1) && (!is_b_jet)){pass_cuts[7] = 1;}
      if ( (pass_cuts[6]==1) && (!is_b_jet)){pass_cuts[8] = 1;}
      // events passing jet pt cut
      if ((pass_cuts[7]==1) && (Jet_leading_vec.Pt() > lead_jet_pt)){pass_cuts[9] = 1;}
      if ((pass_cuts[8]==1) && (Jet_leading_vec.Pt() > lead_jet_pt)){pass_cuts[10] = 1;}
      // events passing jet eta cut
      if ((pass_cuts[9]==1) && (abs(Jet_leading_vec.Eta()) < jet_eta_max) ){pass_cuts[11] =1;}
      if ((pass_cuts[10]==1) && (abs(Jet_leading_vec.Eta()) < jet_eta_max) ){pass_cuts[12] =1;}
      // events passing dphitaujet cut 
      double tau1_jet_dphi = abs(normalizedDphi(Jet_leading_vec.Phi() - Tau1HadTLV.Phi()));
      if ((pass_cuts[11] == 1) && (tau1_jet_dphi > Dphitaujetmax)){pass_cuts[13] = 1;}
      if ((pass_cuts[12] == 1) && (tau1_jet_dphi > Dphitaujetmax)){pass_cuts[14] = 1;}
      // events passing mt cut
      double transmass = TMath::Sqrt(TMath::Abs(2*Tau1HadTLV.Pt()*MET*(1-TMath::Cos(TMath::Abs(Tau1HadTLV.Phi() - MET_phi)))));
      if ((pass_cuts[11] == 1) && (transmass > transversemassmin)){pass_cuts[15] = 1;}
      if ((pass_cuts[12] == 1) && (transmass > transversemassmin)){pass_cuts[16] = 1;}
      // events passing also MET cut
      if ((pass_cuts[11] == 1) && (MET > met_min)){pass_cuts[17] = 1;}
      if ((pass_cuts[12] == 1) && (MET > met_min)){pass_cuts[18] = 1;}
      // events passing dphijetmet cut
      double jet_met_dphi = abs(normalizedDphi(Jet_leading_vec.Phi() - MET_phi));
      if ((pass_cuts[11] == 1) && (jet_met_dphi > deltajetmetphi)){pass_cuts[19]=1;}
      if ((pass_cuts[12] == 1) && (jet_met_dphi > deltajetmetphi)){pass_cuts[20]=1;}
      // events passing effective mass cut
      double efective_mass = TMath::Sqrt(Jet_leading_vec.Pt()*Jet_leading_vec.Pt()+Tau1HadTLV.Pt()*Tau1HadTLV.Pt()+MET*MET);
      if ((pass_cuts[11] == 1) && (efective_mass > efectivemassmin)){pass_cuts[21] = 1;}        
      if ((pass_cuts[12] == 1) && (efective_mass > efectivemassmin)){pass_cuts[22] = 1;}
      
      // events passing previous cuts with exactly 1 tau
      
      // save some important histograms
      // _hmap_lead_jet_pT[0]->Fill(Jet_leading_vec.Pt());  
      
      for (int i = 0; i < nDir; i++){
	_hmap_Nevents[i]->Fill(0.0);
	_hmap_n_jets[i]->Fill(njets_counter);
	_hmap_n_tau[i]->Fill(ntau_counter);
	
	if ( pass_cuts[i] == 1){
	  double jet_met_dphi = abs(normalizedDphi(Jet_leading_vec.Phi() - MET_phi));
	  _hmap_jet_met_Dphi[i]->Fill(abs(jet_met_dphi));
	  _hmap_jet_met_metDphi[i]->Fill(abs(jet_met_dphi),MET);
	  _hmap_Nevents[i]->Fill(1.0);
	  _hmap_met[i]->Fill(MET);
	  _hmap_lead_jet_pT[i]->Fill(Jet_leading_vec.Pt());
	  _hmap_lead_jet_eta[i]->Fill(Jet_leading_vec.Eta());
	  _hmap_lead_jet_phi[i]->Fill(Jet_leading_vec.Phi());
	  _hmap_RM[i]->Fill(MET/Jet_leading_vec.Pt());
	  if(Tau1HadTLV.Pt() > 1.){ 
	    _hmap_tau1_pT[i]->Fill(Tau1HadTLV.Pt());
	    _hmap_tau1_eta[i]->Fill(Tau1HadTLV.Eta());
	    _hmap_tau1_phi[i]->Fill(Tau1HadTLV.Phi());
	    double tau1_met_dphi = abs(normalizedDphi(Tau1HadTLV.Phi() - MET_phi));
	    _hmap_tau1_met_Dphi[i]->Fill(tau1_met_dphi);
	    double transmass = TMath::Sqrt(TMath::Abs(2*Tau1HadTLV.Pt()*MET*(1-TMath::Cos(TMath::Abs(Tau1HadTLV.Phi() - MET_phi)))));
	    _hmap_transverse_mass[i]->Fill(transmass);
            if(Jet_leading_vec.DeltaR(Tau1HadTLV) > .1){
	      double tau1_jet_Dphi = TMath::Abs(Jet_leading_vec.Phi() - Tau1HadTLV.Phi());
	      _hmap_tau1_jet_Dphi[i]->Fill(tau1_jet_Dphi); 
	      double efective_mass = TMath::Sqrt(Jet_leading_vec.Pt()*Jet_leading_vec.Pt()+Tau1HadTLV.Pt()*Tau1HadTLV.Pt()+MET*MET);
	      _hmap_efective_mass[i]->Fill(efective_mass);       
	    }
	  }
	  if(Tau2HadTLV.Pt() > 1.0){
	    _hmap_tau2_pT[i]->Fill(Tau2HadTLV.Pt());
	    _hmap_tau2_eta[i]->Fill(Tau2HadTLV.Eta());
	    _hmap_tau2_phi[i]->Fill(Tau2HadTLV.Phi());
	    double tau2_met_dphi = abs(normalizedDphi(Tau2HadTLV.Phi() - MET_phi));
	    _hmap_tau2_met_Dphi[i]->Fill(tau2_met_dphi);
	  }
	}
      } 
    }//for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  
  theFile->cd();
  Numbertaus_D1->Write();  
  Numbertaus_D2->Write();
  Numbertaus_2->Write();
  Mothers->Write();
  for (int d = 0; d < nDir; d++)
    {
      cdDir[d]->cd();
      
      _hmap_Nevents[d]->Write();
      _hmap_n_jets[d]->Write();
      _hmap_n_tau[d]->Write();
      _hmap_lead_jet_pT[d]->Write();
      _hmap_lead_jet_eta[d]->Write();
      _hmap_lead_jet_phi[d]->Write();
      _hmap_RM[d]->Write();
      _hmap_tau1_pT[d]->Write();
      _hmap_tau1_eta[d]->Write();
      _hmap_tau1_phi[d]->Write();
      _hmap_tau2_pT[d]->Write();
      _hmap_tau2_eta[d]->Write();
      _hmap_tau2_phi[d]->Write();
      _hmap_jet_met_Dphi[d]->Write();
      _hmap_met[d]->Write();
      _hmap_tau1_met_Dphi[d]->Write();
      _hmap_tau2_met_Dphi[d]->Write();
      _hmap_jet_met_metDphi[d]->Write();
      _hmap_tau1_jet_Dphi[d]->Write();
      _hmap_transverse_mass[d]->Write();
      _hmap_efective_mass[d]->Write();
      
    }
  
  //   cdDir[0]->cd();
  theFile->Close();
  
}



PhenoAnalysis::~PhenoAnalysis()
{
  // do anything here that needs to be done at desctruction time
}

bool PhenoAnalysis::TauID(TLorentzVector TauCand) {
  
  bool passedTauID = false;
  if (TauCand.Pt() > 0.){
    double r = drand48();
    if (r <= 0.8){passedTauID = true;}
  }
  return passedTauID;
}

bool PhenoAnalysis::TauID2(TLorentzVector TauCand) {
  
  bool passedTauID = false;
  if (TauCand.Pt() > 0.){
    double r = drand48();
    if (r <= 0.64){passedTauID = true;}
  }
  return passedTauID;
}

bool PhenoAnalysis::TauIDJet(TLorentzVector jet){
  bool passedTauIDJet = false;
  if(jet.Pt() > 0.){
    double x = drand48();
    if(x < 0.03){passedTauIDJet = true;}
  }
  return passedTauIDJet;
}

TLorentzVector PhenoAnalysis::TauSmearing(TLorentzVector TauCand) {
  
  if (TauCand.Pt() > 5.){
    double smeared_pt  = 0.03*TauCand.Pt() + TauCand.Pt();
    double smeared_eta = TauCand.Eta();
    double smeared_phi = TauCand.Phi();
    double smeared_e   = TauCand.E(); 
    TauCand.SetPtEtaPhiE(smeared_pt, smeared_eta, smeared_phi, smeared_e);
  }
  return TauCand;
}


double PhenoAnalysis::calculateE(double eta, double pt, double mass){
  
  double theta = TMath::ATan(TMath::Exp(-eta)); 
  double cos_theta = TMath::Cos(2*theta);
  double p= pt/cos_theta;
  double e = sqrt(pow(p, 2) + pow(mass, 2));
  
  return e;
  
}

double PhenoAnalysis::normalizedDphi(double phi){
  const double PI  = 3.141592653589793238463;
  double twoPI = 2.0*PI;
  if ( phi < -PI ){phi += twoPI;}
  if ( phi > PI ){phi -= twoPI;}
  return phi;
}
void PhenoAnalysis::crateHistoMasps (int directories)
{
  for (int i = 0; i < directories; i++)
    {
      _hmap_Nevents[i]       = new TH1F("Nevents", "Nevents", 3,0.,3); 
      _hmap_lead_jet_pT[i]   = new TH1F("jet_lead_pT",    "Pt leading jet", 100, 0., 1000.);
      _hmap_lead_jet_eta[i]  = new TH1F("jet_lead_eta",   "#eta jet", 50, -5.0, 5.0);
      _hmap_lead_jet_phi[i]  = new TH1F("jet_lead_phi",   "#phi jet", 70, -3.6, 3.6);
      _hmap_RM[i]            = new TH1F("RM",   "rm", 100, 0, 3); 
      _hmap_n_jets[i]        = new TH1F("N_jets",         "N(jet)", 4, 0., 4);
      _hmap_n_tau[i]         = new TH1F("N_taus",          "", 4, 0., 4);
      _hmap_tau1_pT[i]       = new TH1F("tau1_pT",        "p_{T}(#tau_{1})", 200, 0., 600.);
      _hmap_tau1_eta[i]      = new TH1F("tau1_eta",       "#eta(#tau_{1})", 50, -3.5, 3.5);
      _hmap_tau1_phi[i]      = new TH1F("tau1_phi",       "#phi(#tau_{1})", 70, -3.6, 3.6);
      _hmap_tau2_pT[i]       = new TH1F("tau2_pT",        "p_{T}(#tau_{2})", 150, 0., 300.);
      _hmap_tau2_eta[i]      = new TH1F("tau2_eta",       "#eta(#tau_{2})", 50, -3.5, 3.5);
      _hmap_tau2_phi[i]      = new TH1F("tau2_phi",       "#phi(#tau_{2})", 70, -3.6, 3.6);
      _hmap_jet_met_Dphi[i]  = new TH1F("jet_met_Dphi",   "#Delta #phi(jet, MET)", 32, 0, 3.2);
      _hmap_transverse_mass[i] = new TH1F("transverse_mass",   "Transverse Mass", 60, 0, 600);
      _hmap_met[i]           = new TH1F("met",            "Met", 120, 0., 1200.); 
      _hmap_tau1_met_Dphi[i]  = new TH1F("tau_met_Dphi",   "#Delta #phi(#tau_{1}, MET)", 32, 0, 3.2);
      _hmap_tau2_met_Dphi[i]  = new TH1F("tau_met_Dphi",   "#Delta #phi(#tau_{2}, MET)", 32, 0, 3.2);
      _hmap_tau1_jet_Dphi[i]  = new TH1F("tau1_jet_Dphi",   "#Delta #phi(#tau_{1}, jet)", 32, 0, 3.2);
      _hmap_efective_mass[i]  = new TH1F("efective_mass", "Efective_mass", 80, 0, 800); 
      _hmap_jet_met_metDphi[i] = new TH2F("jet_met_metDphi", "#Delta #phi(jet,MET), MET", 32, 0, 3.2, 120, 0, 1200);
    }
}
