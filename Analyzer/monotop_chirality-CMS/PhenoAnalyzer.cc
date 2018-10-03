////////////////////////////////////////////////////////////////
//                                                            //
// Author: Andrés Flórez, Universidad de los Andes, Colombia  //
//                                                            //
////////////////////////////////////////////////////////////////


#include <iostream>
#include "ROOTFunctions.h"
#include "PhenoAnalyzer.h"
#include "DelphesFunctions.h"

int main(int argc, char *argv[]) {
  
  TChain chain("Delphes");
  chain.Add(argv[1]);
  TFile * HistoOutputFile = new TFile(argv[2], "RECREATE");
  int nDir = 9;
  TDirectory *theDirectory[nDir+1];
  theDirectory[0]  = HistoOutputFile->mkdir("No_cuts");
  theDirectory[1]  = HistoOutputFile->mkdir("After_jet1_pt");
  theDirectory[2]  = HistoOutputFile->mkdir("After_bjet_pt");
  theDirectory[3]  = HistoOutputFile->mkdir("After_bjet_eta");
  theDirectory[4]  = HistoOutputFile->mkdir("After_electron");
  theDirectory[5]  = HistoOutputFile->mkdir("After_muon");
  theDirectory[6]  = HistoOutputFile->mkdir("After_bpluslep");
  theDirectory[7]  = HistoOutputFile->mkdir("After_met_only");
  theDirectory[8]  = HistoOutputFile->mkdir("After_mt_only");
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
  
  double DR_jet_elec_max  = params->GetValue ("DR_jet_elec_max", 0.);
  double DR_jet_muon_max  = params->GetValue ("DR_jet_muon_max", 0.); 
  double DR_elec_muon_max = params->GetValue ("DR_elec_muon_max", 0.);
  double jet_pt_min       = params->GetValue ("jet_pt_min", 0.);
  double jet_eta_max      = params->GetValue ("jet_eta_max", 0.);
  double b_pt_min         = params->GetValue ("b_pt_min", 0.);
  double b_eta_max        = params->GetValue ("b_eta_max", 0.);
  double elec_pt_min      = params->GetValue ("elec_pt_min", 0.);
  double muon_pt_min      = params->GetValue ("muon_pt_min", 0.);
  double met_min          = params->GetValue ("met_min", 0.);
  double mt_min           = params->GetValue ("mt_min", 0.);
  double blpt_min         = params->GetValue ("blpt_min", 0.);
 
  crateHistoMasps(nDir);
  
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET"); 

  MissingET *METpointer; 
  Double_t electron_mass = 0.000510998902;
  Double_t muon_mass = 0.105658369;  
  
  //Loop over entries  
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    {
      treeReader->ReadEntry(entry);
      //Missing ET
      METpointer = (MissingET*) branchMissingET->At(0);
      double MET = METpointer->MET;
      double MET_phi = METpointer->Phi;


      TLorentzVector jet_i(0., 0., 0., 0.);
      TLorentzVector elec_1(0., 0., 0., 0.);
      TLorentzVector muon_1(0., 0., 0., 0.);
      TLorentzVector bjet_i(0.,0.,0.,0.);
      TLorentzVector Tau_1(0.,0.,0.,0.);
      TLorentzVector Tau_2(0.,0.,0.,0.);

      TLorentzVector jet_1(0., 0., 0., 0.);
      TLorentzVector jet_2(0., 0., 0., 0.);
       
      int pass_cuts[nDir] = {0};
      
      vector<int> tau_index; 
      //Particle counter   0=njets, 1=nelectrons, 2=nmuons, 3=bjets, 4=taus
      int particle_counter[5] = {0};
      //Particle identificators
      bool is_jet_elec_overlap = false;
      bool is_jet_muon_overlap = false;
      bool is_b_jet = false;
      bool filled_tau1 = false;
      bool filled_tau2 = false;
      bool filled_first_jet = false;
      bool passed_mass_w = false;
      //Soft initial cuts
      double jet_ipt = 10.0;
      double jet_ieta = 5.0;
      double elec_1pt = 10.0;
      double elec_1eta = 2.5;
      double muon_1pt = 10.0;
      double muon_1eta = 2.5;
      double bjet_pt = 10.0;
      double bjet_eta = 5.0;
      double tau_pt = 10.0;
      double tau_eta = 2.5;
      double jet1_pt = 10.0;
      double jet1_eta = 5.0;
      double jet2_pt = 10.0;
      double jet2_eta = 5.0;
      // Transverse mass
      double MT = 0.0;
      //Loop over jets    
      for (int j = 0; j < branchJet->GetEntriesFast(); j++)
	{ 
          is_jet_elec_overlap = false;
	  is_jet_muon_overlap = false;
 	  is_b_jet = false; 
 
          // Find any jet passing pT and eTa cut	  
          Jet *jet = (Jet*) branchJet->At(j);
	  double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
	  if((jet->PT > jet_ipt) && (TMath::Abs(jet->Eta) < jet_ieta)){
            particle_counter[0]++;
	    jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
            
	    // Find if there are electrons overlapping with the jets
	    for (int e = 0; e < branchElectron->GetEntriesFast(); e++){
	      Electron *electron = (Electron*) branchElectron->At(e);
	      double electron_energy = calculateE(electron->Eta, electron->PT, electron_mass);
	      elec_1.SetPtEtaPhiE(electron->PT, electron->Eta, electron->Phi, electron_energy);
	      double DR_jet_electron = jet_i.DeltaR(elec_1);
	      if ((electron->PT > elec_1pt) && (TMath::Abs(electron->Eta) < elec_1eta)){
		if ( DR_jet_electron < DR_jet_elec_max ){
		  is_jet_elec_overlap = true;
		  elec_1.SetPtEtaPhiE(0.0, 0.0, 0.0, 0.0);
		  break;} 
		  else{
		  particle_counter[1]++;
		  elec_1.SetPtEtaPhiE(electron->PT, electron->Eta, electron->Phi, electron_energy);
		  }
	      }
	    }

	    // Find if there are muons overlapping with the jets
	    for (int m = 0; m < branchMuon->GetEntriesFast(); m++){
	      Muon *muon = (Muon*) branchMuon->At(m);
	      double muon_energy = calculateE(muon->Eta, muon->PT, muon_mass);
	      muon_1.SetPtEtaPhiE(muon->PT, muon->Eta, muon->Phi, muon_energy);
	      double DR_jet_muon = jet_i.DeltaR(muon_1);
	      if ((muon->PT > muon_1pt) && (TMath::Abs(muon->Eta) < muon_1eta)){ 
		if (DR_jet_muon < DR_jet_muon_max){ 
		  is_jet_muon_overlap = true;
		  muon_1.SetPtEtaPhiE(0.0, 0.0, 0.0, 0.0);
		  break;} 
		  else{
		  particle_counter[2]++;
		  muon_1.SetPtEtaPhiE(muon->PT, muon->Eta, muon->Phi, muon_energy);
		  }
	      }
	    }
          } //soft jets cuts

	  // If there no jets overlapping with electrons or muons, then fill the info for the b-jets....
	  if((!is_jet_elec_overlap) && (!is_jet_muon_overlap)){
          //--------------------------------------------------------------------
	    if ((jet->PT > bjet_pt) && (TMath::Abs(jet->Eta) < bjet_eta) && (jet->BTag == 1)){
	      bjet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);           
	      particle_counter[3]++;
              is_b_jet = true;
	    }
          //--------------------------------------------------------------------
            if ((jet->TauTag == 1) && (jet->BTag == 0)){
              particle_counter[4]++;
              if ((jet->PT > tau_pt) && (TMath::Abs(jet->Eta) < tau_eta))
                if (!filled_tau1){
                  Tau_1.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
                  filled_tau1 = true;
                  continue;                  
                }
                if (filled_tau1 && !filled_tau2){
                  Tau_2.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
                  filled_tau2 = true;
                }
	    }
	   //-------------------------------------------------------------------- 
	    if ((jet->TauTag == 0) && ((jet->BTag == 0))){
            if ((jet->PT > jet1_pt) && (TMath::Abs(jet->Eta) < jet1_eta))
              double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
	      if(filled_first_jet == false){
		jet_1.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		jet1_pt = jet->PT;
                filled_first_jet = true;
                continue;
	      }
              if (filled_first_jet == true){
              if ((jet->PT > jet2_pt) && (TMath::Abs(jet->Eta) < jet2_eta))
                jet_2.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
                float mass_w = (jet_1 + jet_2).M();
                if ((mass_w < 100.0) && (mass_w > 60.0)){passed_mass_w = true;} 
              }              
	    }
	   //--------------------------------------------------------------------- 
	  } //is not electron or muon overlap
	  
	}    //fisnish loop over branchjets 
    
      //Other parameters
      double bpt=0.;
      double bpluslep=0.;
      double jet1_met_dphi=0.;
      double jet2_met_dphi=0.;
      double dphi_met_lep=0.;
      jet1_met_dphi=normalizedDphi(jet_1.Phi()-MET_phi);
      jet2_met_dphi=normalizedDphi(jet_2.Phi()-MET_phi);
      //Calculating transverse mass
      dphi_met_lep=normalizedDphi(elec_1.Phi()-MET_phi);
      MT = transverseMass(elec_1.Pt(), MET, dphi_met_lep);    
      
      //cuts	
      pass_cuts[0]=1;
      if(jet_1.Pt() > jet_pt_min){pass_cuts[1] = 1;}
      if((bjet_i.Pt() > b_pt_min)&&(pass_cuts[1] == 1)){pass_cuts[2] = 1;}
      if((TMath::Abs(bjet_i.Eta()) < jet_eta_max) && (pass_cuts[2] == 1)){pass_cuts[3] = 1;}
      if((elec_1.Pt() > elec_pt_min) && (particle_counter[2]==0) && (pass_cuts[3] == 1)){
      bpluslep=bjet_i.Pt()+elec_1.Pt();
      bpt = bjet_i.Pt()/bpluslep;
      pass_cuts[4] = 1;}  
      if((muon_1.Pt() > muon_pt_min) && (particle_counter[1]==0) && (pass_cuts[3] == 1)){
      bpluslep=bjet_i.Pt()+muon_1.Pt();
      bpt = bjet_i.Pt()/bpluslep;
      pass_cuts[5] = 1;} 
      if ((pass_cuts[4] == 1) && (bpluslep > blpt_min)){pass_cuts[6] =1;} 
      if ((MET > met_min) && (pass_cuts[5] == 1)){pass_cuts[7] = 1;}
      if ((MT > mt_min) &&  (pass_cuts[4] == 1) ){pass_cuts[8] = 1;} 
      


      // save some important histograms
      for (int i = 0; i < nDir; i++){
        if (pass_cuts[i] == 1){
          _hmap_Nevents[i]->Fill(0.0);  
	  _hmap_N_jets[i]->Fill(particle_counter[0]);
          _hmap_N_elec[i]->Fill(particle_counter[1]);
          _hmap_N_muon[i]->Fill(particle_counter[2]);
	  _hmap_N_bjet[i]->Fill(particle_counter[3]);
          _hmap_N_tau[i]->Fill(particle_counter[4]);
          

	  if(jet_1.Pt()>0.){ 
	    _hmap_jet1_pT[i]->Fill(jet_1.Pt());
	    _hmap_jet1_eta[i]->Fill(jet_1.Eta());
	    _hmap_jet1_phi[i]->Fill(jet_1.Phi());
            _hmap_jet_1_met_dphi[i]->Fill(jet1_met_dphi);
            _hmap_met_jet1_met_dphi[i]->Fill(jet1_met_dphi,MET);
	  }
	  
          if(jet_2.Pt()>0.){
            _hmap_jet2_pT[i]->Fill(jet_2.Pt());
            _hmap_jet2_eta[i]->Fill(jet_2.Eta());
            _hmap_jet2_phi[i]->Fill(jet_2.Phi());
            _hmap_jet_2_met_dphi[i]->Fill(jet2_met_dphi);
          }
	  
	  if(Tau_1.Pt()>0.){
            _hmap_tau1_pT[i]->Fill(Tau_1.Pt());
            _hmap_tau1_eta[i]->Fill(Tau_1.Eta());
            _hmap_tau1_phi[i]->Fill(Tau_1.Phi());
	  }
	  if(Tau_2.Pt()>0.){
            _hmap_tau2_pT[i]->Fill(Tau_2.Pt());
            _hmap_tau2_eta[i]->Fill(Tau_2.Eta());
            _hmap_tau2_phi[i]->Fill(Tau_2.Phi());
          }
	  if(bjet_i.Pt()>0.){
             if ((elec_1.Pt() > 0.0) || (muon_1.Pt() > 0.0))
             {_hmap_bjet_ratio_pT[i]->Fill(bpt);
            _hmap_bjet_plus_lep[i]->Fill(bpluslep);}	
	    _hmap_bjet_pT[i]->Fill(bjet_i.Pt());
	    _hmap_bjet_eta[i]->Fill(bjet_i.Eta());
	    _hmap_bjet_phi[i]->Fill(bjet_i.Phi());
	    _hmap_bjet_energy[i]->Fill(bjet_i.E()*0.002);
	  }
          if ((elec_1.Pt() > 0.0) || (muon_1.Pt() > 0.0))
          { _hmap_mt[i]->Fill(MT);
          } 
	  if (elec_1.Pt() > 0.0) {
            _hmap_met_lep[i]->Fill(elec_1.Pt()+MET); 
	    _hmap_elec_pT[i]->Fill(elec_1.Pt());
	    _hmap_elec_eta[i]->Fill(elec_1.Eta());
	    _hmap_elec_phi[i]->Fill(elec_1.Phi());
	  }
	  if (muon_1.Pt() > 0.0) {
	    _hmap_muon_pT[i]->Fill(muon_1.Pt());
	    _hmap_muon_eta[i]->Fill(muon_1.Eta());
	    _hmap_muon_phi[i]->Fill(muon_1.Phi());
	  }
	  _hmap_met[i]->Fill(MET);
	} 
      }
      
     }//Finish loop over entries
  
  
  theFile->cd();
  for (int d = 0; d < nDir; d++)
    {
      cdDir[d]->cd();
      _hmap_Nevents[d]->Write();
      _hmap_N_jets[d]->Write();
      _hmap_jet1_pT[d]->Write();
      _hmap_jet1_eta[d]->Write();
      _hmap_jet1_phi[d]->Write();
      _hmap_jet2_pT[d]->Write();
      _hmap_jet2_eta[d]->Write();
      _hmap_jet2_phi[d]->Write();
      _hmap_N_tau[d]->Write();
      _hmap_tau1_pT[d]->Write();
      _hmap_tau1_eta[d]->Write();
      _hmap_tau1_phi[d]->Write();
      _hmap_tau2_pT[d]->Write();
      _hmap_tau2_eta[d]->Write();
      _hmap_tau2_phi[d]->Write();
      _hmap_N_bjet[d]->Write();
      _hmap_bjet_plus_lep[d]->Write();
      _hmap_bjet_ratio_pT[d]->Write();
      _hmap_bjet_pT[d]->Write();
      _hmap_bjet_eta[d]->Write();
      _hmap_bjet_phi[d]->Write();
      _hmap_bjet_energy[d]->Write();
      _hmap_N_elec[d]->Write();
      _hmap_elec_pT[d]->Write();
      _hmap_elec_eta[d]->Write();
      _hmap_elec_phi[d]->Write();
      _hmap_N_muon[d]->Write();
      _hmap_muon_pT[d]->Write();
      _hmap_muon_eta[d]->Write();
      _hmap_muon_phi[d]->Write();
      _hmap_met[d]->Write();
      _hmap_met_lep[d]->Write();
      _hmap_mt[d]->Write();
      _hmap_jet_1_met_dphi[d]->Write();
      _hmap_jet_2_met_dphi[d]->Write();
      _hmap_met_jet1_met_dphi[d]->Write();    
    }
  
  theFile->Close();
  
} //finish class

PhenoAnalysis::~PhenoAnalysis()
{
  // do anything here that needs to be done at desctruction time
}

double PhenoAnalysis::calculateE(double eta, double pt, double mass){
  
  double theta = TMath::ATan(TMath::Exp(-eta)); 
  double cos_theta = TMath::Cos(2*theta);
  double p= pt/cos_theta;
  double e = sqrt(pow(p, 2) + pow(mass, 2));
  
  return e;
  
}


double PhenoAnalysis::normalizedDphi(double phi){
  double twoPI = 2.0*TMath::Pi();
  if ( phi < -TMath::Pi() ){phi += twoPI;}
  if ( phi > TMath::Pi() ){phi = twoPI-phi;}
  else phi = TMath::Abs(phi); 
  return phi;
}

double PhenoAnalysis::transverseMass(double ptl, double met, double deltaphi){

  double Mt = TMath::Sqrt(2*ptl*met*(1-cos(deltaphi))); 
  return Mt;
}

void PhenoAnalysis::crateHistoMasps (int directories)
{
  for (int i = 0; i < directories; i++)
    {
      _hmap_Nevents[i]		= new TH1F("Nevents", "Nevents", 5,-0.5,4.5); 
      
      _hmap_N_jets[i]		= new TH1F("N_jets",  "N_jets", 11, -0.5, 10.5);    

      _hmap_jet1_pT[i]		= new TH1F("jet1_pT",  "P_{T} jet1 ", 100, 0., 1000.);
      _hmap_jet1_eta[i]		= new TH1F("jet1_eta", "#eta jet1 ", 50, -5.0, 5.0);
      _hmap_jet1_phi[i]		= new TH1F("jet1_phi", "#phi jet1 ", 70, -3.6, 3.6);

      _hmap_jet2_pT[i]          = new TH1F("jet2_pT",  "P_{T} jet2 ", 100, 0., 1000.);
      _hmap_jet2_eta[i]         = new TH1F("jet2_eta", "#eta jet2 ", 50, -5.0, 5.0);
      _hmap_jet2_phi[i]         = new TH1F("jet2_phi", "#phi jet2 ", 70, -3.6, 3.6);
      
      _hmap_N_tau[i]		= new TH1F("N_tau",  "N_tau", 5, -0.5, 4.5);

      _hmap_tau1_pT[i]		= new TH1F("tau1_pT",  "p_{T}(#tau_{1})", 100, 0., 1000.);
      _hmap_tau1_eta[i]		= new TH1F("tau1_eta", "#eta(#tau_{1})", 50, -5.0, 5.0);
      _hmap_tau1_phi[i]		= new TH1F("tau1_phi", "#phi(#tau_{1})", 72, -3.6, 3.6);
      
      _hmap_tau2_pT[i]          = new TH1F("tau2_pT",  "p_{T}(#tau_{2})", 100, 0., 1000.);
      _hmap_tau2_eta[i]         = new TH1F("tau2_eta", "#eta(#tau_{2})", 50, -5.0, 5.0);
      _hmap_tau2_phi[i]         = new TH1F("tau2_phi", "#phi(#tau_{2})", 72, -3.6, 3.6);

      _hmap_bjet_plus_lep[i]    = new TH1F("bjet_plus_lep", "p_{T} bjet+lep", 50, 0., 500);
      
      _hmap_N_bjet[i]		= new TH1F("N_bjet",    "N_bjet", 5, -0.5, 4.5);
      _hmap_bjet_ratio_pT[i]	= new TH1F("bjet_ratio_pT", "p_{T} bjet", 50, 0, 1);
      _hmap_bjet_pT[i]		= new TH1F("bjet_pT",  "P_{T} bjet ", 100, 0., 1000.);
      _hmap_bjet_eta[i]		= new TH1F("bjet_eta", "#eta bjet ", 50, -5.0, 5.0);
      _hmap_bjet_phi[i]		= new TH1F("bjet_phi", "#phi bjet ", 70, -3.6, 3.6);
      _hmap_bjet_energy[i]	= new TH1F("bjet_energy", "Energy b quark", 50, 0, 2);
      
      _hmap_N_elec[i]           = new TH1F("N_elec",    "N_elec", 11, -0.5, 10.5);
      _hmap_elec_pT[i]          = new TH1F("elec_pT", "p_{T} e ", 100, 0., 1000.);
      _hmap_elec_eta[i]         = new TH1F("elec_eta","#eta e ", 50, -5.0, 5.0);
      _hmap_elec_phi[i]         = new TH1F("elec_phi", "#phi e", 70, -3.6, 3.6);      
      
      _hmap_N_muon[i]           = new TH1F("N_muon",    "N_muon", 11, -0.5, 10.5);
      _hmap_muon_pT[i]          = new TH1F("muon_pT", "p_{T} #mu ", 100, 0., 1000.);
      _hmap_muon_eta[i]         = new TH1F("muon_eta","#eta #mu ", 50, -5.0, 5.0);
      _hmap_muon_phi[i]         = new TH1F("muon_phi", "#phi #mu", 70, -3.6, 3.6); 
      
      _hmap_met[i]		= new TH1F("met", "MET", 200, 0. , 2000.);
      _hmap_met_lep[i]          = new TH1F("met_lep", "P_{T}(l)+Met", 200, 0. , 2000.);
      _hmap_mt[i]               = new TH1F("mt", "Transverse mass", 200, 0. , 2000.);	
      _hmap_jet_1_met_dphi[i]   = new TH1F("jet_1_met_dphi", "#Delta #phi(jet1, MET)", 32, 0, 3.2); 
      _hmap_jet_2_met_dphi[i]   = new TH1F("jet_2_met_dphi", "#Delta #phi(jet2, MET)", 32, 0, 3.2);     
      _hmap_met_jet1_met_dphi[i]= new TH2F("met_jet1_met_dphi", "#Delta #phi(jet1,MET), MET", 32, 0., TMath::Pi(), 50, 0., 500.);
    }
}

