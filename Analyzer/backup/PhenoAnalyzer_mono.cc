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
  int nDir = 6;
  TDirectory *theDirectory[nDir+1];
  theDirectory[0]  = HistoOutputFile->mkdir("No_cuts");
  theDirectory[1]  = HistoOutputFile->mkdir("After_jetpt");
  theDirectory[2]  = HistoOutputFile->mkdir("After_jeteta");
  theDirectory[3]  = HistoOutputFile->mkdir("After_bpt");
  theDirectory[4]  = HistoOutputFile->mkdir("After_beta");
  theDirectory[5]  = HistoOutputFile->mkdir("After_lpt");
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
  double jet_pt_min       = params->GetValue ("jet_pt_min", 0.);
  double jet_eta_max      = params->GetValue ("jet_eta_max", 0.);
  double b_pt_min         = params->GetValue ("b_pt_min", 30.0);
  double b_eta_max        = params->GetValue ("b_eta_max", 0.);
  double elec_pt_min      = params->GetValue ("elec_pt_min", 0.);
  double muon_pt_min      = params->GetValue ("muon_pt_min", 0.);
  double met_min          = params->GetValue ("met_min", 0.);
  crateHistoMasps(nDir);
  
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET"); 
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle"); 
  MissingET *METpointer; 
  
  
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    {
      TLorentzVector muon_i(0., 0., 0., 0.);
      TLorentzVector muon_igen(0., 0., 0., 0.); 
      TLorentzVector jet_i(0., 0., 0., 0.);
      TLorentzVector elec_i(0., 0., 0., 0.);
      TLorentzVector elec_igen(0., 0., 0., 0.);
      TLorentzVector Tau1_vec(0.,0.,0.,0.);
      TLorentzVector b_vec(0.,0.,0.,0.);
      TLorentzVector b_vecgen(0.,0.,0.,0.); 
      int pass_cuts[nDir] = {0};
      
      vector<int> tau_index;
      vector<int> b_index;
      bool is_b_jet = false;
      bool is_jet_elec_overlap = false;
      bool is_jet_muon_overlap = false;
      treeReader->ReadEntry(entry);
      int njets_counter = 0;
      int tau_counter=0;
      int b_counter=0;
      int ntau_counter = 0;
      int nb_counter=0;
      int muon_counter=0;
      int elec_counter=0;
      double jetphi =0; 
      METpointer = (MissingET*) branchMissingET->At(0);
      double MET = METpointer->MET;
      double MET_phi = METpointer->Phi;
      double bpt=0.; 
      double bptgen=0.;

      
      double jet_met_dphi=0.;
      double tau_met_dphi=0.;
      double jet_tau_dphi=0.;
      
      
      pass_cuts[0]=1; 
      if(branchJet->GetEntries() > 0)
	{
     
	  for (int j = 0; j < branchJet->GetEntriesFast(); j++)
	    { 
	      
	      Jet *jet = (Jet*) branchJet->At(j);
	      if((jet->BTag == 0)&&(jet->TauTag == 0))
 		{
		  double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
		  jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		  njets_counter++;
                }
	      
	      if((jet->BTag == 1) &&(jet->TauTag == 0)){
		is_b_jet = true;
		b_index.push_back(j);
		nb_counter++;
	      }
	      if(b_index.size() > 0){
		
		for(int i = 0; i < b_index.size(); i++)
		  {
		    int b_i = b_index.at(i);
		    Jet *jet = (Jet*) branchJet->At(b_i);
		    double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
		    b_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		    b_counter++;    
		  }
	      }
	      //overlap muons
	      for (int k= 0; k < branchMuon->GetEntriesFast(); k++){
		Muon *m = (Muon*) branchMuon->At(k);
		double muon_energy = calculateE(m->Eta, m->PT, 0.105658369);
		muon_i.SetPtEtaPhiE(m->PT, m->Eta, m->Phi, muon_energy);             
		double DR_jet_elec = b_vec.DeltaR(elec_i);
		if (DR_jet_elec < DR_jet_elec_max){
		  is_jet_muon_overlap = true;
		  break;}
		muon_counter++;			 
	      } 
	      
	      //overlap electrons
	      for (int el = 0; el < branchElectron->GetEntriesFast(); el++){
		Electron *elec = (Electron*) branchElectron->At(el);
		double elec_energy = calculateE(elec->Eta, elec->PT, 0.000510998902);
		elec_i.SetPtEtaPhiE(elec->PT, elec->Eta, elec->Phi, elec_energy);
		double DR_jet_elec = b_vec.DeltaR(elec_i);
		if (DR_jet_elec < DR_jet_elec_max){
		  is_jet_elec_overlap = true;
		  break;}
		elec_counter++; 	
	      }
	      
	      
	      if(jet->TauTag == 1){
		tau_index.push_back(j);
		ntau_counter++;
              }
	      
	      if(tau_index.size()>0){
		for(int i = 0; i < tau_index.size(); i++)
		  {
		    int tau_i = tau_index.at(i);
		    Jet *jet = (Jet*) branchJet->At(tau_i);		
		    double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
		    Tau1_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);	
		    tau_counter++;
		  }
		
	      }
              
	}    //finaliza el ciclo sobre las entradas. 

       // Loop on generator level	  
       for(int i = 0; i < branchGenParticle->GetEntriesFast(); i++){	  
       GenParticle *b = (GenParticle*) branchGenParticle->At(i);
       if(TMath::Abs(b->PID) == 5){
       double b_energy = calculateE(b->Eta, b->PT, b->Mass);
       b_vecgen.SetPtEtaPhiE(b->PT, b->Eta, b->Phi, b_energy);}

       GenParticle *e = (GenParticle*) branchGenParticle->At(i); 
       if(TMath::Abs(e->PID) == 11){
       double e_energy = calculateE(e->Eta, e->PT, e->Mass);
       elec_igen.SetPtEtaPhiE(e->PT, e->Eta, e->Phi, e_energy);}

       GenParticle *m = (GenParticle*) branchGenParticle->At(i);
       if(TMath::Abs(m->PID) == 13){
       double m_energy = calculateE(m->Eta, m->PT, m->Mass);
       muon_igen.SetPtEtaPhiE(m->PT, m->Eta, m->Phi, m_energy);}

       }//Loop generator finish
      
       //cuts delphes level	
      if(jet_i.Pt() > jet_pt_min){pass_cuts[1] = 1;} 
      if((TMath::Abs(jet_i.Eta()) < jet_eta_max) && (pass_cuts[1] == 1)){pass_cuts[2] = 1;} 
      if((b_vec.Pt() > b_pt_min) && (pass_cuts[2] == 1)){pass_cuts[3] = 1;}    
      if((TMath::Abs(b_vec.Eta()) < b_eta_max) && (pass_cuts[3] == 1)){pass_cuts[4] = 1;}
      if(pass_cuts[4] == 1){
      if((elec_i.Pt() > 0.) && (muon_i.Pt() > 0.)){pass_cuts[5] = 0;}
      else
      if((elec_i.Pt() > elec_pt_min) || (muon_i.Pt() > muon_pt_min)){pass_cuts[5] = 1;}
      }

      //cuts generation level
      if(jet_i.Pt() > jet_pt_min){pass_cuts[1] = 1;}
      if((TMath::Abs(jet_i.Eta()) < jet_eta_max) && (pass_cuts[1] == 1)){pass_cuts[2] = 1;}
      if((b_vecgen.Pt() > b_pt_min) && (pass_cuts[2] == 1)){pass_cuts[3] = 1;}
      if((TMath::Abs(b_vecgen.Eta()) < b_eta_max) && (pass_cuts[3] == 1)){pass_cuts[4] = 1;}
      if(pass_cuts[4] == 1){
      if((elec_igen.Pt() > 0.) && (muon_igen.Pt() > 0.)){pass_cuts[5] = 0;}
      else
      if((elec_igen.Pt() > elec_pt_min) || (muon_igen.Pt() > muon_pt_min)){pass_cuts[5] = 1;}
      }

      jet_met_dphi = normalizedDphi(jet_i.Phi()-MET_phi);
      jet_tau_dphi = normalizedDphi(jet_i.Phi()-Tau1_vec.Phi());
      tau_met_dphi = normalizedDphi(Tau1_vec.Phi()-MET_phi);
      bpt = b_vec.Pt()/(b_vec.Pt()+elec_i.Pt()+muon_i.Pt());
      bptgen = b_vecgen.Pt()/(b_vecgen.Pt()+elec_igen.Pt()+muon_igen.Pt());  
            
      }// Finaliza el if
    
      
      
      // save some important histograms
      for (int i = 0; i < nDir; i++){
	_hmap_Nevents[i]->Fill(0.0);
        if (pass_cuts[i] == 1){
	  
	  _hmap_n_jets[i]->Fill(njets_counter);	  
 	  if(jet_i.Pt()>0.){
	  _hmap_jet_pT[i]->Fill(jet_i.Pt());
	  _hmap_jet_eta[i]->Fill(jet_i.Eta());
	  _hmap_jet_phi[i]->Fill(jet_i.Phi());
	  }
	  _hmap_n_tau[i]->Fill(ntau_counter);
	  if(Tau1_vec.Pt()>0.){	 
          _hmap_tau1_pT[i]->Fill(Tau1_vec.Pt());
          _hmap_tau1_eta[i]->Fill(Tau1_vec.Eta());
          _hmap_tau1_phi[i]->Fill(Tau1_vec.Phi());
	  }
	  _hmap_n_b[i]->Fill(nb_counter);
	  if(b_vec.Pt()>0.){
          if((elec_i.Pt() > 0.) || (muon_i.Pt() > 0.))
	  _hmap_b_pT[i]->Fill(bpt);
	  _hmap_b_eta[i]->Fill(b_vec.Eta());
	  _hmap_b_phi[i]->Fill(b_vec.Phi());
	  _hmap_b_energy[i]->Fill(b_vec.E()*0.002);
	  }

          if(b_vecgen.Pt()>0.){
          if((elec_igen.Pt() > 0.) || (muon_igen.Pt() > 0.))
          _hmap_b_pTgen[i]->Fill(bptgen);
          _hmap_b_energygen[i]->Fill(b_vecgen.E()*0.002);
          }

	  if(elec_i.Pt()>0.){
	  _hmap_elec_pT[i]->Fill(elec_i.Pt());
	  _hmap_elec_eta[i]->Fill(elec_i.Eta());
	  _hmap_elec_phi[i]->Fill(elec_i.Phi());
	  }
  	  if(muon_i.Pt()>0.){		
	  _hmap_muon_pT[i]->Fill(muon_i.Pt());
          _hmap_muon_eta[i]->Fill(muon_i.Eta());
	  _hmap_muon_phi[i]->Fill(muon_i.Phi());
	  }
	  _hmap_met[i]->Fill(MET);
          if(jet_i.Pt()>0.){
	  _hmap_jet_met_dphi[i]->Fill(jet_met_dphi);	  
	  _hmap_met_jetmetdphi[i]->Fill(jet_met_dphi,MET);
          }
	 if(Tau1_vec.Pt()>0.){
 	 _hmap_met_taumetdphi[i]->Fill(tau_met_dphi,MET);
 	 _hmap_met_jettaudphi[i]->Fill(jet_tau_dphi,MET);
	 }  
	} 
      }
      
    }
  
  
  theFile->cd();
  for (int d = 0; d < nDir; d++)
    {
      cdDir[d]->cd();
      
      _hmap_Nevents[d]->Write();
      _hmap_n_jets[d]->Write();
      _hmap_jet_pT[d]->Write();
      _hmap_jet_eta[d]->Write();
      _hmap_jet_phi[d]->Write();
      _hmap_n_tau[d]->Write();
      _hmap_tau1_pT[d]->Write();
      _hmap_tau1_eta[d]->Write();
      _hmap_tau1_phi[d]->Write();  
      _hmap_n_b[d]->Write();
      _hmap_b_pTgen[d]->Write();
      _hmap_b_pT[d]->Write();
      _hmap_b_eta[d]->Write();
      _hmap_b_phi[d]->Write();
      _hmap_b_energy[d]->Write();
      _hmap_b_energygen[d]->Write(); 
      _hmap_elec_pT[d]->Write();
      _hmap_elec_eta[d]->Write();
      _hmap_elec_phi[d]->Write();
      _hmap_muon_pT[d]->Write();
      _hmap_muon_eta[d]->Write();
      _hmap_muon_phi[d]->Write();
      _hmap_met[d]->Write();
      _hmap_jet_met_dphi[d]->Write();
      _hmap_met_jetmetdphi[d]->Write();    
      _hmap_met_jettaudphi[d]->Write();
      _hmap_met_taumetdphi[d]->Write(); 
    }
  
  theFile->Close();
  
}

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
  else phi = abs(phi); 
  return phi;
}



void PhenoAnalysis::crateHistoMasps (int directories)
{
  for (int i = 0; i < directories; i++)
    {
      _hmap_Nevents[i]       = new TH1F("Nevents", "Nevents", 3,0.,3); 
      _hmap_jet_pT[i]        = new TH1F("jet_pT",  "P_{t} jet ", 100, 0., 1000.);
      _hmap_jet_eta[i]       = new TH1F("jet_eta", "#eta jet ", 50, -5.0, 5.0);
      _hmap_jet_phi[i]       = new TH1F("jet_phi", "#phi jet ", 70, -3.6, 3.6);
      _hmap_n_jets[i]        = new TH1F("n_jets",  "Number of jets", 11, -0.5, 10.5);
      _hmap_tau1_pT[i]       = new TH1F("tau1_pT", "p_{T} #tau_{1}", 100, 0., 1000.);
      _hmap_tau1_eta[i]      = new TH1F("tau1_eta","#eta #tau_{1}", 50, -3.5, 3.5);
      _hmap_tau1_phi[i]      = new TH1F("tau1_phi", "#phi #tau_{1}", 70, -3.6, 3.6);
      _hmap_n_tau[i]         = new TH1F("n_tau",    "Number of taus", 5, -0.5, 4.5);
      _hmap_b_pTgen[i]          = new TH1F("b_pTgen", "p_{T} b quark", 50, 0, 1);
      _hmap_b_energygen[i]      = new TH1F("b_energygen", "Energy b quark", 50, 0, 2);
      _hmap_b_pT[i]          = new TH1F("b_pT", "p_{T} b quark", 50, 0, 1);
      _hmap_b_eta[i]         = new TH1F("b_eta","#eta b quark", 50, -3.5, 3.5);
      _hmap_b_phi[i]         = new TH1F("b_phi", "#phi b quark", 70, -3.6, 3.6);
      _hmap_b_energy[i]      = new TH1F("b_energy", "Energy b quark", 50, 0, 2);
      _hmap_n_b[i]           = new TH1F("n_b",    "Number of b", 5, -0.5, 4.5);
      _hmap_elec_pT[i]          = new TH1F("elec_pT", "p_{T} e ", 50, 0., 500);
      _hmap_elec_eta[i]         = new TH1F("elec_eta","#eta e ", 50, -3.5, 3.5);
      _hmap_elec_phi[i]         = new TH1F("elec_phi", "#phi e", 70, -3.6, 3.6);      
      _hmap_muon_pT[i]          = new TH1F("muon_pT", "p_{T} #mu ", 50, 0., 500);
      _hmap_muon_eta[i]         = new TH1F("muon_eta","#eta #mu ", 50, -3.5, 3.5);
      _hmap_muon_phi[i]         = new TH1F("muon_phi", "#phi #mu", 70, -3.6, 3.6); 
      _hmap_jet_met_dphi[i]  = new TH1F("jet_met_dphi", "#Delta #phi(jet, MET)", 32, 0, 3.2);
      _hmap_met[i]	     = new TH1F("met", "MET", 200, 0. , 2000);	
      _hmap_met_jetmetdphi[i] = new TH2F("met_jetmetdphi", "#Delta #phi(jet,MET), MET", 32, 0, TMath::Pi(), 100, 0, 100);
      _hmap_met_taumetdphi[i] = new TH2F("met_taumetdphi", "#Delta #phi(#tau,MET), MET", 32, 0, TMath::Pi(), 100, 0, 100);   
      _hmap_met_jettaudphi[i] = new TH2F("met_jettaudphi", "#Delta #phi(jet, #tau), MET", 32, 0, TMath::Pi(), 100, 0, 100);
    }
}

