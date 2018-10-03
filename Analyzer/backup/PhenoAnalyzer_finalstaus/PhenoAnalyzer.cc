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
  int nDir = 25;
  TDirectory *theDirectory[nDir];
  theDirectory[0]  = HistoOutputFile->mkdir("Staus_No_cuts");
  theDirectory[1]  = HistoOutputFile->mkdir("Staus_with_one_Tau");
  theDirectory[2]  = HistoOutputFile->mkdir("Staus_with_two_Taus");
  theDirectory[3]  = HistoOutputFile->mkdir("Staus_After_Tau_pt_min_1Tau");
  theDirectory[4]  = HistoOutputFile->mkdir("Staus_After_Tau_pt_min_2Taus");
  theDirectory[5]  = HistoOutputFile->mkdir("Staus_After_Tau_Eta_min_1Tau");
  theDirectory[6]  = HistoOutputFile->mkdir("Staus_After_Tau_Eta_min_2Taus");
  theDirectory[7]  = HistoOutputFile->mkdir("Staus_without_elec_and_muon_1Tau");
  theDirectory[8]  = HistoOutputFile->mkdir("Staus_without_elec_and_muon_2Tau");
  theDirectory[9]  = HistoOutputFile->mkdir("Staus_After_b_jet_veto_1Tau");
  theDirectory[10]  = HistoOutputFile->mkdir("Staus_After_b_jet_veto_2Taus");
  theDirectory[11]  = HistoOutputFile->mkdir("Staus_After_jet_pt_cut_1Tau");
  theDirectory[12]  = HistoOutputFile->mkdir("Staus_After_jet_pt_cut_2Taus");
  theDirectory[13]  = HistoOutputFile->mkdir("Staus_After_jet_Eta_cut_1Taus");
  theDirectory[14]  = HistoOutputFile->mkdir("Staus_After_jet_Eta_cut_2Taus");
  theDirectory[15]  = HistoOutputFile->mkdir("Staus_After_MET_1Tau");
  theDirectory[16]  = HistoOutputFile->mkdir("Staus_After_MET_2Taus");
  theDirectory[17]  = HistoOutputFile->mkdir("Staus_Dphi_tau_jet_cut_1Tau");
  theDirectory[18]  = HistoOutputFile->mkdir("Staus_Dphi_tau_jet_cut_2Taus");
  theDirectory[19]  = HistoOutputFile->mkdir("Staus_Transverse_mass_cut_1Tau");
  theDirectory[20]  = HistoOutputFile->mkdir("Staus_Transverse_mass_cut_2Taus");
  theDirectory[21]  = HistoOutputFile->mkdir("Staus_Dphi_jet_met_1Tau");
  theDirectory[22]  = HistoOutputFile->mkdir("Staus_Dphi_jet_met_2Taus"); 
  theDirectory[23]  = HistoOutputFile->mkdir("Staus_Efective_mass_cut_1Tau");
  theDirectory[24]  = HistoOutputFile->mkdir("Staus_Efective_mass_cut_2Taus");
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
  
  double lead_jet_pt      = params->GetValue ("lead_jet_pt", 100.);
  double jet_eta_max      = params->GetValue ("jet_eta_max", 5.0);
  double tau1_pt_min      = params->GetValue ("tau1_pt_min", 10.);
  double tau2_pt_min      = params->GetValue ("tau2_pt_min", 10.);
  double tau_eta_max      = params->GetValue ("tau_eta_max", 2.1);
  double b_jet_pt_min     = params->GetValue ("b_jet_pt_min", 20.0);
  double DR_jet_lep_max  = params->GetValue ("DR_jet_lep_max", 0.3);
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
  MissingET *METpointer; 
  
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  //for(Int_t entry = 0; entry < 100; ++entry)  
     {
      //cout<<"--------------Evento: "<<entry<<"--------------"<<endl;
      treeReader->ReadEntry(entry);
      int pass_cuts[nDir] = {0};
      TLorentzVector Jet_leading_vec(0., 0., 0., 0.);
      TLorentzVector jet_i(0., 0., 0., 0.);
      TLorentzVector elec_i(0., 0., 0., 0.);
      TLorentzVector muon_i(0., 0., 0., 0.);
      TLorentzVector Tau1HadTLV (0., 0., 0., 0.);
      TLorentzVector Tau2HadTLV (0., 0., 0., 0.);    
      TLorentzVector Tau3HadTLV (0., 0., 0., 0.);
      double p_visible[2]={0., 0.};
      bool is_b_jet = false;
      bool fill_tau1 = false;
      bool fill_tau2 = false;
      bool fill_tau3 = false;
      METpointer = (MissingET*) branchMissingET->At(0);
      double Met = METpointer->MET;
      double MET = 0.;
      double MET_1 = 0.;
      double HT = 0.;
      double HT_lept = 0.;
      double Met_phi = METpointer->Phi;
      double MET_phi = 0.;
      double MET_1_phi = 0.;
      int njets_counter = 0;
      int ntau_counter = 0; 
      int nmuon_counter = 0;
      int nelec_counter = 0;
      int index_jet_lead = 0;
      double jet_min_pt = 20.0;
      
      for (int el = 0; el < branchElectron->GetEntriesFast(); el++){
	bool is_jet_elec_overlap = false;
        Electron *elec = (Electron*) branchElectron->At(el);
	double elec_energy = calculateE(elec->Eta, elec->PT, 0.000510998902);
	elec_i.SetPtEtaPhiE(elec->PT, elec->Eta, elec->Phi, elec_energy);
	if((elec_i.Pt() > 20.) && (elec_i.Eta() < 2.5)){
	  for (int j = 0; j < branchJet->GetEntriesFast(); j++)
	    { 
              Jet *jet = (Jet*) branchJet->At(j);
              double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
              jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
              double DR_jet_elec = jet_i.DeltaR(elec_i);
              if ((DR_jet_elec < DR_jet_lep_max) ){
                is_jet_elec_overlap = true;
                break;
	      } else{
                //cout << "Momento electron x "<<(elec->PT)*(cos(elec->Phi))<<endl;
                p_visible[0]+=(elec->PT)*(cos(elec->Phi));
                p_visible[1]+=(elec->PT)*(sin(elec->Phi));
                HT_lept += elec->PT;
                nelec_counter++; //number of electrons
              }
	    }
	}
      }
      
      for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++){
	bool is_jet_muon_overlap = false;
        Muon *muon = (Muon*) branchMuon->At(muo);
	double muon_energy = calculateE(muon->Eta, muon->PT, 0.1056583715);
	muon_i.SetPtEtaPhiE(muon->PT, muon->Eta, muon->Phi, muon_energy);
	if ((muon->PT > 20.) &&  (abs(muon->Eta) < 2.5)){
	  for (int j = 0; j < branchJet->GetEntriesFast(); j++)
	    {
	      Jet *jet = (Jet*) branchJet->At(j);
	      double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
	      jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);      
	      double DR_jet_muon = jet_i.DeltaR(muon_i);
	      if ((DR_jet_muon < DR_jet_lep_max) ){
		is_jet_muon_overlap = true;
		break;
	      } else{
                //cout << "Momento muon x "<<(muon->PT)*(cos(muon->Phi))<<endl;
                p_visible[0]+=(muon->PT)*(cos(muon->Phi));
                p_visible[1]+=(muon->PT)*(sin(muon->Phi));
                HT_lept += muon->PT;
                nmuon_counter++;
	      }
	    }
	}
      }
      
      //////////////////Search leading jet and taus///////////
      for (int j = 0; j < branchJet->GetEntriesFast(); j++)
	{
	  bool is_jet_elec_overlap = false;
          bool is_jet_muon_overlap = false;
	  Jet *jet = (Jet*) branchJet->At(j);
          double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
	  jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
          //cout <<j<<"  jet Pt"<<jet->PT<<endl; 
          // Remove overlaps of all jets with electrons
	  for (int el = 0; el < branchElectron->GetEntriesFast(); el++){
	    Electron *elec = (Electron*) branchElectron->At(el);
	    double elec_energy = calculateE(elec->Eta, elec->PT, 0.000510998902);
	    elec_i.SetPtEtaPhiE(elec->PT, elec->Eta, elec->Phi, elec_energy);
            if((elec_i.Pt() > 20.) && (elec_i.Eta() < 2.5)){
	      double DR_jet_elec = jet_i.DeltaR(elec_i);
	      if ((DR_jet_elec < DR_jet_lep_max) ){ 
		is_jet_elec_overlap = true;
		break;
	      }
            }
	  } 
	  
          for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++){
            Muon *muon = (Muon*) branchMuon->At(muo);
            double muon_energy = calculateE(muon->Eta, muon->PT, 0.1056583715);
            muon_i.SetPtEtaPhiE(muon->PT, muon->Eta, muon->Phi, muon_energy);
            double DR_jet_muon = jet_i.DeltaR(muon_i); 
            if ((muon->PT > 20.) &&  (abs(muon->Eta) < 2.5)){
	      if ((DR_jet_muon < DR_jet_lep_max) ){                                             
		is_jet_muon_overlap = true;
		break;
	      }
            }
	  }
	  
	  // cout << "suma momento electron y muon componente x "<<p_visible[0]<<endl;
	  if((!is_jet_elec_overlap) && (!is_jet_muon_overlap)){
	    if ((jet->PT > b_jet_pt_min) && (jet->BTag == 1)){is_b_jet = true;}
	    if (jet->TauTag == 0){
	      if(jet->PT > jet_min_pt){
                Jet_leading_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		jet_min_pt = jet->PT;
		njets_counter++;
                index_jet_lead = j;
	      }
	    }
	    //cout<<"index jet lead "<<index_jet_lead<<endl;
	  }
          //add moment of jets without leading jet, here are also the taus
          if(jet->TauTag == 1){

            ntau_counter++;
            p_visible[0] += (jet->PT)*(cos(jet->Phi));
            p_visible[1] += (jet->PT)*(sin(jet->Phi));
            HT_lept += jet->PT;

            double tau_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
            if((fill_tau1 == false) && (fill_tau2 == false) && (fill_tau3 == false)){
	      Tau1HadTLV.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, tau_energy);
	      fill_tau1 = true;
	      continue;
	    }
	    
	    if((fill_tau1 == true) && (fill_tau2 == false) && (fill_tau3 == false)){
	      Tau2HadTLV.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, tau_energy);
              fill_tau1 = true;
	      continue;
	    }
	    
	    if((fill_tau1 == true) && (fill_tau2 == true) && (fill_tau3 == false)){
	      Tau3HadTLV.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, tau_energy);
	      fill_tau3 = true;
	      continue;
	    }
          }
	}
      //cout << "suma momento electron y muon componente x "<<p_visible[0]<<endl;
      for (int j = 0; j < branchJet->GetEntriesFast(); j++)
        {
          Jet *jet = (Jet*) branchJet->At(j);
          if((j != index_jet_lead) && (jet->PT > 20.0) && (jet->TauTag == 0)){
            //cout << "componente x del momento de Jet "<<(jet->PT)*(cos(jet->Phi))<<endl;
            p_visible[0]+=(jet->PT)*(cos(jet->Phi));
            p_visible[1]+=(jet->PT)*(sin(jet->Phi));
            HT += jet->PT;
            HT_lept += jet->PT;
	  }
	  //cout<< "suma momento electron, muon componente y jets x "<<p_visible[0]<<endl;
	  //cout<<"-------------------------------------------------------"<<endl;
        }
      ///Variable MET = Met + pT_lead, where Met is the missing energy transverse standard.
      MET = Met + Jet_leading_vec.Pt();
      MET_phi = Met_phi + Jet_leading_vec.Phi();
      ////Variable MET_1: suma vectorial del momento de las particulas visibles. p_visible[0] y [1]  
      MET_1 = sqrt(p_visible[0]*p_visible[0]+p_visible[1]*p_visible[1]);  
      MET_1_phi = acos(p_visible[0]/sqrt(p_visible[0]*p_visible[0]+p_visible[1]*p_visible[1]));
      // Events with no cuts
      
      if (Jet_leading_vec.Pt() > 20.0){
	pass_cuts[0] = 1;
      } 
      // for events with exactly 1 Tau
      if ((pass_cuts[0] == 1) && (ntau_counter == 1)){
        if(Jet_leading_vec.DeltaR(Tau1HadTLV) > 0.3){
          pass_cuts[1] = 1; 
        }
      }
      ///Exactly 2 Taus
      if((pass_cuts[0] == 1) && (ntau_counter == 2)){
        if ((Jet_leading_vec.DeltaR(Tau1HadTLV) > 0.3) && (Jet_leading_vec.DeltaR(Tau2HadTLV) > 0.3)){
	  pass_cuts[2] = 1;
        }
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
      // events without electrones and muons
      if ((pass_cuts[5]==1) && (nelec_counter == 0) && (nmuon_counter == 0)){pass_cuts[7] = 1;}
      if ((pass_cuts[6]==1) && (nelec_counter == 0) && (nmuon_counter == 0)){pass_cuts[8] = 1;}
      // events passing also b-jet veto requirement
      if ((pass_cuts[7]==1) && (!is_b_jet)){pass_cuts[9] = 1;}
      if ((pass_cuts[8]==1) && (!is_b_jet)){pass_cuts[10] = 1;}
      // events passing jet pt cut
      if ((pass_cuts[9]==1) && (Jet_leading_vec.Pt() > lead_jet_pt)){pass_cuts[11] = 1;}
      if ((pass_cuts[10]==1) && (Jet_leading_vec.Pt() > lead_jet_pt)){pass_cuts[12] = 1;}
      // events passing jet eta cut
      if ((pass_cuts[11]==1) && (abs(Jet_leading_vec.Eta()) < jet_eta_max) ){pass_cuts[13] =1;}
      if ((pass_cuts[12]==1) && (abs(Jet_leading_vec.Eta()) < jet_eta_max) ){pass_cuts[14] =1;}
      // events passing met cut 
      if ((pass_cuts[13] == 1) && (MET > met_min)){pass_cuts[15] = 1;}
      if ((pass_cuts[14] == 1) && (MET > met_min)){pass_cuts[16] = 1;}
      // events passing dphitaujet cut 
      double tau1_jet_dphi = abs(normalizedDphi(Jet_leading_vec.Phi() - Tau1HadTLV.Phi()));
      if ((pass_cuts[15] == 1) && (tau1_jet_dphi > Dphitaujetmax)){pass_cuts[17] = 1;}
      if ((pass_cuts[16] == 1) && (tau1_jet_dphi > Dphitaujetmax)){pass_cuts[18] = 1;}
      // events passing mt cut
      double transmass = TMath::Sqrt(TMath::Abs(2*Tau1HadTLV.Pt()*MET*(1-TMath::Cos(TMath::Abs(Tau1HadTLV.Phi() - Met_phi)))));
      if ((pass_cuts[15] == 1) && (transmass > transversemassmin)){pass_cuts[19] = 1;}
      if ((pass_cuts[16] == 1) && (transmass > transversemassmin)){pass_cuts[20] = 1;}
      // events passing dphijetmet cut
      double jet_met_dphi = abs(normalizedDphi(Jet_leading_vec.Phi() - Met_phi));
      if ((pass_cuts[15] == 1) && (jet_met_dphi > deltajetmetphi)){pass_cuts[21]=1;}
      if ((pass_cuts[16] == 1) && (jet_met_dphi > deltajetmetphi)){pass_cuts[22]=1;}
      // events passing effective mass cut
      double efective_mass = TMath::Sqrt(Jet_leading_vec.Pt()*Jet_leading_vec.Pt()+Tau1HadTLV.Pt()*Tau1HadTLV.Pt()+MET*MET);
      if ((pass_cuts[15] == 1) && (efective_mass > efectivemassmin)){pass_cuts[23] = 1;}        
      if ((pass_cuts[16] == 1) && (efective_mass > efectivemassmin)){pass_cuts[24] = 1;}
      
      
      for (int i = 0; i < nDir; i++){
	_hmap_Nevents[i]->Fill(0.0);
	_hmap_n_jets[i]->Fill(njets_counter);
	_hmap_n_tau[i]->Fill(ntau_counter);
	
	if ( pass_cuts[i] == 1){
	  double jet_met_dphi = abs(normalizedDphi(Jet_leading_vec.Phi() - MET_phi));
	  _hmap_jet_met_Dphi[i]->Fill(abs(jet_met_dphi));
	  _hmap_jet_met_metDphi[i]->Fill(abs(jet_met_dphi),MET);
	  _hmap_Nevents[i]->Fill(1.0);
          _hmap_met[i]->Fill(Met);
	  _hmap_MET[i]->Fill(MET);
          _hmap_MET_1[i]->Fill(MET_1);           
	  _hmap_lead_jet_pT[i]->Fill(Jet_leading_vec.Pt());
	  _hmap_lead_jet_eta[i]->Fill(Jet_leading_vec.Eta());
	  _hmap_lead_jet_phi[i]->Fill(Jet_leading_vec.Phi());
	  _hmap_rm[i]->Fill(Met/Jet_leading_vec.Pt());
          _hmap_RM[i]->Fill(MET/Jet_leading_vec.Pt());
          _hmap_RM_1[i]->Fill(MET_1/Jet_leading_vec.Pt());
          _hmap_met_HT[i]->Fill(Met/HT);
          _hmap_MET_HT[i]->Fill(MET/HT);
          _hmap_MET_1_HT[i]->Fill(MET_1/HT);
          _hmap_met_sqrtHT[i]->Fill(Met/TMath::Sqrt(HT));
          _hmap_MET_sqrtHT[i]->Fill(MET/TMath::Sqrt(HT));
          _hmap_MET_1_sqrtHT[i]->Fill(MET_1/TMath::Sqrt(HT));         
          _hmap_jpT_MET[i]->Fill(Jet_leading_vec.Pt(), MET);
	  if(Tau1HadTLV.Pt() > 1.){ 
	    _hmap_tau1_pT[i]->Fill(Tau1HadTLV.Pt());
	    _hmap_tau1_eta[i]->Fill(Tau1HadTLV.Eta());
	    _hmap_tau1_phi[i]->Fill(Tau1HadTLV.Phi());
	    double tau1_met_dphi = abs(normalizedDphi(Tau1HadTLV.Phi() - MET_phi));
	    _hmap_tau1_met_Dphi[i]->Fill(tau1_met_dphi);
	    _hmap_Dphi_tau_met_with_met[i]->Fill(tau1_met_dphi, MET);
            double transmass = TMath::Sqrt(TMath::Abs(2*Tau1HadTLV.Pt()*Met*(1-TMath::Cos(TMath::Abs(Tau1HadTLV.Phi() - Met_phi)))));
	    double TransMass = TMath::Sqrt(TMath::Abs(2*Tau1HadTLV.Pt()*MET*(1-TMath::Cos(TMath::Abs(Tau1HadTLV.Phi() - MET_phi)))));
            double TransMass_1 = TMath::Sqrt(TMath::Abs(2*Tau1HadTLV.Pt()*MET*(1-TMath::Cos(TMath::Abs(Tau1HadTLV.Phi() - MET_1_phi)))));
            _hmap_transverse_mass[i]->Fill(transmass);
            _hmap_Transverse_Mass[i]->Fill(TransMass);  
            _hmap_Transverse_Mass_1[i]->Fill(TransMass_1);    
            _hmap_MT_with_MET[i]->Fill(transmass, MET);
            double MT_on_pT_lead = transmass/Jet_leading_vec.Pt();
            _hmap_MT_on_pT_lead[i]->Fill(MT_on_pT_lead);
            _hmap_MT_with_Dphi_tau_met[i]->Fill(tau1_met_dphi, transmass);
            _hmap_MT_with_Dphi_jet_met[i]->Fill(jet_met_dphi, transmass);
	    double tau1_jet_Dphi = abs(normalizedDphi(Jet_leading_vec.Phi() - Tau1HadTLV.Phi()));
	    _hmap_tau1_jet_Dphi[i]->Fill(tau1_jet_Dphi);
            _hmap_Dphi_tau_jet_with_jPt[i]->Fill(tau1_jet_Dphi, Jet_leading_vec.Pt()); 
	    double efective_mass = TMath::Sqrt(Jet_leading_vec.Pt()*Jet_leading_vec.Pt()+Tau1HadTLV.Pt()*Tau1HadTLV.Pt()+MET*MET);
	    _hmap_efective_mass[i]->Fill(efective_mass); 
            _hmap_Meff_with_MET[i]->Fill(efective_mass, MET);      
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
  for (int d = 0; d < nDir; d++)
    {
      cdDir[d]->cd();
      _hmap_Nevents[d]->Write();
      _hmap_n_jets[d]->Write();
      _hmap_n_tau[d]->Write();
      _hmap_lead_jet_pT[d]->Write();
      _hmap_lead_jet_eta[d]->Write();
      _hmap_lead_jet_phi[d]->Write();
      _hmap_rm[d]->Write();
      _hmap_RM[d]->Write();
      _hmap_RM_1[d]->Write();
      _hmap_tau1_pT[d]->Write();
      _hmap_tau1_eta[d]->Write();
      _hmap_tau1_phi[d]->Write();
      _hmap_tau2_pT[d]->Write();
      _hmap_tau2_eta[d]->Write();
      _hmap_tau2_phi[d]->Write();
      _hmap_jpT_MET[d]->Write();
      _hmap_jet_met_Dphi[d]->Write();
      _hmap_met[d]->Write();
      _hmap_MET[d]->Write();
      _hmap_MET_1[d]->Write();
      _hmap_met_HT[d]->Write();
      _hmap_MET_HT[d]->Write(); 
      _hmap_MET_1_HT[d]->Write();
      _hmap_met_sqrtHT[d]->Write();
      _hmap_MET_sqrtHT[d]->Write();
      _hmap_MET_1_sqrtHT[d]->Write();
      _hmap_tau1_met_Dphi[d]->Write();
      _hmap_tau2_met_Dphi[d]->Write();
      _hmap_jet_met_metDphi[d]->Write();
      _hmap_tau1_jet_Dphi[d]->Write();
      _hmap_transverse_mass[d]->Write();
      _hmap_Transverse_Mass[d]->Write();
      _hmap_Transverse_Mass_1[d]->Write();
      _hmap_MT_on_pT_lead[d]->Write();
      _hmap_efective_mass[d]->Write();
      _hmap_Dphi_tau_met_with_met[d]->Write();
      _hmap_Dphi_tau_jet_with_jPt[d]->Write();     
      _hmap_MT_with_MET[d]->Write();
      _hmap_MT_with_Dphi_tau_met[d]->Write();
      _hmap_MT_with_Dphi_jet_met[d]->Write();         
      _hmap_Meff_with_MET[d]->Write(); 
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
      _hmap_rm[i]            = new TH1F("rm",   "rm", 100, 0, 10.);
      _hmap_RM[i]            = new TH1F("RM",   "RM", 100, 0, 10.); 
      _hmap_RM_1[i]          = new TH1F("RM_1",   "RM_1", 100, 0, 10.);
      _hmap_n_jets[i]        = new TH1F("N_jets",         "N(jet)", 4, 0., 4);
      _hmap_n_tau[i]         = new TH1F("N_taus",          "", 4, 0., 4);
      _hmap_tau1_pT[i]       = new TH1F("tau1_pT",        "p_{T}(#tau_{1})", 200, 0., 600.);
      _hmap_tau1_eta[i]      = new TH1F("tau1_eta",       "#eta(#tau_{1})", 50, -3.5, 3.5);
      _hmap_tau1_phi[i]      = new TH1F("tau1_phi",       "#phi(#tau_{1})", 70, -3.6, 3.6);
      _hmap_tau2_pT[i]       = new TH1F("tau2_pT",        "p_{T}(#tau_{2})", 150, 0., 300.);
      _hmap_tau2_eta[i]      = new TH1F("tau2_eta",       "#eta(#tau_{2})", 50, -3.5, 3.5);
      _hmap_tau2_phi[i]      = new TH1F("tau2_phi",       "#phi(#tau_{2})", 70, -3.6, 3.6);
      _hmap_jet_met_Dphi[i]  = new TH1F("jet_met_Dphi",   "#Delta #phi(jet, MET)", 32, 0, 3.2);
      _hmap_transverse_mass[i] = new TH1F("transverse_mass",   "transverse mass", 60, 0, 600);
      _hmap_Transverse_Mass[i] = new TH1F("Transverse_Mass",   "Transverse Mass", 60, 0, 600);
      _hmap_Transverse_Mass_1[i] = new TH1F("Transverse_Mass_1",   "Transverse Mass_1", 60, 0, 600);
      _hmap_met[i]           = new TH1F("met",            "met", 120, 0., 1200.); 
      _hmap_MET[i]           = new TH1F("MET",            "MET", 120, 0., 1200.);
      _hmap_MET_1[i]           = new TH1F("MET_1",            "MET_1", 1200, 0., 1200.);
      _hmap_met_HT[i]        = new TH1F("met_HT",   "met_HT", 100, 0, 10.);
      _hmap_MET_HT[i]        = new TH1F("MET_HT",   "MET_HT", 100, 0, 10.);
      _hmap_MET_1_HT[i]        = new TH1F("MET_1_HT",   "MET_1_HT", 100, 0, 30.);
      _hmap_met_sqrtHT[i]        = new TH1F("met_sqrtHT",   "met_sqrtHT", 100, 0, 20.);
      _hmap_MET_sqrtHT[i]        = new TH1F("MET_sqrtHT",   "MET_sqrtHT", 100, 0, 50.);
      _hmap_MET_1_sqrtHT[i]        = new TH1F("MET_1_sqrtHT",   "MET_1_sqrtHT", 100, 0, 20.);
      _hmap_MT_on_pT_lead[i] = new TH1F("MT_on_pT_lead",  "MT_on_pT_lead", 64, 0., 8.);
      _hmap_tau1_met_Dphi[i]  = new TH1F("tau_met_Dphi",   "#Delta #phi(#tau_{1}, MET)", 32, 0, 3.2);
      _hmap_tau2_met_Dphi[i]  = new TH1F("tau2_met_Dphi",   "#Delta #phi(#tau_{2}, MET)", 32, 0, 3.2);
      _hmap_tau1_jet_Dphi[i]  = new TH1F("tau1_jet_Dphi",   "#Delta #phi(#tau_{1}, jet)", 64, 0, 3.2);
      _hmap_efective_mass[i]  = new TH1F("efective_mass", "Efective_mass", 80, 0, 800); 
      _hmap_jet_met_metDphi[i] = new TH2F("jet_met_metDphi", "#Delta #phi(jet,MET), MET", 128, 0, 3.2, 150, 0, 300);
      _hmap_Dphi_tau_met_with_met[i] = new TH2F("Dphi_taumet_met", "#Delta #phi(#tau,MET), MET", 128, 0, 3.2, 150, 0, 300);
      _hmap_Dphi_tau_jet_with_jPt[i] = new TH2F("Dphi_taujet_jetPt", "#Delta #phi(#tau,j_{1}), PT(j_{1})", 128, 0, 3.2, 150, 0, 300);
      _hmap_MT_with_MET[i] = new TH2F("MT_MET", "M_{T}, MET", 150, 0, 300, 150, 0, 300);
      _hmap_MT_with_Dphi_tau_met[i] = new TH2F("MT_Dphi_tau_met", "#Delta #phi(#tau,MET), M_{T}", 128, 0, 3.2, 150, 0, 300);
      _hmap_MT_with_Dphi_jet_met[i] = new TH2F("MT_Dphi_tau_jet", "#Delta #phi(#tau,j_{1}) M_{T}", 128, 0, 3.2, 150, 0, 300);
      _hmap_Meff_with_MET[i] = new TH2F("Meff_MET", "Meff, MET", 150, 0, 300, 150, 0, 300); 
      _hmap_jpT_MET[i] = new TH2F("jpT_MET", "jpT, MET", 200, 0, 400, 200, 0, 400);    
    }
}
