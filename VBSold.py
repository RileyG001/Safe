#include "SimpleAnalysisFramework/AnalysisClass.h"

DefineAnalysis(VBS);


void VBS::Init(){

  // UCI add your own variable 
  addRegions({"loose"});
  addRegions({"UCI_test_region"});

  // counting variables 
  addHistogram("jets_n", 20, -0.5, 19.5);
  addHistogram("bjets_n",20, -0.5, 19.5);
  addHistogram("signal_electrons_n", 20, -0.5, 19.5);
  addHistogram("signal_muons_n", 20, -0.5, 19.5);
  addHistogram("signal_taus_n", 20, -0.5, 19.5);
  addHistogram("signal_leptons_n", 20, -0.5, 19.5);
  addHistogram("baseline_electrons_n", 20, -0.5, 19.5);
  addHistogram("baseline_muons_n", 20, -0.5, 19.5);
  addHistogram("baseline_taus_n", 20, -0.5, 19.5);
  addHistogram("baseline_leptons_n", 20, -0.5, 19.5);

  // weight variables 
  addHistogram("gen_filt_met",200,0,2000);
  addHistogram("gen_filt_ht",200,0,2000);
  addHistogram("mc_weight", 1, 0, 1);

  addHistogram("meff_incl",200,0,2000);
  addHistogram("meff_4j",200,0,2000);
  addHistogram("met",100,0,1000);
  addHistogram("met_phi",80,-4,4);
  addHistogram("mTb_min",100,0,1000);
  addHistogram("mCT_bb",100,0,1000);
  addHistogram("dphi_min",40,0,4);
  addHistogram("dphi_1jet",40,0,4);
  addHistogram("m_bb",200,0,1000);
  addHistogram("m_non_bb",200,0,1000);
  addHistogram("ZCR_meff_4j",200,0,2000);
  addHistogram("ZCR_met",100,0,1000);
  addHistogram("Z_mass",200,0,1000);

  // lepton variables 
  addHistogram("pt_lep_1",200,0,1000);
  addHistogram("pt_lep_2",200,0,1000);
  addHistogram("eta_lep_1",80,-4,4);
  addHistogram("eta_lep_2",80,-4,4);
  addHistogram("phi_lep_1",80,-4,4);
  addHistogram("phi_lep_2",80,-4,4);

  // jet variables 
  addHistogram("pt_jet_1",200,0,2000);
  addHistogram("pt_jet_2",200,0,2000);
  addHistogram("pt_jet_3",200,0,2000);
  addHistogram("pt_jet_4",200,0,2000);
  addHistogram("pt_jet_5",200,0,2000);
  addHistogram("pt_jet_6",200,0,2000);
  addHistogram("eta_jet_1",80,-5,5);
  addHistogram("eta_jet_2",80,-5,5);
  addHistogram("eta_jet_3",80,-5,5);
  addHistogram("eta_jet_4",80,-5,5);
  addHistogram("eta_jet_5",80,-5,5);
  addHistogram("eta_jet_6",80,-5,5);

  addHistogram("phi_jet_1",80,-4,4);
  addHistogram("phi_jet_2",80,-4,4);

  // b-tagged jet variables 
  addHistogram("pt_bjet_1",200,0,2000);
  addHistogram("pt_bjet_2",200,0,2000);
  addHistogram("pt_bjet_3",200,0,2000);
  addHistogram("pt_bjet_4",200,0,2000);


  // define your histograms here 
  addHistogram("UCI_variable1",200,0,1000);

}

static int jetFlavor(AnalysisObject &jet) {
 
  if (jet.pass(TrueBJet)) return 5;
  if (jet.pass(TrueCJet)) return 4;
  if (jet.pass(TrueTau)) return 15;
  return 0;
}

void VBS::ProcessEvent(AnalysisEvent *event){

  // baseline electrons are requested to pass the loose likelihood identification criteria and have pT > 8 GeV and |eta| < 5
  // baseline muons are required to pass the Medium selections and to have pT > 20 GeV, |eta| < 5
  // baseline taus are required to pass TauRNNLoose selections and to have pt> 8 GeV |eta| < 5 
  // baseline jets are required to have pt > 20 GeV and |eta| < 5 

  // retreive your baseline objects 
  auto electrons  = event->getElectrons(8, 5, ELooseLH);
  auto muons      = event->getMuons(8, 5, MuMedium);
  auto taus       = event->getTaus(8, 5, TauRNNLoose);
  auto candJets   = event->getJets(20., 5);

  // retreive the missing transverse energy vector - which contains the magnitude and the direction
  auto metVec     = event->getMET();
  double met      = metVec.Et(); // magnitude 
  double met_phi  = metVec.Phi(); // direction

  // If there is a bad jet, then the event is thrown away and we don't process it. 
  bool there_is_a_bad_jet = countObjects(candJets, 20, 5, NOT(LooseBadJet)) !=0 ;
  if ( there_is_a_bad_jet ) return;

  // calculate some radius quantities for the jets, muon and electons - these use magic numbers but its either 0.4 or 0.04 + 10/lep.Pt().
  auto radiusCalcJet  = [] (const AnalysisObject& , const AnalysisObject& muon) { return std::min(0.4, 0.04 + 10/muon.Pt()); };
  auto radiusCalcMuon = [] (const AnalysisObject& muon, const AnalysisObject& ) { return std::min(0.4, 0.04 + 10/muon.Pt()); };
  auto radiusCalcElec = [] (const AnalysisObject& elec, const AnalysisObject& ) { return std::min(0.4, 0.04 + 10/elec.Pt()); };
  //  auto radiusCalcTau = [] (const AnalysisObject& tau, const AnalysisObject& ) { return std::min(0.4, 0.04 + 10/tau.Pt()); };

  // apply JVT
  // apply a jet-vertex-tracking requirement to the candidate jets. 
  // This means that the jet must pass some criteria of coming from the primary vertex. 
  candJets = filterObjects(candJets, 20, 5, JVT50Jet);

  // Overlap removal
  // We want our jets, electrons and muons to not be overlapping. In order to do this we do something call overlap removal. 
  // To see my information go to "SimpleAnalysisFramework/src/AnalysisClass.cxx" where the procedure is defined. 

  electrons  = overlapRemoval(electrons, muons, 0.01);
  candJets   = overlapRemoval(candJets, electrons, 0.2, NOT(BTag77MV2c10));
  electrons  = overlapRemoval(electrons, candJets, radiusCalcElec);
  candJets   = overlapRemoval(candJets, muons, radiusCalcJet, LessThan3Tracks);
  muons      = overlapRemoval(muons, candJets, radiusCalcMuon);
 
  // require signal jets to be 30 GeV
  // require signal electrons pass Medium Likelihood criteria and isolated using LooseTrackOnly
  // signal muons are required to be isolated using LooseTrackOnly

  auto signalJets      = filterObjects(candJets, 5);
  auto signalElectrons = filterObjects(electrons, 8, 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoBoosted);
  auto signalMuons     = filterObjects(muons, 8, 2.5, MuD0Sigma3|MuZ05mm|MuIsoBoosted|MuNotCosmic);
  auto signalTaus = filterObjects(taus, 20, 2.5, TauRNNLoose); 
  
  // combine into signalLeptons for easy counting
  auto signalLeptons   = signalElectrons + signalMuons + signalTaus;

  // get the b-tagged jets 
  auto candBJets    = filterObjects(signalJets, 30., 5);
  auto bjets        = filterObjects(signalJets, 30., 5, BTag77MV2c10);
  auto nonbjets     = filterObjects(signalJets, 30., 5, NOT(BTag77MV2c10));

  // count the objects that we are interested in 
  int bjets_n       = bjets.size();
  int n_nonbjets    = nonbjets.size();
  int jets_n        = signalJets.size();

  // count the baseline leptons 
  int baseline_electrons_n   = electrons.size();
  int baseline_muons_n       = muons.size();
  int baseline_taus_n       = taus.size();
  int baseline_leptons_n     = electrons.size()+muons.size()+taus.size();

  // count the signal leptons 
  int signal_electrons_n   = signalElectrons.size();
  int signal_muons_n       = signalMuons.size();
  int signal_taus_n       = signalTaus.size();
  int signal_leptons_n     = signalLeptons.size();

  // We need at least 4 jets - if there are less, then we don't process the event. 
  if(jets_n < 4) return;


  // From here we begin calculating all the variables that we want to put in our histograms. 

  // inclusive meff - all jets + leptons
  float meff_incl = met + sumObjectsPt(signalJets) + sumObjectsPt(signalLeptons);
  // meff from 4 leading jets only
  float meff_4j = met + sumObjectsPt(signalJets, 4);
  // min mT of leading 3 bjets
  float mTb_min = -999;
  if (bjets.size()>0) mTb_min = calcMTmin(bjets, metVec);
  // min mCT_bb of leading 2 bjets
  float mCT_bb = -999;
  if (bjets.size()>1) mCT_bb = calcMCT(bjets[0], bjets[1]);
  // dphimin between leading 4 signal jets and met
  float dphi_min  = minDphi(metVec, signalJets, 4);
  // dPhi(j1, MET) for Gbb
  float dphi_1jet  = minDphi(metVec, signalJets, 1);
  // invariant mass of two leading b-tagged jets
  float m_bb = -999.;
  if (bjets.size()>1) m_bb = (bjets[0] + bjets[1]).M();
  // invariant mass of two leading non-b-tagged jets
  float m_non_bb = -999.;
  if (nonbjets.size()>1) m_non_bb = (nonbjets[0] + nonbjets[1]).M();

  // leading and subleading lepton pt
  float pt_lep_1 = -999;
  if (signalLeptons.size()>0) pt_lep_1 = signalLeptons[0].Pt();
  float pt_lep_2 = -999;
  if (signalLeptons.size()>1) pt_lep_2 = signalLeptons[1].Pt();

  // leading and subleading lepton eta
  float eta_lep_1 = -999;
  if (signalLeptons.size()>0) eta_lep_1 = signalLeptons[0].Eta();
  float eta_lep_2 = -999;
  if (signalLeptons.size()>1) eta_lep_2 = signalLeptons[1].Eta();

  // leading and subleading lepton phi
  float phi_lep_1 = -999;
  if (signalLeptons.size()>0) phi_lep_1 = signalLeptons[0].Phi();
  float phi_lep_2 = -999;
  if (signalLeptons.size()>1) phi_lep_2 = signalLeptons[1].Phi();

  // leading and subleading jet pt
  float pt_jet_1 = -999;
  if (signalJets.size()>0) pt_jet_1 = signalJets[0].Pt();
  float pt_jet_2 = -999;
  if (signalJets.size()>1) pt_jet_2 = signalJets[1].Pt();
  float pt_jet_3 = -999;
  if (signalJets.size()>2) pt_jet_3 = signalJets[2].Pt();
  float pt_jet_4 = -999;
  if (signalJets.size()>3) pt_jet_4 = signalJets[3].Pt();
  float pt_jet_5 = -999;
  if (signalJets.size()>4) pt_jet_5 = signalJets[4].Pt();
  float pt_jet_6 = -999;
  if (signalJets.size()>5) pt_jet_6 = signalJets[5].Pt();

  // leading and subleading bjet pt
  float pt_bjet_1 = -999;
  if (bjets.size()>0) pt_bjet_1 = bjets[0].Pt();
  float pt_bjet_2 = -999;
  if (bjets.size()>1) pt_bjet_2 = bjets[1].Pt();
  float pt_bjet_3 = -999;
  if (bjets.size()>2) pt_bjet_3 = bjets[2].Pt();
  float pt_bjet_4 = -999;
  if (bjets.size()>3) pt_bjet_4 = bjets[3].Pt();

  // leading and subleading jet eta
  float eta_jet_1 = -999;
  if (signalJets.size()>0) eta_jet_1 = signalJets[0].Eta();
  float eta_jet_2 = -999;
  if (signalJets.size()>1) eta_jet_2 = signalJets[1].Eta();

  float eta_jet_3 = -999;
  if (signalJets.size()>2) eta_jet_3 = signalJets[2].Eta();
  float eta_jet_4 = -999;
  if (signalJets.size()>3) eta_jet_4 = signalJets[3].Eta();

  float eta_jet_5 = -999;
  if (signalJets.size()>4) eta_jet_5 = signalJets[4].Eta();
  float eta_jet_6 = -999;
  if (signalJets.size()>5) eta_jet_6 = signalJets[5].Eta();



  // leading and subleading jet phi
  float phi_jet_1 = -999;
  if (signalJets.size()>0) phi_jet_1 = signalJets[0].Phi();
  float phi_jet_2 = -999;
  if (signalJets.size()>1) phi_jet_2 = signalJets[1].Phi();

  // Z+jets CR info
  TLorentzVector tlv_Z;
  if((signal_electrons_n == 2) && (signal_muons_n < 1)) {
    tlv_Z = signalElectrons[0]+signalElectrons[1];
  }
  if((signal_muons_n == 2) && (signal_electrons_n < 1)) {
    tlv_Z = signalMuons[0]+signalMuons[1];
  }

  float Z_pt = tlv_Z.Pt();
  float Z_phi = tlv_Z.Phi();
  float Z_mass = tlv_Z.M();

  TVector2 tv2_met; tv2_met.SetMagPhi(met, metVec.Phi());
  TVector2 tv2_Z; tv2_Z.SetMagPhi(Z_pt, Z_phi);
  TVector2 total = tv2_met+tv2_Z;

  float ZCR_met = total.Mod();
  float ZCR_meff_4j = ZCR_met + sumObjectsPt(signalJets, 4);

  float gen_filt_met = event->getGenMET();
  float gen_filt_ht  = event->getGenHT();
  int channel_number = event->getMCNumber();

  // UCI variable 
  float UCI_variable1 = Z_mass;

  // After we have calculated everything, we fill the corresponding histogram with the variable value. 

  // Fill the counting variables 
  fill("jets_n", jets_n);
  fill("bjets_n", bjets_n);
  fill("signal_electrons_n", signal_electrons_n);
  fill("signal_muons_n", signal_muons_n);
  fill("signal_taus_n", signal_taus_n);
  fill("signal_leptons_n", signal_leptons_n);
  fill("baseline_electrons_n", baseline_electrons_n);
  fill("baseline_electrons_n", baseline_taus_n);
  fill("baseline_muons_n", baseline_muons_n);
  fill("baseline_leptons_n", baseline_leptons_n);

 // Fill the jet momenta
  fill("pt_jet_1", pt_jet_1);
  fill("pt_jet_2", pt_jet_2);
  fill("pt_jet_3", pt_jet_3);
  fill("pt_jet_4", pt_jet_4);
  fill("pt_jet_5", pt_jet_5);
  fill("pt_jet_6", pt_jet_6);

  // Fill the jet eta 
  fill("eta_jet_1", eta_jet_1);
  fill("eta_jet_2", eta_jet_2);
  fill("eta_jet_3", eta_jet_3);
  fill("eta_jet_4", eta_jet_4);
  fill("eta_jet_5", eta_jet_5);
  fill("eta_jet_6", eta_jet_6);

  // fill the jet phi 
  fill("phi_jet_1", phi_jet_1);
  fill("phi_jet_2", phi_jet_2);

  // Fill the bjet momenta 
  fill("pt_bjet_1", pt_bjet_1);
  fill("pt_bjet_2", pt_bjet_2);
  fill("pt_bjet_3", pt_bjet_1);
  fill("pt_bjet_4", pt_bjet_1);

  // fill the lepton pt 
  fill("pt_lep_1",pt_lep_1);
  fill("pt_lep_2",pt_lep_2);

  // fill the lepton eta 
  fill("eta_lep_1",eta_lep_1);
  fill("eta_lep_2",eta_lep_2);

  // fill the lepton phi 
  fill("phi_lep_1",phi_lep_1);
  fill("phi_lep_2",phi_lep_2);

  // fill the weight 
  fill("mc_weight", event->getMCWeights()[0]);

  // Fill the derived variables 
  fill("meff_incl", meff_incl);
  fill("meff_4j", meff_4j);
  fill("met", met);
  fill("met_phi", met_phi);
  fill("mTb_min", mTb_min);
  fill("mCT_bb", mCT_bb);
  fill("dphi_min", dphi_min);
  fill("dphi_1jet", dphi_1jet);
  fill("m_bb", m_bb);
  fill("m_non_bb", m_non_bb);
  fill("ZCR_meff_4j", ZCR_meff_4j);
  fill("ZCR_met", ZCR_met);
  fill("Z_mass", Z_mass);

  // UCI variable 
  fill("UCI_variable1", UCI_variable1);


  // categorise the event as in the region "loose"
  accept("loose");

  // UCI region
  bool your_condition = true;
  if (your_condition){
    accept("UCI_test_region");
  }


  // once we have filled the histograms, then we specify which of the histograms are added to the output root file 
  // We call output root files "n-tuples". 
  ntupVar("jets_n", jets_n);
  ntupVar("bjets_n", bjets_n);
  ntupVar("baseline_electrons_n", static_cast<int>(electrons.size()));
  ntupVar("baseline_muons_n", static_cast<int>(muons.size()));
  ntupVar("baseline_leptons_n", static_cast<int>(electrons.size()+muons.size()));
  ntupVar("signal_electrons_n", static_cast<int>(signalElectrons.size()));
  ntupVar("signal_muons_n", static_cast<int>(signalMuons.size()));
  ntupVar("signal_leptons_n", signal_leptons_n);

  ntupVar("signalJe:s0_pass_BTag77MV2c10", signalJets[0].pass(BTag77MV2c10));
  ntupVar("signalJets0_Pt", signalJets[0].Pt());
  ntupVar("signalJets3_Pt", signalJets[3].Pt());

  ntupVar("truth_id0", jetFlavor(signalJets[0]));
  ntupVar("truth_id1", jetFlavor(signalJets[1]));
  ntupVar("truth_id2", jetFlavor(signalJets[2]));
  ntupVar("truth_id3", jetFlavor(signalJets[3]));
  ntupVar("mc_weight", event->getMCWeights()[0]);

  ntupVar("gen_filt_met", gen_filt_met);
  ntupVar("gen_filt_ht", gen_filt_ht);
  ntupVar("channel_number", channel_number);

  ntupVar("meff_incl", meff_incl);
  ntupVar("meff_4j", meff_4j);
  ntupVar("met", met);
  ntupVar("met_phi", met_phi);
  ntupVar("mTb_min", mTb_min);
  ntupVar("mCT_bb", mCT_bb);
  ntupVar("dphi_min", dphi_min);
  ntupVar("dphi_1jet", dphi_1jet);
  ntupVar("m_bb", m_bb);
  ntupVar("m_non_bb", m_non_bb);
  ntupVar("pt_lep_1",pt_lep_1);
  ntupVar("pt_lep_2",pt_lep_2);
  ntupVar("eta_lep_1",eta_lep_1);
  ntupVar("eta_lep_2",eta_lep_2);
  ntupVar("phi_lep_1",phi_lep_1);
  ntupVar("phi_lep_2",phi_lep_2);
  ntupVar("ZCR_meff_4j", ZCR_meff_4j);
  ntupVar("ZCR_met", ZCR_met);
  ntupVar("Z_mass", Z_mass);

  ntupVar("pt_jet_1", pt_jet_1);
  ntupVar("pt_jet_2", pt_jet_2);
  ntupVar("pt_jet_3", pt_jet_3);
  ntupVar("pt_jet_4", pt_jet_4);
  ntupVar("pt_jet_5", pt_jet_5);
  ntupVar("pt_jet_6", pt_jet_6);

  ntupVar("eta_jet_1",eta_jet_1);
  ntupVar("eta_jet_2",eta_jet_2);
  ntupVar("phi_jet_1",phi_jet_1);
  ntupVar("phi_jet_2",phi_jet_2);

  ntupVar("pt_bjet_1", pt_bjet_1);
  ntupVar("pt_bjet_2", pt_bjet_2);
  ntupVar("pt_bjet_3", pt_bjet_3);
  ntupVar("pt_bjet_4", pt_bjet_4);

  // UCI variable 
  ntupVar("UCI_variable1", UCI_variable1);


  return;
}
