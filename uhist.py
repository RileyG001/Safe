import ROOT
from ROOT import TCanvas, TFile, TH1D, TLegend
import os


class hist:

    def __init__(self, h1, name):
            self.h1=h1
            self.name=name


    def axistitles(self):
        self.h1.GetYaxis().SetTitle("Weighted Number of Entries")#not all of them are, waiting till thrusday lecture when the prof/jason go over them
        if self.name in xname.keys():
            self.h1.GetXaxis().SetTitle(xname[self.name])
        else:
            self.h1.GetXaxis().SetTitle("No X-axis title")
        #self.h1.SetOption("hist C"); # self.h1.SetOption("hist E"); 
        #self.h1.Draw("hist C")
        return self.h1


    def pdata(self):
        binmax = self.h1.GetMaximumBin()
        x_max = self.h1.GetXaxis().GetBinCenter(binmax)
        y_max=self.h1.GetMaximum()
        rv="{0}{1}{2},{3}{4}{5}{6}".format(str(name),": Max: (",str(x_max),str(y_max),"), Mean: ",str(h1.GetMean(1)),"\n")
        return rv

#    def logit(self,clist,c):
 #           clist.Add(c)
            
  #          c_log = ROOT.TCanvas(f"log_{self.name}")
   #         c_log.SetLogx()
    #        c_log.SetTitle("Log")
     #       self.h1.Draw()
      #      clist.Add(c_log)
       #     return clist



xname = {
    'meff_incl' : "Met+sum of signal Jets and Leptons, (GeV)",
    'meff_4j' : "Met+1st 4 Jets, (GeV)",
    'met' : "Missing Transverse Energy, (GeV})",
    'met_phi' : "Phi of MET",
    'mTb_min' : "Min Merge and Dump",
    'mCT_bb' : "Contransverse Mass Variable, (GeV)",
    'dphi_min' : "Minimum Phi_{f}-Phi_{i}, leading 4 signal jets and met",
    'dphi_1jet' : "1st Jet's Phi_{f}-Phi_{i}",
    'm_bb' : "Invariant M of two leading b-tagged Jets, (GeV)",
    'm_non_bb' : "Invariant M of two leading non-b-tagged Jets, (GeV)",
    'ZCR_meff_4j' : "ZCR_{met}+the sume of the 1st 4 Jets",
    'ZCR_met' : "Zero Contanct Resolution for MET",
    'Z_mass' : "M of Z, (GeV)",
    'pt_lep_1' : "P_{t}, 1st Lepton, (GeV)",
    'pt_lep_2' : "P_{t}, 2nd Lepton, (GeV)",
    'eta_lep_1' : "Eta of the 1st Lepton",
    'eta_lep_2' : "Eta of the 2nd Lepton",
    'phi_lep_1' : "Phi of the 1st Lepton",
    'phi_lep_2' : "Phi of the 2nf Lepton",
    'pt_jet_1' : "P_{t}, 1st Jet, (GeV)",
    'pt_jet_2' : "P_{t}, 2nd Jet, (GeV)",
    'pt_jet_3' : "P_{t}, 3rd Jet, (GeV)",
    'pt_jet_4' : "P_{t}, 4th Jet, (GeV)",
    'pt_jet_5' : "P_{t}, 5th Jet, (GeV)",
    'pt_jet_6' : "P_{t}, 6th Jet, (GeV)",
    'eta_jet_1' : "Eta of the 1st Jet",
    'eta_jet_2' : "Eta of the 2nd Jet",
    'eta_jet_3' : "Eta of the 3rd Jet",
    'eta_jet_4' : "Eta of the 4th Jet",
    'eta_jet_5' : "Eta of the 5th Jet",
    'eta_jet_6' : "Eta of the 6th Jet",
    'phi_jet_1' : "Phi of the 1st Jet",
    'phi_jet_2' : "Phi of the 2nd Jet",
    'pt_bjet_1' : "P_{t}, 1st Jet from the bottom quark",
    'pt_bjet_2' : "P_{t}, 2nd Jet from the bottom quark",
    'pt_bjet_3' : "P_{t}, 3rd Jet from the bottom quark",
    'pt_bjet_4' : "P_{t}, 4th Jet from the bottom quark",
    'jets_n' : "Number of Jets",
    'bjets_n' : "Number of Jets from the bottom quark",
    'signal_electrons_n' : "Number of signal",
    'signal_muons_n' : "Number of signal Muons",
    'signal_taus_n' : "Number of signal Taus",
    'signal_leptons_n' : "Number of Leptons",
    'baseline_electrons_n' : "Number of baseline Electrons",
    'baseline_muons_n' : "Number of baseline Muons",
    'baseline_taus_n' : "Number of baseline Taus",
    'baseline_leptons_n' : "Number of baseline Leptons",
    'gen_filt_met' : "NA",
    'gen_filt_ht' : "NA",
    'mc_weight' : "NA",
    'dR_j1_j2' : "dR of jets 1 and 2",
    'dR_j1_j3' : "dR of jets 1 and 3",
    'dR_j1_j4' : "dR of jets 1 and 4",
    'dR_j1_j5' : "dR of jets 1 and 5",
    'dR_j1_j6' : "dR of jets 1 and 6",
    'dR_j2_j3' : "dR of jets 2 and 3",
    'dR_j2_j4' : "dR of jets 2 and 4",
    'dR_j2_j5' : "dR of jets 2 and 5",
    'dR_j2_j6' : "dR of jets 2 and 6",
    'dR_j3_j4' : "dR of jets 3 and 4",
    'dR_j3_j5' : "dR of jets 3 and 5",
    'dR_j3_j6' : "dR of jets 3 and 6",
    'dR_j4_j5' : "dR of jets 4 and 5",
    'dR_j4_j6' : "dR of jets 4 and 6",
    'dR_j5_j6' : "dR of jets 5 and 6",
    'dR_bj1_bj2' : "dR of bjets 1 and 2",
    'dR_bj1_bj3' : "dR of bjets 1 and 3",
    'dR_bj1_bj4' : "dR of bjets 1 and 4",
    'dR_bj2_bj3' : "dR of bjets 2 and 3",
    'dR_bj2_bj4' : "dR of bjets 2 and 4",
    'dR_bj3_bj4' : "dR of bjets 3 and 4",
    'dR_l1_l2' : "dR of lep 1 and 2",

}
# x=[]
# y=[]

if __name__=="__main__":
    
    input_file= "/export/nfs0home/rmgleaso/data/uclhc/uci/rmgleaso/atlas/VBS_MCTruth/run/VBS.root"
    f = ROOT.TFile(input_file)
    print("Open VBS.root")
    
    reset = True
    output = "/export/nfs0home/rmgleaso/data/uclhc/uci/rmgleaso/atlas/VBS_MCTruth/run/VBSnew.root"
    routput = output.split("/")[-1] 
    print("Open output {}".format(routput))
    if(os.path.exists(output) and reset): 
        os.remove(output)
        outp = ROOT.TFile(output, "RECREATE")
    else:
        outp = ROOT.TFile(output, "RECREATE")


    file1 = open("histoinfo.txt", "w")
    print("Open histoinfo.txt")
    # contents = file1.read().split("\n")
    # file1.seek(0)                        # <- This is the missing piece
    # file1.truncate()
    file1.write('New contents\n')
    
    print("Grabbing Keys")
    lnames = [] 
    for i in f.GetListOfKeys() :
        lnames.append(i.GetName())

    print(lnames)

    print("Loop")
    clist = ROOT.TList()
    for i in range(len(lnames)):
       
        name = lnames[i]
        print(name)
        h1=f.Get(name)
        histo=hist(h1,name)
        h1=histo.axistitles()

        file1.write(histo.pdata())

        c = ROOT.TCanvas(name)
        c.cd()
        h1.Draw("hist C")
        #c.Update()
        #clist.Add(c)
        c.Write()
       
        #clist=histo.logit(clist,c)
        c_name = "log_" + name 
        c_log = ROOT.TCanvas(c_name)
        c_log.cd()
        #print(name)
        h1.Draw("hist C")
        c_log.SetLogx()
        c_log.SetTitle("Log")
        #c_log.Update()
        c_log.Write()
        #clist.Add(c_log)
        
        #break
        #c2_log = ROOT.TCanvas(f"dumb_{name}")
        #h1.Draw()



    #clist.Write()
    #clist.Delete()
    print("Done, now go to X2Go and check the graphs")
    outp.Close()
