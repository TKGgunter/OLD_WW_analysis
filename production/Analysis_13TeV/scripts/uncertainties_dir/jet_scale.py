#Thoth Gunter

import os, sys, pickle, copy
import datetime

cdir = os.getcwd()
cdir = "/home/gunter/WW_analysis/production/Analysis_13TeV/scripts/uncertainties_dir"

sys.path.append(cdir + "/../../")
from prep_ana_II import *
sys.path.append(cdir + "/../../tools/")
sys.path.append(os.getcwd() + "/../../tools/fit/")
sys.path.append(cdir + "/../../tools/jecsys/")
from JECUncertainty import JECUncertainty, jecUncertainties
from lepton_eff import muonEff, electronEff
from pile_up import pileUpFunction, Min_Bias
from cross_section_calc import calc_cross_stuff, cross_calc, pseudo_data_yield_sum 
import fit


def create_pseudo(df, df_da, scales):
  pseudo = {}
  pseudo["tot"] = pseudo_data_yield_sum(df, df_da, scales=scales)
  pseudo["j0"] =  pseudo_data_yield_sum(df[df.numb_jets == 0], df_da[df_da.numb_jets == 0], scales=scales)
  pseudo["j1"] =  pseudo_data_yield_sum(df[df.numb_jets == 1], df_da[df_da.numb_jets == 1], scales=scales)
  pseudo["j2"] =  pseudo_data_yield_sum(df[df.numb_jets == 2], df_da[df_da.numb_jets == 2], scales=scales)
  return pseudo

def create_Xs( ana_obj, pseudo):
  Xs = {"tot":0, "j0":0, "j1":0, "j2":0}
  scales = ana_obj.scales
  df = ana_obj.df
  df_da = ana_obj.df_da
  df_ww = ana_obj.df_ww
  df_ggww = ana_obj.df_ggww

  Xs["tot"] = cross_calc(df, df_da, df_ww, df_ggww, scales,  fiducial=True, pseudo=pseudo["tot"])  
  Xs["j0"] =  cross_calc(df[df.numb_jets == 0], df_da[df_da.numb_jets == 0], df_ww[df_ww.numb_jets == 0], df_ggww[df_ggww.numb_jets == 0], scales,  fiducial=True, pseudo=pseudo["j0"])  
  Xs["j1"] =  cross_calc(df[df.numb_jets == 1], df_da[df_da.numb_jets == 1], df_ww[df_ww.numb_jets == 1], df_ggww[df_ggww.numb_jets == 1], scales,  fiducial=True, pseudo=pseudo["j1"])  
  Xs["j2"] =  cross_calc(df[df.numb_jets == 2], df_da[df_da.numb_jets == 2], df_ww[df_ww.numb_jets == 2], df_ggww[df_ggww.numb_jets == 2], scales,  fiducial=True, pseudo=pseudo["j2"])  
  return Xs


"""
ToDo: ?
+?Add unc term to prep_ana
+?jec_setup should be moved to tools
+need more documentation
"""
jet_pt_cut_global = 30 # 31.2


class jec_setup():
  """
  jec_setup:
  This object loads the data, mc, and jet energy correction
  module and data files neded to complete the jet scale measurement. 
   
  """
  def __init__(self, unc="jet", flavor=""):

    self.flavor = flavor
    self.unc_type = unc

    ana_obj = analysis_setup( unc )
    self.scales = ana_obj.scales
    self.df     = ana_obj.df
    self.df_da  = ana_obj.df_da
    self.df_ww  = ana_obj.df_ww 
    self.df_ggww  = ana_obj.df_ggww 

    self.columns = columns_jets_unc
    if unc == "lhe":
      self.columns = columns_lhe

    if flavor == "diff":
      df = self.df
      df_da = self.df_da
      df_ww = self.df_ww
      
      df = df[df.lep1_type != df.lep2_type]
      df_da = df_da[df_da.lep1_type != df_da.lep2_type]
      df_ww = df_ww[df_ww.lep1_type != df_ww.lep2_type]

      self.df = df
      self.df_da = df_da
      self.df_ww = df_ww

    if flavor == "same":
      df = self.df
      df_da = self.df_da
      df_ww = self.df_ww
      
      df = df[df.lep1_type == df.lep2_type]
      df_da = df_da[df_da.lep1_type == df_da.lep2_type]
      df_ww = df_ww[df_ww.lep1_type == df_ww.lep2_type]

      self.df = df
      self.df_da = df_da
      self.df_ww = df_ww

    self.df["org_numb_jets"]    = 1 * self.df.numb_jets.values
    self.df_da["org_numb_jets"] = 1 * self.df_da.numb_jets.values

    self.rfs            = load_randomForest()
    self.jec_obj        = JECUncertainty()
    
    self.verbose        = False
    self.JECUnc         = None

    self.path = "/home/gunter/WW_analysis/production/Analysis_13TeV/scripts/uncertainties_dir/"


  def reset(self):
    ana_obj = analysis_setup( self.unc_type )
    self.df     = ana_obj.df
    self.df_da  = ana_obj.df_da
    self.df_ww  = ana_obj.df_ww 
  
    if self.flavor == "diff":
      df = self.df
      df_da = self.df_da
      df_ww = self.df_ww
      
      df = df[df.lep1_type != df.lep2_type]
      df_da = df_da[df_da.lep1_type != df_da.lep2_type]
      df_ww = df_ww[df_ww.lep1_type != df_ww.lep2_type]

      self.df = df
      self.df_da = df_da
      self.df_ww = df_ww

    if self.flavor == "same":
      df = self.df
      df_da = self.df_da
      df_ww = self.df_ww
      
      df = df[df.lep1_type == df.lep2_type]
      df_da = df_da[df_da.lep1_type == df_da.lep2_type]
      df_ww = df_ww[df_ww.lep1_type == df_ww.lep2_type]

      self.df = df
      self.df_da = df_da
      self.df_ww = df_ww

  def loadJECUnc(self):
    self.JECUnc = pickle.load(open(self.path+"data/jecDump.pkl", "r")) 

  def calcJECUnc(self):
    df_ww    = self.df_ww
    pt_data  = df_ww[df_ww.jet1_pt > 30].jet1_pt.values[:5000]
    eta_data = df_ww[df_ww.jet1_pt > 30].jet1_eta.values[:5000] 
    pt_list  = list(pt_data) 
    eta_list = list(eta_data) 
    jec_unc  = np.array(jecUncertainties( pt_list, eta_list, self.jec_obj))
    
    #Migrate through 
    pt_bins = [10] * 7 + [50] * 6 + [10000]
    eta_bins = [1] * 10

    bins = np.zeros((len(pt_bins), len(eta_bins)))
    #Iterate over bins and average ....
    pt_ = 30.
    for it, pt in enumerate(pt_bins):
      #print "pt:", pt_, pt_ + pt
      eta_ = -5.
      for jt, eta in enumerate(eta_bins):
        #print "eta:", eta_, eta_ + eta
         
        bins[it][jt] = jec_unc[(pt_data >= pt_) & (pt_data < pt_ + pt) & (eta_data >= eta_) & (eta_data < eta_ + eta)].mean()
        eta_ += eta
      pt_ += pt
        
    bins[np.isnan(bins)] = bins[~np.isnan(bins)].max()
    return bins, pt_bins, eta_bins, pt_, eta_



  def saveJECUnc(self):
    """
    Saves the binned jet energy corrections from a subset of the ww data frame.
    Format: dict([bins, pt_bins, eta_bins, inital pt, initial eta])
    """
    bins, pt_bins, eta_bins, pt_, eta_ = self.calcJECUnc()
    pickle.dump( {"bins": bins, "pt_bins": pt_bins, "eta_bins": eta_bins, "pt_": pt_, "eta_": eta_}, open(self.path+"data/jecDump.pkl", "w"))


  def rf_ana(self, df, fDY, fTT):
    """
    Custom random forest selector that allows user to switch between different rf tags. 
    """
    fDY = df[fDY] > .9
    fTT = df[fTT] > .6
    
    return df[fDY & fTT]


  def plot_bins(self):
    if self.JECUnc == None:
      print "JECUnc has not been initialized"
      return 
    bins = self.JECUnc["bins"]
    pt_bins = self.JECUnc["pt_bins"]
    eta_bins = self.JECUnc["eta_bins"]

    plt.figure(1, figsize=(25,7))
    plt.imshow(np.rot90(bins), interpolation='nearest')


    x = [i for i in range(14)]
    x_l = []
    x_ = 30
    for i in range(14):
        x_ = x_ + pt_bins[i]
        x_l.append( x_ )

    plt.xticks(x, x_l)

    y = [i for i in range(10)]
    plt.yticks(y, np.array(eta_bins) + np.array(y) - 5.5)

    plt.ylabel("eta")
    plt.xlabel("pt")

    plt.colorbar() 
    plt.savefig("data/2d_jec_unc_binned.png")
    plt.show()


  def resultsJEC(self):
    """
    Official JEC scale uncertainties
    """
    df_ww = self.df_ww
    jet_scale_shift_new(df_ww, rf=rfs)

    print df_ww[["jet1_pt", "jet1_unc", "jet2_pt", "jet2_unc", "numb_jets"]][:10]
    print np.histogram(df_ww.numb_jets, range=(-.5, 6.5), bins =7)[0]
    print np.histogram(rf_ana(df_ww, "fDY", "fTT").numb_jets, range=(-.5, 6.5), bins =7)[0]
###############################################
###############################################

def kill_jets( df, pt_cut= jet_pt_cut_global ):
  """
  Edit number of jets in event based on pt cut
  """
  #?Not quite what I should be doing..
  n_jets = np.zeros(df.shape[0]) #np.maximum(np.zeros(df.shape[0]), df.numb_jets.values - 6)
  for k in df.keys():
    if "jet" in k[:3] and "pt" in k:
      cut = (df[k] > pt_cut)
      n_jets[cut.values] = n_jets[cut.values] + 1 

  df["numb_jets"] = n_jets



def jet_scale_shift_flat(data, jet_pt_shift=1., pt_cut=jet_pt_cut_global, rf=None):
  """
  Jet scale flat shift 
  """  
  #?Is this working as advertized
  vec_ht_values_sin_orig = copy.copy(data.HT.values) * 0
  vec_ht_values_cos_orig = copy.copy(data.HT.values) * 0

  vec_ht_values_sin_post = copy.copy(data.HT.values) * 0
  vec_ht_values_cos_post = copy.copy(data.HT.values) * 0


  data.HT = data.HT * 0
  #Scale pt of each jet
  for k in data.keys():
    if "jet" in k[:3] and "pt" in k:

      #Original vec_ht values Energy scale  
      vec_ht_values_sin_orig += data[k] * np.sin(data[k[:4] + "_phi"]) 
      vec_ht_values_cos_orig += data[k] * np.cos(data[k[:4] + "_phi"]) 
      

      data[k] = data[k] * jet_pt_shift
      #NEW TO CORRECT FOR EVENTS WITH LOST JETS
      ht_lost_jet = data[k] >= pt_cut
      data.HT.values[ht_lost_jet] = data[ht_lost_jet].HT + data[ht_lost_jet][k]

      #Post uncertainty corrected  vec_ht values 
      vec_ht_values_sin_post += data[k] * np.sin(data[k[:4] + "_phi"]) 
      vec_ht_values_cos_post += data[k] * np.cos(data[k[:4] + "_phi"]) 


  #MET results
  #?Doesn't exist yet 
  previous_met = copy.copy(data.metMod.values)
  data.metMod  = np.sqrt( (data.metMod * np.sin(data.met_phi) - vec_ht_values_sin_orig + vec_ht_values_sin_post )**2 +\
                          (data.metMod * np.cos(data.met_phi) - vec_ht_values_cos_orig + vec_ht_values_cos_post )**2 ) 
  data.METProj = data.metMod.values / previous_met * data.METProj
  data.recoil  = np.sqrt( (data.metMod * np.sin(data.met_phi) + data.lep1_pt * np.sin(data.lep1_phi) + data.lep2_pt * np.sin(data.lep2_phi))**2 +\
                          (data.metMod * np.cos(data.met_phi) + data.lep1_pt * np.cos(data.lep1_phi) + data.lep2_pt * np.cos(data.lep2_phi))**2 ) 

  
  #Update number of jets
  kill_jets( data, pt_cut )

  data.dPhiMETJet.values[data.numb_jets.values == 0] = -1
  data.dPhiLLJet.values[data.numb_jets.values == 0] = -1

  if rf != None:
  #print "Recreating random forest scores."
    pred_fTT = rf["clf_fTT"].predict_proba(np.float32(data[rf["features_fTT"]].values))
    data["pred_fTT_WW"] = pred_fTT[:,0]

    temp = data[rf["features_fDY"]]
    temp = temp.replace([np.inf,-np.inf], 0)
    pred_fDY = rf["clf_fDY"].predict_proba(np.float32(temp.values))
    data["pred_fDY_WW"] = pred_fDY[:,0]




def jet_scale_shift_official(data, pt_cut=jet_pt_cut_global, up_down="up", rf=None):
  #?Is this working as advertized
  vec_ht_values_sin_orig = copy.copy(data.HT.values) * 0
  vec_ht_values_cos_orig = copy.copy(data.HT.values) * 0

  vec_ht_values_sin_post = copy.copy(data.HT.values) * 0
  vec_ht_values_cos_post = copy.copy(data.HT.values) * 0


  if "jet1_unc" not in data.keys():
    print data.keys()
    print "\033[1mJet1 unc is not available...\033[0m"
  ht_values = copy.copy(data.HT.values)
  data.HT = data.HT * 0
  #Scale pt of each jet
  if up_down == "up":
    thingoo = 1
  else:
    thingoo = -1
  for k in data.keys():
    if "bjet" in k: continue
    if "jet" in k[:3] and "unc" in k:

      #Original vec_ht values Energy scale  
      vec_ht_values_sin_orig += data[k[:4] + "_pt"] * np.sin(data[k[:4] + "_phi"]) 
      vec_ht_values_cos_orig += data[k[:4] + "_pt"] * np.cos(data[k[:4] + "_phi"]) 
      

      data[k[:-3]+"pt"] = data[k[:-3]+"pt"] * np.abs(thingoo + data[k])
  #######?Same loop maybe
  for k in data.keys():
    if "jet" in k[:3] and "pt" in k:
      ht_lost_jet = data[k] >= pt_cut
      data.HT.values[ht_lost_jet] = data[ht_lost_jet].HT.values + data[ht_lost_jet][k].values

      #Post uncertainty corrected  vec_ht values 
      vec_ht_values_sin_post += data[k] * np.sin(data[k[:4] + "_phi"]) 
      vec_ht_values_cos_post += data[k] * np.cos(data[k[:4] + "_phi"]) 
      #print "SIN STUFF", k,  np.sin(data[data[k[:4] + "_pt"] > 0][k[:4] + "_phi"].values[:5])


  #MET results
  #?met projected isn't right need to change
  previous_met = copy.copy(data.metMod.values )
  data.metMod  = np.sqrt( (data.metMod * np.sin(data.met_phi) - vec_ht_values_sin_orig + vec_ht_values_sin_post )**2 +\
                          (data.metMod * np.cos(data.met_phi) - vec_ht_values_cos_orig + vec_ht_values_cos_post )**2 ) 
  data.METProj = data.metMod.values / previous_met * data.METProj # np.abs(data.METProj - (ht_values - data.HT))
  data.recoil  = np.sqrt( (data.metMod * np.sin(data.met_phi) + data.lep1_pt * np.sin(data.lep1_phi) + data.lep2_pt * np.sin(data.lep2_phi))**2 +\
                          (data.metMod * np.cos(data.met_phi) + data.lep1_pt * np.cos(data.lep1_phi) + data.lep2_pt * np.cos(data.lep2_phi))**2 ) 


  
  #Update number of jets
  kill_jets( data, pt_cut )

  #Add dPhiXJet stuff 
  data.dPhiMETJet.values[data.numb_jets.values == 0] = -1
  data.dPhiLLJet.values[data.numb_jets.values == 0] = -1


  if rf != None:
  #print "Recreating random forest scores."
    pred_fTT = rf["clf_fTT"].predict_proba(np.float32(data[rf["features_fTT"]].values))
    data["fTT"] = pred_fTT[:,0]

    temp = data[rf["features_fDY"]]
    temp = temp.replace([np.inf,-np.inf], 0)
    pred_fDY = rf["clf_fDY"].predict_proba(np.float32(temp.values))
    data["fDY"] = pred_fDY[:,0]



def apply_binned_jet_uncertainty( df, bins, pt_bins, eta_bins):
  #Apply binned uncertainty to WW and see how jet bin hist change
  #bins, pt_bins, eta_bins

  for k in df.keys():
    if "jet" in k[:3] and "eta" in k:
      df[ k[:-3]+"unc"] = 0.

  pt_ = 30
  for it, pt in enumerate(pt_bins):
    #print "pt:", pt_, pt_ + pt
    eta_ = -5
    for jt, eta in enumerate(eta_bins):
      #print "eta:", eta_, eta_ + eta
      for k in df.keys():
        if "jet" in k[:3] and "eta" in k:
          if pt_ == 30:
            selection = (df[k[:-3] + "pt"] >= 20.) &\
                        (df[k[:-3] +"pt"] < pt_ + pt) &\
                        (df[k] >= eta_) &\
                        (df[k] < eta_ + eta)
          else: 
            selection = (df[k[:-3] + "pt"] >= pt_) &\
                        (df[k[:-3] +"pt"] < pt_ + pt) &\
                        (df[k] >= eta_) &\
                        (df[k] < eta_ + eta)
          df[k[:-3]+"unc"].values[selection] = bins[it][jt]

      eta_ += eta
    pt_ += pt
    







def jet_process_njets(df, process, pre_rf="rf"):
  np.set_printoptions(suppress=True)
  if pre_rf == "rf":
    tmp_df = df[(df.process == process) & ( df.type.str.contains("rf"))][["process", "type", "njets"]].reset_index(drop=True)
  else:
    tmp_df = df[(df.process == process) & ( ~df.type.str.contains("rf"))][["process", "type", "njets"]].reset_index(drop=True) 
    
  print tmp_df 
  delta = (np.abs( tmp_df.njets[0] - tmp_df.njets[1]) + np.abs(tmp_df.njets[0] - tmp_df.njets[2]) ) / 2.
  print "delta", delta
  print "% delta", delta / tmp_df.njets[0] * 100. 






def preselection_JES( jec_obj, file_name="", verbose=True ):


  #Different flavor cut
  #Mll cut 30
  #lep 1 > 25
  #lep 2 > 20
  jec_obj.df = jec_obj.df[(jec_obj.df.mll > 30) & (jec_obj.df.lep1_pt > 25) & (jec_obj.df.lep2_pt > 20) & (jec_obj.df.lep_Type > 0)]


  def yield_string(df, string):
    raw_numbers = {}
    for process in jec_obj.df.process.unique():
      j0 = 0
      j1 = 0
      j2 = 0
      for process_decay in jec_obj.df[jec_obj.df.process == process].process_decay.unique():
        if process_decay == "GluGluWWTo2L2Nu":
          continue
        #j0 += jec_obj.df[(jec_obj.df.process_decay == process_decay) & (jec_obj.df.numb_jets == 0) & (jec_obj.df.lep1_Charge != jec_obj.df.lep2_Charge)].shape[0] * 1.000000  
        #j1 += jec_obj.df[(jec_obj.df.process_decay == process_decay) & (jec_obj.df.numb_jets == 1) & (jec_obj.df.lep1_Charge != jec_obj.df.lep2_Charge)].shape[0] * 1.000000 
        #j2 += jec_obj.df[(jec_obj.df.process_decay == process_decay) & (jec_obj.df.numb_jets == 2) & (jec_obj.df.lep1_Charge != jec_obj.df.lep2_Charge)].shape[0] * 1.000000 
        j0 += jec_obj.df[(jec_obj.df.process_decay == process_decay) & (jec_obj.df.numb_jets == 0) & (jec_obj.df.lep1_Charge != jec_obj.df.lep2_Charge)].weight.sum() * scales[process_decay] 
        j1 += jec_obj.df[(jec_obj.df.process_decay == process_decay) & (jec_obj.df.numb_jets == 1) & (jec_obj.df.lep1_Charge != jec_obj.df.lep2_Charge)].weight.sum() * scales[process_decay] 
        j2 += jec_obj.df[(jec_obj.df.process_decay == process_decay) & (jec_obj.df.numb_jets == 2) & (jec_obj.df.lep1_Charge != jec_obj.df.lep2_Charge)].weight.sum() * scales[process_decay] 


      raw_numbers[process] = np.array([j0,j1,j2], dtype=float)
      string += process + "\t" +  str(j0) + "\t" + \
                                  str(j1) + "\t" + \
                                  str(j2) + "\n"  
    return string, raw_numbers

  print "Original Jets Value counts", jec_obj.df["org_numb_jets"].value_counts()
  jet_scale_shift_flat(jec_obj.df, jet_pt_shift=1.0, pt_cut=jet_pt_cut_global, rf=jec_obj.rfs)
  nominal_string, nominal_raw = yield_string(jec_obj.df, "NOMINAL\n")

  print jec_obj.df[jec_obj.df.numb_jets > 2][["jet1_pt", "jet2_pt", "jet3_pt", "numb_jets", "metMod"]].head()

  #UP
  jec_obj.reset()
  jec_obj.apply_pre_cuts()

  apply_binned_jet_uncertainty(jec_obj.df, jec_obj.JECUnc["bins"], jec_obj.JECUnc["pt_bins"], jec_obj.JECUnc["eta_bins"])  
  jet_scale_shift_official(jec_obj.df, rf=jec_obj.rfs, up_down="up")
  up_string, up_raw = yield_string(jec_obj.df, "UP\n")
  print "UP"
  print jec_obj.df[jec_obj.df.numb_jets > 2][["jet1_pt", "jet2_pt","jet3_pt", "numb_jets", "metMod"]].head()

  #DOWN
  jec_obj.reset()
  jec_obj.apply_pre_cuts()

  apply_binned_jet_uncertainty(jec_obj.df, jec_obj.JECUnc["bins"], jec_obj.JECUnc["pt_bins"], jec_obj.JECUnc["eta_bins"])  
  jet_scale_shift_official(jec_obj.df, rf=jec_obj.rfs, up_down="down")
  down_string, down_raw = yield_string(jec_obj.df, "DOWN\n")
  print "DOWN"
  print jec_obj.df[jec_obj.df.numb_jets > 2][["jet1_pt", "jet2_pt", "jet3_pt", "numb_jets", "metMod"]].head()


  bkg_raw = np.array([0.]*3)
  tot_raw = np.array([0.]*3)
  for k in nominal_raw:
    if "WW" not in k:
      bkg_raw += up_raw[k]
      tot_raw += nominal_raw[k]

  bkg_raw = bkg_raw / tot_raw


  #?true date
  #date = "17oct"
  date = datetime.date.today() 
  f = open( "results/preselection_jes_"+str(date.month)+"_"+str(date.day)+file_name+".txt", "w")
  f.write("Process \t 0 Jet \t 1 Jet \t 2 Jet\n")
  f.write(nominal_string + up_string + down_string)
  f.write( "% delta"+ "".join([ str((k, (abs(up_raw[k] - nominal_raw[k]) + abs(down_raw[k] - nominal_raw[k])) / ( 2 * nominal_raw[k]) * 100 )) + "\n" for k in nominal_raw]))
  f.write("\nfrac UP"+ "".join([ str((k, (up_raw[k] / nominal_raw[k] * 100 ))) + "\n" for k in nominal_raw]))
  f.write("\nfrac DOWN"+ "".join([ str((k, (down_raw[k] / nominal_raw[k] * 100)))+ "\n"  for k in nominal_raw]))
  f.write("\nup bkg unc " + str(bkg_raw) )
  f.write("\nfrac of WW" + str(nominal_raw["WW"] / (jec_obj.df[jec_obj.df.process_decay == "WW"].shape[0] * 1.000)))#.weight.sum() * scales["WW"])))
  f.close()




def make_plots(df, feature, prefix):
  df[feature].hist(alpha=0.5, bins=50 )

  plt.title(prefix + " " +  " ".join(feature.split("_")))
  plt.legend(["orig", "scaled"])
  plt.savefig("data/"+" ".join(prefix.split("_"))+"_"+feature+".png")
  #plt.show()
  plt.close()



def cross_section_JES( jec_obj, file_name, verbose=False, flavor=""):

  scales = jec_obj.scales
  def yield_string(df, string):
    raw_numbers = {}
    for process in jec_obj.df.process.unique():
      j0 = 0
      j1 = 0
      j2 = 0
      for process_decay in jec_obj.df[jec_obj.df.process == process].process_decay.unique():
        j0 += jec_obj.rf_ana(jec_obj.df[(jec_obj.df.process_decay == process_decay) & (jec_obj.df.numb_jets == 0) & (jec_obj.df.lep1_Charge != jec_obj.df.lep2_Charge)]).weight.sum() * scales[process_decay] 
        j1 += jec_obj.rf_ana(jec_obj.df[(jec_obj.df.process_decay == process_decay) & (jec_obj.df.numb_jets == 1) & (jec_obj.df.lep1_Charge != jec_obj.df.lep2_Charge)]).weight.sum() * scales[process_decay] 
        j2 += jec_obj.rf_ana(jec_obj.df[(jec_obj.df.process_decay == process_decay) & (jec_obj.df.numb_jets == 2) & (jec_obj.df.lep1_Charge != jec_obj.df.lep2_Charge)]).weight.sum() * scales[process_decay] 


      raw_numbers[process] = np.array([j0,j1,j2], dtype=float)
      string += process + "\t" +  str(j0) + "\t" + \
                                  str(j1) + "\t" + \
                                  str(j2) + "\n"  
    return string, raw_numbers



  # Cross-sections #
  # Nominal cross-sections and yields #
  pseudo = pseudo_data_yield_sum(rf_ana(jec_obj.df), rf_ana(jec_obj.df_da), scales=scales)
  #jec_obj.df_ww["gen_weight"] = jec_obj.df_ww.gen_weight.values/ np.abs(jec_obj.df_ww.gen_weight.values)
  #jec_obj.df_ww["weight"] = jec_obj.df_ww.weight.values * jec_obj.df_ww.gen_weight.values
  print "Test of Xs before jet scale flat", cross_calc(jec_obj.df, jec_obj.df_da, jec_obj.df_ww, jec_obj.df_ggww, scales,  fiducial=True, pseudo=pseudo), scales["WW"]


  jet_scale_shift_flat(jec_obj.df, jet_pt_shift=1.0, pt_cut=jet_pt_cut_global, rf=jec_obj.rfs)

  pseudo = create_pseudo(rf_ana(jec_obj.df), rf_ana(jec_obj.df_da), scales=scales)
  orig = create_Xs(jec_obj, pseudo) #{"tot":0, "j0":0, "j1":0, "j2":0}
  print "Post nominal jet scale shift: ", orig


  #########################
  #Fits
  nominal_fit_result = fit.comprehensive_fit(jec_obj.df, jec_obj.df_da, "metMod", scales)
  #########################


  eff_orig = jec_obj.rf_ana(jec_obj.df[jec_obj.df.process == "WW"]).shape[0] / float(jec_obj.df_ww[jec_obj.df_ww.process == "WW"].shape[0])


  # Yields #
  nominal_string = "nominal:\n"
  nominal_string, nominal_raw = yield_string(jec_obj.df, nominal_string)
  for feature in jec_obj.rfs["features_fTT"] + ["pred_fTT_WW", "pred_fDY_WW"]:
    make_plots(jec_obj.df, feature, "nominalCross")
                   
                   
  jec_obj.reset()
  jec_obj.apply_pre_cuts()
  # Up JES cross-sections and yields #
  apply_binned_jet_uncertainty(jec_obj.df, jec_obj.JECUnc["bins"], jec_obj.JECUnc["pt_bins"], jec_obj.JECUnc["eta_bins"])  
  jet_scale_shift_official(jec_obj.df, rf=jec_obj.rfs)

  jec_obj.df["pred_fTT_WW"] = jec_obj.df.fTT 
  jec_obj.df["pred_fDY_WW"] = jec_obj.df.fDY

  #Yields#
  up_string = "up:\n"
  up_string, up_raw = yield_string(jec_obj.df, up_string)
  #########################
  #Fits
  up_fit_result = fit.comprehensive_fit(jec_obj.df, jec_obj.df_da, "metMod", scales)
  #########################

  up = create_Xs(jec_obj, pseudo) #{"tot":0, "j0":0, "j1":0, "j2":0}
  print up

  eff_up = jec_obj.rf_ana(jec_obj.df[jec_obj.df.process == "WW"]).shape[0] / float(jec_obj.df_ww[jec_obj.df_ww.process == "WW"].shape[0])
  for feature in jec_obj.rfs["features_fTT"] + ["pred_fTT_WW", "pred_fDY_WW"]:
    make_plots(jec_obj.df, feature, "upCross")


  ####################################
  jec_obj.reset()
  jec_obj.apply_pre_cuts()
  # Down JES cross-sections and yields #
  apply_binned_jet_uncertainty(jec_obj.df, jec_obj.JECUnc["bins"], jec_obj.JECUnc["pt_bins"], jec_obj.JECUnc["eta_bins"])  
  jet_scale_shift_official(jec_obj.df, rf=jec_obj.rfs, up_down="down")

  jec_obj.df["pred_fTT_WW"] = jec_obj.df.fTT 
  jec_obj.df["pred_fDY_WW"] = jec_obj.df.fDY

  # Yields #
  down_string = "down:\n"
  down_string, down_raw = yield_string(jec_obj, down_string)

  down = create_Xs(jec_obj, pseudo) #{"tot":0, "j0":0, "j1":0, "j2":0}
  print down
  #########################
  #Fits
  down_fit_result = fit.comprehensive_fit(jec_obj.df, jec_obj.df_da, "metMod", scales)
  #########################

  eff_down = jec_obj.rf_ana(jec_obj.df[jec_obj.df.process == "WW"] ).shape[0] / float(jec_obj.df_ww[jec_obj.df_ww.process == "WW"].shape[0])
  for feature in jec_obj.rfs["features_fTT"] + ["pred_fTT_WW", "pred_fDY_WW"]:
    make_plots(jec_obj.df, feature, "downCross")
  ########################################
  if verbose == True:
    print "Cross section JES"
    print "Orig:", orig 
    print "up", up 
    print "down", down


    print ""
    print "UP delta  ",    [(k, (abs(up[k]   - orig[k]) ) ) for k in orig]
    print "UP % delta",    [(k, (abs(up[k]   - orig[k]) ) / (orig[k]) * 100.) for k in orig]

    print ""
    print "DOWN delta  ",  [(k, (abs(down[k] - orig[k]))) for k in orig]
    print "DOWN % delta",  [(k, (abs(down[k] - orig[k]))  / (orig[k]) * 100.) for k in orig]

    print ""
    print "Tot delta  ",  [(k, (abs(up[k] - orig[k]) + abs(down[k] - orig[k])) / 2.) for k in orig]
    print "Tot % delta",  [(k, (abs(up[k] - orig[k]) + abs(down[k] - orig[k])) / (2 * orig[k]) * 100.) for k in orig]


    print "Process \t 0 Jet \t 1 Jet \t 2 Jet"
    print nominal_string
    print up_string
    print down_string

    print "Tot Signal % delta", abs(sum([up_raw[k].sum() for k in up_raw if k == "WW"])  - sum([nominal_raw[k].sum() for k in nominal_raw if k == "WW"])) / sum([nominal_raw[k].sum() for k in nominal_raw if k == "WW"]) 
    print "Tot Bkg by bin : Up ", sum([up_raw[k] for k in up_raw if k != "WW"]),  "Nominal: ", sum([nominal_raw[k] for k in nominal_raw if k != "WW"])
    print "Tot Bkg: Up ", sum([up_raw[k].sum() for k in up_raw if k != "WW"]),  "Nominal: ", sum([nominal_raw[k].sum() for k in nominal_raw if k != "WW"])
    print "Tot Bkg % delta", abs(sum([up_raw[k].sum() for k in up_raw if k != "WW"])  - sum([nominal_raw[k].sum() for k in nominal_raw if k != "WW"])) / sum([nominal_raw[k].sum() for k in nominal_raw if k != "WW"])
  #########################
  print "Fit results"
  fit_processes = ["WW", "Top", "DY"] 
  print fit_processes
  print "nominal: ", nominal_fit_result.x
  print "up: ", up_fit_result.x
  print "down: ", down_fit_result.x
  print "Diffs(%)"
  for it, ele in enumerate(nominal_fit_result.x):
    print  fit_processes[it], (abs(ele - up_fit_result.x[it]) + abs(ele - down_fit_result.x[it]))/2. * 100.
  #########################
  

  bkg_raw = np.array([0.]*3)
  tot_raw = np.array([0.]*3)
  for k in nominal_raw:
    if "WW" not in k:
      bkg_raw += up_raw[k]
      tot_raw += nominal_raw[k]

  bkg_raw = bkg_raw / tot_raw

  #?true date
  if flavor != "":
    flavor = "_" + flavor
  date = datetime.date.today() 
  f = open( "results/mar/jes_"+str(date.year) + "_" + str(date.month) + "_" +flavor+".txt", "w")
  f.write("Cross section JES\n")
  f.write("X-section:"+str(orig) + "\n")
  f.write("Orig: "+     str([(k, (orig[k] - orig[k]) / orig[k] * 100)     for k in orig]) + "\t Efficiency: " + str(eff_orig) +"\n")
  f.write("up: "+       str([(k,abs(up[k] - orig[k]) / orig[k] * 100)   for k in orig]) + "\t Efficiency: " + str(eff_up) +"\n") 
  f.write("down: "+     str([(k,abs(down[k] - orig[k]) / orig[k] * 100) for k in orig]) + "\t Efficiency: " + str(eff_down) +"\n")
  f.write("% delta: "+  str([(k,(abs(up[k] - orig[k]) + abs(down[k] - orig[k])) / (2 * orig[k]) * 100.) for k in orig]) + "\n")
  f.write("Process \t 0 Jet \t 1 Jet \t 2 Jet\n")
  f.write(nominal_string + up_string + down_string)
  f.write( "% delta"+ "".join([ str((k, (abs(up_raw[k] - nominal_raw[k]) + abs(down_raw[k] - nominal_raw[k])) / ( 2 * nominal_raw[k]) * 100 )) + "\n" for k in nominal_raw]))
  f.write("\nfrac UP"+ "".join([ str((k, (up_raw[k] / nominal_raw[k] * 100 ))) + "\n" for k in nominal_raw]))
  f.write("\nfrac DOWN"+ "".join([ str((k, (down_raw[k] / nominal_raw[k] * 100)))+ "\n"  for k in nominal_raw]))
  f.write("\nup bkg unc " + str(bkg_raw) )
  f.write("\nfrac of WW" + str(nominal_raw["WW"] / (jec_obj.df[jec_obj.df.process_decay == "WW"].weight.sum() * scales["WW"])))
  #########################
  f.write("\n\n\nFit results\n")
  f.write(str(fit_processes) + "\n")
  f.write("nominal: "+ str(nominal_fit_result.x) + "\n")
  f.write("up: "+ str(up_fit_result.x) + "\n")
  f.write("down: "+ str(down_fit_result.x) + "\n")
  f.write("Diffs(%)" + "\n")
  for it, ele in enumerate(nominal_fit_result.x):
    f.write(str(fit_processes[it]) + " " + str((abs(ele - up_fit_result.x[it]) + abs(ele - down_fit_result.x[it]))/2. * 100.) + "\n")
  #########################
  f.close()









if __name__ == "__main__":
  import warnings
  warnings.filterwarnings('ignore')

  for flavor in ["", "same", "diff"]:
    jec_obj = analysis_setup(unc="jet", flavor=flavor)
    jec_obj.loadJECUnc()
    jec_obj.apply_pre_cuts()

    cross_section_JES( jec_obj, "data/jes_unc.txt", verbose=True, flavor=flavor)

  exit()
  jec_obj = jec_setup()
  #jec_obj.saveJECUnc() 
  jec_obj.loadJECUnc()


#  print jec_obj.rfs

  ################ 
  # Yields #
  #Run if jec or flat has not been saved
#  print "\033[31mOfficial\033[0m"
#  jec_obj.reset()
#  print "\033[31mFlat\033[0m"
#


  #jec_df = pickle.load(open("data/official_njets.pkl", "r")) 
  #flat_df = pickle.load(open("data/flat_njets.pkl", "r")) 

  #Official JEC
  #for process in jec_df.process.unique():
  #  jet_process_njets(jec_df, process, pre_rf="pre")
  #for process in jec_df.process.unique():
  #  jet_process_njets(jec_df, process, pre_rf="rf")



  #preselection_JES(jec_obj)
  cross_section_JES( jec_obj, "data/jes_unc.txt", verbose=True )

  #////////////////////////////////////////////////////////////////////////
  #////////////////////////////////////////////////////////////////////////
  def make_plots(df, feature, shift_func, prefix, function_inputs):
    df_ = copy.copy(df)

    df_[feature].hist(alpha=0.5, bins=50 )
    #jet_scale_shift_official(jec_obj.df, rf=jec_obj.rfs, up_down="down")
    shift_func(df_, up_down=function_inputs["jet_pt_shift"], rf=function_inputs["rf"])
    if  "pred" in feature:
      df_["pred_fTT_WW"] = df_.fTT
      df_["pred_fDY_WW"] = df_.fDY
    df_[feature].hist(alpha=0.5, bins=50 )

    plt.title(prefix + " " +  " ".join(feature.split("_")))
    plt.legend(["orig", "scaled"])
    plt.savefig("data/"+prefix+"_"+feature+".png")
    #plt.show()
    plt.close()
    #print feature, "orig", df[feature].min(), df[feature].max()
    #print feature, "post", df_[feature].min(), df_[feature].max()

  def make_matrix(jec_obj):
    jet_scale_shift_flat(jec_obj.df, jet_pt_shift=1.0, pt_cut=jet_pt_cut_global, rf=jec_obj.rfs)
    apply_binned_jet_uncertainty(jec_obj.df, jec_obj.JECUnc["bins"], jec_obj.JECUnc["pt_bins"], jec_obj.JECUnc["eta_bins"])  


    orig_jets_dic = {}
    up_jets_dic = {}
    down_jets_dic = {}
    orig_numb_jets = copy.copy(jec_obj.df.numb_jets.values)

    jet_scale_shift_official(jec_obj.df, up_down="up", rf=jec_obj.rfs)
    print "UP"
    for jet in sorted(np.unique(orig_numb_jets)):
      for jet_ in sorted(jec_obj.df.numb_jets.unique()):
        up_jets_dic[(jet, jet_)] = jec_obj.df[(orig_numb_jets == jet) & (jec_obj.df.numb_jets == jet_)].shape[0]
        print round(float(up_jets_dic[(jet, jet_)]) / max(jec_obj.df[(orig_numb_jets == jet)].shape[0], 1.), 4) , "\t\t\t", 
      print ""

    #?Broke need to fix
    jec_obj.reset()
    apply_binned_jet_uncertainty(jec_obj.df, jec_obj.JECUnc["bins"], jec_obj.JECUnc["pt_bins"], jec_obj.JECUnc["eta_bins"])  
    jet_scale_shift_official(jec_obj.df, up_down="down", rf=jec_obj.rfs)
    print "DOWN"
    for jet in sorted(np.unique(orig_numb_jets)):
      for jet_ in sorted(jec_obj.df.numb_jets.unique()):
        down_jets_dic[(jet, jet_)] = jec_obj.df[(orig_numb_jets == jet) & (jec_obj.df.numb_jets == jet_)].shape[0]
        print round(float(down_jets_dic[(jet, jet_)]) / max(jec_obj.df[(orig_numb_jets == jet)].shape[0], 1.), 4) , "\t\t\t", 
      print ""



