import root_pandas as rp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import copy, string 
import pickle, sys, socket

from physics_selections import *
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator
from ast import literal_eval as make_tuple
plt.rc('text', usetex=True)

computername = socket.gethostname().lower()
if computername == "Ifrit".lower():
  sys.path.append("/home/gunter/WW_analysis/production/Analysis_13TeV/scripts/uncertainties_dir")
  #Set analysis home directory
  home_dir = "/home/gunter/WW_analysis"
  post_data_dir = home_dir + "/production/Analysis_13TeV/data/"
  data_path = home_dir + "/data"
  production_path = home_dir + "/production/Analysis_13TeV/"

elif computername == "swordfish":
  print("not implemented")
  sys.path.append("/home/gunter/WW_analysis/production/Analysis_13TeV/scripts/uncertainties_dir")
  #Set analysis home directory
  home_dir = "/home/gunter/WW_analysis"
  post_data_dir = home_dir + "/production/Analysis_13TeV/data/"
  data_path = home_dir + "/data"
  production_path = home_dir + "/production/Analysis_13TeV/"


print "home", home_dir


print("Loading binning options...")
binning_options = pd.read_csv( production_path + "/binning_options.txt", index_col=None, sep=";" )


#Set energy and lumi defaults
page = "13TeV"
lumi_amount="35.9"

page_temp = "13TeV"

print("Loading plotting options")
plotting_options = pd.read_csv(production_path + "/plotting_options_"+page+".csv", index_col=None, sep=";")


#NOTE: 
#Defunct we instantiation at the top of file

#home_dir_ = "/home/gunter/WW_analysis" 
#data_path = home_dir_ + "/data_alpha"#"/data_old_alpha" 



##################
# Process Scalings
scales = {}
for key in plotting_options.process_decay.unique():
  scales[key] = eval(plotting_options[plotting_options.process_decay == key]["scale"].values[0]) 

# unc  process
unc_mc_process = {}
for key in plotting_options.process_decay.unique():
  unc_mc_process[key] = plotting_options[plotting_options.process_decay == key]["unc"].values[0]

print("unc_mc_process and scales are parameter dictionaries")


##################
# Palettes 
palettes = {}
for key in plotting_options.keys():
  if "color" in key:
    if key == "color":
      palettes["default"] = {}
    else:
      palettes[key] = {}
    for process in plotting_options.process.unique():
      if key == "color": palettes["default"][process] = eval( plotting_options[ plotting_options.process == process ][key].values[0] )
      else: palettes[key][process] = eval( plotting_options[ plotting_options.process == process ][key].values[0] )
####################      
# Features
features = ['numb_BJet', 'HT', 'lep_Type', 'numb_jets', 'numb_jets_24', 'numb_gen_jets', 'numb_gen_jets_24', 'lep1_pt',
            'jet1_pt', 'lep2_pt', 'jet2_pt', 'metMod', 'dPhiLL', 
            'METProj', 'qT', 'dPhiLLJet', 'met_phi', 'lep3_pt',
            'tot_npv', 'mll', 'METProj_sin', 'met_over_sET',
            'recoil', 'jet1_pt', 'dPhiLLMET', 'dPhiMETJet', 'mllMET']
columns = list(set(['process', 'process_decay', 'lep1_type', 'lep2_type',
                    'gen_weight',  'lep1_Charge', 'lep2_Charge', 'lep_Type', 'numbExtraLep',
                    'tot_npv', 'gen_pu', 'pileup_weight', 'pileup_up_weight', 'pileup_down_weight',
                    'pileup_weight_mf', 'pileup_up_weight_mf','pileup_down_weight_mf',
                    'mll', 'numb_jets', 'metMod', 'numb_BJet',
                    'lep1_pt', 'lep2_pt', 'lep3_pt', 'dPhiLL', 'lep1_phi', 'lep2_phi', 'lep1_eta', 'lep2_eta',
                    'jet1_csv', 'jet1_pt', 'jet2_pt', 'HT', 'numb_BJet_gen', 'jet1_phi', 'jet2_phi', 'jet1_eta', 'jet2_eta', 
                    'METProj_sin', 'met_over_sET', 'METProj', 'met_phi', 
                    'weight', 'bjet_weight', 'id_weight', 'gen_weight', 'el1_id_weight', 'el2_id_weight', 'mu_id_weight',
                    'dPhiLLMET',  'mllMET', 'qT', 'recoil', 'dPhiLLJet', 'dPhiMETJet', 'ww_pt', "metFilter_flag"] + features)) 

columns_lep_unc   = list(set(columns +['lep1_id_error_high', 'lep2_id_error_high', 'lep1_id_error_low', 'lep2_id_error_low']))
columns_jets_unc  = list(set(columns +['jet1_csv', 'jet2_csv','jet3_csv', 'jet4_csv','jet5_csv', 'jet6_csv', 
                                       'jet1_flv', 'jet2_flv','jet3_flv', 'jet4_flv','jet5_flv', 'jet6_flv', 
                                       'jet1_pt', 'jet2_pt', 'jet3_pt', 'jet4_pt','jet5_pt', 'jet6_pt','HT', 'jet1_csv',
                                       'jet1_phi', 'jet2_phi', 'jet3_phi', 'jet4_phi', 'jet5_phi', 'jet6_phi', 'jet1_eta',
                                       'jet2_eta', 'jet3_eta', 'jet4_eta', 'jet5_eta', 'jet6_eta', 'numb_BJet_gen', 'bjet_weight', 'bjet_unc'])) 
columns_lhe       = list(set(columns + ['pdf_weight'] + ["qcd_weight"+str(i) for i in range(1, 9)] + ['qcd_weight']))

def load_origMC( add_columns=None, columns=columns):
  ###################
  if add_columns != None:
    columns = list(set(columns + add_columns))

  df_list = []
  #MC
  df_dy0 = rp.read_root(data_path+"/dyjetstoll_m-50_complete.root", columns=columns)
  df_dy_m_10 = rp.read_root(data_path+"/dyjetstoll_m-10to50_complete.root", columns=columns)
  df_ww = rp.read_root(data_path+"/ww_complete.root", columns=columns)


  ### HIGHER ORDER WW PT CORRECTIONS
  import ROOT
  f = ROOT.TFile.Open("/home/gunter/WW_analysis/production/Analysis_13TeV/tools/WW_pTresummation_13-master/MyRatioWWpTHistogramAll.root", "read")
  wwpt_hist = f.Get("wwpt")
  for i in range(1000):
      weight =  wwpt_hist.GetBinContent(i)
      df_ww.weight.values[(df_ww.ww_pt >= i*0.5 ).values & (df_ww.ww_pt < i*0.5 + 0.5).values] *= weight
  #####################333

  df_list += [df_dy0, df_dy_m_10, df_ww]

  #TTbar and single Top
  df_tt_l = rp.read_root(data_path+"/ttbar_leptonic_complete.root", columns=columns)
  df_tt_sl = rp.read_root(data_path+"/ttbar_semileptonic_complete.root", columns=columns)
  df_tbar_tw = rp.read_root(data_path+"/tbar_tw-_complete.root", columns=columns)
  #df_tbar_s = rp.read_root(data_path+"/tbar_s-_complete.root", columns=columns)
  df_tbar_t = rp.read_root(data_path+"/tbar_t-_complete.root", columns=columns)
  df_t_tw = rp.read_root(data_path+"/t_tw-_complete.root", columns=columns)
  df_t_s = rp.read_root(data_path+"/t_s-_complete.root", columns=columns)
  df_t_t = rp.read_root(data_path+"/t_t-_complete.root", columns=columns)

  df_list += [df_tt_l, df_tt_sl, df_tbar_tw, df_tbar_t, df_t_tw, df_t_s, df_t_t]

  #ZZ and WZ
  df_zz_lq = rp.read_root(data_path+"/zzjetsto2l2q_complete.root", columns=columns)
  df_zz_ln = rp.read_root(data_path+"/zzjetsto2l2nu_complete.root", columns=columns)
  df_wz_lq = rp.read_root(data_path+"/wzjetsto2l2q_complete.root", columns=columns)
  df_wz_ln = rp.read_root(data_path+"/wzjetsto3lnu_complete.root", columns=columns)

  df_list += [df_zz_lq, df_zz_ln, df_wz_lq, df_wz_ln]

  #WJets
  df_wj1 = rp.read_root(data_path+"/w1jetstolnu_complete.root", columns=columns)
  df_wj2 = rp.read_root(data_path+"/w2jetstolnu_complete.root", columns=columns)
  df_wj3 = rp.read_root(data_path+"/w3jetstolnu_complete.root", columns=columns)
  df_wj4 = rp.read_root(data_path+"/w4jetstolnu_complete.root", columns=columns)

  df_list += [df_wj1, df_wj2, df_wj3, df_wj4]

  #Other stuff 
  df_wg = rp.read_root(data_path+"/wg_complete.root", columns=columns)
  df_ggh = rp.read_root(data_path+"/glugluh_complete.root", columns=columns)  
  df_ggww = rp.read_root(data_path+"/glugluww_complete.root", columns=columns)  

  df_list += [df_wg, df_ggh, df_ggww]

  #Some smaller higgs stuff
  df_ggh_tautau = rp.read_root(data_path+"/glugluhtotautau_complete.root", columns=columns)
  df_hzj_ww = rp.read_root(data_path+"/hzj_htowwto2l2nu_complete.root", columns=columns)
  df_vbfh_ww = rp.read_root(data_path+"/vbfhtowwto2l2nu_complete.root", columns=columns)
  df_hwplus = rp.read_root(data_path+"/hwplusj_htoww_complete.root", columns=columns)
  df_hwminus = rp.read_root(data_path+"/hwminusj_htoww_complete.root", columns=columns)

  df_list += [df_ggh_tautau, df_hzj_ww, df_vbfh_ww, df_hwplus, df_hwminus]


  #Finishing up
  for i in df_list:
    if "gen_weight" not in i.keys():
      print i.head()
  df = pd.concat(df_list)

  df["gen_weight"] = df.gen_weight.values/ np.abs(df.gen_weight.values)
  df["weight"] = df.weight.values * df.gen_weight.values

  df = df[df.lep2_pt > 20]
  df = df.reset_index()

  return df



def load_origDATA(data_path=data_path, columns=columns):

  df_da_b = rp.read_root(data_path+"/muon_2016B_complete.root", columns=columns+["runNum", "eventNumb"])
  df_da_c = rp.read_root(data_path+"/muon_2016C_complete.root", columns=columns+["runNum", "eventNumb"])
  df_da_d = rp.read_root(data_path+"/muon_2016D_complete.root", columns=columns+["runNum", "eventNumb"])
  df_da_e = rp.read_root(data_path+"/muon_2016E_complete.root", columns=columns+["runNum", "eventNumb"])
  df_da_f = rp.read_root(data_path+"/muon_2016F_complete.root", columns=columns+["runNum", "eventNumb"])
  df_da_g = rp.read_root(data_path+"/muon_2016G_complete.root", columns=columns+["runNum", "eventNumb"])
  df_da_h = rp.read_root(data_path+"/muon_2016H_complete.root", columns=columns+["runNum", "eventNumb"])

  df_da = pd.concat([df_da_b, df_da_c, df_da_d, df_da_e, df_da_f, df_da_g, df_da_h])

  df_da_b_e = rp.read_root(data_path+"/electron_2016B_complete.root", columns=columns+["runNum", "eventNumb"])
  df_da_c_e = rp.read_root(data_path+"/electron_2016C_complete.root", columns=columns+["runNum", "eventNumb"])
  df_da_d_e = rp.read_root(data_path+"/electron_2016D_complete.root", columns=columns+["runNum", "eventNumb"])
  df_da_e_e = rp.read_root(data_path+"/electron_2016E_complete.root", columns=columns+["runNum", "eventNumb"])
  df_da_f_e = rp.read_root(data_path+"/electron_2016F_complete.root", columns=columns+["runNum", "eventNumb"])
  df_da_g_e = rp.read_root(data_path+"/electron_2016G_complete.root", columns=columns+["runNum", "eventNumb"])
  df_da_h_e = rp.read_root(data_path+"/electron_2016H_complete.root", columns=columns+["runNum", "eventNumb"])

  df_da_e = pd.concat([df_da_b_e, df_da_c_e, df_da_d_e, df_da_e_e, df_da_f_e, df_da_g_e, df_da_h_e])

  df_da = pd.concat( [df_da, df_da_e] ).reset_index()
  df_da = df_da.drop(np.where(df_da[["runNum", "eventNumb"]].duplicated())[0])
  df_da = df_da[df_da.lep2_pt > 25]
  return df_da


def load_preselMC(unc=""):
  if unc != "":
    unc += "_"
  return pd.read_hdf( post_data_dir + unc + "preselMC.hdf", "table")

def load_preselDATA(unc=""):
  if unc != "":
    unc += "_"
  return pd.read_hdf( post_data_dir + unc + "preselDATA.hdf", "table")

def load_anti_preselMC():
  return pd.read_hdf( post_data_dir + "anti_preselMC.hdf", "table")

def load_anti_preselDATA():
  return pd.read_hdf( post_data_dir + "anti_preselDATA.hdf", "table")

def load_dy_trainset():
  return pd.read_hdf( post_data_dir + "train_dy.hdf","table") 

def load_tt_trainset():
  return pd.read_hdf( post_data_dir + "train_top.hdf","table")

def load_testset(unc=""):
  scales = pickle.load(open(post_data_dir + "scales_test.pkl", "r"))
  if unc != "":
    unc += "_"
  return scales, pd.read_hdf( post_data_dir + unc + "test.hdf","table")


def load_presel_w_fDY_fTT_MC(unc=""):
  if unc != "":
    unc += "_"
  return pd.read_hdf( post_data_dir + unc + "preselMC_rf.hdf", "table")

def load_presel_w_fDY_fTT_DATA(unc=""):
  if unc != "":
    unc += "_"
  return pd.read_hdf( post_data_dir + unc + "preselDATA_rf.hdf", "table")

def load_rfselMC():
  print "Coming Soon"

def load_rfselDATA():
  print "Coming Soon"

def load_randomForest():
  from sklearn.externals import joblib 
  return joblib.load( post_data_dir + "RF/mar_fDY_fTT.jbl")


def create_pseudoDATA(columns=None):
  df = load_origMC()
  df = df.sample(frac=0.99, replace=True)
  return df

def loose_pre_cuts( df ): 

  dif_lep = df.lep_Type > 0
  sam_lep = df.lep_Type < 0
  z_mass = (df.mll < 81) | (df.mll > 101 ) 
  nBjet = df.numb_BJet == 0
  extra_lep = df.numbExtraLep == 0
  quality_cuts = df.mll > 28 

  lep_pt = df.lep2_pt > 20

  s_df = df[ sam_lep & z_mass & nBjet & extra_lep & quality_cuts & lep_pt]
  d_df = df[ dif_lep & nBjet & extra_lep & quality_cuts & lep_pt]

  return pd.concat([ s_df, d_df ])


def write_preselMC(pre_unc="", columns=columns):
  if pre_unc == "lhe":
    columns = columns_lhe
  elif pre_unc == "jet":
    columns = columns_jets_unc  
  elif pre_unc == "lep":
    columns = columns_lep_unc
  print columns 
  if pre_unc != "":
    pre_unc += "_"
  df = load_origMC( columns=columns )
  loose_pre_cuts(df).to_hdf( post_data_dir+pre_unc+"preselMC.hdf", 'table', complevel=3)

def write_preselDATA(pre_unc="", columns=columns):
  print columns 
  if pre_unc == "lhe":
    columns = columns_lhe
  elif pre_unc == "jet":
    columns = columns_jets_unc  
  elif pre_unc == "lep":
    columns = columns_lep_unc
  if pre_unc != "":
    pre_unc += "_"
  df = load_origDATA( columns=columns )
  loose_pre_cuts(df).to_hdf( post_data_dir+pre_unc+"preselDATA.hdf", 'table', complevel=3)

def write_preselRF(pre_unc="", data_mc=""):

  if pre_unc == "lhe":
    columns = columns_lhe
  elif pre_unc == "jet":
    columns = columns_jets_unc  
  elif pre_unc == "lep":
    columns = columns_lep_unc

  if data_mc == "MC":
    df = load_preselMC(pre_unc)
  else:
    df = load_preselDATA(pre_unc)
  from sklearn.externals import joblib

  random_forests = load_randomForest()
  features_fDY = random_forests["features_fDY"]
  clf_fDY = random_forests["clf_fDY"]
  features_fTT = random_forests["features_fTT"]
  clf_fTT = random_forests["clf_fTT"]

  #Predict MC
  temp = df[features_fTT] 
  temp = temp.replace([np.inf, -np.inf, np.nan], 0)
  pred_fTT = clf_fTT.predict_proba(np.float32(temp.values))
  df["pred_fTT_WW"] = pred_fTT[:,0]

  temp = df[features_fDY] 
  temp = temp.replace([np.inf, -np.inf, np.nan], 0)
  pred_fDY = clf_fDY.predict_proba(np.float32(temp.values))
  df["pred_fDY_WW"] = pred_fDY[:,0]

  if pre_unc != "":
    pre_unc += "_"
  df.to_hdf(post_data_dir+pre_unc+"presel"+data_mc+"_rf.hdf", "table", complevel=3)


def write_preselRFMC(pre_unc=""):
  write_preselRF(pre_unc, data_mc="MC")


def write_preselRFDATA(pre_unc=""):
  write_preselRF(pre_unc, data_mc="DATA")



class analysis_setup():
  """
  analysis_setup:
  This object loads the data, mc, and jet energy correction
  module and data files neded to complete the jet scale measurement. 
   
  """
  def __init__(self, unc="jet", flavor=""):
    self.path = "/home/gunter/WW_analysis/production/Analysis_13TeV/scripts/uncertainties_dir/"
    self.flavor = flavor
    self.unc_type = unc

    self.scales, self.df    = load_testset(unc=unc)
    self.df_da = load_presel_w_fDY_fTT_DATA(unc)
    print process_yields(self.df, self.df_da, scales=self.scales)
    self.columns = columns_jets_unc

    self.df = self.df[self.df.metFilter_flag == 0]
    self.df_da = self.df_da[self.df_da.metFilter_flag == 0]
    self.df_da["weight"] = np.array([1.0] * self.df_da.shape[0]) 

    if unc == "lhe":
      self.columns = columns_lhe
    self.df_ww = rp.read_root(data_path+"/ww_complete.root", columns=self.columns)
    self.df_ggww = rp.read_root(data_path+"/glugluww_complete.root", columns=self.columns)  

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


  def reset(self):
    self.__init__(self.unc_type, self.flavor)

  def apply_pre_cuts(self, diff_charge=False):
    self.df = pre_cuts(self.df, diff_charge=diff_charge)
    self.df_da = pre_cuts(self.df_da, diff_charge=diff_charge)

  def apply_flat_jet_correction(self):
    jet_scale_shift_flat(self.df, jet_pt_shift=1.0, pt_cut=30, rf=self.rfs)
    jet_scale_shift_flat(self.df_da, jet_pt_shift=1.0, pt_cut=30, rf=self.rfs)

  def rf_ana(self, df):
    return df[(df.pred_fDY_WW > .9)&(df.pred_fTT_WW > .6) ]

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



def kill_jets( df, pt_cut=30 ):
  """
  Edit number of jets in event based on pt cut
  """
  #?Not quite what I should be doing..
  n_jets = np.zeros(df.shape[0]) #np.maximum(np.zeros(df.shape[0]), df.numb_jets.values - 6)
  for k in df.keys():
    if "jet" in k and "pt" in k:
      cut = (df[k] > pt_cut)
      n_jets[cut.values] = n_jets[cut.values] + 1 

  df["numb_jets"] = n_jets


def jet_scale_shift_flat(data, jet_pt_shift=1., pt_cut=30, rf=None):
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
    if "jet" in k and "pt" in k:

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
    temp = data[rf["features_fTT"]]
    temp = temp.replace([np.inf,-np.inf, np.nan], 0)
    pred_fTT = rf["clf_fTT"].predict_proba(np.float32(temp.values))
    data["pred_fTT_WW"] = pred_fTT[:,0]

    temp = data[rf["features_fDY"]]
    temp = temp.replace([np.inf,-np.inf, np.nan], 0)
    pred_fDY = rf["clf_fDY"].predict_proba(np.float32(temp.values))
    data["pred_fDY_WW"] = pred_fDY[:,0]










#########End Of What I Want to Include in This File##########


################## Extra rando shit
# set-up tables and 

def create_table( data, round_digit ):
    for flavor in data.keys():
        print flavor
        print "\t",{ process : round(data[flavor][process], round_digit) for process in data[flavor].keys() }

def combine_unc( data ):
    comb_unc = {}
    for process in data[data.keys()[0]].keys():
        comb_unc[process] = 0
    for flavor in data.keys():
        for process in data[flavor].keys():
            comb_unc[process] += pow( data[flavor][process], 2 )
    print { process : round( pow( comb_unc[process], .5) , 2) for process in comb_unc.keys()}


###############################################
#### Create kinematic histograms
def make_control_plots( df_mc, df_data, date_tag, selection_tag, energy_dir= "13TeV", control_regions= None, scales=scales ):
  """
  Creates all control region plots.
  """
  def full_ana( df, diff_charge=True ):
    return df

  def same_ana( df, diff_charge= True ):
    return df[df.lep_Type < 0]

  def diff_ana( df, diff_charge= True ):
    return df[ df.lep_Type > 0]

  if control_regions == None:
    #"Z_tt": Z_tt_ana
    control_regions = {"WW": WW_ana, "TT": TT_ana, "DY": DY_ana, "full": full_ana, "same": same_ana, "diff": diff_ana, "Z_tt": Z_tt_ana}  
  else:
    control_regions["same"] = same_ana
    control_regions["diff"] = diff_ana


  for key in control_regions.keys():
    print key
    create_kinematic_hist( control_regions[key](df_mc, diff_charge= False), control_regions[key](df_data, diff_charge= False), prefix= energy_dir + "/" + selection_tag + "/" + key + "/" + date_tag , scales=scales)

def create_kinematic_hist(df_mc, df_data, prefix="", scales=scales):
  """
  Creates all the basic histograms you'll ever need:
  create_kinematic_hist(df)
  """

  range = (0,250)
  bins  = 100

  features = [ 'numb_BJet', 'HT', 'lep_Type',
       'numb_jets', 'lep1_pt',
       'jet1_pt', 'lep2_pt',
       'jet2_pt', 'metMod', 'dPhiLL', 'METProj',
       'qT', 'dPhiLLJet', 'met_phi',
       'lep3_pt', 'tot_npv',
       'mll', 'METProj_sin', 'met_over_sET',
        'recoil', 'jet1_pt', 'dPhiLLMET', 'dPhiMETJet', 'mllMET']
  if "pred_fDY_WW" in df_mc.keys() and "pred_fDY_WW" in df_data.keys():
    features.append("pred_fDY_WW")
  if "pred_fTT_WW" in df_mc.keys() and "pred_fTT_WW" in df_data.keys():
    features.append("pred_fTT_WW")

  for feature in features:
    if feature not in df_mc.keys(): continue
    print feature, df_mc.shape
    a, b, figs, ax = full_bin_plot(df_mc, df_data, feature, scales=scales,)
    figs.savefig(production_path + '/plots/' + prefix + "_"+ feature + "_prod.png")
    figs.savefig(production_path + '/plots/' + prefix + "_"+ feature + "_prod.pdf")
  print "\nFeature plots done.", prefix
  return
##############################################

def two_tree_process_map( df, pred_names, bins=10, scales=scales):
    bins_i = bins
    bins_j = bins
    if type(bins) == tuple:
        bins_i = bins[0]
        bins_j = bins[1]
    results = {}

    for decay in df.process_decay.unique():
        ax_i = np.array([ float(i) / float(bins_i) for i in xrange(bins_i)])
        ax_j = np.array([ float(j) / float(bins_j) for j in xrange(bins_j)])
        a= df[df.process_decay == decay][pred_names[0]].values.reshape( (df[df.process_decay == decay].shape[0], 1) ) > ax_i
        b= df[df.process_decay == decay][pred_names[1]].values.reshape( (df[df.process_decay == decay].shape[0], 1) ) > ax_j
        ones = np.ones((ax_i.shape[0],df[df.process_decay == decay].shape[0],1), dtype=np.bool)
        a_ = ones == a
        results_ = a_ & (b.transpose().reshape((ax_i.shape[0],df[df.process_decay == decay].shape[0],1)) == np.ones((df[df.process_decay == decay].shape[0],ax_i.shape[0])))
        if decay in scales.keys():
            if process in results.keys(): 
                print "if ", process,  decay
                results[decay] += results_.sum(axis=1) * scales[decay]
            else:
                print decay
                results[decay] = results_.sum(axis=1) * scales[decay]
    return results, [ax_i,ax_j]

def yield_asymetry( process_map, df):
    results = {}
    process_names = ["WW", "DY", "Top"]
    for process in process_names:
        for decay in process_map[0].keys():
            if process in df[ df.process_decay == decay].process.unique():
                if process not in results:
                    results[process] = process_map[0][decay]
                else:
                    results[process] += process_map[0][decay]
    #results = (process_map[0]["WW"] - (process_map[0]["DY"] + process_map[0]["Top"])) / (process_map[0]["WW"] + process_map[0]["DY"] + process_map[0]["Top"])
    results_Numerator = results["WW"] - ( results["DY"] + results["Top"])
    results_Denominator = results["WW"] + results["DY"] + results["Top"]
    return results_Numerator / results_Denominator, process_map[1]

def calc_norm_unc( data ):
    coeff = 1./(19.7e3*.122*(3*.108)**2*(data["WW"] /  np.max(data["WW"])))
    norm = [ unc_mc_process[process]*(data[process])**0.5 for process in data.keys() if "WW" not in process]
    sum_norm = np.zeros(data[process].shape)
    for ele in norm:
        sum_norm += ele**2
    return coeff*(sum_norm + unc_mc_process["WW"]/data["WW"])**.5

def calc_stat_unc( data ):
    coeff = 1./(19.7e3*.122*(3*.108)**2*(data["WW"] /  np.max(data["WW"])))
    norm = [ scales[process]*(data[process]/scales[process])**0.5 for process in data.keys() if "WW" not in process]
    
    sum_norm = np.zeros(data[process].shape)
    for ele in norm:
        sum_norm += ele**2 
    return coeff*(sum_norm + scales["WW"]/data["WW"] )**.5

def full_stat( data ):
    coeff = 1./(19.7e3*.122*(3*.108)**2*(data["WW"] /  np.max(data["WW"])))
    stat = np.zeros(data["WW"].shape)
    for i in data.keys():
        stat += data[i]        
    return coeff*(stat)**.5

def unc_map( process_map ):

    ax_i = process_map[1][0]
    ax_j = process_map[1][1]
    
    unc_sum = np.power(calc_norm_unc( process_map[0] )**2 + calc_stat_unc( process_map[0] )**2,.5 )#+ full_stat( process_map[0] )**2,.5)
    print unc_sum.min()
    unc_sum = unc_sum / unc_sum.min()
    return unc_sum, [ax_i,ax_j]



###############################################
#### Histogramming and Binning
def bin_df( df, binned_feature, binning_options=binning_options, plotting_options=plotting_options, scales=scales, range=None, bins=None, lumi_scale=1, density=False, weights=True, weights_arr_=None):
  """
  bin_df( df, binned_feature, binning_options=binning_options, plotting_options=plotting_options scales=None, range=None, bins=None, lumi_scale=1, density=False)
  """
  binned_results = {}
  defaults = ["pt", "met", "eta", "phi", "numb", "dphi"]
  bins_default = binning_options[ binning_options.feature == "default"].binning.values[0]
  range_default = make_tuple(binning_options[ binning_options.feature == "default"].range.values[0])
  y_label_default =  binning_options[ binning_options.feature == "default"].y_label.values[0]
  title_default = binned_feature

  count_defaults = 0
  for default in defaults:
    if default in binned_feature.lower():
      count_defaults += 1
      range_default = make_tuple(binning_options[ binning_options.feature == "default_"+default ].range.values[0])
      bins_default = binning_options[ binning_options.feature == "default_"+default ].binning.values[0]
      y_label_default = binning_options[ binning_options.feature == "default_"+default ].y_label.values[0]
      title_default = binning_options[ binning_options.feature == "default_"+default ].title.values[0]
  if count_defaults > 1:
    if "_" in title_default:
      title_default = " ".join(binned_feature.split("_"))
      y_label_default = "Entries"



  if binned_feature in binning_options.feature.values:
      bins_ = binning_options[ binning_options.feature == binned_feature ].binning.values[0]
      range_ = make_tuple( binning_options[ binning_options.feature == binned_feature ].range.values[0] )
      y_label = binning_options[ binning_options.feature == binned_feature ].y_label.values[0]
      title = binning_options[ binning_options.feature == binned_feature].title.values[0]

  else:
    bins_   = bins_default 
    range_  = range_default 
    y_label = y_label_default 
    title = title_default

  if bins == None:
    bins = bins_
  if range == None:
    range = range_

  if "???" in title:
    title  = string.replace(title, "???", binned_feature) 
    if "_" in title:
      title = " ".join(title.split("_"))

  binned_results["plotting"] = {"y_label": y_label, "title": title}


  unique_df_processes = df.process_decay.unique()

  for process in plotting_options.process_decay.unique():
    if process in unique_df_processes:
      df_process = df[df.process_decay == process]
      #print process, scales[process], range, bins
      if weights == True:
        if type(weights_arr_) == type(None):
          weights_arr = df_process.weight.values
        else:
          #print (df.process_decay == process).shape, weights_arr.shape
          weights_arr = weights_arr_[(df.process_decay == process).values]
        sq_weights_arr = weights_arr**2
      else:
        weights_arr = None
        sq_weights_arr = None
      binned_results[process] = list( np.histogram( df_process[binned_feature], bins=bins, range=range, weights=weights_arr, density=density ) )
      binned_results[process][0] = binned_results[process][0] * lumi_scale * scales[process]
      binned_results[process].append( (binned_results[process][1][1:]  - binned_results[process][1][:-1]) / 2. + binned_results[process][1][:-1] )
      binned_results[process].append( np.histogram( df_process[binned_feature], bins=bins, range=range, weights=sq_weights_arr, density=density )[0] )
      binned_results[process][3] = binned_results[process][3] * lumi_scale**2. * scales[process]**2.

  return binned_results



def plot_hist( bins, plotting_options=plotting_options, processes=None, x_range=None, y_range=None, title=None, y_label=None, colors=None, logy=True, x_minor_ticks=True, lumi_amount="19.7", ax=None):
  """
  Histogramming stuffs
  plot_hist( bins, processes=[ "WW", "TT", "WZ", "ZZ", "DY"], x_range=None, y_range=None, title=None, y_label=None, color=colors, logy=True, x_minor_ticks=True)
  """
  if ax == None: fig, ax = plt.subplots(figsize=(11, 9))
  if colors == None:
    colors = palettes["default"]
  elif type(colors) == str:
    colors = palettes[colors]

  if "plotting" in bins:
    if y_label == None:
      y_label = bins["plotting"]["y_label"]    
    if title == None:
      title = bins["plotting"]["title"] 
    
  if processes == None:
    processes = plotting_options.process.unique()
#  if "_" in title:
#    title = " ".join(title.split("_"))

  minorLocator = AutoMinorLocator()

  tot_bins = {}
  sum_yerr = np.zeros( len( bins[ bins.keys()[0] ][3] ) )
  for process in processes:#plotting_options.process.unique():
    for process_decay in plotting_options[ plotting_options.process == process ].process_decay.unique():
      if process_decay in bins.keys():
        sum_yerr += bins[process_decay][3]
        if process not in tot_bins.keys():
          tot_bins[process] = copy.deepcopy(bins[process_decay])
        else: 
          tot_bins[process][0] += bins[process_decay][0]
          

#Plotting rects
  rect = []
  sum_bins = np.zeros( len( tot_bins[ tot_bins.keys()[0] ][0] ) )
  last_color = None
  for process in processes:
#########
    if process in tot_bins.keys() and process in colors.keys():
########
      bottom = sum_bins
      if int(matplotlib.__version__.split(".")[0]) == 1:  
        rect.append(ax.bar( tot_bins[process][1][:-1], tot_bins[process][0],
                     tot_bins[process][1][1] - tot_bins[process][1][0] , color = colors[process],
                     edgecolor = colors[process], bottom=bottom ))
      if int(matplotlib.__version__.split(".")[0]) >= 2:  
        rect.append(ax.bar( tot_bins[process][2], tot_bins[process][0],
                      tot_bins[process][1][1] - tot_bins[process][1][0] , color = colors[process],
                      edgecolor = colors[process], bottom=bottom ))
      sum_bins +=tot_bins[process][0]
      last_color = colors[process]


  #Yerror
  process_ = tot_bins.keys()[0]
  sum_yerr = np.sqrt(sum_yerr)
  for i, yerr in enumerate(sum_yerr): 
    ax.fill( [tot_bins[process_][1][i], tot_bins[process_][1][i+1], tot_bins[process_][1][i+1], tot_bins[process_][1][i] ],\
              [sum_bins[i] - yerr, sum_bins[i] - yerr, sum_bins[i] + yerr, sum_bins[i] + yerr], fill=False, hatch='//', edgecolor='0.45' )


  #Configurables
  if logy == True: 
    if sum_bins.sum() > 0: 
      ax.set_yscale("log", nonposy='clip')
  if type(x_range)==tuple: ax.set_xlim(x_range)
  if y_range==None: 
    bottom = 1 

    bottoms = [100, 10] 
    for _bottom in bottoms:
      if _bottom * 5 < sum_bins[~np.isinf(sum_bins)].min():
        continue
      else:
        if bottom < _bottom:
          bottom = _bottom

    if logy == True: ax.set_ylim( bottom=bottom,  top= sum_bins[~np.isinf(sum_bins)].max()*5.)
    else: ax.set_ylim( bottom=0,  top= sum_bins.max()*2.)
  elif type(y_range)==tuple: ax.set_ylim( y_range )


  ax.xaxis.labelpad = 20
  ax.yaxis.labelpad = 15


  if y_label != None:
      ax.set_ylabel(y_label, fontsize=22, fontname='Bitstream Vera Sans', )
  if title != None:
      plt.xlabel( title, fontname='Bitstream Vera Sans', fontsize=24)#position=(1., 0.), va='bottom', ha='right',)

  #plt.rc('text', usetex=True)
  page_ = [page, ""]
  if page == "13TeV":
    page_ = ["13", "TeV"]
  else:
    page_ = ["8", "TeV"]
  ax.set_title(r"\textbf{CMS} Work in Progress \hspace{8cm} $"+ lumi_amount +" fb^{-1}$ $\sqrt{s}="+page_[0]+" \mathrm{"+page_[1]+"}$", fontname='Bitstream Vera Sans', fontsize=24)

  ####################################
  #Add minor tick marks to the x-axis
  if x_minor_ticks == False:
      loc = matplotlib.ticker.MultipleLocator(base=1) # this locator puts ticks at regular intervals
      ax.xaxis.set_major_locator(loc)
  else:
      ax.xaxis.set_minor_locator(minorLocator)
  
###################################
  ax.yaxis.set_tick_params(length=10, labelsize=22)
  ax.yaxis.set_tick_params(which='minor',length=5)
  ax.xaxis.set_tick_params(length=10, labelsize=22)
  ax.xaxis.set_tick_params(which='minor',length=5)

  ax.yaxis.grid(color='gray', linestyle='dashed')
  
  plt.xticks()
  plt.tight_layout()


  processes_return = [ process for process in processes if process in tot_bins.keys() ]
  return ax, rect, processes_return
   
def plot_errbar( bins, process='Da',label="data",  ax=None ):
  """
  Plot errbar
  plot_errbar( bins, process='Da' )
  bins: dictionary of list containing bin contents, edges and middles
  """
  if ax == None: plotline, caplines, barlinecols = plt.errorbar( bins[process][2], bins[process][0], yerr=np.sqrt(bins[process][0]), ecolor='black',color="black",fmt="o", label=label )
  else: plotline, caplines, barlinecols = ax.errorbar( bins[process][2], bins[process][0], yerr=np.sqrt(bins[process][0]), ecolor='black',color="black",fmt="o", label=label )

  return plotline, caplines, barlinecols 

def plot_ratio( bins_1, bins_2, y_label=None, x_label=None, ax=None):
  """
  Plot ratio plots
  plot_ratio( bins_1, bins_2, y_label=None, x_label=None, ax=None):
  bins_(1/2): list of three numpy arrays [ bins contents, bin edges, bin centers ]
  """
  process_list_1 = [ k for k in  bins_1.keys() if type(bins_1[k]) == list ] 
  tot_1 = np.zeros( bins_1[ process_list_1[0] ][0].shape[0] )
  unc_1 = np.zeros( bins_1[ process_list_1[0] ][0].shape[0] )
  for process in process_list_1:#bins_1.keys():
    tot_1 += bins_1[process][0]
    unc_1 += bins_1[process][3]

  key = ''
  process_list_2 = [ k for k in  bins_2.keys() if type(bins_2[k]) == list ] 
  tot_2 = np.zeros( bins_2[ process_list_2[0] ][0].shape[0] )
  unc_2 = np.zeros( bins_2[ process_list_2[0] ][0].shape[0] )
  for process in process_list_2:#bins_2.keys():
      tot_2 += bins_2[process][0]
      unc_2 += bins_2[process][3]
      key = process


  tot_1_OVER_tot_2 = tot_1 / tot_2 
  yerr = tot_1_OVER_tot_2 * np.sqrt( unc_1 / tot_1**2. +  unc_2 / tot_2**2.  )
  if ax==None : plotline, caplines, barlinecols =  plt.errorbar( bins_2[key][2], tot_1_OVER_tot_2, yerr= yerr, ecolor='black',color="black",fmt="o" )
  else: 
    plotline, caplines, barlinecols =  ax.errorbar( bins_2[key][2], tot_1_OVER_tot_2, yerr= yerr, ecolor='black',color="black",fmt="o" )
    ax.set_ylim( bottom=0.25,  top= 1.75)
    ax.set_ylabel('DATA/MC', fontname='Bitstream Vera Sans', fontsize=20)
    ax.locator_params(axis='y', nbins=4)
    plt.tight_layout()

  ax.xaxis.set_tick_params(labelsize=22)
  ax.yaxis.set_tick_params(labelsize=22)
  ax.yaxis.grid(color='gray', linestyle='dashed')
  return plotline, caplines, barlinecols




def full_plot(bins_mc, bins_data,  processes=[ "WW","Higgs", "WG", "WJ", "Top", "WZ", "ZZ", "DY"], x_range=None, y_range=None, title=None, y_label=None, color=None, logy=True, x_minor_ticks=True, lumi_amount=lumi_amount, ax=None):


  if type(color) == str:
    color = palettes[color]
  elif type(color) != dict:
    color = palettes["default"]
  if ax==None: fig, ax = plt.subplots(2, figsize=(11,8), sharex=True, gridspec_kw = {'height_ratios':[3, 1]})

  hist_stuff = plot_hist( bins_mc, processes=processes, x_range=x_range,\
                          y_range=y_range, title=title, y_label=y_label,\
                          colors=color, logy=logy, x_minor_ticks=x_minor_ticks,\
                          ax=ax[0], lumi_amount=lumi_amount)

  err_stuff  = plot_errbar( bins_data, ax=ax[0])

  if page == "8TeV": 
    if len( [i for i, ele in enumerate(hist_stuff[2]) if ele == "WG"]) > 0:
        hist_stuff[2][[i for i, ele in enumerate(hist_stuff[2]) if ele == "WG"][0]] = "WG(*)"

  ax[0].legend(hist_stuff[1] + [err_stuff], hist_stuff[2] + ['Da'], frameon=False, fontsize="x-large", numpoints=1)

  plot_ratio( bins_data, bins_mc, ax=ax[1])

  fig.subplots_adjust(hspace=0.1)

  return fig, ax

def plot_two_rf( map_arr ):
  fig, ax = plt.subplots(figsize=(11,9))
  pcolor_plot = ax.pcolor(map_arr[0])
  plt.colorbar(pcolor_plot)
  plt.xticks([i for i in range( len(map_arr[1][0])) if i%10==0], [i for e, i in enumerate(map_arr[1][0]) if e%10==0])
  plt.yticks([i for i in range( len(map_arr[1][0])) if i%10==0], [i for e, i in enumerate(map_arr[1][0]) if e%10==0])
  plt.xlabel("TT RF")
  plt.ylabel("DY RF")

def process_yields( df, df_da=None, query=None, processes= ['WW', 'DY', 'Top', 'WZ', 'ZZ', 'WG', 'WJ', 'Higgs'], scales=scales ):
  """
  returns a dataframe
  """

  if type(query) == type(""):
    df = df.query(query)
    if type(df_da) == type(df):
      df_da = df_da.query(query)

  tot_same = 0
  tot_diff = 0

  ds_ov_ss = 2.05
  
  yield_dic = {}
  yield_dic["Both Flavor"] = []
  yield_dic["Diff Flavor"] = []
  yield_dic["Same Flavor"] = []
  yield_dic["Process"] = []
  
  for process in processes:
    if "WJ" in process and type(df_da) != type(None): 
      continue
    if process not in df.process.unique(): 
      print "Skipping", process, "not in data frame"#, processes, df.process.unique()
      continue
    sum_process_same = 0
    sum_process_diff = 0
    for decay in df[df.process == process].process_decay.unique():
      if decay not in scales.keys():
        print "Skipping", process, "process ", decay, " decay  not in scale keys"
        continue
      if process == "WW":
        sum_process_same_ = df[(df.process_decay==decay) & (df.lep_Type < 0) & (df.lep1_Charge != df.lep2_Charge)].weight.sum() * scales[decay]
        sum_process_diff_ = df[(df.process_decay==decay) & (df.lep_Type > 0) & (df.lep1_Charge != df.lep2_Charge)].weight.sum() * scales[decay]
        yield_dic["Same Flavor"].append( int(round(sum_process_same_)) )
        yield_dic["Diff Flavor"].append( int(round(sum_process_diff_)) )
        yield_dic["Both Flavor"].append( int(round(sum_process_diff_ + sum_process_same_)) )
        yield_dic["Process"].append( decay )

      sum_process_same += df[(df.process_decay==decay) & (df.lep_Type < 0) & (df.lep1_Charge != df.lep2_Charge)].weight.sum() * scales[decay]
      sum_process_diff += df[(df.process_decay==decay) & (df.lep_Type > 0) & (df.lep1_Charge != df.lep2_Charge)].weight.sum() * scales[decay]

    yield_dic["Same Flavor"].append( int(round(sum_process_same)) )
    yield_dic["Diff Flavor"].append( int(round(sum_process_diff)) )
    yield_dic["Both Flavor"].append( int(round(sum_process_diff + sum_process_same)) )
    if process == "WG":
      yield_dic["Process"].append( process + "(*)" )
    elif process == "WW":
      yield_dic["Process"].append( "Total: " + process )
    else: yield_dic["Process"].append( process )

    #print process, sum_process_same, sum_process_diff
    tot_same += sum_process_same
    tot_diff += sum_process_diff

  if type(df_da) != type(None):
    sum_process_same = 0
    sum_process_diff = 0
    for process in processes:#df.process.unique():
      if process in df.process.unique():
        if "WJ" == process and type(df_da) != type(None): continue
        for decay in df[df.process == process].process_decay.unique():
          if decay in scales.keys():
            sum_process_same += df[(df.process_decay==decay) & (df.lep_Type < 0) & (df.lep1_Charge == df.lep2_Charge)].weight.sum() * scales[decay]
            sum_process_diff += df[(df.process_decay==decay) & (df.lep_Type > 0) & (df.lep1_Charge == df.lep2_Charge)].weight.sum() * scales[decay]
            #print decay, df[(df.process_decay==decay) & (df.lep1_Charge == df.lep2_Charge)].weight.sum() * scales[decay]
    #print df_da[(df_da.lep1_Charge == df_da.lep2_Charge) ].shape[0], sum_process_diff + sum_process_same
    sum_process_same =  max([df_da[(df_da.lep_Type < 0) & (df_da.lep1_Charge == df_da.lep2_Charge) ].shape[0] - sum_process_same, 0]) * ds_ov_ss
    sum_process_diff =  max([df_da[(df_da.lep_Type > 0) & (df_da.lep1_Charge == df_da.lep2_Charge) ].shape[0] - sum_process_diff, 0]) * ds_ov_ss  
    yield_dic["Same Flavor"].append( int(round(sum_process_same)) )
    yield_dic["Diff Flavor"].append( int(round(sum_process_diff)) )
    yield_dic["Both Flavor"].append( int(round(sum_process_diff + sum_process_same)) )
    yield_dic["Process"].append( 'WJ' )    

    tot_same += sum_process_same
    tot_diff += sum_process_diff


  yield_dic["Same Flavor"].append( int(round(tot_same)) )
  yield_dic["Diff Flavor"].append( int(round(tot_diff)) )
  yield_dic["Both Flavor"].append( int(round(tot_diff + tot_same)) )
  yield_dic["Process"].append( "Total" )


  if type(df_da) != type(None):
    data_same = df_da[(df_da.lep_Type < 0) & (df_da.lep1_Charge != df_da.lep2_Charge)].shape[0]
    data_diff = df_da[(df_da.lep_Type > 0) & (df_da.lep1_Charge != df_da.lep2_Charge)].shape[0]
    yield_dic["Same Flavor"].append( int(round(data_same)) )
    yield_dic["Diff Flavor"].append( int(round(data_diff)) )
    yield_dic["Both Flavor"].append( int(round(data_diff + data_same)) )
    yield_dic["Process"].append( "DATA" )

  yield_df = pd.DataFrame( yield_dic )
  
  return yield_df

def save_df_to_html( df, file_name, columns=["Process", "Same Flavor", "Diff Flavor"], header="<h3><b>Yields</b></h3>\n"):
  """
  Save file to tables directory
  """
  f = open(production_path+"/tables/"+file_name, "w")
  f.write(header)
  f.write('<div style="float:left; width:60%">\n')
  f.write( df.to_html(columns=columns, index=False) )
  f.write('</div>\n')
  f.write('<div style="float:right; width:40%">\n')
  f.write('<p style="position:absolute; right:0; bottom:80px; background-color:#62364C; color:white"><b>Purity:</b><br>Same Flavor: ' +\
     str(round(float(df[df.Process == 'WW']["Same Flavor"].max()) / df[df.Process == 'Total']["Same Flavor"].values[0], 2) ) + '<br>Different Flavor: '+\
     str(round(float(df[df.Process == 'WW']["Diff Flavor"].max()) / df[df.Process == 'Total']["Diff Flavor"].values[0], 2) )+'</p>')
  f.write('</div>')
  f.close()





#Ok this is a temp
def full_bin_plot(df, df_da, binned_feature, query=None,scales=scales, logy=True, x_range=None, y_range=None):
    """
    Returns a, b, c, d
    """
    df_diff    = df[df.lep1_Charge != df.lep2_Charge]
    df_da_diff = df_da[df_da.lep1_Charge != df_da.lep2_Charge]
    
    df_same    = df[df.lep1_Charge == df.lep2_Charge]
    df_da_same = df_da[df_da.lep1_Charge == df_da.lep2_Charge]
    
    if query != None:
        df_diff = df_diff.query(query)
        df_da_diff = df_da_diff.query(query)
        
        df_same = df_same.query(query)
        df_da_same = df_da_same.query(query)
    
    df_bin    = bin_df(df_diff, scales=scales, binned_feature=binned_feature, range=x_range)
    df_da_bin = bin_df(df_da_diff, binned_feature=binned_feature, range=x_range)

    #Same charge
    df_bin_same    = bin_df(df_same, scales=scales, binned_feature=binned_feature, range=x_range)
    df_da_bin_same = bin_df(df_da_same, binned_feature=binned_feature, range=x_range)

    #Loop through WJets set to zero and replace
    sum_mc = np.zeros(df_bin[df.process_decay.values[0]][0].shape[0])
    sum_wj = np.zeros(df_bin[df.process_decay.values[0]][0].shape[0])
    for k in df_bin_same:
        if k in ['W1JetsToLNu', 'W2JetsToLNu', 'W3JetsToLNu', 'W4JetsToLNu']: 
          sum_wj += df_bin_same[k][0]
          continue 
        if k in ['plotting']: continue 
        sum_mc += df_bin_same[k][0]

    if "Da" in df_da_bin_same:
        sum_mc = np.maximum(df_da_bin_same["Da"][0] - sum_mc, 0 * sum_mc)
        sf_wj = sum_mc #/ sum_wj
        sf_wj[np.isinf(sf_wj)] = 0.
        sf_wj[np.isnan(sf_wj)] = 0.

    ##Added a tab 
        ds_ov_ss = 2.05

        process_temp = None
        sum_wj = sum_wj * 0
        for process in ['W1JetsToLNu', 'W2JetsToLNu', 'W3JetsToLNu', 'W4JetsToLNu']:
            if process in df_bin.keys():
                process_temp = process
                ###
                sum_wj += df_bin[process][0] 
                ###
                df_bin[process][0] = df_bin[process][0] * 0
                df_bin[process][3] = df_bin[process][3] * 0


        if process_temp != None:
            #print process_temp, sum_mc
            df_bin[process_temp][0] = ds_ov_ss * sf_wj
            df_bin[process_temp][3] = sf_wj**0.5
    ##Tab was added above 

    if logy == None and y_range == None:
      logy = True
      y_range = None
      if df_da_bin["Da"][0].max() < 1000: 
        logy = False
        y_range = (0, df_da_bin["Da"][0].max() * 1.3)

    figs, ax = full_plot(df_bin, df_da_bin, color="color_1", logy=logy, x_range=x_range, y_range=y_range)
    return df_bin, df_da_bin, figs, ax 
###################

