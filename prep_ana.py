import root_pandas as rp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import copy, string 
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator
from ast import literal_eval as make_tuple
plt.rc('text', usetex=True)



home_dir = "/home/gunter/WW_analysis"


print("Loading plotting specs...")
binning_options = pd.read_csv( home_dir+"/binning_options.txt", index_col=None, sep=";" )

page = "8TeV"
lumi_amount="19.7"
page_temp = raw_input("8 or 13 TeV:")
if page_temp != "8TeV" and page_temp != "13TeV":
  print("Loading 8TeV")
else:
  page = page_temp
  lumi_amount="35.9"
plotting_options = pd.read_csv(home_dir+"/plotting_options_"+page+".csv", index_col=None, sep=";")

print("unc_mc_process and scales as dictionaries")

data_path = home_dir + "/data_" + page #+ "_old"
if input("Load MC and Data?"):
  columns = None
  #if page == "13TeV" :
  columns = [ 'process', 'process_decay',
      'weight', 'lep1_Charge', 'lep2_Charge', 'lep_Type', 'numbExtraLep',
      'tot_npv', 'mll', 'numb_jets', 'metMod', 'numb_BJet',
      'lep1_pt', 'lep2_pt', 'lep3_pt', 'dPhiLL',
      'jet1_pt', 'jet2_pt', 'HT', 'jet1_csv',
      'METProj_sin', 'met_over_sET','METProj_trk_sin', 'METProj',
      'dPhiLLMET',  'mllMET', 'qT', 'recoil', 'dPhiLLJet', 'dPhiMETJet'] #'met_phi'

  df_dy0 = rp.read_root(data_path+"/dyjetstoll_m-50_complete.root", columns=columns)
  df_dy_m_10 = rp.read_root(data_path+"/dyjetstoll_m-10to50_complete.root", columns=columns)
  df_ww = rp.read_root(data_path+"/ww_complete.root", columns=columns)
  df_tt_l = rp.read_root(data_path+"/ttbar_leptonic_complete.root", columns=columns)
  df_tt_sl = rp.read_root(data_path+"/ttbar_semileptonic_complete.root", columns=columns)
  df_tbar_tw = rp.read_root(data_path+"/tbar_tw-_complete.root", columns=columns)
  df_tbar_s = rp.read_root(data_path+"/tbar_s-_complete.root", columns=columns)
  df_tbar_t = rp.read_root(data_path+"/tbar_t-_complete.root", columns=columns)
  df_t_tw = rp.read_root(data_path+"/t_tw-_complete.root", columns=columns)
  df_t_s = rp.read_root(data_path+"/t_s-_complete.root", columns=columns)
  df_t_t = rp.read_root(data_path+"/t_t-_complete.root", columns=columns)
  df_zz_lq = rp.read_root(data_path+"/zzjetsto2l2q_complete.root", columns=columns)
  df_zz_ln = rp.read_root(data_path+"/zzjetsto2l2nu_complete.root", columns=columns)
  df_wz_lq = rp.read_root(data_path+"/wzjetsto2l2q_complete.root", columns=columns)
  df_wz_ln = rp.read_root(data_path+"/wzjetsto3lnu_complete.root", columns=columns)

  #May want to turn this off
  df_wj1 = rp.read_root(data_path+"/w1jetstolnu_complete.root", columns=columns)
  df_wj2 = rp.read_root(data_path+"/w2jetstolnu_complete.root", columns=columns)
  df_wj3 = rp.read_root(data_path+"/w3jetstolnu_complete.root", columns=columns)
  df_wj4 = rp.read_root(data_path+"/w4jetstolnu_complete.root", columns=columns)

  #Other stuff that I might want to turn off
  df_wg = rp.read_root(data_path+"/wg_complete.root", columns=columns)

  df = pd.concat([df_dy0, df_dy_m_10, df_ww, df_tt_l, df_tt_sl, df_tbar_s, df_tbar_t, df_tbar_tw, df_t_s, df_t_t, df_t_tw, df_zz_lq, df_zz_ln, df_wz_ln, df_wz_lq,
                  df_wj1, df_wj2, df_wj3, df_wj4, df_wg])#, df_hg, df_wgs, df_ggWW ])

  df = df.reset_index()

  if input("Correct MET?"):
    def apply_met_correction( array_values, array_labels):
      a = array_values
      if page == "8TeV":
        a[ array_labels == "DY"] = array_values[ array_labels == "DY"] * 1.059
      else:
        a[ array_labels == "DY"] = array_values[ array_labels == "DY"] * 1.059
      return a 
    #a = 1. * df.metMod.values
    #a[df.process.values == "DY"] = df.metMod.values[df.process.values == "DY"] * 1.059
    df["metMod"] = apply_met_correction(df.metMod.values, df.process.values)
    df["met_over_sET"] = apply_met_correction(df.met_over_sET.values, df.process.values)
    df["METProj"] = apply_met_correction(df.METProj.values, df.process.values)




print( "df = pd.concat([df_dy0, df_dy1, df_dy2, df_dy3, df_dy4, df_dy_m_10, df_ww, df_tt_l, df_tt_sl, df_zz_ln, df_wz_ln, df_wz_lq ])")
##################
# Process Scalings
scales = {}
for key in plotting_options.process_decay.unique():
  scales[key] = eval(plotting_options[plotting_options.process_decay == key]["scale"].values[0]) 

# unc  process
unc_mc_process = {}#{ "WW": .05, "DY": .03, "TT": 0.05, "ZZ": 0.1, "WZ":.1 }
for key in plotting_options.process_decay.unique():
  unc_mc_process[key] = plotting_options[plotting_options.process_decay == key]["unc"].values[0]



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
      
    

##################
#:def rm_duplicates_df(df_mu, df_el):
#    df_el_mu = pd.concat([df_mu[["runNum", "eventNumb", "lep1_pt", "process_decay"]], df_el[["runNum", "eventNumb", "lep1_pt", "process_decay"]]])    
#
#    df_el_mu = df_el_mu[ df_el_mu[["runNum", "eventNumb", "lep1_pt",]].duplicated() ]
#    
#    tot_run_event = df_el_mu.runNum.values * 10**12 + df_el_mu.eventNumb.values
#    el_run_event = df_el.runNum.values * 10**12 + df_el.eventNumb.values
#    
#    mask_index = np.where( np.in1d(el_run_event, tot_run_event) )
#
#    return df_el.drop(mask_index)


##################
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

def deltaPhi( dPhi ):
  """
  """
  if dPhi == None:
    dPhi = 0
  while dPhi >= np.pi: dPhi -= 2.*np.pi
  while dPhi < -np.pi: dPhi += 2.*np.pi
  return dPhi

def proj_calc( met, met_phi, phi_l1, phi_2):
  """
  """
  met_proj = 0
  dPhi = abs( deltaPhi( phi_l1 - met_phi) )
  if dPhi > abs( deltaPhi( phi_l2 - met_phi) ):
    dPhi = abs( deltaPhi( phi_l2 - met_phi) )
  

  if dPhi < np.pi / 2.:
    met_proj = abs( sin(dPhi) * met )
  else:
    met_proj = met

  return met_proj

def met_over_sum_et_calc( met, obj1, obj2):
  """
  """
  electron_mass = 0.511 * 10**-3
  muon_mass = 105. * 10**-3

  def set_mass( obj ):
    if obj["flavor"] == 13:
      obj["mass"] = muon_mass
  set_mass( obj1 )
  set_mass( obj2 )
  sum_et = ( obj1["pt"]**2 + obj1["mass"]**2 + obj1["pt"]**2 + obj1["mass"]**2 )**.5

  return met / sum_et

def pre_cuts( df, diff_charge= True):
  dif_lep = df.lep_Type > 0
  sam_lep = df.lep_Type < 0
  z_mass = (df.mll < 76) | (df.mll > 106 )
  nBjet = df.numb_BJet == 0
  extra_lep = df.numbExtraLep == 0
  quality_cuts = df.mll > 30 
  charge = df.lep1_Charge != df.lep2_Charge

  if diff_charge == False:
    s_df = df[ sam_lep & z_mass & nBjet & extra_lep & quality_cuts ]
    d_df = df[ dif_lep & nBjet & extra_lep & quality_cuts ]
  
  if diff_charge == True:
    s_df = df[ sam_lep & z_mass & nBjet & extra_lep & quality_cuts & charge ]
    d_df = df[ dif_lep & nBjet & extra_lep & quality_cuts & charge ]  
  
  return pd.concat([ s_df, d_df])


def cuts_ana( df, flavor="both" ):
  """
    Basic WW cuts analysis !!!CLEAN THIS SHIT UP!!!
    cuts_ana( df ) 
    df: data frame
  """
  initial_cuts = (df.mll > 30 ) & (df.METProj > 50) & (df.numbExtraLep == 0) & (df.HT < 50.) & (df.numb_jets < 2)
  same_flavor_cut = (df.lep_Type < 0)
  dif_flavor_cut = (df.lep_Type > 0)
  z_peak_cut = (df.mll > 106) | (df.mll < 76)
  met_proj_cut = (df.METProj > 50)
  met_diff_flavor = (df.metMod > 50) & (df.metMod < 120)
  mllmet = df.mllMET > 100

  basic_sf_cuts = initial_cuts & same_flavor_cut & met_proj_cut & z_peak_cut & (df.qT > 35)
  basic_df_cuts = initial_cuts & dif_flavor_cut & met_proj_cut & met_diff_flavor & mllmet 
  if flavor=="both": return pd.concat( [df[basic_sf_cuts], df[basic_df_cuts]] )
  elif flavor=="same": return df[basic_sf_cuts]
  elif flavor=="diff": return df[basic_df_cuts]
  else: 
    print "Broke"
    return -999

def WW_ana( df , diff_charge= True): 
  initial_cuts = (df.mll > 50 )  & (df.numbExtraLep == 0) & (df.numb_jets ==  0) & ( df.metMod > 50 )
  dif_flavor_cut = (df.lep_Type > 0)

  if diff_charge == True:
    initial_cuts = initial_cuts & (df.lep1_Charge != df.lep2_Charge)

  basic_df_0j_cuts = initial_cuts & dif_flavor_cut 
  return df[basic_df_0j_cuts]


def TT_ana( df, diff_charge= True ): 
  initial_cuts = (df.mll > 50 ) & (df.qT > 30) & (df.metMod > 60) & (df.numbExtraLep == 0) & (df.numb_jets > 0)
  dif_flavor_cut = (df.lep_Type > 0)
  
  basic_df_0j_cuts = initial_cuts & dif_flavor_cut 
  return df[basic_df_0j_cuts] 

def DY_ana( df, diff_charge = True ): 
  initial_cuts = (df.mll > 50 )  & (df.numbExtraLep == 0) & (df.numb_jets == 0) & (df.metMod < 80)

  return df[ initial_cuts ]

def Z_tt_ana( df, diff_charge = True ): 
    initial_cuts = (df.mll > 50 )  & (df.numbExtraLep == 0) & (df.numb_jets ==  0) & ( df.qT < 10 ) & (df.lep1_pt > 30) & (df.lep2_pt > 20)
    
    dif_flavor_cut = (df.lep_Type > 0)
    
    m_1 = (df.mll < 120) & (df.mll > 60)
    met_1 = (df.metMod > 30) & (df.metMod < 50)

    return df[ initial_cuts & met_1 & m_1 & dif_flavor_cut ] 

def wjets_selection(df_mc, df_data, feature):
    bins_mc_ = bin_df( df_mc[df_mc.lep1_Charge == df_mc.lep2_Charge], feature, )#weights=False) 
    bins_data_ = bin_df( df_data[df_data.lep1_Charge == df_data.lep2_Charge], feature, )

    sum_mc = np.zeros(bins_mc_[bins_mc_.keys()[0]][0].shape[0])
    for k in bins_mc_.keys():
        if k in ['ttbar_semileptonic', 'ttbar_leptonic', 'DYJetsToLL_M-50', 'WZJetsTo3LNu', 'WW',
     'Tbar_tW-channel','T_t-channel','ZZJetsTo2L2Nu', 'WGToLNuG']:
            sum_mc += bins_mc_[k][0]

    sum_mc = bins_data_["Da"][0] - sum_mc
    return sum_mc

###############################################
#### Create kinematic histograms
def make_control_plots( df_mc, df_data, date_tag, selection_tag, energy_dir= "8TeV" ):
  """
  Creates all control region plots.
  """
  def full_ana( df, diff_charge=True ):
    return df

  def same_ana( df, diff_charge= True ):
    return df[df.lep_Type < 0]

  def diff_ana( df, diff_charge= True ):
    return df[ df.lep_Type > 0]

  control_regions = {"WW": WW_ana, "TT": TT_ana, "DY": DY_ana, "Z_tt": Z_tt_ana, "full": full_ana, "same": same_ana, "diff": diff_ana}  

  for key in control_regions.keys():
    print key
    create_kinematic_hist( control_regions[key](df_mc, diff_charge= False), control_regions[key](df_data, diff_charge= False), prefix= energy_dir + "/" + selection_tag + "/" + key + "/" + date_tag )

def create_kinematic_hist(df_mc, df_data, prefix="", scales=scales):
  """
  Creates all the basic histograms you'll ever need:
  create_kinematic_hist(df)
  """

  range = (0,250)
  bins  = 100

  features = [ 'numb_BJet', 'HT',
       'numb_jets', 'lep1_pt',
       'jet1_pt', 'lep2_pt',
       'jet2_pt', 'metMod', 'dPhiLL', 'METProj',
       'qT', 'dPhiLLJet', 'met_phi',
       'lep3_pt', 'tot_npv',
       'mll', 'METProj_sin', 'met_over_sET','METProj_trk_sin',
        'recoil', 'jet1_pt', 'dPhiLLMET', 'dPhiMETJet', 'mllMET']

  for feature in features:
    if feature not in df_mc.keys(): continue
    print feature
    bins_mc = bin_df( df_mc[df_mc.lep1_Charge != df_mc.lep2_Charge], feature, scales=scales) 
    bins_data = bin_df( df_data[df_data.lep1_Charge != df_data.lep2_Charge], feature, ) 

    #########################
    #WJets stuff
    bins_mc_ = bin_df( df_mc[df_mc.lep1_Charge == df_mc.lep2_Charge], feature, scales=scales) 
    bins_data_ = bin_df( df_data[df_data.lep1_Charge == df_data.lep2_Charge], feature, ) 

    sum_mc = np.zeros(bins_mc["WW"][0].shape[0])
    for k in bins_mc_.keys():
      if k in ['ttbar_semileptonic', 'ttbar_leptonic', 'DYJetsToLL_M-50', 'WZJetsTo3LNu', 'WW',
   'Tbar_tW-channel','T_t-channel','ZZJetsTo2L2Nu', 'WGToLNuG']:
        sum_mc += bins_mc_[k][0]
    
    if "Da" in bins_data_.keys():
      sum_mc = bins_data_["Da"][0] - sum_mc

    process_temp = None
    for process in ['W1JetsToLNu', 'W2JetsToLNu', 'W3JetsToLNu', 'W4JetsToLNu']:
      if process in bins_mc.keys():
        process_temp = process
        bins_mc[process][0] = bins_mc[process][0] * 0
        bins_mc[process][3] = bins_mc[process][3] * 0

    if process_temp != None:
      bins_mc[process_temp][0] = sum_mc
      bins_mc[process_temp][3] = sum_mc
    ##############################


    logy = True
    y_range = None
    if bins_data["Da"][0].max() < 1000: 
      logy = False
      y_range = (0, bins_data["Da"][0].max() * 1.3)

    figs, ax = full_plot(bins_mc, bins_data, color="color_1", logy=logy, y_range=y_range)
    figs.savefig(home_dir + '/plots/' + prefix + "_"+ feature + ".png")

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
def bin_df( df, binned_feature, binning_options=binning_options, plotting_options=plotting_options, scales=scales, range=None, bins=None, lumi_scale=1, density=False, weights=True):
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
      #print default, binning_options[ binning_options.feature == "default_"+default ].range
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
    #print process
    if process in unique_df_processes:
      df_process = df[df.process_decay == process]
      #print process, scales[process], range, bins
      if weights == True:
        weights_arr = df_process.weight.values
        sq_weights_arr = df_process.weight.values**2
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
        sum_yerr += bins[process_decay][3] #+ float(lumi_amount)**2 * scales[process_decay]**2 * unc_mc_process[process_decay]**2
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
      rect.append(ax.bar( tot_bins[process][1][:-1], tot_bins[process][0],
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
  #ax.bar( tot_bins[process_][1][:-1], 2*sum_yerr, tot_bins[process_][1][1] - tot_bins[process_][1][0], bottom= sum_bins - sum_yerr , alpha= .05, color="none", edgecolor='black', hatch='//')


  #Configurables
  if logy == True: ax.set_yscale("log", nonposy='clip')
  if x_range!=None: ax.set_xlim(x_range)
  if y_range==None: 
    if logy == True: ax.set_ylim( bottom=1,  top= sum_bins.max()*30.)
    else: ax.set_ylim( bottom=0,  top= sum_bins.max()*2.)
  elif type(y_range)==tuple: ax.set_ylim( y_range )


  ax.xaxis.labelpad = 20
  ax.yaxis.labelpad = 15


  if y_label != None:
      ax.set_ylabel(y_label, fontsize=22, fontname='Bitstream Vera Sans', )
  if title != None:
      plt.xlabel( title, fontname='Bitstream Vera Sans', fontsize=24)#position=(1., 0.), va='bottom', ha='right',)

  #plt.rc('text', usetex=True)
  ax.set_title(r"\textbf{CMS} Preliminary \hspace{8cm} $"+ lumi_amount +" fb^{-1}$ $\sqrt{s}="+page+"$", fontname='Bitstream Vera Sans', fontsize=24)

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
  #test = None
  #for k in  tot_bins.keys():
  #  if test == None:
  #    test = copy.copy(tot_bins[k][0])
  #  else :
  #    test += tot_bins[k][0]
  #print test
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
    ax.set_ylabel('MC/DATA', fontname='Bitstream Vera Sans', fontsize=20)
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

  hist_stuff = plot_hist( bins_mc, processes=processes, x_range=x_range, y_range=y_range, title=title, y_label=y_label, colors=color, logy=logy, x_minor_ticks=x_minor_ticks, ax=ax[0], lumi_amount=lumi_amount)
  err_stuff  = plot_errbar( bins_data, ax=ax[0])
  #wg_loc = [for i in hist_stuff[2] if i == "WG"][0]
  if page == "8TeV": 
    if len( [i for i, ele in enumerate(hist_stuff[2]) if ele == "WG"]) > 0:
        hist_stuff[2][[i for i, ele in enumerate(hist_stuff[2]) if ele == "WG"][0]] = "WG(*)"
  ax[0].legend(hist_stuff[1] + [err_stuff], hist_stuff[2] + ['Da'], frameon=False, fontsize="x-large", numpoints=1)

  plot_ratio( bins_mc, bins_data, ax=ax[1])

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

def process_yields( df, df_da=None, processes= ['WW', 'DY', 'Top', 'WZ', 'ZZ', 'WG', 'WJ', 'Higgs'], scales=scales ):
  """
  return a dataframe
  """
  tot_same = 0
  tot_diff = 0
  
  yield_dic = {}
  yield_dic["Diff Flavor"] = []
  yield_dic["Same Flavor"] = []
  yield_dic["Process"] = []
  
  for process in processes:
    if "WJ" in process and type(df_da) != type(None): 
      continue
    if process in df.process.unique():
      sum_process_same = 0
      sum_process_diff = 0
      for decay in df[df.process == process].process_decay.unique():
        if decay in scales.keys():
          if process == "WW":
            sum_process_same_ = df[(df.process_decay==decay) & (df.lep_Type < 0) & (df.lep1_Charge != df.lep2_Charge)].weight.sum() * scales[decay]
            sum_process_diff_ = df[(df.process_decay==decay) & (df.lep_Type > 0) & (df.lep1_Charge != df.lep2_Charge)].weight.sum() * scales[decay]
            yield_dic["Same Flavor"].append( int(round(sum_process_same_)) )
            yield_dic["Diff Flavor"].append( int(round(sum_process_diff_)) )
            yield_dic["Process"].append( decay )

          sum_process_same += df[(df.process_decay==decay) & (df.lep_Type < 0) & (df.lep1_Charge != df.lep2_Charge)].weight.sum() * scales[decay]
          sum_process_diff += df[(df.process_decay==decay) & (df.lep_Type > 0) & (df.lep1_Charge != df.lep2_Charge)].weight.sum() * scales[decay]
      print process, sum_process_same, sum_process_diff

      yield_dic["Same Flavor"].append( int(round(sum_process_same)) )
      yield_dic["Diff Flavor"].append( int(round(sum_process_diff)) )
      if process == "WG":
        yield_dic["Process"].append( process + "(*)" )
      else: yield_dic["Process"].append( process )

      tot_same += sum_process_same
      tot_diff += sum_process_diff

  if type(df_da) != type(None):
    for process in processes:#df.process.unique():
      if process in df.process.unique():
        sum_process_same = 0
        sum_process_diff = 0
        for decay in df[df.process == process].process_decay.unique():
          if "WJ" in decay and df_da != None: continue
          if decay in scales.keys():
            sum_process_same += df[(df.process_decay==decay) & (df.lep_Type < 0) & (df.lep1_Charge == df.lep2_Charge)].weight.sum() * scales[decay]
            sum_process_diff += df[(df.process_decay==decay) & (df.lep_Type > 0) & (df.lep1_Charge == df.lep2_Charge)].weight.sum() * scales[decay]
    sum_process_same =  df_da[ (df_da.lep_Type < 0) & (df_da.lep1_Charge == df_da.lep2_Charge) ].weight.sum() - sum_process_same
    sum_process_diff =  df_da[ (df_da.lep_Type > 0) & (df_da.lep1_Charge == df_da.lep2_Charge) ].weight.sum() - sum_process_diff
    yield_dic["Same Flavor"].append( int(round(sum_process_same)) )
    yield_dic["Diff Flavor"].append( int(round(sum_process_diff)) )
    yield_dic["Process"].append( 'WJ' )    

    tot_same += sum_process_same
    tot_diff += sum_process_diff


  yield_dic["Same Flavor"].append( int(round(tot_same)) )
  yield_dic["Diff Flavor"].append( int(round(tot_diff)) )
  yield_dic["Process"].append( "Total" )


  if type(df_da) != type(None):
    data_same = df_da[(df_da.lep_Type < 0) & (df_da.lep1_Charge != df_da.lep2_Charge)].shape[0]
    data_diff = df_da[(df_da.lep_Type > 0) & (df_da.lep1_Charge != df_da.lep2_Charge)].shape[0]
    yield_dic["Same Flavor"].append( int(round(data_same)) )
    yield_dic["Diff Flavor"].append( int(round(data_diff)) )
    yield_dic["Process"].append( "DATA" )

  print yield_dic
  yield_df = pd.DataFrame( yield_dic )
  
  return yield_df

def save_df_to_html( df, file_name, columns=["Process", "Same Flavor", "Diff Flavor"], header="<h3><b>Yields</b></h3>\n"):
  """
  Save file to tables directory
  """
  f = open(home_dir+"/plots/tables/"+file_name, "w")
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
