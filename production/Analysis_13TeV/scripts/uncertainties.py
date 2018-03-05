#Thoth Gunter

import os, sys, time
sys.path.append(os.getcwd() + "/../")
sys.path.append(os.getcwd() + "/../tools/")

from prep_ana_II import * 
from cross_section_calc import calc_cross_stuff, cross_calc, stat_unc_calc, normalization_unc_calc 
import pile_up
np.random.seed(108)



###############################################################
###############################################################
def toy_setup( df, lep_id= 13):
  if lep_id == 13:
    lep_id_weight = "mu_id_weight"

  if lep_id == 11:
    lep_id_weight = "elX_id_weight"#note X

  id_eff  = 1 * df[(df.lep1_type == lep_id) & (df.lep2_type == lep_id)][lep_id_weight].values
  low_1   = 1 * df[(df.lep1_type == lep_id) & (df.lep2_type == lep_id)].lep1_id_error_low.values
  high_1  = 1 * df[(df.lep1_type == lep_id) & (df.lep2_type == lep_id)].lep1_id_error_high.values
  low_2   = 1 * df[(df.lep1_type == lep_id) & (df.lep2_type == lep_id)].lep2_id_error_low.values
  high_2  = 1 * df[(df.lep1_type == lep_id) & (df.lep2_type == lep_id)].lep2_id_error_high.values

  return id_eff, low_1, low_2, high_1, high_2



def toy(id_eff, low_1, low_2, high_1, high_2):
  #print "Doing a toy"
  rt_eff = 1 * id_eff #1 * np.array to make sure we return a copy

  unique_errors = np.unique(low_1)
  for i in unique_errors:
    std = ((low_1[low_1 == i] +  high_1[low_1 == i]) / 2)[0] #standard deviation
    delta = np.random.normal(0, std,) #delta 1
    rt_eff[low_1 == i] = id_eff[low_1 == i] + delta 
    rt_eff[low_2 == i] = id_eff[low_2 == i] + delta 

  return rt_eff


def toy_el_temp( el1_id_eff, el2_id_eff, df_,  error_dic):
    el1_eff = 1 * el1_id_eff #1 * np.array to make sure we return a copy
    el2_eff = 1 * el2_id_eff
    for i in error_dic:
        for j in error_dic[i]:
            b_1 = (i[0] < df_.lep1_eta.abs()) &  (df_.lep1_eta.abs() < i[1]) & (df_.el1_id_weight == j)
            b_2 = (i[0] < df_.lep2_eta.abs()) &  (df_.lep2_eta.abs() < i[1]) & (df_.el2_id_weight == j)
            c = np.random.normal(0, error_dic[i][j])
            el1_eff[b_1] = el1_id_eff[b_1] + c
            el2_eff[b_2] = el2_id_eff[b_2] + c

    return el1_eff * el2_eff

def toy_both(mu_id_eff, el1_id_eff, el2_id_eff, lep1_type, lep2_type, mu_unique, el_unique, low_1, high_1, low_2, high_2):
  mu_eff = 1 * mu_id_eff
  el1_eff = 1 * el1_id_eff #1 * np.array to make sure we return a copy
  el2_eff = 1 * el2_id_eff

  mu1_mask = lep1_type == 13
  mu2_mask = lep2_type == 13
  el1_mask = lep1_type == 11
  el2_mask = lep2_type == 11


  #Electrons
  for i in el_unique:
    low_1_mask = low_1 == i
    low_2_mask = low_2 == i
    b_1 = (( low_1[low_1_mask] +  high_1[ low_1_mask ]) / 2)[0] #standard deviation
    b_2 = (( low_2[low_2_mask] +  high_2[ low_2_mask ]) / 2)[0] #standard deviation
    c_1 = np.random.normal(0, b_1,) #delta 1
    c_2 = np.random.normal(0, b_2,) #delta 2
    el1_eff[low_1_mask & el1_mask] = el1_eff[low_1_mask & el1_mask] + c_1
    el2_eff[low_2_mask & el2_mask] = el2_eff[low_2_mask & el2_mask] + c_2
  """Restructor. Should follow the example set by Muons
    for i in el_unc:
      for j in el_unc[i]:
        b_1 = (i[0] < np.abs(df.lep1_eta)) &  (np.abs(df.lep1_eta) < i[1]) & (el1_eff == j)
        b_2 = (i[0] < np.abs(df.lep2_eta)) &  (np.abs(df.lep2_eta) < i[1]) & (el2_eff == j)
        c = np.random.normal(0, el_unc[i][j])
        el1_eff[b_1 & el1_mask ] = el1_id_eff[b_1 & el1_mask] +  c
        el2_eff[b_2 & el2_mask] = el2_id_eff[b_2 & el2_mask] + c
  """

  #Muons
  for i in mu_unique:
    low_1_mask = low_1 == i
    low_2_mask = low_2 == i


    b_1 = ((low_1[low_1_mask] +  high_1[ low_1_mask ]) / 2)[0] #standard deviation
    #b_2 = ((low_2[low_2_mask] +  high_2[ low_2_mask ]) / 2)[0] #standard deviation
    c_1 = np.random.normal(0, b_1,) #delta 1
    #c_2 = np.random.normal(0, b_2,) #delta 2
    mu_eff[low_1_mask & mu1_mask] = mu_eff[low_1_mask & mu1_mask] + 2*c_1
    #mu_eff[low_2_mask & mu2_mask] = mu_eff[low_2_mask & mu2_mask] + c_2

  #Remove Zeros
  mu_eff[mu_eff == 0] = 1.
  el1_eff[el1_eff == 0] = 1.
  el2_eff[el2_eff == 0] = 1.

  return mu_eff * el1_eff * el2_eff

###############################################################
###############################################################

def lepton_scale_factor(X_section, flavor, ww_eff):
  ###Prob works but takes too long
  #########################
  #print("loading data")
  columns = ["numbExtraLep", "lep_Type", "mll", "numb_BJet", "lep1_Charge", "lep2_Charge", "process_decay", "mu_id_weight", "lep1_id_error_low", "lep1_id_error_high", "lep2_id_error_low", "lep2_id_error_high", "lep1_eta", "lep2_eta", "lep1_pt", "lep2_pt", "lep1_type", "lep2_type", "el1_id_weight", "el2_id_weight", "gen_weight", "weight"]
  df_ww = rp.read_root(data_path+"/ww_complete.root", columns=columns)
  #df_dy = rp.read_root(data_path+"/dyjetstoll_m-50_complete.root", columns=columns)
  #df_tt = df_tt_l = rp.read_root(data_path+"/ttbar_leptonic_complete.root", columns=columns)
  df = pd.concat([df_ww, ])#df_dy, df_tt])
  df = pre_cuts(df)

  if flavor == 'diff':
    flavor_cut = df.lep_Type > 0
    df = df[flavor_cut]
  elif flavor == 'same':
    flavor_cut = df.lep_Type < 0 
    df = df[flavor_cut]
  
  df["scales"] = 1.
  df["gen_weight"] = df.gen_weight.values/ np.abs(df.gen_weight.values)
  df["weight"] = df.weight.values * df.gen_weight.values

  #print("Finished loading data")

  for process in scales:
    if process in df.process_decay.unique():
      df.loc[df.process_decay == process, "scales"] = scales[process]
  #########################
  toy_set = []
  ###############################################################
  """
  Re work so that electron and muon are together
  """
  ###############################################################
  ###############################################################
  #MUON
  ###############################################################
  id_eff, low_1, low_2, high_1, high_2 = toy_setup(df, lep_id= 13)
  #  for i in xrange(10):
  #  toy_set.append(toy(id_eff, low_1, low_2, high_1, high_2))
  #toy_set = np.array(toy_set) 


  weights = df[(df.lep1_type == 13) & (df.lep2_type == 13)].scales.values
  #print("toy set shape: ", toy_set.shape, weights.shape)
  #muon_unc = np.average(toy_set, axis=1, weights=weights).std()#toy_set.mean(1).std()
  #print("Muon unc: ", muon_unc)

  ###############################################################
  #ELECTRON
  ###############################################################
  #!I don't need this for 13TeV I took care of it properly in BLT
  #errors = {(0.0, 0.8): {0.827 : 0.021, 0.924 :0.010, 0.960 : 0.003 , 0.978: 0.001, 0.981: 0.001, 0.982: 0.002 },\
  #(0.8, 1.442): {0.948: (0.024+ 0.023) / 2, 0.932: 0.012, 0.936 : 0.004, 0.958: 0.002, 0.969: 0.001, 0.969: 0.002},\
  #(1.442, 1.556): {1.073: (0.117+ 0.107) / 2, 0.808:  (0.045+ 0.042) / 2, 0.933: (0.015+ 0.017) / 2, 0.907: 0.008, 0.904: 0.004, 0.926: 0.011},\
  #(1.556, 2.0): {0.854: (0.048 + 0.047) / 2, 0.853: 0.022, 0.879: 0.007, 0.909: 0.003, 0.942: 0.002, 0.957: 0.004 },\
  #(2.0, 2.5):{1.007: (0.047+ 0.046)/2, 0.903: 0.029, 0.974: 0.004, 0.987: 0.004, 0.991: 0.003, 0.999: 0.005 }}





  el1_id_eff = 1 * df[(df.lep1_type == 11) & (df.lep2_type == 11)].el1_id_weight.values
  el2_id_eff = 1 * df[(df.lep1_type == 11) & (df.lep2_type == 11)].el2_id_weight.values

  #toy_set_el = []
  #for i in xrange(10):
  #    toy_set_el.append(toy_el_temp( el1_id_eff, el2_id_eff, df[(df.lep1_type == 11) & (df.lep2_type == 11)],  errors))
  #toy_set_el = np.array(toy_set_el)

  #print "Electron stuff", toy_set_el.mean(1)[:3], "Mean: ", toy_set_el.mean(), "Uncertainty: ", toy_set_el.mean(1).std()

  ###############################################################
  #BOTH MUON AND ELECTRON
  ###############################################################
  mu_id_eff = 1 * df.mu_id_weight.values 
  el1_id_eff = 1 * df.el1_id_weight.values
  el2_id_eff = 1 * df.el2_id_weight.values

  lep1_type = 1 * df.lep1_type.values
  lep2_type = 1 * df.lep2_type.values

  low_1 = 1 * df.lep1_id_error_low.values
  high_1 = 1 * df.lep1_id_error_high.values
  low_2 = 1 * df.lep2_id_error_low.values
  high_2 = 1 * df.lep2_id_error_high.values

  lep1_eta = 1 * df.lep1_eta.values
  lep2_eta = 1 * df.lep2_eta.values

  #el_unc = errors

  mu_unique = np.unique(low_1[df.lep1_type == 13])
  el_unique = np.unique(low_1[df.lep1_type == 11])
  toy_set_both = []
  for i in xrange(10):
    #print "toy", i
    toy_set_both.append(toy_both(mu_id_eff, el1_id_eff, el2_id_eff, lep1_type, lep2_type, mu_unique, el_unique, low_1, high_1, low_2, high_2 ))
  toy_set_both = np.array(toy_set_both)

  #print "Both stuff", toy_set_both.mean(1)[:3], "Mean: ", toy_set_both.mean(), "Uncertainty: ", toy_set_both.mean(1).std()


  ###############################################################
  #CLEAN UP
  ###############################################################
  del df_ww
  #del df_dy
  #del df_tt
  return toy_set_both.mean(1).std() * (ww_eff)**-1 * X_section

###############################################################
###############################################################

###############################################################
###############################################################
def bjet_unc_toy(df):
  #Input df[(df.numb_BJet_gen > 0) ]
  jet_csv_mask = [ df.jet1_csv.values <= .846,
               df.jet2_csv.values <= .846,
               df.jet3_csv.values <= .846,
               df.jet1_csv.values > .846,
               df.jet2_csv.values > .846,
               df.jet3_csv.values > .846]
  #Needs to be updated

  mc_eff = np.array([
   0.4442664407048948,
   0.5102205372678276,
   0.5542576310070667,
   0.5823890951134911,
   0.6001000347947112,
   0.6196382554577722,
   0.6298785946506259,
   0.6329302867842538,
   0.5800544905623328,
   0.5224910394265233,
   0.47150500116306115,
   0.40217391304347827,
   0.33751743375174337,
   0.2570093457943925,
   0.2661290322580645,
   0.21739130434782608])


  bjet_unc = [0.0264008, 0.0272757, 0.0275565, 0.0248745, 0.0218456, 0.0253845, 0.0239588, 0.0271791, 0.0273912,
    0.0379822, 0.0411624, 0.0786307, 0.0866832, 0.0942053, 0.102403, 0.1]
  bjet_mc_eff_unc = [
  0.002044870533737016,
  0.0016060593437708893,
  0.0014327658052221096,
  0.0013855675110511305,
  0.0014172562032056706,
  0.0011053255585410577,
  0.001330940168220748,
  0.0013147241396643421,
  0.001985532350816485,
  0.003536128885636905,
  0.005681122677667734,
  0.009330720680568763,
  0.016278904251102247,
  0.030427507469274185,
  0.04294916994144168,
  0.09408167988273519]#np.sqrt(csv_bjet) / true_bjet

  ptmin = [30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800]
  ptmax = [40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800, 1000]


  bjet_iw = np.ones(df.shape[0])
  for i in xrange(len(ptmin)):
      #roll the dice
      offset_unc = np.random.normal(0,bjet_unc[i])
      offset_mc_eff = np.random.normal(0,bjet_mc_eff_unc[i])
      #loop over jets
      for jet_iter in range(1,len(jet_csv_mask)+1):
          jet_interest= 'jet'+str(jet_iter)




          jet_pt_mask = ( df[jet_interest+'_pt'] > ptmin[i] ) & ( df[jet_interest+'_pt'] < ptmax[i] )
          mc_eff_ = mc_eff[i] + offset_mc_eff



          csvt = lambda x: 0.884016*((1.+(0.0331508*x))/(1.+(0.0285096*x)))
          w = csvt(df[ jet_csv_mask[jet_iter-1] & jet_pt_mask ][jet_interest+'_pt'].values) + offset_unc

          if jet_iter < 3:
            bjet_iw[ jet_csv_mask[jet_iter-1].astype(bool) & jet_pt_mask.values ] *= (1 - w * mc_eff_) / (1 - mc_eff_)    #THIS IS WEIRD I'M CALCING INEFFS
          else:
            bjet_iw[ jet_csv_mask[jet_iter-1].astype(bool) & jet_pt_mask.values ] *= (w * mc_eff_) / (mc_eff_)    #THIS IS WEIRD I'M CALCING INEFFS
          



  return bjet_iw

def bjet_scale_factor( X_section, flavor ):

  ############################################
  ############################################
  def calc_unc_bjet(df, eps=.63, D= 1922276, s=1):
    """
    eps: efficiency of bjet
    D: ?number of data events? <= Should change per flavor 
    s: scale
    """

    B = (1 - eps)**2 * df[df.numb_BJet_gen == 2].shape[0] + (1 - eps) * df[df.numb_BJet_gen == 1].shape[0] + df[df.numb_BJet_gen == 1].shape[0]
    dB = 2*(1 - eps) * df[df.numb_BJet_gen == 2].shape[0] + df[df.numb_BJet_gen == 1].shape[0]

    B *= scales["ttbar_leptonic"] * s
    dB *= scales['ttbar_leptonic'] * s

    #D = 1922276 #Post preselection number

    return dB / (D - B) * 0.03 #0.03 is a conservative estimate taken from the bjet unc under Test scale factors header 
  ############################################
  ############################################

  jet_columns = ["jet1_csv", "jet2_csv", "jet3_csv","jet1_flv", "jet2_flv", "jet3_flv", "numb_BJet_gen", "numb_jets", "numb_BJet", "jet1_pt", "jet2_pt", "jet3_pt", "lep_Type"]

  df_tt_l = rp.read_root("/home/gunter/WW_analysis/data_alpha/ttbar_leptonic_complete.root",
                         columns=jet_columns)
  df_tt_sl = rp.read_root(data_path+"/ttbar_semileptonic_complete.root",
                          columns=jet_columns)
  df_top = pd.concat([df_tt_l, df_tt_sl])

  if flavor == "diff": 
    df_top = df_top[df_top.lep_Type > 0]
  elif flavor == "same":
    df_top = df_top[df_top.lep_Type < 0]
  bjet_unc = calc_unc_bjet(df_top[df_top.numb_BJet == 0], eps=.63, D= 5383, s= 2015./ (df_top[df_top.numb_BJet == 0].shape[0] * scales['ttbar_leptonic']) )
  return bjet_unc * X_section

  """
  def bjet...( ratio )
  #Suspect
  columns=['jet1_pt', 'jet2_pt', 'jet3_pt', 'jet4_pt', 'jet5_pt', 'jet6_pt',
           'jet1_csv', 'jet2_csv', 'jet3_csv', #'jet4_csv', 'jet5_csv', 'jet6_csv',
           'jet1_flv', 'jet2_flv', 'jet3_flv', #'jet4_flv', 'jet5_flv', 'jet6_flv',
           'numb_jets', 'numb_BJet_gen', 'numb_BJet', 'bjet_unc', 'bjet_weight']


  df_tt_l = rp.read_root(data_path+"/ttbar_leptonic_complete.root", columns=columns)
  df_tt_sl = rp.read_root(data_path+"/ttbar_semileptonic_complete.root", columns=columns)
  #df_tbar_tw = rp.read_root(data_path+"/tbar_tw-_complete.root", columns=columns)
  #df_tbar_s = rp.read_root(data_path+"/tbar_s-_complete.root", columns=columns)
  #df_tbar_t = rp.read_root(data_path+"/tbar_t-_complete.root", columns=columns)
  #df_t_tw = rp.read_root(data_path+"/t_tw-_complete.root", columns=columns)
  #df_t_s = rp.read_root(data_path+"/t_s-_complete.root", columns=columns)
  #df_t_t = rp.read_root(data_path+"/t_t-_complete.root", columns=columns)



  df = pd.concat([df_tt_l, df_tt_sl,])#df_tbar_s, df_tbar_t, df_tbar_tw, df_t_s, df_t_t, df_t_tw])
  df = df[(df.numb_BJet_gen > 0)] #There must be a taggable jet



  toy_set = []
  for i in xrange(10):
      toy_set.append(bjet_unc_toy(df))
  toy_set = np.array(toy_set) 

  del df_tt_l
  del df_tt_sl
  #del df_tbar_tw 
  #del df_tbar_s
  #del df_tbar_t
  #del df_t_tw
  #del df_t_s
  #del df_t_t
  print "Bjet uncertainty ", toy_set.mean(1).std()
  return toy_set.mean(1).std() * 2015 / ( float(lumi_amount) * 10**3 * ratio ) #2015 is the current, sept 15, number of ttbar events after rf
  """
###############################################################
###############################################################

###############################################################
###############################################################
def lep_scale_shift(data, lep_pt=1.003, rf=None): #0.2% muon scale and 0.3% electron scale
  #leptons
  data.lep1_pt = data.lep1_pt * lep_pt
  data.lep2_pt = data.lep2_pt * lep_pt
  data.mll = data.mll * lep_pt**.5
  data.qT = data.qT * lep_pt
  #ADD recoil

  #MET results
  data.metMod = data.metMod - (data.lep1_pt + data.lep2_pt)/ lep_pt * (lep_pt - 1) 
  data.metProj = data.METProj - (data.lep1_pt + data.lep2_pt) / lep_pt * (lep_pt - 1) 

  #rerun rf score calculation
  if rf != None:
   # print "Recreating random forest scores."
    pred_fTT = rf["clf_fTT"].predict_proba(np.float32(data[rf["features_fTT"]].values))
    data["pred_fTT_WW"] = pred_fTT[:,0]

    temp = data[rf["features_fDY"]]
    temp = temp.replace([np.inf,-np.inf], 0)
    pred_fDY = rf["clf_fDY"].predict_proba(np.float32(temp.values))
    data["pred_fDY_WW"] = pred_fDY[:,0]

def lep_pt_scale(flavor=None):
  df = load_presel_w_fDY_fTT_MC()
  df_da = load_presel_w_fDY_fTT_DATA()
  df_ww = rp.read_root(data_path+"/ww_complete.root", columns=columns)

  rfs = load_randomForest()

  #All the mask we will need
  mc_lepton_Z =  (df.mll < 76) |  (df.mll > 106)
  data_lepton_Z =  (df_da.mll < 76) |  (df_da.mll > 106)
  mc_leptons_pt = (df.lep1_pt > 27) & (df.lep2_pt > 25) & (df.mll > 30) & mc_lepton_Z
  data_leptons_pt = (df_da.lep1_pt > 27) & (df_da.lep2_pt > 25) & (df_da.mll > 30) & data_lepton_Z

  mc_same =  df.lep_Type < 0
  data_same =  df_da.lep_Type < 0
  mc_diff = df.lep_Type > 0
  data_diff = df_da.lep_Type > 0
  
  if flavor == "same":
    mc_leptons_pt =  mc_leptons_pt & mc_same
    data_leptons_pt =  data_leptons_pt & data_same
  if flavor == "diff":
    mc_leptons_pt =  mc_leptons_pt & mc_diff
    data_leptons_pt =  data_leptons_pt & data_diff

  #it works df is a pointer to an object, there was no implicit copy. If the values pointed to change the function will see it
  def masks(mc_leptons_pt):
    mc_lepton_Z =  (df.mll < 76) |  (df.mll > 106)

    mc_same =  df.lep_Type < 0
    data_same =  df_da.lep_Type < 0
    mc_diff = df.lep_Type > 0
    data_diff = df_da.lep_Type > 0

    if flavor == "same":
      mc_lepton_Z =  mc_lepton_Z & mc_same
    if flavor == "diff":
      mc_lepton_Z =  mc_lepton_Z & mc_diff

    mc_leptons_pt = (df.lep1_pt > 27) & (df.lep2_pt > 25) & (df.mll > 30) & mc_lepton_Z 
    return mc_leptons_pt 

  no_change = cross_calc(df[mc_leptons_pt ], df_da[data_leptons_pt], df_ww, scales,  flavor=flavor, fiducial=True)

  #Scale shift and redo pre calculation
  lep_scale_shift(df, rf=rfs)
  mc_leptons_pt = masks( mc_leptons_pt )

  change = cross_calc(df[mc_leptons_pt ], df_da[data_leptons_pt ], df_ww, scales, flavor=flavor, fiducial=True)

  #print "Lepton scale up: ",  no_change, change, no_change - change
  up = abs(no_change - change)

  #Shift leptons back to starting value and redo pre calc
  lep_scale_shift(df, lep_pt=1. / 1.003, rf=rfs)
  mc_leptons_pt = masks( mc_leptons_pt )
  #First calc
  no_change = cross_calc(df[mc_leptons_pt ], df_da[data_leptons_pt ], df_ww, scales, flavor=flavor, fiducial=True)

  #Scale shift to low end
  lep_scale_shift(df, lep_pt=1. / 1.003, rf=rfs)
  mc_leptons_pt = masks( mc_leptons_pt )

  #Pre calc and redo cross calc
  change = cross_calc(df[mc_leptons_pt], df_da[data_leptons_pt], df_ww, scales, flavor=flavor, fiducial=True)
  #print "Lepton scale down: ",  no_change, change, no_change - change
  down = abs(no_change - change)

  #del df
  #del df_da
  return (up + down) / 2.
###############################################################
###############################################################

###############################################################
#Jet pt Scale
###############################################################
def jet_scale_shift(data, jet_pt=1.025, rf=None): #2.5% jet scale
  #jets
  data.HT = data.HT * jet_pt
  data.jet1_pt = data.jet1_pt * jet_pt

  #MET results
  data.metMod  = data.metMod - data.HT / jet_pt * (jet_pt - 1)
  data.metProj = data.METProj - data.HT / jet_pt * (jet_pt - 1)
  data.recoil  = data.recoil - data.HT / jet_pt * (jet_pt - 1) 
  if rf != None:
    #print "Recreating random forest scores."
    pred_fTT = rf["clf_fTT"].predict_proba(np.float32(data[rf["features_fTT"]].values))
    data["pred_fTT_WW"] = pred_fTT[:,0]

    temp = data[rf["features_fDY"]]
    temp = temp.replace([np.inf,-np.inf], 0)
    pred_fDY = rf["clf_fDY"].predict_proba(np.float32(temp.values))
    data["pred_fDY_WW"] = pred_fDY[:,0]



def jet_pt_scale( X_section, flavor, verbose=False):
  """
  We define effect of jet momentum scale uncertainties for WW as:

  $$\bar{\epsilon} = \frac{\epsilon_{0} N_{0} + \epsilon_{1} N_{1} + \epsilon_{2} N_{2} + \ldots} { N_{0} + N_{1} + N_{2} + \ldots}$$

  Where $\epsilon_{i}$ is the efficiency of the WW after we apply all analysis cuts.
  """
  #########################################################
  #########################################################
  def kill_jets( df, pt_cut= 31 ):
      #Edit number of jets per event
      n_jets = np.zeros(df.shape[0])
      for k in df.keys():
          if "jet" in k and "pt" in k:
              cut = (df[k] > pt_cut)
              n_jets[cut.values] = n_jets[cut.values] + 1 

      df["numb_jets"] = n_jets


  def jet_scale_shift(data, jet_pt=1.025, pt_cut= 31, rf=None): #2.5% jet scale
      #jets
      data.HT = data.HT * 0
      #Scale pt of each jet
      for k in df.keys():
          if "jet" in k and "pt" in k:
              data[k] = data[k] * jet_pt
              #NEW TO CORRECT FOR EVENTS WITH LOST JETS
              ht_lost_jet = data[k] >= pt_cut
              data.HT.values[ht_lost_jet] = data[ht_lost_jet].HT + data[ht_lost_jet][k]

      #MET results
      data.metMod  = data.metMod - data.HT / jet_pt * (jet_pt - 1)
      data.metProj = data.METProj - data.HT / jet_pt * (jet_pt - 1)
      data.recoil  = data.recoil - data.HT / jet_pt * (jet_pt - 1)
      
      #Update number of jets
      kill_jets( data, pt_cut )
      if rf != None:
      #print "Recreating random forest scores."
          pred_fTT = rf["clf_fTT"].predict_proba(np.float32(data[rf["features_fTT"]].values))
          data["fTT"] = pred_fTT[:,0]

          temp = data[rf["features_fDY"]]
          temp = temp.replace([np.inf,-np.inf], 0)
          pred_fDY = rf["clf_fDY"].predict_proba(np.float32(temp.values))
          data["fDY"] = pred_fDY[:,0]

  def rf_ana(df, fDY, fTT):
      fDY = df[fDY] > .9
      fTT = df[fTT] > .6
      
      return df[fDY & fTT]
  #########################################################
  #########################################################


  df = load_presel_w_fDY_fTT_MC()
  df_da = load_presel_w_fDY_fTT_DATA()
  df_ww = rp.read_root(data_path+"/ww_complete.root", columns=columns)

  if flavor == "diff":
    df    = df[df.lep_Type > 0] 
    df_da = df_da[df_da.lep_Type > 0]
    df_ww = df_ww[df_ww.lep_Type > 0] 
  elif flavor == "same":
    df    = df[df.lep_Type < 0] 
    df_da = df_da[df_da.lep_Type < 0]
    df_ww = df_ww[df_ww.lep_Type < 0] 
    


  rfs = load_randomForest()

  df["org_numb_jets"] = 1 * df.numb_jets.values
  df_da["org_numb_jets"] = 1 * df_da.numb_jets.values


  #With preselection cuts, org rf
  pre_sel_ww = np.histogram(df_ww[df_ww.process == "WW"].numb_jets.values, bins=7, range=(-0.5, 6.5), )[0].astype(float)
  pre_rfsel_ww = np.histogram(rf_ana(df[df.process == "WW"], "pred_fDY_WW", "pred_fTT_WW").numb_jets.values, bins=7, range=(-0.5, 6.5), )[0].astype(float)

  #Apply shift With shift preselection, with rf
  jet_scale_shift(df, jet_pt=0.975, pt_cut= 30, rf=rfs)
  pre_shift_rfsel_low = np.histogram(rf_ana(df[df.process == "WW"], "fDY", "fTT").numb_jets.values, bins=7, range=(-0.5, 6.5), )[0].astype(float)

  #Apply shift With shift preselection, with rf
  jet_scale_shift(df, jet_pt=1.025, pt_cut= 31, rf=rfs)
  jet_scale_shift(df, jet_pt=1.025, pt_cut= 31, rf=rfs)
  pre_shift_rfsel_high = np.histogram(rf_ana(df[df.process == "WW"], "fDY", "fTT").numb_jets.values, bins=7, range=(-0.5, 6.5), )[0].astype(float)

  #print pre_rfsel_ww.sum(), pre_sel_ww.sum(), pre_shift_rfsel_low.sum(),  pre_shift_rfsel_low.sum()
  epsilon_low = ( ( pre_rfsel_ww / pre_sel_ww ) * pre_shift_rfsel_low ).sum() / pre_shift_rfsel_low.sum()
  epsilon_high = ( ( pre_rfsel_ww / pre_sel_ww ) * pre_shift_rfsel_high ).sum() / pre_shift_rfsel_high.sum()
  epsilon_orig = ( ( pre_rfsel_ww / pre_sel_ww ) * pre_shift_rfsel_high ).sum() / pre_rfsel_ww.sum()


  if verbose == True:
    print "\t epsilon \t bar epsilon \t N jets"
    print "Low:  ", epsilon_low, abs(epsilon_low - epsilon_orig) / epsilon_orig, [i for i in pre_shift_rfsel_low[:4]]
    print "High: ", epsilon_high, abs(epsilon_high - epsilon_orig) / epsilon_orig, [i for i in pre_shift_rfsel_high[:4]]
    print "Orig: ", epsilon_orig
    print "Tot jet scale unc: ", (abs(epsilon_low - epsilon_orig) + abs(epsilon_high - epsilon_orig) )/ (2. * epsilon_orig )

  tot_jet_scale_unc = (abs(epsilon_low - epsilon_orig) + abs(epsilon_high - epsilon_orig) )/ (2. * epsilon_orig )
  return tot_jet_scale_unc * X_section   



###############################################################
###############################################################
#Just changed from .017 to .037 Aug 2 changed back
def lep_resolution(data, lep_resolution=.017, rf=None): #.006 is the muon resolution, 1.7 - 4.5% for electrons (I can't delieve the electrons are that big)
    #Leptons 1 & 2 
    lep1_delta = np.random.normal(loc=[0]*data.shape[0], scale=lep_resolution * data.lep1_pt.values)
    lep2_delta = np.random.normal(loc=[0]*data.shape[0], scale=lep_resolution * data.lep2_pt.values)

    #print lep1_delta[:10]
    #print data.lep1_pt.values[:10]
    data.lep1_pt = data.lep1_pt + lep1_delta
    data.lep2_pt = data.lep2_pt + lep2_delta
    #data.mll = 

    #MET
    data.metMod = data.metMod - lep1_delta - lep2_delta
    #data.recoil = data.recoil + 2*lep1_delta + 2*lep1_delta
    if rf != None:
      #print "Recreating random forest scores."
      pred_fTT = rf["clf_fTT"].predict_proba(np.float32(data[rf["features_fTT"]].values))
      data["pred_fTT_WW"] = pred_fTT[:,0]

      temp = data[rf["features_fDY"]]
      temp = temp.replace([np.inf,-np.inf], 0)
      pred_fDY = rf["clf_fDY"].predict_proba(np.float32(temp.values))
      data["pred_fDY_WW"] = pred_fDY[:,0]
    return lep1_delta, lep2_delta


def lep_pt_resolution( flavor = 'both'):
  df = load_presel_w_fDY_fTT_MC()
  df_da = load_presel_w_fDY_fTT_DATA()
  df_ww = rp.read_root(data_path+"/ww_complete.root", columns=columns)
  rfs = load_randomForest()

  data_lepton_Z =  (df_da.mll < 76) |  (df_da.mll > 106)
  data_leptons_pt = (df_da.lep1_pt > 27) & (df_da.lep2_pt > 17) & (df_da.mll > 30) & data_lepton_Z

  mc_same =  df.lep_Type < 0
  data_same =  df_da.lep_Type < 0
  mc_diff = df.lep_Type > 0
  data_diff = df_da.lep_Type > 0

  if flavor == "same":
    data_leptons_pt =  data_leptons_pt & data_same
  elif flavor == "diff":
    data_leptons_pt =  data_leptons_pt & data_diff
  else:
    print "No flavor match. Will proceed with all flavors."


  def lep_masks( df ):
      mc_lepton_Z =  (df.mll < 76) |  (df.mll > 106)
      if flavor == "same":
        mc_lepton_Z =  mc_lepton_Z & mc_same
      if flavor == "diff":
        mc_lepton_Z =  mc_lepton_Z & mc_diff

      mc_leptons_pt = (df.lep1_pt > 27) & (df.lep2_pt > 17) & (df.mll > 30) & mc_lepton_Z
      return mc_leptons_pt

  mc_leptons_pt = lep_masks( df )
  no_change = cross_calc(df[mc_leptons_pt ], df_da[data_leptons_pt],  df_ww, scales, fiducial=True)

  #Scale shift and redo pre calculation
  a = lep_resolution(df, .003, rf=rfs)
  mc_leptons_pt = lep_masks( df )

  change = cross_calc(df[mc_leptons_pt], df_da[data_leptons_pt],  df_ww, scales, fiducial=True)

  print "Lepton resolution up: ",  no_change, change, no_change - change
  del df_ww
  del rfs
  return abs(no_change - change)

###############################################################
###############################################################


###############################################################
###############################################################
#?Should be using this
#?https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
def jet_resolution(data, jet_scale=.1, rf = None):
    #Jets
    ht_delta =np.random.normal(loc=[0]*data.shape[0], scale=jet_scale * data.jet1_pt) + np.random.normal(loc=[0]*data.shape[0], scale=jet_scale * data.jet2_pt)#np.random.normal(loc=[0]*data.shape[0], scale=jet_scale ) * np.sqrt(data.numb_jets)
    data.HT = data.HT + ht_delta 
    data.jet1_pt = data.jet1_pt + np.random.normal(loc=[0]*data.shape[0], scale=jet_scale * data.jet1_pt)
    data.jet2_pt = data.jet2_pt + np.random.normal(loc=[0]*data.shape[0], scale=jet_scale * data.jet1_pt)
    #MET
    data.metMod = data.metMod - ht_delta
    #Update leptons and jets
    #rerun rf score calculation
    if rf != None:
      #print "Recreating random forest scores."
      pred_fTT = rf["clf_fTT"].predict_proba(np.float32(data[rf["features_fTT"]].values))
      data["pred_fTT_WW"] = pred_fTT[:,0]
      pred_fDY = rf["clf_fDY"].predict_proba(np.float32(data[rf["features_fDY"]].values))
      data["pred_fDY_WW"] = pred_fDY[:,0]


def jet_pt_resolution( flavor ):
  df = load_presel_w_fDY_fTT_MC()
  df_da = load_presel_w_fDY_fTT_DATA()
  df_ww = rp.read_root(data_path+"/ww_complete.root", columns=columns)
  rfs = load_randomForest()

  if flavor == "diff":
    df    = df[df.lep_Type > 0] 
    df_da = df_da[df_da.lep_Type > 0]
    df_ww = df_ww[df_ww.lep_Type > 0] 
  elif flavor == "same":
    df    = df[df.lep_Type < 0] 
    df_da = df_da[df_da.lep_Type < 0]
    df_ww = df_ww[df_ww.lep_Type < 0] 


  mc_jets_pt = ((df.jet1_pt > 31) | (df.jet1_pt < 5)) & ((df.jet2_pt > 31 )| (df.jet2_pt < 5)) 
  data_jets_pt = ((df_da.jet1_pt > 31) | (df_da.jet1_pt < 5)) & ((df_da.jet2_pt > 31 )| (df_da.jet2_pt < 5)) 
  def masks_jet( mc_leptons_jet):
      mc_jets_pt = ((df.jet1_pt > 31) | (df.jet1_pt < 5)) & ((df.jet2_pt > 31 )| (df.jet2_pt < 5)) 
      return mc_jets_pt

  flavor = 'both'
  mc_jets_pt = masks_jet(mc_jets_pt)
  no_change = cross_calc(df[mc_jets_pt], df_da[data_jets_pt], df_ww, scales, flavor=flavor, fiducial=True) 
  jet_resolution(df_da, rf=rfs)
  mc_jets_pt = masks_jet(mc_jets_pt)
  change = cross_calc(df[mc_jets_pt], df_da[data_jets_pt], df_ww, scales, flavor=flavor, fiducial=True)

  #print "Jet resolution: ",  no_change, change, no_change - change
  return abs(no_change - change)

###############################################################
###############################################################

###############################################################
###############################################################

def pileup_unc( flavor='both', flag=1):
  df = load_presel_w_fDY_fTT_MC()
  df_da = load_presel_w_fDY_fTT_DATA()
  df_ww = rp.read_root(data_path+"/ww_complete.root", columns=columns)

  orig = cross_calc(df, df_da, df_ww, scales, fiducial=True, flavor=flavor) #REMOVE
  difs = []

  if flag == 0:
    df = df[df.metFilter_flag ==0]
    df_da = df_da[df_da.metFilter_flag ==0]
    df_ww = df_ww[df_ww.metFilter_flag ==0]
    

  #print "orig", orig
  for i in [ pile_up.Min_Bias.HIGH, pile_up.Min_Bias.MID, pile_up.Min_Bias.LOW]:

    gw = np.ones(df_ww.shape[0])
    gw *= df_ww.gen_weight.values
    l = list(df_ww.gen_pu.values)
    w = np.array(pile_up.pileUpFunction(l, i, flag))
    df_ww["w_alt"] = w * gw


    gw = np.ones(df.shape[0])
    gw *= df.gen_weight.values
    l = list(df.gen_pu.values)
    w = np.array(pile_up.pileUpFunction(l, i, flag))
    df["w_alt"] = w * gw


    kwargs = calc_cross_stuff(df, df_da, df_ww, scales, weights="w_alt", flavor=flavor)
    #print i, cross_calc(df, df_da, df_ww, scales, fiducial=True, **kwargs)

    if i != pile_up.Min_Bias.MID: difs.append( cross_calc(df, df_da, df_ww, scales, fiducial=True, **kwargs))
    else: orig = cross_calc(df, df_da, df_ww, scales, fiducial=True, **kwargs)

#  print orig, difs
  return sum([ (i - orig)**2 for i in difs])**.5

###############################################################
###############################################################
def background_uncertainty( flavor = 'both'):
  df = load_presel_w_fDY_fTT_MC()
  df_da = load_presel_w_fDY_fTT_DATA()
  df_ww = rp.read_root(data_path+"/ww_complete.root", columns=columns)


  normalization_unc = normalization_unc_calc(df, df_da, df_ww, scales, unc_mc_process, fiducial=True, flavor=flavor)
  #print "Systematics ", normalization_unc, normalization_unc/cross_calc(df, df_da, df_ww, scales, fiducial=True, flavor=flavor), cross_calc(df, df_da, df_ww, scales, fiducial=True, flavor=flavor)  
  return normalization_unc


###############################################################
###############################################################


def statistical( flavor = 'both'):
  df    = load_presel_w_fDY_fTT_MC()
  df_da = load_presel_w_fDY_fTT_DATA()
  df_ww = rp.read_root(data_path+"/ww_complete.root", columns=columns)

  stat_unc = stat_unc_calc(rf_ana(df), rf_ana(df_da), df_ww, scales, unc_mc_process, fiducial=True, flavor=flavor)
  #print "Systematics ", stat_unc, stat_unc/cross_calc(df, df_da, df_ww, scales, fiducial=True, flavor=flavor), cross_calc(df, df_da, df_ww, scales, fiducial=True, flavor=flavor)  
  return stat_unc

###############################################################
###############################################################
def ww_unc(flavor):
  from sklearn.externals import joblib
  df    = load_presel_w_fDY_fTT_MC()
  df_da = load_presel_w_fDY_fTT_DATA()

  random_forests = joblib.load("../data/RF/aug_05_fDY_fTT.jbl")
  features_fDY = random_forests["features_fDY"] 
  clf_fDY = random_forests["clf_fDY"]
  features_fTT = random_forests["features_fTT"] 
  clf_fTT = random_forests["clf_fTT"]

  df_ww = rp.read_root(data_path+"/ww_complete.root", columns=columns)
  
  pred_fTT = clf_fTT.predict_proba(np.float32(df_ww[features_fTT].values))
  df_ww["pred_fTT_WW"] = pred_fTT[:,0]

  pred_fDY = clf_fDY.predict_proba(np.float32(df_ww[features_fDY].values))
  df_ww["pred_fDY_WW"] = pred_fDY[:,0]

  #that_dic = pkl.load(open("../data/WWpt/8TeV_nnnlo.pkl", "r"))
  resumm_weights = {}
#  resumm_weights["central"] = pd.read_csv("../data/WWpt/resummation_central_weights.csv", index_col=False)
#  resumm_weights["resum_up"] = pd.read_csv("../data/WWpt/resummation_resum_up_weights.csv", index_col=False )
#  resumm_weights["resum_down"] = pd.read_csv("../data/WWpt/resummation_resum_down_weights.csv", index_col=False)
#  resumm_weights["scale_up"] = pd.read_csv("../data/WWpt/resummation_scale_up_weights.csv", index_col=False )
#  resumm_weights["scale_down"] = pd.read_csv("../data/WWpt/resummation_scale_down_weights.csv", index_col=False)
  resumm_weights["central"] = pd.read_csv("../data/WWpt/resummation_central_weights.csv", index_col=False)
  resumm_weights["resum_up"] = pd.read_csv("../data/WWpt/resummation_resum_up_weights.csv", index_col=False )
  resumm_weights["resum_down"] = pd.read_csv("../data/WWpt/resummation_resum_down_weights.csv", index_col=False)
  resumm_weights["scale_up"] = pd.read_csv("../data/WWpt/resummation_scale_up_weights.csv", index_col=False )
  resumm_weights["scale_down"] = pd.read_csv("../data/WWpt/resummation_scale_down_weights.csv", index_col=False)
  temp_df = df_ww[df_ww.ww_pt > -1]

  #print "tot shape", temp_df.shape
  temp_df["ww_weights"] = [1.] * temp_df.shape[0]

  orig = 0
  q_ = []
  r_ = []


  for k in resumm_weights:
    for it, row in resumm_weights[k].iterrows():
      temp_df.ww_weights.values[ (temp_df.ww_pt >= row["low"]).values & (temp_df.ww_pt < row["high"]).values]  = row["weight"]

      #if it == 5:
      #  print k, row["low"], row["high"],row["weight"], temp_df.ww_weights.values[ (temp_df.ww_pt >= row["low"]).values & (temp_df.ww_pt < row["high"]).values].shape

    rf_df = rf_ana(temp_df)
    if flavor == "diff":
      flavor_cut = rf_df.lep_Type > 0 
    elif flavor == "same":
      flavor_cut = rf_df.lep_Type < 0 
    else:
      flavor_cut = rf_df.lep_Type > -99

    if k == "central":
      orig = rf_df[flavor_cut].ww_weights.sum()/temp_df.ww_weights.sum()
      #print "Weights diff correction", rf_ana(temp_df).weight.sum(), (rf_ana(temp_df).weight.values * rf_ana(temp_df).ww_weights.values).sum()
    if "resum" in k:
      r_.append(rf_df[flavor_cut].ww_weights.sum()/temp_df.ww_weights.sum())
    if "scale" in k :
      q_.append(rf_df[flavor_cut].ww_weights.sum()/temp_df.ww_weights.sum())



  """
  orig = 0
  q_ = []
  r_ = []

  ###
  for two in range(5):
      temp_df["ww_weights"] = [that_dic[8][two][3][1][-1]] * temp_df.shape[0]
      for it, i in enumerate(that_dic[8][two][3][0]):
          if it+1 < len(that_dic[8][two][3][0]): 
              temp_df.ww_weights.values[(temp_df.ww_pt >= i) & (temp_df.ww_pt < that_dic[8][two][3][0][it+1])] = that_dic[8][two][3][1][it]
      #print two, rf_ana(temp_df).ww_weights.sum()/temp_df.ww_weights.sum()
      if two == 0:
          orig = rf_ana(temp_df).ww_weights.sum()/temp_df.ww_weights.sum()
      elif two > 0 and two < 3:
          q_.append(rf_ana(temp_df).ww_weights.sum()/temp_df.ww_weights.sum())
      elif two >= 3 and two < 5:
          r_.append(rf_ana(temp_df).ww_weights.sum()/temp_df.ww_weights.sum())
      else:
          print two, "We got a problem"
  ####
  """

  #print "orig", orig
  #print "r_", r_
  #print "q_", q_
  r_ = np.array([abs(orig-i) for i in r_]).mean()
  q_ = np.array([abs(orig-i) for i in q_]).mean()
  ww_pt_unc = (q_**2 + r_**2)**.5

#  print temp_df.weight.sum()/rf_ana(temp_df).weight.sum(), orig, cross_calc(df, df_da, df_ww, scales, fiducial=True), q_, r_ 
#
#  print "ww pt r: cross unc", cross_calc(df, df_da, df_ww, scales, fiducial=True) * temp_df.weight.sum()/rf_ana(temp_df).weight.sum() * r_, "%: ", temp_df.weight.sum()/rf_ana(temp_df).weight.sum() * r_ * 100.
#
#  print "ww_pt q: cross unc", cross_calc(df, df_da, df_ww, scales, fiducial=True) * temp_df.weight.sum()/rf_ana(temp_df).weight.sum() * q_, "%: ", temp_df.weight.sum()/rf_ana(temp_df).weight.sum() * q_ * 100. 
#
#  print "ww pt total: cross unc", cross_calc(df, df_da, df_ww, scales, fiducial=True) * temp_df.weight.sum()/rf_ana(temp_df).weight.sum() * ww_pt_unc, "%: ", temp_df.weight.sum()/rf_ana(temp_df).weight.sum() * ww_pt_unc * 100.
  return cross_calc(df, df_da, df_ww, scales, fiducial=True) * temp_df.weight.sum()/rf_ana(temp_df).weight.sum() * ww_pt_unc



###############################################################
###############################################################
def main():

  flavor = "both"

  df = load_presel_w_fDY_fTT_MC()
  df_da = load_presel_w_fDY_fTT_DATA()
  df_ww = rp.read_root(data_path+"/ww_complete.root", columns=columns)

  kwargs = calc_cross_stuff(df, df_da, df_ww, scales, flavor=flavor)
  X_section = cross_calc(df, df_da, df_ww, scales, fiducial=True, **kwargs)
  ratio_s_t = kwargs["ratio_s_t"]
  print "Cross-section", X_section, "WW Efficiency", ratio_s_t

  del df
  del df_da
  del df_ww

  def print_unc( label, eval_func ):
    print label, eval_func, eval_func / X_section, eval_func / X_section * 100 

  print "\tunc\t\tunc/x_sec\t\t%"
#  for flavor in ["both", "diff", "same"]:
  print_unc("background uncertainty",  background_uncertainty( flavor = flavor))#Same Flavor is large 
#  print_unc("lepton scale factor",  lepton_scale_factor(X_section, flavor=flavor, ww_eff=ratio_s_t))
#  print_unc("lepton pt scale",      lep_pt_scale(flavor=flavor))
#  print_unc("lep pt resolution",    lep_pt_resolution(flavor=flavor))
#
#  print_unc("bjet scale factor",    bjet_scale_factor(X_section, flavor))
#  print_unc("jet pt scale",         jet_pt_scale( X_section, flavor ))
#  print_unc("jet pt resolution",    jet_pt_resolution( flavor ))
#
#  print_unc("pileup unc",           pileup_unc(flavor, 0))
#  print_unc("Theory ww",            ww_unc( flavor )) #Many problems for 13 TeV flavors not right among others.
#  print_unc("Statistical",          statistical( flavor=flavor))
#




if __name__ == '__main__':
  main()



