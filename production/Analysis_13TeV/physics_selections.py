# Thoth Gunter
import numpy as np
import pandas as pd


def pre_cuts( df, diff_charge= True):
  """
    Initial selection apply before preforming analysis. Initial selection includes a same sign Z peak cut, 3 or more lepton cut, b jet cut and invariant mass cut of 30.  
    pre_cuts( df, diff_charge)
    df: data frame
    diff_charge: bool (Apply opposite sign charge cut)
  """
  dif_lep = df.lep_Type > 0
  sam_lep = df.lep_Type < 0
  z_mass = (df.mll < 74) | (df.mll > 106 )
  nBjet = df.numb_BJet == 0
  extra_lep = df.numbExtraLep == 0
  quality_cuts = df.mll > 30 
  charge = df.lep1_Charge != df.lep2_Charge
  lep_pt = (df.lep1_pt > 25) & (df.lep2_pt > 25)

  if diff_charge == False:
    s_df = df[ sam_lep & z_mass & nBjet & extra_lep & quality_cuts & lep_pt]
    d_df = df[ dif_lep & nBjet & extra_lep & quality_cuts & lep_pt]
  
  if diff_charge == True:
    s_df = df[ sam_lep & z_mass & nBjet & extra_lep & quality_cuts & charge & lep_pt]
    d_df = df[ dif_lep & nBjet & extra_lep & quality_cuts & charge & lep_pt]  
  
  return pd.concat([ s_df, d_df])


def cuts_ana( df, flavor="both" ):
  """
    Indepth WW cut selection analysis.
    cuts_ana( df ) 
    df: data frame
    flavor: "both", "same", "diff"
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


def full_ana( df, diff_charge=True ):
  return df

def same_ana( df, diff_charge= True ):
  return df[df.lep_Type < 0]

def diff_ana( df, diff_charge= True ):
  return df[ df.lep_Type > 0]


def WW_ana( df , diff_charge= True): 
  """
    Basic cuts anaylsis
    WW_ana( df , diff_charge= True)
    df: data frame
    diff_charge: bool (Apply opposite sign charge cut)
  """
  initial_cuts = (df.mll > 50 )  & (df.numbExtraLep == 0) & (df.numb_jets ==  0) & ( df.metMod > 50 )
  dif_flavor_cut = (df.lep_Type > 0)

  if diff_charge == True:
    initial_cuts = initial_cuts & (df.lep1_Charge != df.lep2_Charge)

  basic_df_0j_cuts = initial_cuts & dif_flavor_cut 
  return df[basic_df_0j_cuts]


def TT_ana( df, diff_charge= True ):
  """
    Basic TTbar/Single Top analysis
    TT_ana( df , diff_charge= True)
    df: data frame
    diff_charge: bool (Apply opposite sign charge cut)
  """ 
  initial_cuts = (df.mll > 50 ) & (df.qT > 30) & (df.metMod > 60) & (df.numbExtraLep == 0) & (df.numb_jets > 0)
  dif_flavor_cut = (df.lep_Type > 0)
  
  basic_df_0j_cuts = initial_cuts & dif_flavor_cut 
  return df[basic_df_0j_cuts] 

def DY_ana( df, diff_charge = True ): 
  """
    Basic Drell Yan analysis
    DY_ana( df , diff_charge= True)
    df: data frame
    diff_charge: bool (Apply opposite sign charge cut)
  """ 
  initial_cuts = (df.mll > 50 )  & (df.numbExtraLep == 0) & (df.numb_jets == 0) & (df.metMod < 80)

  return df[ initial_cuts ]

def Z_tt_ana( df, diff_charge = True ): 
  """
    Basic Z to tau tau analysis
    Z_tt_ana( df , diff_charge= True)
    df: data frame
    diff_charge: bool (Apply opposite sign charge cut)
  """ 
  initial_cuts = (df.mll > 50 )  & (df.numbExtraLep == 0) & (df.numb_jets ==  0) & ( df.qT < 10 ) & (df.lep1_pt > 30) & (df.lep2_pt > 20)
  
  dif_flavor_cut = (df.lep_Type > 0)
  
  m_1 = (df.mll < 120) & (df.mll > 60)
  met_1 = (df.metMod > 30) & (df.metMod < 50)

  return df[ initial_cuts & met_1 & m_1 & dif_flavor_cut ] 

def wjets_selection(df_mc, df_data, feature):
  """
    Basic wjets analysis
    wjets_selection( df , diff_charge= True)
    df: data frame
    diff_charge: bool (Apply opposite sign charge cut)
  """ 
  bins_mc_ = bin_df( df_mc[df_mc.lep1_Charge == df_mc.lep2_Charge], feature, )#weights=False) 
  bins_data_ = bin_df( df_data[df_data.lep1_Charge == df_data.lep2_Charge], feature, )

  sum_mc = np.zeros(bins_mc_[bins_mc_.keys()[0]][0].shape[0])
  for k in bins_mc_.keys():
      if k in ['ttbar_semileptonic', 'ttbar_leptonic', 'DYJetsToLL_M-50', 'WZJetsTo3LNu', 'WW',
   'Tbar_tW-channel','T_t-channel','ZZJetsTo2L2Nu', 'WGToLNuG']:
          sum_mc += bins_mc_[k][0]

  sum_mc = bins_data_["Da"][0] - sum_mc
  return sum_mc




def rf_ana(df, flavor="both"):
  """
    cuts:
    df.pred_fDY_WW > .9
    df.pred_fTT_WW > .6
  """
  fDY = df.pred_fDY_WW > .9
  fTT = df.pred_fTT_WW > .6

  if flavor == "same":
    flv = df.lep_Type < 0
  elif flavor == "diff":
    flv = df.lep_Type > 0
  else:
    flv = df.lep_Type > -99

  return df[fDY & fTT & flv]





def same_sign_ana(df):
  
  return df[df.lep1_Charge == df.lep2_Charge]


