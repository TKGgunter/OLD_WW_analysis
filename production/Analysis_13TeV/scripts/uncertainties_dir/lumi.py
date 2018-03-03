#Thoth Gunter

import os, sys, pickle, copy
import datetime
sys.path.append(os.getcwd() + "/../../")
from prep_ana_II import *
sys.path.append(os.getcwd() + "/../../tools/")
sys.path.append(os.getcwd() + "/../../tools/fit/")
sys.path.append(os.getcwd() + "/../../tools/jecsys/")
from jet_scale import jec_setup, jet_scale_shift_flat, pseudo_data_yield_sum 
from cross_section_calc import calc_cross_stuff, cross_calc, stat_unc_calc, normalization_unc_calc
import fit
import warnings
warnings.filterwarnings("ignore")



lumi_unc = 0.025


def computate_lumi(flavor):
  """
  Compute luminosity uncertainty.

  This is how we compute the lumi uncertainty:
      df weights * (1 +/- 0.025) 
  """
  #Setup analysis frame work
  ana_obj = analysis_setup()
  ana_obj.apply_pre_cuts() 
  #ana_obj.apply_flat_jet_correction() 

  df = ana_obj.df
  df_da = ana_obj.df_da
  df_ww = ana_obj.df_ww
  df_ggww = ana_obj.df_ggww
  scales = ana_obj.scales

  if flavor == "diff":
    df =    df[   df.lep1_type    != df.lep2_type]
    df_da = df_da[df_da.lep1_type != df_da.lep2_type]
    df_ww = df_ww[df_ww.lep1_type != df_ww.lep2_type]
    df_ggww = df_ggww[df_ggww.lep1_type != df_ggww.lep2_type]
  if flavor == "same":
    df =    df[   df.lep1_type    == df.lep2_type]
    df_da = df_da[df_da.lep1_type == df_da.lep2_type]
    df_ww = df_ww[df_ww.lep1_type == df_ww.lep2_type]
    df_ggww = df_ggww[df_ggww.lep1_type == df_ggww.lep2_type]


  pseudo = {}
  pseudo["tot"] = pseudo_data_yield_sum(rf_ana(df), rf_ana(df_da), scales=scales)

  orig = {"tot":0, "j0":0, "j1":0, "j2":0}
  orig["tot"] =  cross_calc(df, df_da, df_ww, df_ggww, scales,  fiducial=True, pseudo=pseudo["tot"])  


  #Up Variation
  up_down = 1 
  df["weight"] = df.weight.values *  (1 + up_down * lumi_unc)
  df_ww["weight"] = df_ww.weight.values *  (1 + up_down * lumi_unc)
  df_ggww["weight"] = df_ggww.weight.values *  (1 + up_down * lumi_unc)
  up = {"tot":0, "j0":0, "j1":0, "j2":0}
  up["tot"] =  cross_calc(df, df_da, df_ww, df_ggww, scales,  fiducial=True, pseudo=pseudo["tot"]) * (1 + up_down * lumi_unc)**-1



  #Down Variation
  ana_obj.reset() 
  ana_obj.apply_pre_cuts() 
  #ana_obj.apply_flat_jet_correction() 

  df = ana_obj.df
  df_da = ana_obj.df_da
  df_ww = ana_obj.df_ww


  if flavor == "diff":
    df =    df[   df.lep1_type    != df.lep2_type]
    df_da = df_da[df_da.lep1_type != df_da.lep2_type]
    df_ww = df_ww[df_ww.lep1_type != df_ww.lep2_type]
  if flavor == "same":
    df =    df[   df.lep1_type    == df.lep2_type]
    df_da = df_da[df_da.lep1_type == df_da.lep2_type]
    df_ww = df_ww[df_ww.lep1_type == df_ww.lep2_type]

  up_down = -1 
  df["weight"] = df.weight.values * (1 + up_down * lumi_unc)
  df_ww["weight"] = df_ww.weight.values * (1 + up_down * lumi_unc)
  df_ggww["weight"] = df_ggww.weight.values * (1 + up_down * lumi_unc)

  down = {"tot":0, "j0":0, "j1":0, "j2":0}
  down["tot"] =  cross_calc(df, df_da, df_ww, df_ggww, scales,  fiducial=True, pseudo=pseudo["tot"]) * (1 + up_down * lumi_unc)**-1 


  print orig
  print up
  print down
  print "unc: ", (abs(up["tot"] - orig["tot"]) + abs(down["tot"] - orig["tot"])) / (2 * orig["tot"]) * 100.

  if flavor != "":
    flavor = "_" + flavor 
  f = open("results/jan/lumi" + flavor + ".txt", "w")
  f.write("Orig: " + str(orig)+"\n") 
  f.write("Up: "   + str(up)+"\n") 
  f.write("Down: " + str(down)+"\n") 
  f.write("unc: "  + str((abs(up["tot"] - orig["tot"]) + abs(down["tot"] - orig["tot"])) / (2 * orig["tot"]) * 100.))


if __name__ == "__main__":
  for flavor in ["", "same", "diff"]:
    print flavor
    computate_lumi(flavor)
