#Thoth Gunter

import os, sys, pickle, copy
import datetime
sys.path.append(os.getcwd() + "/../../")
sys.path.append(os.getcwd() + "/../../tools/")
sys.path.append(os.getcwd() + "/../../tools/fit/")
sys.path.append(os.getcwd() + "/../../tools/jecsys/")

from prep_ana_II import *
from jet_scale import jec_setup, jet_scale_shift_flat, pseudo_data_yield_sum 
from cross_section_calc import calc_cross_stuff, cross_calc, stat_unc_calc, normalization_unc_calc
from electron_idiso import yield_string
import fit
import warnings
warnings.filterwarnings("ignore")




"""
ToDo:
+?Review whatacturally gets done
+?Add comments to guide user
+?figure out what is the issue, if there is one
"""

def plots(df):
  """
  Maybe soon
  """
  return None


def muon_idiso_unc(df, up_down="up"):
  """
  muon_scale should scale the weights of samples with atleast one muon by the uncertainty on that muons scale factor
  So: w -> w * (1 +/- unc)
  """
  weights = np.ones(df.shape[0])
  weights *= df.weight.values

  cut_lep_type = df.lep_Type != -2  

  if up_down == "up":
    pm = 1
    error_string1 = "lep1_id_error_high"
    error_string2 = "lep2_id_error_high"
  else:
    pm = -1
    error_string1 = "lep1_id_error_low"
    error_string2 = "lep2_id_error_low"

  
  lep1_cuts = cut_lep_type & (df.lep1_type == 13)
  lep2_cuts = cut_lep_type & (df.lep2_type == 13)



  weights[lep1_cuts] =  weights[lep1_cuts] * (1 + pm * df[lep1_cuts][error_string1].values) 
  weights[lep2_cuts] =  weights[lep2_cuts] * (1 + pm * df[lep2_cuts][error_string2].values) 


  return weights
   
   

#Muon IsoID uncertainties
def compute_muon( ana_obj, flavor):
  rfs = ana_obj.rfs
  scales = ana_obj.scales
  df = pre_cuts(ana_obj.df, diff_charge= False)

  pseudo = {}
  pseudo["tot"] = pseudo_data_yield_sum(rf_ana(ana_obj.df), rf_ana(ana_obj.df_da), scales=scales)

  orig = {"tot":0, "j0":0, "j1":0, "j2":0}
  orig["tot"] =  cross_calc(ana_obj.df, ana_obj.df_da, ana_obj.df_ww, ana_obj.df_ggww, scales,  fiducial=True, pseudo=pseudo["tot"])  

  #########################
  #Fits
  nominal_fit_result = fit.comprehensive_fit(df, ana_obj.df_da, "metMod", scales)
  #########################


  print "orig", df.weight[df.lep_Type == -1].mean(), df[df.lep_Type == -1].weight.max(), df[df.lep_Type == -1].weight.min(),  df[df.lep_Type == -1].weight.std(), "why are these weights so small...Was expecting 90s avg"
  idiso_unc_weights = muon_idiso_unc(df)

  print "up shift", np.mean(idiso_unc_weights[df.lep_Type == -1]), np.std(idiso_unc_weights)  
  print "up errors", df.lep1_id_error_high.mean(), df.lep2_id_error_high.std()
  up_weights = idiso_unc_weights 




  idiso_unc_weights = muon_idiso_unc(df, up_down="down")
  print "\ndown shift", np.mean(idiso_unc_weights[df.lep_Type == -1]), np.std(idiso_unc_weights)  
  print "down errors", df.lep1_id_error_low.mean(), df.lep2_id_error_low.std()
  down_weights = idiso_unc_weights 



  print "\n\n"
  nominal_string, nominal_raw = yield_string(df[df.lep_Type != -2], "nominal\n", scales, rf_ana)
  df.weight = up_weights 
  up_string, up_raw = yield_string(df[df.lep_Type != -2], "up\n", scales, rf_ana)
  up = {"tot":0, "j0":0, "j1":0, "j2":0}
  up["tot"] =  cross_calc(df, ana_obj.df_da, ana_obj.df_ww, ana_obj.df_ggww, scales,  fiducial=True, pseudo=pseudo["tot"])  

  #########################
  #Fits
  up_fit_result = fit.comprehensive_fit(df, ana_obj.df_da, "metMod", scales)
  #########################



  #Calc down
  df.weight = down_weights 
  down_string, down_raw = yield_string(df[df.lep_Type != -2], "down\n", scales, rf_ana)
  down_string, down_raw = yield_string(df[df.lep_Type != -1], "down\n", scales, rf_ana)
  down = {"tot":0, "j0":0, "j1":0, "j2":0}
  down["tot"] =  cross_calc(df, ana_obj.df_da, ana_obj.df_ww, ana_obj.df_ggww, scales,  fiducial=True, pseudo=pseudo["tot"])  
  #########################
  #Fits
  down_fit_result = fit.comprehensive_fit(df, ana_obj.df_da, "metMod", scales)
  #########################

  print "Nominal", orig
  print "Up", up
  print "Down", down
  

  
  #########################
  print "Fit results"
  fit_processes = ["WW", "Top", "DY"] 
  print fit_processes
  print "nominal: ", nominal_fit_result.x
  print "up: ", up_fit_result.x
  print "down: ", down_fit_result.x
  print "Diffs(%)"
  for it, ele in enumerate(nominal_fit_result.x):
    print  fit_processes[it], (abs(ele - up_fit_result.x[it]) + abs(ele - down_fit_result.x[it]))/2. * 100
  #########################


  if flavor != "":
    flavor = "_" + flavor

  #Print results
  date = datetime.date.today() 
  f = open("results/jan/muon_isoid" + str(date.year) + "_" + str(date.month) + "_" + str(date.day) +  flavor + ".txt", "w")
  f.write(nominal_string)
  f.write("\n\n\n"+up_string)
  f.write(down_string)
  f.write("\n\n\nRatio between orignal process and altered\n")
  f.write("\nUP\n" + "\n".join([ process + ": " + str(up_raw[process] / nominal_raw[process])  for process in up_raw]) + "\n")
  f.write("\nDOWN\n" + "\n".join([  process + ": " + str(down_raw[process] / nominal_raw[process])  for process in down_raw]) + "\n") 
  f.write("\nCross Sections\n")
  f.write(str(orig) + "\n")
  f.write(str(up) + "\n")
  f.write(str(down) + "\n")
  f.write("Cross Sections Unc (%)\n")
  f.write("Up: "    + str([jet + " " + str((orig[jet] - up[jet]) / (orig[jet]) * 100.)  for jet in ["tot"]] ) + "\n")
  f.write("Down: "  + str([jet + " " + str((orig[jet] - down[jet]) / (orig[jet]) * 100.)  for jet in ["tot"]]) + "\n")
  f.write("Total: " + str([jet + " " + str((abs(orig[jet] - down[jet]) + abs(orig[jet] - up[jet])) / (2 * orig[jet]) * 100.)  for jet in ["tot"]]) + "\n")

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





if __name__ == "__main__":
  for flavor in ["",]:# "diff", "same"]:
    if flavor == "":
      print "total"
      ana_obj = analysis_setup("lep")
      ana_obj.apply_pre_cuts() 
      ana_obj.apply_flat_jet_correction() 
      compute_muon(ana_obj, flavor)
    elif flavor == "diff":
      print "diff"
      ana_obj = analysis_setup("lep")
      ana_obj.apply_pre_cuts() 
      ana_obj.apply_flat_jet_correction() 
      ana_obj.df = ana_obj.df[ana_obj.df.lep1_type != ana_obj.df.lep2_type]
      ana_obj.df_da = ana_obj.df_da[ana_obj.df_da.lep1_type != ana_obj.df_da.lep2_type]
      ana_obj.df_ww = ana_obj.df_ww[ana_obj.df_ww.lep1_type != ana_obj.df_ww.lep2_type]
      compute_muon(ana_obj, flavor)
    elif flavor == "same":
      print "same"
      ana_obj = analysis_setup("lep")
      ana_obj.apply_pre_cuts() 
      ana_obj.apply_flat_jet_correction() 
      ana_obj.df = ana_obj.df[ana_obj.df.lep1_type == ana_obj.df.lep2_type]
      ana_obj.df_da = ana_obj.df_da[ana_obj.df_da.lep1_type == ana_obj.df_da.lep2_type]
      ana_obj.df_ww = ana_obj.df_ww[ana_obj.df_ww.lep1_type == ana_obj.df_ww.lep2_type]
      compute_muon(ana_obj, flavor)










