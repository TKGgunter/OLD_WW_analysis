#Thoth Gunter

import os, sys, pickle, copy
import datetime
sys.path.append(os.getcwd() + "/../../")
from prep_ana_II import *
sys.path.append(os.getcwd() + "/../../tools/")
sys.path.append(os.getcwd() + "/../../tools/fit/")
sys.path.append(os.getcwd() + "/../../tools/jecsys/")
from jet_scale import jet_scale_shift_flat, pseudo_data_yield_sum, create_pseudo, create_Xs
from cross_section_calc import calc_cross_stuff, cross_calc, stat_unc_calc, normalization_unc_calc
import fit
import warnings
warnings.filterwarnings("ignore")



"""
ToDo:
+?Use functions from tool.py
"""

def plots(df):
  """
  Maybe soon
  """
  return None



#
#compute_unc(ana_obj, nominal_up_down)




def electron_idiso_unc(df, up_down="up"):#?should default to both...but both is complicated
  """
  electron_scale should scale the weights of samples with atleast one electron by the uncertainty on that electrons scale factor
  So: w -> w * (1 +/- unc)
  """
  weights = np.ones(df.shape[0])
  weights *= df.weight.values

  cut_lep_type = df.lep_Type != -1
  if up_down == "up":
    pm = 1
    error_string1 = "lep1_id_error_high"
    error_string2 = "lep2_id_error_high"
  else:
    pm = -1
    error_string1 = "lep1_id_error_low"
    error_string2 = "lep2_id_error_low"

  #?what is the electron type number
  lep1_cuts = cut_lep_type & (df.lep1_type == 11)
  lep2_cuts = cut_lep_type & (df.lep2_type == 11)
  weights[lep1_cuts] =  weights[lep1_cuts] * (1 + pm * df[lep1_cuts][error_string1].values) 
  weights[lep2_cuts] =  weights[lep2_cuts] * (1 + pm * df[lep2_cuts][error_string2].values) 


  #value_plot = (1 + pm * df[lep1_cuts & lep2_cuts][error_string1].values)*\
  #              (1 + pm * df[lep1_cuts & lep2_cuts][error_string2].values)
  #plt.hist(value_plot, bins=50, range=(0,2))
  #plt.show()

  return weights
   
   
def yield_string(df, string, scales, rf_ana):
  raw_numbers = {}
  string += "process \t jet0\t jet1\t jet2\n"
  for process in df.process.unique():
    j0 = 0
    j1 = 0
    j2 = 0
    
    for process_decay in df[df.process == process].process_decay.unique():
      j0 += rf_ana(df[(df.process_decay == process_decay) & (df.numb_jets == 0) & (df.lep1_Charge != df.lep2_Charge)]).weight.sum() * scales[process_decay] 
      j1 += rf_ana(df[(df.process_decay == process_decay) & (df.numb_jets == 1) & (df.lep1_Charge != df.lep2_Charge)]).weight.sum() * scales[process_decay] 
      j2 += rf_ana(df[(df.process_decay == process_decay) & (df.numb_jets == 2) & (df.lep1_Charge != df.lep2_Charge)]).weight.sum() * scales[process_decay] 


    raw_numbers[process] = np.array([j0,j1,j2], dtype=float)
    string += process + "\t" +  str(j0) + "\t" + \
                                str(j1) + "\t" + \
                                str(j2) + "\n"  
  return string, raw_numbers
  




  Xs["tot"] = cross_calc(df, df_da, df_ww, df_ggww, scales,  fiducial=True, pseudo=pseudo["tot"])  
  Xs["j0"] =  cross_calc(df[df.numb_jets == 0], df_da[df_da.numb_jets == 0], df_ww[df_ww.numb_jets == 0], df_ggww[df_ggww.numb_jets == 0], scales,  fiducial=True, pseudo=pseudo["j0"])  
  Xs["j1"] =  cross_calc(df[df.numb_jets == 1], df_da[df_da.numb_jets == 1], df_ww[df_ww.numb_jets == 1], df_ggww[df_ggww.numb_jets == 1], scales,  fiducial=True, pseudo=pseudo["j1"])  
  Xs["j2"] =  cross_calc(df[df.numb_jets == 2], df_da[df_da.numb_jets == 2], df_ww[df_ww.numb_jets == 2], df_ggww[df_ggww.numb_jets == 2], scales,  fiducial=True, pseudo=pseudo["j2"])  
  return Xs

def compute_electron(ana_obj, flavor):
  rfs    = ana_obj.rfs
  scales = ana_obj.scales

  pseudo = create_pseudo(rf_ana(ana_obj.df), rf_ana(ana_obj.df_da), scales)
  orig = create_Xs(ana_obj, pseudo)

  #########################
  #Fits
  nominal_fit_result = fit.comprehensive_fit(ana_obj.df, ana_obj.df_da, "metMod", scales)
  #########################


  print "orig", ana_obj.df.weight[ana_obj.df.lep_Type == -2].mean(), ana_obj.df[ana_obj.df.lep_Type == -2].weight.max(), ana_obj.df[ana_obj.df.lep_Type == -2].weight.min(),  ana_obj.df[ana_obj.df.lep_Type == -2].weight.std(), "why are these weights so small...Was expecting 90s avg"
  idiso_unc_weights = electron_idiso_unc(ana_obj.df)
  print "up shift", np.mean(idiso_unc_weights[ana_obj.df.lep_Type == -2]), np.std(idiso_unc_weights)  
  print "up errors", ana_obj.df.lep1_id_error_high.mean(), ana_obj.df.lep2_id_error_high.std()
  up_weights = idiso_unc_weights 



  idiso_unc_weights = electron_idiso_unc(ana_obj.df, up_down="down")
  print "\ndown shift", np.mean(idiso_unc_weights[ana_obj.df.lep_Type == -2]), np.std(idiso_unc_weights)  
  print "down errors", ana_obj.df.lep1_id_error_low.mean(), ana_obj.df.lep2_id_error_low.std()
  down_weights = idiso_unc_weights 


  print "\n\n"
  nominal_string, nominal_raw = yield_string(ana_obj.df[ana_obj.df.lep_Type != -1], "nominal\n", scales, ana_obj.rf_ana)
  print nominal_raw

  ana_obj.df.weight = up_weights 
  up_string, up_raw = yield_string(ana_obj.df[ana_obj.df.lep_Type != -1], "up\n", scales, ana_obj.rf_ana)
  print "UP", "\n".join([ str(up_raw[process] / nominal_raw[process])  for process in up_raw])
  up = create_Xs(ana_obj, pseudo)

  #########################
  #Fits
  up_fit_result = fit.comprehensive_fit(ana_obj.df, ana_obj.df_da, "metMod", scales)
  #########################


  ana_obj.df.weight = down_weights 
  down_string, down_raw = yield_string(ana_obj.df[ana_obj.df.lep_Type != -1], "down\n", scales, ana_obj.rf_ana)
  print "DOWN", "\n".join([ str(down_raw[process] / nominal_raw[process])  for process in up_raw]) 
  down = create_Xs(ana_obj, pseudo)
  #########################
  #Fits
  down_fit_result = fit.comprehensive_fit(ana_obj.df, ana_obj.df_da, "metMod", scales)
  #########################

  print orig
  print down

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



  #df.weight.hist()
  #plt.show()
  
  if flavor != "":
    flavor = "_" + flavor

  date = datetime.date.today() 
  print date.year, date.month
  f = open("results/mar/electron_isoid" + str(date.year) + "_" + str(date.month) + "_" + str(date.day) + flavor + ".txt", "w")
  f.write("Yields\n")
  f.write(nominal_string)
  f.write(up_string)
  f.write(down_string)
  f.write("\nUP\n" + "\n".join([ process + ": " + str(up_raw[process] / nominal_raw[process])  for process in up_raw]))
  f.write( "\nDOWN\n" + "\n".join([ process + ": " + str(down_raw[process] / nominal_raw[process])  for process in up_raw])) 
  f.write("\nCross Sections\n")
  f.write(str(orig) + "\n")
  f.write(str(up) + "\n")
  f.write(str(down) + "\n")
  f.write("\n\nCross Sections Unc (%)\n")
  f.write("Up: "    + str([jet + " " + str((orig[jet] - up[jet]) / (orig[jet]) * 100.)  for jet in up]) + "\n")
  f.write("Down: "  + str([jet + " " + str((orig[jet] - down[jet]) / (orig[jet]) * 100.)  for jet in up]) + "\n")
  f.write("Total: " + str([jet + " " + str((abs(orig[jet] - down[jet]) + abs(orig[jet] - up[jet])) / (2 * orig[jet]) * 100.)  for jet in up]) + "\n")
  #########################
  f.write("\n\n\nFit results\n")
  f.write(str(fit_processes) + "\n")
  f.write("nominal: "+ str(nominal_fit_result.x) + "\n")
  f.write("up: "+ str(up_fit_result.x) + "\n")
  f.write("down: "+ str(down_fit_result.x) + "\n")
  f.write("Diffs(%)" + "\n")
  for it, ele in enumerate(nominal_fit_result.x):
    f.write(str(fit_processes[it]) + str((abs(ele - up_fit_result.x[it]) + abs(ele - down_fit_result.x[it]))/2. * 100.))
  #########################



if __name__ == "__main__":
  for flavor in ["", "diff", "same"]:
    if flavor == "":
      print "total"
      ana_obj = analysis_setup("lep")
      ana_obj.apply_pre_cuts() 
      #ana_obj.apply_flat_jet_correction() 
      compute_electron(ana_obj, flavor)
    elif flavor == "diff":
      print "diff"
      ana_obj = analysis_setup("lep")
      ana_obj.apply_pre_cuts() 
      #ana_obj.apply_flat_jet_correction() 
      ana_obj.df = ana_obj.df[ana_obj.df.lep1_type != ana_obj.df.lep2_type]
      ana_obj.df_da = ana_obj.df_da[ana_obj.df_da.lep1_type != ana_obj.df_da.lep2_type]
      ana_obj.df_ww = ana_obj.df_ww[ana_obj.df_ww.lep1_type != ana_obj.df_ww.lep2_type]
      compute_electron(ana_obj, flavor)
    elif flavor == "same":
      print "same"
      ana_obj = analysis_setup("lep")
      ana_obj.apply_pre_cuts() 
      #ana_obj.apply_flat_jet_correction() 
      ana_obj.df = ana_obj.df[ana_obj.df.lep1_type == ana_obj.df.lep2_type]
      ana_obj.df_da = ana_obj.df_da[ana_obj.df_da.lep1_type == ana_obj.df_da.lep2_type]
      ana_obj.df_ww = ana_obj.df_ww[ana_obj.df_ww.lep1_type == ana_obj.df_ww.lep2_type]
      compute_electron(ana_obj, flavor)




