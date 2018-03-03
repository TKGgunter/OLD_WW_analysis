#Thoth Gunter

import os, sys, pickle, copy
sys.path.append(os.getcwd() + "/../../")
from prep_ana_II import *
sys.path.append(os.getcwd() + "/../../tools/")
sys.path.append(os.getcwd() + "/../../tools/fit/")
sys.path.append(os.getcwd() + "/../../tools/jecsys/")
from jet_scale import jec_setup, jet_scale_shift_flat, pseudo_data_yield_sum 
from cross_section_calc import calc_cross_stuff, cross_calc, stat_unc_calc, normalization_unc_calc
import fit


"""
ToDo:
+?preselection yields
+?send to file
"""


def plot_jet_eta( df ):
  df[df.numb_jets > 0].jet1_eta.hist( bins=50)
  plt.show()
 




def jer_calc( df,  ):
  #https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
  etalist  = [(0.0, 0.5), (0.5, 0.8), (0.8, 1.1), (1.1,1.3), (1.3, 1.7), (1.7, 1.9), (1.9,2.1), (2.1, 2.3), (2.3, 2.5), (2.5, 2.8), (2.8, 3.0), (3.0, 3.2), (3.2, 5.0)]
  sflist   = [1.109, 1.138, 1.114, 1.123, 1.084, 1.082, 1.140, 1.067, 1.177, 1.364, 1.857, 1.328, 1.16]                
  unclist  = [0.008, 0.013, 0.013, 0.024, 0.011, 0.035, 0.047, 0.053, 0.041, 0.039, 0.071, 0.022, 0.029]

  for k in df.keys():
    if "jet" in k and "pt" in k:
      cut_jet_pt =  df[k] > 20
  
      for it, eta_i in enumerate(etalist):
        cut_jet_eta = (df[k[:4] + "_eta"].abs() > eta_i[0]) & (df[k[:4] + "_eta"].abs() < eta_i[1])     
        
        #print sf_jer, ( 1 + np.random.normal(scale=unclist[it]) * max(sflist[it]**2 - 1, 0)**.5 )
        sf_jer = ( 1 + np.random.normal(loc= 0.0, scale=unclist[it], size=df[cut_jet_pt & cut_jet_eta].shape[0]) * max(sflist[it]**2 - 1, 0)**.5 ) 
        #if k == "jet1_pt":
        #  print "sf",  sf_jer[:5], sf_jer.min(), sf_jer.max(), sf_jer.mean()
        #  print  df[k].values[cut_jet_pt & cut_jet_eta][:5]
        #  print  (df[k].values[cut_jet_pt & cut_jet_eta] * sf_jer)[:5]
        df[k].values[cut_jet_pt & cut_jet_eta] *= sf_jer 
  





def yield_string(df, string, scales):
  raw_numbers = {}
  for process in df.process.unique():
    j0 = 0
    j1 = 0
    j2 = 0
    for process_decay in df[df.process == process].process_decay.unique():
      j0 += ana_obj.rf_ana(df[(df.process_decay == process_decay) & (df.numb_jets == 0) & (df.lep1_Charge != df.lep2_Charge)]).weight.sum() * scales[process_decay] 
      j1 += ana_obj.rf_ana(df[(df.process_decay == process_decay) & (df.numb_jets == 1) & (df.lep1_Charge != df.lep2_Charge)]).weight.sum() * scales[process_decay] 
      j2 += ana_obj.rf_ana(df[(df.process_decay == process_decay) & (df.numb_jets == 2) & (df.lep1_Charge != df.lep2_Charge)]).weight.sum() * scales[process_decay] 


    raw_numbers[process] = np.array([j0,j1,j2], dtype=float)
    string += process + "\t" +  str(j0) + "\t" + \
                                str(j1) + "\t" + \
                                str(j2) + "\n"  
  return string, raw_numbers


#Nominal

def compute_resolution( ana_obj, flavor ):

  scales = ana_obj.scales
  #jet_scale_shift_flat(ana_obj.df, jet_pt_shift=1., rf=ana_obj.rfs)


  ana_obj.df.metMod.hist(alpha=0.5, bins=50, range=(0, 100))

  pseudo = {}
  pseudo["tot"] = pseudo_data_yield_sum(rf_ana(ana_obj.df), rf_ana(ana_obj.df_da), scales=scales)
  pseudo["j0"] =  pseudo_data_yield_sum(rf_ana(ana_obj.df[ana_obj.df.numb_jets == 0]), rf_ana(ana_obj.df_da[ana_obj.df_da.numb_jets == 0]))
  pseudo["j1"] =  pseudo_data_yield_sum(rf_ana(ana_obj.df[ana_obj.df.numb_jets == 1]), rf_ana(ana_obj.df_da[ana_obj.df_da.numb_jets == 1]))
  pseudo["j2"] =  pseudo_data_yield_sum(rf_ana(ana_obj.df[ana_obj.df.numb_jets == 2]), rf_ana(ana_obj.df_da[ana_obj.df_da.numb_jets == 2]))
  pseudo["j3"] =  pseudo_data_yield_sum(rf_ana(ana_obj.df[ana_obj.df.numb_jets == 3]), rf_ana(ana_obj.df_da[ana_obj.df_da.numb_jets == 3]))

  orig = {"tot":0, "j0":0, "j1":0, "j2":0}
  orig["tot"] =  cross_calc(ana_obj.df, ana_obj.df_da, ana_obj.df_ww, ana_obj.df_ggww, scales,  fiducial=True, pseudo=pseudo["tot"])  
  orig["j0"] =  cross_calc(ana_obj.df[ana_obj.df.numb_jets == 0],\
                           ana_obj.df_da[ana_obj.df_da.numb_jets == 0],\
                           ana_obj.df_ww[ana_obj.df_ww.numb_jets == 0],\
                           ana_obj.df_ggww[ana_obj.df_ggww.numb_jets == 0], scales,  fiducial=True, pseudo=pseudo["j0"])  
  orig["j1"] =  cross_calc(ana_obj.df[ana_obj.df.numb_jets == 1],\
                           ana_obj.df_da[ana_obj.df_da.numb_jets == 1],\
                           ana_obj.df_ww[ana_obj.df_ww.numb_jets == 1],\
                           ana_obj.df_ggww[ana_obj.df_ggww.numb_jets == 1], scales,  fiducial=True, pseudo=pseudo["j1"])  
  orig["j2"] =  cross_calc(ana_obj.df[ana_obj.df.numb_jets == 2],\
                           ana_obj.df_da[ana_obj.df_da.numb_jets == 2],\
                           ana_obj.df_ww[ana_obj.df_ww.numb_jets == 2],\
                           ana_obj.df_ggww[ana_obj.df_ggww.numb_jets == 2], scales,  fiducial=True, pseudo=pseudo["j2"])  
  orig["j3"] =  cross_calc(ana_obj.df[ana_obj.df.numb_jets == 3],\
                           ana_obj.df_da[ana_obj.df_da.numb_jets == 3],\
                           ana_obj.df_ww[ana_obj.df_ww.numb_jets == 3],\
                           ana_obj.df_ggww[ana_obj.df_ggww.numb_jets == 3], scales,  fiducial=True, pseudo=pseudo["j3"])  



  nominal_string, nominal_raw = yield_string(ana_obj.df, "Nominal\n", scales)
  #########################
  #Fits
  nominal_fit_result = fit.comprehensive_fit(ana_obj.df, ana_obj.df_da, "metMod", scales)
  #########################

  #Post JER
  ana_obj.reset()
  ana_obj.apply_pre_cuts()

  if flavor == "same":
    ana_obj.df = ana_obj.df[ana_obj.df.lep1_type == ana_obj.df.lep2_type]
    ana_obj.df_da = ana_obj.df_da[ana_obj.df_da.lep1_type == ana_obj.df_da.lep2_type]
    ana_obj.df_ww = ana_obj.df_ww[ana_obj.df_ww.lep1_type == ana_obj.df_ww.lep2_type]

  if flavor == "diff":
    ana_obj.df = ana_obj.df[ana_obj.df.lep1_type != ana_obj.df.lep2_type]
    ana_obj.df_da = ana_obj.df_da[ana_obj.df_da.lep1_type != ana_obj.df_da.lep2_type]
    ana_obj.df_ww = ana_obj.df_ww[ana_obj.df_ww.lep1_type != ana_obj.df_ww.lep2_type]


  jer_calc(ana_obj.df)
  ana_obj.apply_flat_jet_correction()
  #jet_scale_shift_flat(ana_obj.df, jet_pt_shift=1., rf=ana_obj.rfs)
  ana_obj.df.metMod.hist(alpha=0.5, bins=50, range=(0, 100))
  plt.legend(["orig", "jer"])
  plt.show()




  post = {"tot":0, "j0":0, "j1":0, "j2":0}
  post["tot"] =  cross_calc(ana_obj.df, ana_obj.df_da, ana_obj.df_ww, ana_obj.df_ggww, scales,  fiducial=True, pseudo=pseudo["tot"])  
  post["j0"] =  cross_calc(ana_obj.df[ana_obj.df.numb_jets == 0],\
                           ana_obj.df_da[ana_obj.df_da.numb_jets == 0],\
                           ana_obj.df_ww[ana_obj.df_ww.numb_jets == 0],\
                           ana_obj.df_ggww[ana_obj.df_ggww.numb_jets == 0], scales,  fiducial=True, pseudo=pseudo["j0"])  
  post["j1"] =  cross_calc(ana_obj.df[ana_obj.df.numb_jets == 1],\
                           ana_obj.df_da[ana_obj.df_da.numb_jets == 1],\
                           ana_obj.df_ww[ana_obj.df_ww.numb_jets == 1],\
                           ana_obj.df_ggww[ana_obj.df_ggww.numb_jets == 1], scales,  fiducial=True, pseudo=pseudo["j1"])  
  post["j2"] =  cross_calc(ana_obj.df[ana_obj.df.numb_jets == 2],\
                           ana_obj.df_da[ana_obj.df_da.numb_jets == 2],\
                           ana_obj.df_ww[ana_obj.df_ww.numb_jets == 2],\
                           ana_obj.df_ggww[ana_obj.df_ggww.numb_jets == 2], scales,  fiducial=True, pseudo=pseudo["j2"])  
  post["j3"] =  cross_calc(ana_obj.df[ana_obj.df.numb_jets == 3],\
                           ana_obj.df_da[ana_obj.df_da.numb_jets == 3],\
                           ana_obj.df_ww[ana_obj.df_ww.numb_jets == 3],\
                           ana_obj.df_ggww[ana_obj.df_ggww.numb_jets == 3], scales,  fiducial=True, pseudo=pseudo["j3"])  


  #########################
  #Fits
  post_fit_result = fit.comprehensive_fit(ana_obj.df, ana_obj.df_da, "metMod", scales)
  #########################

  jer_string, jer_raw = yield_string(ana_obj.df, "JER\n", scales)

  print nominal_string 
  print orig
  print jer_string 
  print post
  print "ratio", [(k, post[k] / orig[k]) for k in orig]
  print "unc", [ (k, abs(post[k] - orig[k])/orig[k]) for k in orig ]


  #########################
  print "Fit results"
  fit_processes = ["WW", "Top", "DY"] 
  print fit_processes
  print "nominal: ", nominal_fit_result.x
  print "post: ", post_fit_result.x
  print "Diffs(%)"
  for it, ele in enumerate(nominal_fit_result.x):
    print  fit_processes[it], abs(ele - post_fit_result.x[it]) * 100.
  #########################



  if flavor != "":
    flavor = "_"+flavor

  f = open("results/jan/jet_resolution" + flavor + ".txt", "w")
  f.write(nominal_string)
  f.write(jer_string)
  f.write("\n\nUncertianty(%): " + str([ (k, abs(post[k] - orig[k])/orig[k] * 100) for k in orig ]))
  #########################
  f.write("\n\nFit results\n")
  f.write(str(fit_processes) + "\n")
  f.write("nominal: "+ str(nominal_fit_result.x) + "\n")
  f.write("post: "+ str(post_fit_result.x) + "\n")
  f.write("Diffs(%)" + "\n")
  for it, ele in enumerate(nominal_fit_result.x):
    f.write(str(fit_processes[it]) + " " + str(abs(ele - post_fit_result.x[it]) * 100.) + "\n")
  #########################



if __name__ == "__main__": 
  for flavor in ["",]:# "diff", "same"]:
    if flavor == "":
      ana_obj = analysis_setup( "jet")
      ana_obj.apply_pre_cuts() 
      ana_obj.apply_flat_jet_correction() 
      compute_resolution(ana_obj, flavor)
    if flavor == "same":
      ana_obj = analysis_setup("jet")
      ana_obj.apply_pre_cuts() 
      ana_obj.apply_flat_jet_correction() 
      ana_obj.df = ana_obj.df[ana_obj.df.lep1_type == ana_obj.df.lep2_type]
      ana_obj.df_da = ana_obj.df_da[ana_obj.df_da.lep1_type == ana_obj.df_da.lep2_type]
      ana_obj.df_ww = ana_obj.df_ww[ana_obj.df_ww.lep1_type == ana_obj.df_ww.lep2_type]

      jet_scale_shift_flat(ana_obj.df, jet_pt_shift=1., rf=ana_obj.rfs)
      compute_resolution(ana_obj, flavor)
    if flavor == "diff":
      ana_obj = analysis_setup("jet")
      ana_obj.apply_pre_cuts() 
      ana_obj.apply_flat_jet_correction() 
      ana_obj.df = ana_obj.df[ana_obj.df.lep1_type != ana_obj.df.lep2_type]
      ana_obj.df_da = ana_obj.df_da[ana_obj.df_da.lep1_type != ana_obj.df_da.lep2_type]
      ana_obj.df_ww = ana_obj.df_ww[ana_obj.df_ww.lep1_type != ana_obj.df_ww.lep2_type]

      jet_scale_shift_flat(ana_obj.df, jet_pt_shift=1., rf=ana_obj.rfs)
      compute_resolution(ana_obj, flavor)




