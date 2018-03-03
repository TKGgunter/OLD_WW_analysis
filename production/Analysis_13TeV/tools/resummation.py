#Thoth Gunter

import os, sys, time
sys.path.append(os.getcwd() + "/../")

from prep_ana_II import * 
from cross_section_calc import calc_cross_stuff, cross_calc  
import root_pandas as rp
import numpy as np
import pickle as pkl
import pandas as pd

def load_resummation_stuff( file_type ):
  df = load_presel_w_fDY_fTT_MC()
  df_da = load_presel_w_fDY_fTT_DATA()


  #resum_down.dat     resum_up.dat       scale_down.dat   scale_up.dat
  resumm_stuff = {}
  tot_xs = {}
  if file_type == "powheg":
    resumm_stuff["central"] = pd.read_csv("WW_pTresummation_13-master/powheg_2l2nu_nlo.dat", delim_whitespace=True, index_col=False, names=["WWpt", "dSigma"])
    resumm_stuff["resum_down"] = pd.read_csv("WW_pTresummation_13-master/powheg_2l2nu_qdown_nlo.dat", delim_whitespace=True, index_col=False, names=["WWpt", "dSigma"])
    resumm_stuff["resum_up"] = pd.read_csv("WW_pTresummation_13-master/powheg_2l2nu_qup_nlo.dat", delim_whitespace=True, index_col=False, names=["WWpt", "dSigma"])
    resumm_stuff["scale_down"] = pd.read_csv("WW_pTresummation_13-master/powheg_2l2nu_sdown_nlo.dat", delim_whitespace=True, index_col=False, names=["WWpt", "dSigma"])
    resumm_stuff["scale_up"] = pd.read_csv("WW_pTresummation_13-master/powheg_2l2nu_sup_nlo.dat", delim_whitespace=True, index_col=False, names=["WWpt", "dSigma"])

    tot_xs["central"]    = 106.546
    tot_xs["resum_down"] = 106.546
    tot_xs["resum_up"]   = 106.546 
    tot_xs["scale_down"] = 105.742
    tot_xs["scale_up"]   = 108.262
    for i in resumm_stuff:
      print "Resumm: ", i, "intergal: ", resumm_stuff[i].dSigma.sum()
      resumm_stuff[i]["dSigma"] = resumm_stuff[i].dSigma.values * (tot_xs[i] / resumm_stuff[i].dSigma.sum())
      print "Renormalizing bin contents intergal: ", resumm_stuff[i].dSigma.sum()


  elif file_type == "normal":
    resumm_stuff["central"] =     pd.read_csv("WW_pTresummation_13-master/central.dat", delim_whitespace=True, index_col=False, names=["WWpt", "dSigma"])
    resumm_stuff["resum_down"] =  pd.read_csv("WW_pTresummation_13-master/resum_down.dat", delim_whitespace=True, index_col=False, names=["WWpt", "dSigma"])
    resumm_stuff["resum_up"] =    pd.read_csv("WW_pTresummation_13-master/resum_up.dat", delim_whitespace=True, index_col=False, names=["WWpt", "dSigma"])
    resumm_stuff["scale_down"] =  pd.read_csv("WW_pTresummation_13-master/scale_down.dat", delim_whitespace=True, index_col=False, names=["WWpt", "dSigma"])
    resumm_stuff["scale_up"] =    pd.read_csv("WW_pTresummation_13-master/scale_up.dat", delim_whitespace=True, index_col=False, names=["WWpt", "dSigma"])

    tot_xs["central"]    = 106.546
    tot_xs["resum_down"] = 106.546
    tot_xs["resum_up"]   = 106.546 
    tot_xs["scale_down"] = 105.742
    tot_xs["scale_up"]   = 108.262
    for i in resumm_stuff:
      print "Resumm: ", i, "intergal: ", resumm_stuff[i].dSigma.sum()
      resumm_stuff[i]["dSigma"] = resumm_stuff[i].dSigma.values * (tot_xs[i] / resumm_stuff[i].dSigma.sum())
      print "Renormalizing bin contents intergal: ", resumm_stuff[i].dSigma.sum()



  else:
    resumm_stuff["central"] =     pd.read_csv("WW_pTresummation_13-master/central_np.dat", delim_whitespace=True, index_col=False, names=["WWpt", "dSigma"])
    resumm_stuff["resum_down"] =  pd.read_csv("WW_pTresummation_13-master/resum_down_np.dat", delim_whitespace=True, index_col=False, names=["WWpt", "dSigma"])
    resumm_stuff["resum_up"] =    pd.read_csv("WW_pTresummation_13-master/resum_up_np.dat", delim_whitespace=True, index_col=False, names=["WWpt", "dSigma"])
    resumm_stuff["scale_down"] =  pd.read_csv("WW_pTresummation_13-master/scale_down_np.dat", delim_whitespace=True, index_col=False, names=["WWpt", "dSigma"])
    resumm_stuff["scale_up"] =    pd.read_csv("WW_pTresummation_13-master/scale_up_np.dat", delim_whitespace=True, index_col=False, names=["WWpt", "dSigma"])

    tot_xs["central"]    = 106.546
    tot_xs["resum_down"] = 106.546
    tot_xs["resum_up"]   = 106.546 
    tot_xs["scale_down"] = 105.742
    tot_xs["scale_up"]   = 108.262
    for i in resumm_stuff:
      print "Resumm: ", i, "intergal: ", resumm_stuff[i].dSigma.sum()
      resumm_stuff[i]["dSigma"] = resumm_stuff[i].dSigma.values * (tot_xs[i] / resumm_stuff[i].dSigma.sum())
      print "Renormalizing bin contents ", resumm_stuff[i].dSigma.sum()

 
  return resumm_stuff, tot_xs 


def save_ww_weights(df_ww = rp.read_root("../data/WWpt/ww_complete.root", columns=["ww_pt", "weight"]), file_type = ""):
  resummation_dic, xs_dic = load_resummation_stuff(file_type)
  

  home_dir = "/home/gunter/WW_analysis/production/Analysis_13TeV/" 
  page = "13TeV"
  plotting_options = pd.read_csv(home_dir+"/plotting_options_"+page+".csv", index_col=None, sep=";")
  scales = {}
  for key in plotting_options.process_decay.unique():
    scales[key] = eval(plotting_options[plotting_options.process_decay == key]["scale"].values[0]) 

  print df_ww.keys(), df_ww[df_ww.ww_pt < 0 ].shape
  temp_df = df_ww[df_ww.ww_pt > -1]

  j = 0

  for k in resummation_dic:
    save = []
    for it, pt in enumerate(resummation_dic[k].WWpt.values):
      
      if it == 0: save.append([pt, it*.25 + .25, 1, 1, 1])
      elif temp_df[(temp_df.ww_pt >= it*.25) & (temp_df.ww_pt < it*.25+.25)].weight.sum() == 0: 
        save.append([pt, it*.25 + .25, 1, resummation_dic[k].dSigma.values[it], 1])
      else:
        save.append( [it*.25 ,\
                      it*.25 + .25 ,\
                      resummation_dic[k].dSigma.values[it] /(temp_df[(temp_df.ww_pt >= it*.25) & (temp_df.ww_pt < it*.25+.25)].weight.sum() * 118.7 / temp_df.weight.sum()),\
                      resummation_dic[k].dSigma.values[it],\
                      temp_df[(temp_df.ww_pt >= it*.25) & (temp_df.ww_pt < it*.25+.25)].weight.sum() * 118.7 / temp_df.weight.sum() ])

    print save[0]
    save = np.array(save)

    print "Is every thing the same shape: ", save.shape
    save = pd.DataFrame(save, columns=["low", "high", "weight", "theory", "MC"])
    #save.to_csv("test_temp_"+k+".csv", index=False)
    save.to_csv("../data/WWpt/resummation_"+k+"_"+file_type+"_weights.csv", index=False)



if __name__ == "__main__":
  #print generate_ww_weights().shape
  #save_ww_weights()
  #generate_ww_weights(np.array([10.]), np.array([12.]))
  #resummation_stuff = load_resummation_stuff()
  #for k in resummation_stuff:
  #  print k, resummation_stuff[k].head()

  for ft in ["normal", "np", "powheg"]:
    print "saving", ft
    save_ww_weights(file_type= ft)






