import os, sys, pickle, copy
import datetime
sys.path.append(os.getcwd() + "/../../")
sys.path.append(os.getcwd() + "/../../tools/fit/")
from prep_ana_II import *
from cross_section_calc import calc_cross_stuff, cross_calc, stat_unc_calc, normalization_unc_calc
import fit
import warnings
warnings.filterwarnings("ignore")




def compute_resumm_alt(flavor):
  fig, ax1 = plt.subplots()

  ana_obj = analysis_setup("lep")
  ana_obj.apply_pre_cuts() 
  ana_obj.apply_flat_jet_correction() 

  df     = ana_obj.df
  scales = ana_obj.scales
  df_da  = ana_obj.df_da

  if flavor == "same":
    df    = df[df.lep1_type == df.lep2_type]
    df_da = df_da[df_da.lep1_type == df_da.lep2_type]
  elif flavor == "diff":
    df    = df[df.lep1_type != df.lep2_type]
    df_da = df_da[df_da.lep1_type != df_da.lep2_type]


  random_forests = ana_obj.rfs 
  features_fDY = random_forests["features_fDY"]
  clf_fDY = random_forests["clf_fDY"]
  features_fTT = random_forests["features_fTT"]
  clf_fTT = random_forests["clf_fTT"]

  temp_df = df[(df.process_decay == "WW") & (df.ww_pt > -1)]
  temp_df["ww_weights"] = [1.] * temp_df.shape[0]


  import ROOT
  f = ROOT.TFile.Open("/home/gunter/WW_analysis/production/Analysis_13TeV/tools/WW_pTresummation_13-master/MyRatioWWpTHistogramAll.root", "read")

  name_list = ["wwpt", "wwpt_scaleup", "wwpt_scaledown", "wwpt_resumup", "wwpt_resumdown"]
  wwpt_dic = {}
  for name in name_list:
    wwpt_dic[name] = f.Get(name)
  
  orig = 0
  orig0 = 0
  orig1 = 0
  q_ =  {}
  r_ =  {}
  q_0 = {}
  r_0 = {}
  q_1 = {}
  r_1 = {}


  for name in wwpt_dic:
    y_arr = []
    for i in range(1000):
        weight =  wwpt_dic[name].GetBinContent(i)
        if name != "wwpt":
          weight *=  wwpt_dic["wwpt"].GetBinContent(i)
        temp_df.ww_weights.values[ (temp_df.ww_pt >= i*0.5 ).values & (temp_df.ww_pt < i*0.5 + 0.5).values]  = weight
        y_arr.append(temp_df[ (temp_df.ww_pt >= i*0.5) & (temp_df.ww_pt < i*0.5 + 0.5) ].ww_weights.sum())

    tot_temp_value  = temp_df.ww_weights.sum() * scales["WW"]
    tot_temp_value0 = temp_df.query("numb_jets == 0").ww_weights.sum() * scales["WW"]
    tot_temp_value1 = temp_df.query("numb_jets == 1").ww_weights.sum() * scales["WW"]
    temp_value = rf_ana(temp_df).ww_weights.sum() * scales["WW"]
    temp_value0 = rf_ana(temp_df).query("numb_jets == 0").ww_weights.sum() * scales["WW"]
    temp_value1 = rf_ana(temp_df).query("numb_jets == 1").ww_weights.sum() * scales["WW"]
    if "scale" in name:
      q_[name]  = (temp_value , tot_temp_value)
      q_0[name] = (temp_value0, tot_temp_value0)  
      q_1[name] = (temp_value1, tot_temp_value1)
    elif "resum" in name:
      r_[name]  = (temp_value , tot_temp_value)
      r_0[name] = (temp_value0, tot_temp_value0)  
      r_1[name] = (temp_value1, tot_temp_value1)  
    else:
      orig  = (temp_value , tot_temp_value)
      orig0 = (temp_value0, tot_temp_value0)  
      orig1 = (temp_value1, tot_temp_value1)  
      ax1.scatter([i*0.5 for i in range(1000)], y_arr, color='r')
      ax1.set_ylabel('hist of WW', color='r') 
       

  print "Orig", orig
  print "Orig0", orig0
  print "Orig1", orig1
  print ""
  print "scale", q_
  print "scale0", q_0
  print "scale1", q_1
  print ""
  print "resum", r_
  print "resum0", r_0
  print "resum1", r_1

  print "Results:"
  _values  = []
  _0values = []
  _1values = []
  for name in q_:
    print name, "Tot: ", (q_[name][0]/ q_[name][1]) / (orig[0]/orig[1])
    print name, "j0: ", (q_0[name][0]/q_0[name][1]) / (orig0[0]/orig0[1])
    print name, "j1: ", (q_1[name][0]/q_1[name][1]) / (orig1[0]/orig1[1])
    _values.append((q_[name][0]/ q_[name][1]) / (orig[0]/orig[1]))
    _0values.append((q_0[name][0]/ q_0[name][1]) / (orig0[0]/orig0[1]))
    _1values.append((q_1[name][0]/ q_1[name][1]) / (orig1[0]/orig1[1]))

  tot_scale_unc = sum([abs(1 - i) for i in _values])/2. * 100.
  j0_scale_unc = sum([abs(1 - i) for i in _0values])/2. * 100.
  j1_scale_unc = sum([abs(1 - i) for i in _1values])/2. * 100.

  print "Scale tot unc", tot_scale_unc
  print "Scale 0j  unc", j0_scale_unc
  print "Scale 1j  unc", j1_scale_unc

  _values  = []
  _0values = []
  _1values = []
  for name in r_:
    print name, "Tot: ", (r_[name][0] /r_[name][1])  / (orig[0]/orig[1])
    print name, "j0:  ", (r_0[name][0]/r_0[name][1]) / (orig0[0]/orig0[1])
    print name, "j1:  ", (r_1[name][0]/r_1[name][1]) / (orig1[0]/orig1[1])
    _values.append((r_[name][0]/ r_[name][1]) / (orig[0]/orig[1]))
    _0values.append((r_0[name][0]/ r_0[name][1]) / (orig0[0]/orig0[1]))
    _1values.append((r_1[name][0]/ r_1[name][1]) / (orig1[0]/orig1[1]))
  tot_res_unc = sum([abs(1 - i) for i in _values])/2. * 100.
  j0_res_unc = sum([abs(1 - i) for i in _0values])/2. * 100.
  j1_res_unc = sum([abs(1 - i) for i in _1values])/2. * 100.

  print "Res tot unc", tot_res_unc
  print "Res 0j  unc", j0_res_unc
  print "Res 1j  unc", j1_res_unc

  plt.show()

if __name__ == "__main__":
  for flavor in ["",]:# "same", "diff"]:
    compute_resumm_alt(flavor)
    print "\n\n"









