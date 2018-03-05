#Thoth Gunter
import numpy as np
from physics_selections import rf_ana

import os, sys, time
from prep_ana_II import * 




def pseudo_data_yield_sum(df_mc, df_data, flavor="both", weights="weight", scales=scales):

    table = process_yields(df_mc, df_data, scales=scales) 
    return table[ table.Process == "Total" ]["Same Flavor"].values[0] + table[ table.Process == "Total" ]["Diff Flavor"].values[0] 


scale_ww_tot = scales["WW"]
scale_ggww = scales["GluGluWWTo2L2Nu"]


def calc_cross_stuff(df_mc, df_data, df_ww, df_ggww, scales, flavor="both", weights='weight', rf_ana=rf_ana, pseudo=-99, scale_ww_tot=scale_ww_tot, scale_ggww=scale_ggww ):
    lumi = float(lumi_amount) * 10**3
    acc  = 0.188
    Br   = (3*.108)**2.
    
    df_mc = df_mc[df_mc.pred_fTT_WW > .6]
    df_data = df_data[df_data.pred_fTT_WW > .6]
    table = process_yields(rf_ana(df_mc, flavor=flavor), rf_ana(df_data, flavor=flavor), scales=scales) 
    #print table
 
    N_mc = 0.
    N_Wjets = 0.
    for process in table.Process.unique():
        if "Total" in process: continue
        if "DATA" in process: continue
        if process == "Total: WW": continue 
        if process == "WW": continue 
        if process == "GluGluWWTo2L2Nu": continue 
        N_mc += table[table.Process == process]["Same Flavor"].values[0] + table[table.Process == process]["Diff Flavor"].values[0] 
        #print process, table[table.Process == process]["Same Flavor"].values[0] + table[table.Process == process]["Diff Flavor"].values[0] 

        if process == "WJ": N_Wjets += table[table.Process == process]["Same Flavor"].values[0] + table[table.Process == process]["Diff Flavor"].values[0] 


    if pseudo == -99:
      N_data = table[table.Process == "DATA"]["Same Flavor"].values[0] + table[table.Process == "DATA"]["Diff Flavor"].values[0]
    else:
      N_data = pseudo

 
    
    if flavor == 'both':
      query_str = "lep_Type > -3"
    elif flavor == 'same':
      query_str = "lep_Type < 0"
    elif flavor == 'diff':
      query_str = "lep_Type > 0"
    else:
      print "flavor is trash", flavor
      print "Acceptable flavors are: both, same, diff"
      exit()



    N_ww_select = table[table.Process == "WW"]["Same Flavor"].values[0] + table[table.Process == "WW"]["Diff Flavor"].values[0]
    N_ww_tot = df_ww[df_ww.process_decay == "WW"].query(query_str)[weights].sum() * scale_ww_tot

    if type(df_ggww) == type(df_ww):
      N_ww_tot += df_ggww[df_ggww.process_decay == "GluGluWWTo2L2Nu"].query(query_str)[weights].sum() * scale_ggww 
      N_ww_select = table[table.Process == "Total: WW"]["Same Flavor"].values[0] + table[table.Process == "Total: WW"]["Diff Flavor"].values[0]


    #print "N_ww_tot", N_ww_tot, N_ww_select, N_ww_select / N_ww_tot, scale_ww_tot
    ratio_s_t = N_ww_select / N_ww_tot

    return {"lumi": lumi, "acc": acc, "Br": Br, "N_mc": N_mc, "N_data": N_data, "ratio_s_t": ratio_s_t, "N_ww_select":N_ww_select, "N_Wjets": N_Wjets}



def cross_calc(df_mc, df_data, df_ww, df_ggww, scales, flavor="both", fiducial=False, pseudo=-99, scale_ww_tot=35.9e3 * 118.7 * (3*.108)**2 / 1998956, scale_ggww=35.9e3 * 0.84365 / 500000.,  query=None, **kwargs):
    if kwargs:
        var = kwargs
    else:
        if type(query) == type("asdf"):
          df_mc = df_mc.query(query)
          df_data = df_data.query(query)
          df_ww = df_ww.query(query)
          df_ggww = df_ggww.query(query)
        var = calc_cross_stuff(df_mc, df_data, df_ww, df_ggww, scales=scales, flavor=flavor, pseudo=pseudo, scale_ww_tot=scale_ww_tot, scale_ggww=scale_ggww)

    lumi = var["lumi"]
    acc = var["acc"]
    Br = var["Br"]
    N_mc = var["N_mc"]
    N_data = var["N_data"]
    ratio_s_t = var["ratio_s_t"]
    N_ww_select = var["N_ww_select"]
    N_Wjets = var["N_Wjets"]

    #print "var from function", var    
    if fiducial == False:
        return  (N_ww_select) / (lumi * acc * Br *ratio_s_t)
        #return  (N_data - N_mc) / (lumi * acc * Br * ratio_s_t)
    else:
        return (N_ww_select) / (lumi * ratio_s_t ) 
        #return (N_data - N_mc) / (lumi * ratio_s_t ) 



def normalization_unc_calc(df_mc, df_data, df_ww, df_ggww, scales, unc_mc_process, flavor="both", fiducial=False, rf_ana=rf_ana, **kwargs):
    """
    X = (D - MC) / (L STUFF),
    MC = SUM: sf_i * yields_i,
    unc = SUM_QUAD( d/dsf_i X  * UNC_sf)**.5
    UNC_sf = %unc * sf
    """
    if kwargs:
        var = kwargs
    else:
        var = calc_cross_stuff(df_mc, df_data, df_ww, df_ggww, scales=scales, flavor=flavor)
    lumi = var["lumi"]
    Br = var["Br"]
    acc = var["acc"]
    N_mc = var["N_mc"]
    N_data = var["N_data"]
    ratio_s_t = var["ratio_s_t"]
    N_ww_select = var["N_ww_select"]
    N_Wjets = var["N_Wjets"]
    
    cuts_mc = {process: rf_ana(df_mc[(df_mc.process_decay == process) & (df_mc.lep1_Charge != df_mc.lep2_Charge)], flavor=flavor) for process in scales.keys()} 
    process_sys_unc = []
    temp_diboson = 0

    for process in cuts_mc.keys():
      if process not in ['WW', 'GluGluWWTo2L2Nu', 'W1JetsToLNu','W2JetsToLNu','W3JetsToLNu','W4JetsToLNu']:
        print process, "unc on process", unc_mc_process[process], "\traw yields", cuts_mc[process].shape[0], "\tyields", scales[process] * cuts_mc[process].shape[0],\
              "\tunc * yields", unc_mc_process[process] * cuts_mc[process].shape[0],\
              "\t1/(L r) * unc * yields", str(1. / (lumi * ratio_s_t) * unc_mc_process[process] * cuts_mc[process].shape[0])[:4]
        if process in ["WZ", "ZZ",]:
          temp_diboson = temp_diboson + (unc_mc_process[process] * scales[process]) * cuts_mc[process].shape[0]
        else:
          process_sys_unc.append((unc_mc_process[process] * scales[process])**2 * cuts_mc[process].shape[0]**2)

    process_sys_unc.append(temp_diboson**2)    


    WW_sys_unc = (N_data - N_mc) * unc_mc_process["WW"]**2 * ratio_s_t**2. * df_ww.shape[0] * scales["WW"] / N_ww_select**2. #UNDER CONSTRUCTION


    print "Wjets", "unc on process", unc_mc_process["W4JetsToLNu"], "\traw yields", N_Wjets, "\tyields",  N_Wjets,\
          "\tunc * yields", unc_mc_process["W4JetsToLNu"] * N_Wjets,\
          "\t1/(L r) * unc * yields", str(1. / (lumi * ratio_s_t) * unc_mc_process["W4JetsToLNu"] * N_Wjets)[:4]
    Wjets_sys_unc = N_Wjets * unc_mc_process["W4JetsToLNu"] 
   

    print "w/o wj", 1. / (lumi * ratio_s_t) * sum(process_sys_unc)**.5,  "w/ wj", 1. / (lumi * ratio_s_t) * Wjets_sys_unc**.5
    if fiducial == False:
        return 1. / (lumi * acc * Br * ratio_s_t)  * ( sum(process_sys_unc) +  WW_sys_unc + Wjets_sys_unc)**.5
    else:
        return 1. / (lumi * ratio_s_t )  * ( sum(process_sys_unc) + Wjets_sys_unc )**.5





def sum_quad(array):
  #I don't think your doing what you think your doing
  #result = []
  #for i in array:
  #  result.append(i**2)
  #return sum(result)
  return np.sum(array**2)



def stat_unc_calc(df_mc, df_data, df_ww, df_ggww, scales, unc_mc_process, flavor="both", fiducial=False, rf_ana=rf_ana, **kwargs):
    if kwargs:
        var = kwargs
    else:
        var = calc_cross_stuff(df_mc, df_data, df_ww, df_ggww, scales=scales, flavor=flavor)
    lumi = var["lumi"]
    Br = var["Br"]
    acc = var["acc"]
    N_mc = var["N_mc"]
    N_data = var["N_data"]
    ratio_s_t = var["ratio_s_t"]
    N_ww_select = var["N_ww_select"]
    N_Wjets = var["N_Wjets"]

    #This is just an estimate I need to determine the uncertainties from monte carlo properly
    if flavor == "diff":
      df_mc = df_mc[df_mc.lep_Type > 0]
      df_data = df_data[df_data.lep_Type > 0]
      df_ww = df_ww[df_ww.lep_Type > 0]
    elif flavor == "same":
      #print ratio_s_t
      df_mc = df_mc[df_mc.lep_Type < 0]
      df_data = df_data[df_data.lep_Type < 0]
      f_ww = df_ww[df_ww.lep_Type < 0]
    if fiducial == False:
        return 1. / (lumi * acc * Br * ratio_s_t)  * (N_data + sum([scales[process]**2. * df_mc[df_mc.process_decay == process].shape for process in df_mc.process_decay.unique()]))**.5
    else:
      #  print "ASD", N_mc + N_Wjets + N_ww_select,\
      #sum([scales[process] * df_mc[df_mc.process_decay == process].weight.values.sum() for process in df_mc.process_decay.unique()]),\
      #[(process, scales[process] * df_mc[df_mc.process_decay == process].weight.values.sum()) for process in df_mc.process_decay.unique()]
        #return (sum([scales[process]**2. * sum_quad(df_mc[df_mc.process_decay == process].weight.values) for process in df_mc.process_decay.unique()]))**.5 /\
        #       (sum([scales[process] * df_mc[df_mc.process_decay == process].weight.values.sum() for process in df_mc.process_decay.unique()]))
        return (sum([scales[process]**2. * sum_quad(df_mc[df_mc.process_decay == process].weight.values) for process in df_mc.process_decay.unique()]))**.5 /\
                (N_mc + N_Wjets + N_ww_select)





if __name__ == "__main__":
#  df = pd.read_hdf("data/preselMC_rf.hdf", "table")   
#  df_da = pd.read_hdf("data/preselDATA_rf.hdf", "table")  
#  df_ww = rp.read_root(data_path+"/ww_complete.root", columns=columns)
#  df_ggww = rp.read_root(data_path+"/glugluww_complete.root", columns=columns)  

  ana_obj = analysis_setup(unc="lep")
  scales = ana_obj.scales
  ana_obj.apply_pre_cuts()
  df = ana_obj.df
  df_da = ana_obj.df_da
  df_ww = ana_obj.df_ww
  df_ggww = ana_obj.df_ggww

  #print "Calc cross section stuff", calc_cross_stuff(df, df_da, df_ww, df_ggww, scales)
  print "Total Cross section calc: Fiducial:", cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=True), "Total:", cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=False)
  print "0j: Cross section calc: Fiducial:", cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=True, query="numb_jets == 0"), "Total:", cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=False, query="numb_jets == 0")
  print "1j: Cross section calc: Fiducial:", cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=True, query="numb_jets == 1"), "Total:", cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=False, query="numb_jets == 1")
  print "2j: Cross section calc: Fiducial:", cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=True, query="numb_jets == 2"), "Total:", cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=False, query="numb_jets == 2")
  print "3j: Cross section calc: Fiducial:", cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=True, query="numb_jets == 3"), "Total:", cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=False, query="numb_jets == 3")
  print "Stat unc", stat_unc_calc(rf_ana(df), rf_ana(df_da), df_ww, df_ggww, scales, unc_mc_process, fiducial=True)

#  for p in df.process_decay.unique():
#    print p, scales[p]**2, rf_ana(df[df.process_decay == p]).shape[0]




