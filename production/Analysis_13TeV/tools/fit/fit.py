#Thoth Gunter
import numpy as np
import os, sys
sys.path.append(os.getcwd() + "/../../")
from prep_ana_II import bin_df, full_bin_plot, scales
from scipy.optimize import curve_fit, minimize
import matplotlib.pyplot as plt 
 
def c_dy(df, charge_state=""):
    dy_cuts = (df.pred_fDY_WW < .6) & (df.pred_fTT_WW > .6)
    if charge_state == "diff":
        dy_cuts = dy_cuts & (df.lep1_Charge != df.lep2_Charge)
    if charge_state == "same":
        dy_cuts = dy_cuts & (df.lep1_Charge == df.lep2_Charge)
    
    return df[dy_cuts]

def c_tt(df, charge_state=""):
    tt_cuts = (df.pred_fDY_WW > .6) & (df.pred_fTT_WW < .6)
    if charge_state == "diff":
        tt_cuts = tt_cuts & (df.lep1_Charge != df.lep2_Charge)
    if charge_state == "same":
        tt_cuts = tt_cuts & (df.lep1_Charge == df.lep2_Charge)
    
    return df[tt_cuts]

def c_ww(df, charge_state=""):
    ww_cuts = (df.pred_fDY_WW > .9) & (df.pred_fTT_WW > .6)
    if charge_state == "diff":
        ww_cuts = ww_cuts & (df.lep1_Charge != df.lep2_Charge)
    if charge_state == "same":
        ww_cuts = ww_cuts & (df.lep1_Charge == df.lep2_Charge)
    
    return df[ww_cuts]



def pack_data_dic(data_dic, df_da, binned_feature= None):
    charge = "diff"
    if binned_feature == None:
        data_dic["dy"] = float(c_dy(df_da, charge).shape[0])
        data_dic["tt"] = float(c_tt(df_da, charge).shape[0])
        data_dic["sig"]= float(c_ww(df_da, charge).shape[0])
    else:
        data_dic["dy"] = bin_df(c_dy(df_da, charge), binned_feature=binned_feature)['Da'][0].astype(np.float64)
        data_dic["tt"] = bin_df(c_tt(df_da, charge), binned_feature=binned_feature)['Da'][0].astype(np.float64)
        data_dic["sig"]= bin_df(c_ww(df_da, charge), binned_feature=binned_feature)['Da'][0].astype(np.float64)
        
    
def pack_mc_dic(mc_dic, df, df_da, scales, binned_feature=None):
    ds_ov_ss = 2.04

    if binned_feature == None:
        mc_c_dy = {'W1JetsToLNu':0}
        mc_c_tt = {'W1JetsToLNu':0}
        mc_s    = {'W1JetsToLNu':0}
        
        sf_wj_dy = 0
        sf_wj_tt = 0
        sf_wj_s  = 0

        for p in df.process_decay.unique():
            if p in df[df.process == "WJ"].process_decay.unique(): continue
            mc_c_dy[p] = c_dy(df[(df.process_decay == p)], "diff").weight.sum() * scales[p]
            mc_c_tt[p] = c_tt(df[(df.process_decay == p)], "diff").weight.sum() * scales[p] 
            mc_s[p]    = c_ww(df[(df.process_decay == p)], "diff").weight.sum() * scales[p]


            mc_c_dy['W1JetsToLNu'] += c_dy(df[(df.process_decay == p)], "same").weight.sum() * scales[p]
            mc_c_tt['W1JetsToLNu'] += c_tt(df[(df.process_decay == p)], "same").weight.sum() * scales[p] 
            mc_s['W1JetsToLNu']    += c_ww(df[(df.process_decay == p)], "same").weight.sum() * scales[p]


        sf_wj_dy = c_dy(df_da, "same").shape[0] - mc_c_dy['W1JetsToLNu'] 
        sf_wj_tt = c_tt(df_da, "same").shape[0] - mc_c_tt['W1JetsToLNu'] 
        sf_wj_s  = c_ww(df_da, "same").shape[0] - mc_s['W1JetsToLNu']   

        sf_wj_dy = max(0, sf_wj_dy)
        sf_wj_tt = max(0, sf_wj_tt)
        sf_wj_s  = max(0, sf_wj_s)

        sf_wj_dy *= ds_ov_ss 
        sf_wj_tt *= ds_ov_ss  
        sf_wj_s  *= ds_ov_ss  



        mc_c_dy['W1JetsToLNu'] = sf_wj_dy
        mc_c_tt['W1JetsToLNu'] = sf_wj_tt
        mc_s['W1JetsToLNu']    = sf_wj_s
        

        mc_dic["dy"] = mc_c_dy
        mc_dic["tt"] = mc_c_tt
        mc_dic["sig"]= mc_s 
    else:
        func_dic = {"dy": c_dy, "tt": c_tt, "sig": c_ww}
        for control_region in func_dic:
            mc_bins, da_bins, c, d = full_bin_plot(func_dic[control_region](df), func_dic[control_region](df_da), binned_feature, query=None,scales=scales) 
            plt.close(c)

            mc_c = {}
            for process in mc_bins:
                if "plot" in process: continue
                mc_c[process] = mc_bins[process][0].astype(np.float64)
            mc_dic[control_region] = mc_c
                
    
def func_c_dy(mc_dic, aWW, aTT, aDY):
    every_thing_else = 0
    mc_c_dy = mc_dic["dy"]
    
    for i in mc_c_dy:
        if i == "WW" or i == 'ttbar_leptonic' or i == 'DYJetsToLL_M-50' or i == "GluGluWWTo2L2Nu":
            continue
        every_thing_else += mc_c_dy[i]
    return aWW * (mc_c_dy["WW"] + mc_c_dy["GluGluWWTo2L2Nu"]) + aTT * mc_c_dy['ttbar_leptonic'] +  aDY * mc_c_dy['DYJetsToLL_M-50'] + every_thing_else 

def func_c_tt(mc_dic, aWW, aTT, aDY):
    every_thing_else = 0
    mc_c_tt = mc_dic["tt"]
    
    for i in mc_c_tt:
        if i == "WW" or i == 'ttbar_leptonic' or i == 'DYJetsToLL_M-50' or i == "GluGluWWTo2L2Nu":
            continue
        every_thing_else += mc_c_tt[i]
    return aWW * (mc_c_tt["WW"] + mc_c_tt["GluGluWWTo2L2Nu"]) + aTT * mc_c_tt['ttbar_leptonic'] +  aDY * mc_c_tt['DYJetsToLL_M-50'] + every_thing_else 

def func_s(mc_dic, aWW, aTT, aDY):
    every_thing_else = 0
    mc_s = mc_dic["sig"]

    #print "MC_S", mc_s["WW"].sum()
    for i in mc_s:
        if i == "WW" or i == 'ttbar_leptonic' or i == 'DYJetsToLL_M-50' or i == "GluGluWWTo2L2Nu":
            continue
        every_thing_else += mc_s[i]
    return aWW * (mc_s["WW"] + mc_s["GluGluWWTo2L2Nu"]) + aTT * mc_s['ttbar_leptonic'] +  aDY * mc_s['DYJetsToLL_M-50'] + every_thing_else




class MinFunction:
    def __init__(self, mc_control_dic, data_control_dic, verbose=False): 
      self.data_control_dic = data_control_dic
      self.mc_control_dic = mc_control_dic
      self.verbose = verbose
      self.fit_result = None
      self.uncertainty_for = None
      self.order = ["WW", "Top", "DY"]

    def sum(self, v):
        if type(v) == type([]):
            return sum(v)
        elif type(v) == type(np.array([])):
            #Temp solution for is element in data_c_?? is zero
            v[np.isnan(v)] = 0
            v[np.isinf(v)] = 0
            return np.sum(v)
        else:
            return v
    def min_func(self, x):
        data_c_dy = self.data_control_dic["dy"]
        data_c_tt = self.data_control_dic["tt"]
        data_s    = self.data_control_dic["sig"]
        
        fit_c_dy = func_c_dy(self.mc_control_dic, x[0], x[1], x[2]) 
        fit_c_tt = func_c_tt(self.mc_control_dic, x[0], x[1], x[2]) 
        fit_s =    func_s(self.mc_control_dic,    x[0], x[1], x[2])    

        #print x, fit_c_dy , fit_c_tt, fit_s

        fit_c_dy = (fit_c_dy - data_c_dy) / data_c_dy**0.5
        fit_c_tt = (fit_c_tt - data_c_tt) / data_c_tt**0.5
        fit_s    = (fit_s - data_s)     / data_s**0.5
        
        if type(fit_c_dy) == type(np.array([])):
            if self.verbose == True: print x, fit_c_dy[:5]**2 , fit_c_tt[:5]**2, fit_s[:5]**2
        else:
            if self.verbose == True: print x, fit_c_dy , fit_c_tt, fit_s
        constraint_fit = (x[1] - 1)**2/0.05**2 + (x[2] - 1)**2/0.03**2

        if self.verbose == True: print x, self.sum(fit_s**2), self.sum(fit_c_tt**2), self.sum(fit_c_dy**2), constraint_fit
        return   self.sum(fit_s**2) + self.sum(fit_c_tt**2) + self.sum(fit_c_dy**2) + constraint_fit

    def unc_func(self, x):
        
        uncertainty_for = self.uncertainty_for
        fit_params = [1.,1.,1.]
        fit_result = [1.,1.,1.]
        if type(self.fit_result) == type(np.array([])):
            fit_result = self.fit_result
            aWW, aTT, aDY = self.fit_result
        else:
            aWW, aTT, aDY = [1.,1.,1.]
            
        if uncertainty_for == "ww":
           fit_params[0] = x[0] 
           fit_params[1] = aTT 
           fit_params[2] = aDY
        if uncertainty_for == "tt":
           fit_params[0] = aWW 
           fit_params[1] = x[0]
           fit_params[2] = aDY
        if uncertainty_for == "dy":
           fit_params[0] = aWW 
           fit_params[1] = aTT 
           fit_params[2] = x[0]
        if abs(self.min_func(fit_params) - self.min_func(fit_result) - 1) < 0.001:
            print abs(self.min_func(fit_params) - self.min_func(fit_result) - 1), fit_params, fit_result
        return abs(self.min_func(fit_params) - self.min_func(fit_result) - 1)
      

def comprehensive_fit(df, df_da, feature="numb_Bjet", scales=scales):
  data_control_dic = {}
  mc_control_dic = {}
  pack_data_dic(data_control_dic, df_da, feature)
  pack_mc_dic(mc_control_dic, df, df_da, scales, feature)

  func_obj = MinFunction(mc_control_dic, data_control_dic)
  fit_result = minimize(func_obj.min_func, [1.,1., 1.], method='SLSQP')
  return fit_result 


if __name__ == "__main__":
    """
    An example
    """
    import os, sys
    sys.path.append(os.getcwd() + "/../../")
    from prep_ana_II import *
    sys.path.append(os.getcwd() + "/../../tools/")
    sys.path.append(os.getcwd() + "/../../tools/bjetsys/")
    from lepton_eff import muonEff, electronEff
    from pile_up import pileUpFunction, Min_Bias
    from cross_section_calc import calc_cross_stuff, cross_calc, stat_unc_calc, normalization_unc_calc
    import pile_up
    import BTagCalibrationStandalone as BT


    import matplotlib
    import warnings



    scales, df    = load_testset("jet")
    df = pre_cuts(df[(df.lep2_pt > 25) & (df.mll > 30)], diff_charge=False)
    df_da = pre_cuts(load_presel_w_fDY_fTT_DATA(), diff_charge=False)
    df_da = df_da[(df_da.lep2_pt > 25) & (df_da.mll > 30)]

    df = df[df.metFilter_flag == 0]
    df_da = df_da[df_da.metFilter_flag == 0]

    #df_da["weight"] = np.array([1.0] * df_da.shape[0]) 
    warnings.filterwarnings('ignore')





    #############
    #Fit shit
    data_control_dic = {}
    mc_control_dic = {}
    feature = "numb_BJet"
    pack_data_dic(data_control_dic, df_da, feature)
    pack_mc_dic(mc_control_dic, df, df_da, scales, feature)

    func_obj = MinFunction(mc_control_dic, data_control_dic)

    print func_obj.min_func([0,0,0])
    temp1 = minimize(func_obj.min_func, [1.,1., 1.], method='SLSQP')

    print temp1
    func_obj.fit_result = temp1.x
    for i in ["ww", "tt", "dy"]:
      func_obj.uncertainty_for = i 
      if i == "ww":
          t = func_obj.fit_result[0] + 0.01
      if i == "tt":
          t = func_obj.fit_result[1] - 0.001
      if i == "dy":
          t = func_obj.fit_result[2] - 0.001
      print i, t
      temp2 = minimize(func_obj.unc_func, [t], method='SLSQP')
      #print temp2







