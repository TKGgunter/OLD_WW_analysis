#Thoth Gunter
import os, sys, pickle, copy
sys.path.append(os.getcwd() + "/../../")
from prep_ana_II import *





############################################################
############################################################
############################################################
############################################################
#
def yield_string( jer_obj, df, string):
  raw_numbers = {}
  for process in df.process.unique():
    j0 = 0
    j1 = 0
    j2 = 0
    for process_decay in df[df.process == process].process_decay.unique():
      j0 += jer_obj.rf_ana(df[(df.process_decay == process_decay) & (df.numb_jets == 0) & (df.lep1_Charge != df.lep2_Charge)], "pred_fDY_WW", "pred_fTT_WW").weight.sum() * scales[process_decay] 
      j1 += jer_obj.rf_ana(df[(df.process_decay == process_decay) & (df.numb_jets == 1) & (df.lep1_Charge != df.lep2_Charge)], "pred_fDY_WW", "pred_fTT_WW").weight.sum() * scales[process_decay] 
      j2 += jer_obj.rf_ana(df[(df.process_decay == process_decay) & (df.numb_jets == 2) & (df.lep1_Charge != df.lep2_Charge)],"pred_fDY_WW", "pred_fTT_WW" ).weight.sum() * scales[process_decay] 


    raw_numbers[process] = np.array([j0,j1,j2], dtype=float)
    string += process + "\t" +  str(j0) + "\t" +\
                                str(j1) + "\t" +\
                                str(j2) + "\n"  
  return string, raw_numbers



#
def yield_ratio_string(df, string, nominal_numbers, raw_numbers= None):
  if raw_numbers == None:
    yield_temp, raw_numbers = yield_string(df, string)

  raw_ratios = {}

  for process in raw_numbers:
    raw_ratios[process] = {}
    raw_ratios[process]["ratio"] = raw_numbers[process] / nominal_numbers[process] 
    raw_ratios[process]["stat"] = np.sqrt( raw_numbers[process]**2. / nominal_numbers[process]**3. + raw_numbers[process] / nominal_numbers[process]**2. )
     
    string += process + "\t" + str(raw_ratios[process]["ratio"][0]) + "\t+/-  " + str(raw_ratios[process]["stat"][0]) +\
                               str(raw_ratios[process]["ratio"][1]) + "\t+/-  " + str(raw_ratios[process]["stat"][1]) +\
                               str(raw_ratios[process]["ratio"][2]) + "\t+/-  " + str(raw_ratios[process]["stat"][2]) 
  return string, raw_ratios 



#
def make_Xsections_by_jet(df, df_da, df_ww, numb_jets=None, pseudo=None):
  if numb_jets == None:
    cut_numb_jets_df = df.numb_jets > -1 
    cut_numb_jets_da = df_da.numb_jets > -1  
    cut_numb_jets_ww = df_ww.numb_jets > -1  
  else:
    cut_numb_jets_df = df.numb_jets == numb_jets 
    cut_numb_jets_da = df_da.numb_jets == numb_jets 
    cut_numb_jets_ww = df_ww.numb_jets == numb_jets
  

  Xsections = cross_calc(df[cut_numb_jets_df], jer_obj.df_da[cut_numb_jets_da], jer_obj.df_ww[cut_numb_jets_ww], scales,  fiducial=True, pseudo=pseudo) 



def x_per_jet_bin( dic, df, df_da, df_ww, df_ggww, calc_func, scales, args= None): 
  if args == None:
    dic["tot"] = calc_func(df, df_da, scales=scales) 
    dic["j0"] =  calc_func(df[df.numb_jets == 0], df_da[df_da.numb_jets == 0], scales=scales)
    dic["j1"] =  calc_func(df[df.numb_jets == 1], df_da[df_da.numb_jets == 1], scales=scales)
    dic["j2"] =  calc_func(df[df.numb_jets == 2], df_da[df_da.numb_jets == 2], scales=scales)
  else:
    dic["tot"] = calc_func(df, df_da, df_ww, df_ggww, scales=scales, fiducial=args["fiducial"], pseudo=args["pseudo"]["tot"]) 
    dic["j0"] =  calc_func(df[df.numb_jets == 0], df_da[df_da.numb_jets == 0], df_ww[df_ww.numb_jets == 0], df_ggww[df_ggww.numb_jets == 0], scales=scales, fiducial=args["fiducial"], pseudo=args["pseudo"]["j0"] )
    dic["j1"] =  calc_func(df[df.numb_jets == 1], df_da[df_da.numb_jets == 1], df_ww[df_ww.numb_jets == 1], df_ggww[df_ggww.numb_jets == 1], scales=scales, fiducial=args["fiducial"], pseudo=args["pseudo"]["j1"] )
    dic["j2"] =  calc_func(df[df.numb_jets == 2], df_da[df_da.numb_jets == 2], df_ww[df_ww.numb_jets == 2], df_ggww[df_ggww.numb_jets == 2], scales=scales, fiducial=args["fiducial"], pseudo=args["pseudo"]["j2"] )
  #  dic["j3"] =  calc_func(df[df.numb_jets == 3], df_da[df_da.numb_jets == 3], df_ww[df_ww.numb_jets == 3], scales=args["scales"], fiducial=args["fiducial"], pseudo=args["pseudo"]["j3"] )
  #  dic["j4"] =  calc_func(df[df.numb_jets == 4], df_da[df_da.numb_jets == 4], df_ww[df_ww.numb_jets == 4], scales=args["scales"], fiducial=args["fiducial"], pseudo=args["pseudo"]["j4"] )





