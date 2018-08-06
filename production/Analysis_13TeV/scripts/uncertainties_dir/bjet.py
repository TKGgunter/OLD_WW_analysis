# coding: utf-8
import os, sys, pickle, copy
import datetime
sys.path.append(os.getcwd() + "/../../")
from prep_ana_II import *
from physics_selections import pre_cuts
sys.path.append(os.getcwd() + "/../../tools/")
sys.path.append(os.getcwd() + "/../../tools/fit/")
sys.path.append(os.getcwd() + "/../../tools/bjetsys/")
from cross_section_calc import calc_cross_stuff, cross_calc, pseudo_data_yield_sum 

import BTagCalibrationStandalone as BT
from jet_scale import jet_scale_shift_flat
import fit



#ToDo:
#+?Fix how things are saved.  Its nasty and hard to understand



  #LOAD DATA AND BTAGGING CSV FILES
###############################################
print "loading csv data test"
f = open("data/CSVv2_Moriond17_B_H.csv")
strings = f.read()
f.close()
calibration = BT.BTagCalibration("csv")
calibration.readCSV(strings)

reader = BT.BTagCalibrationReader(BT.OperatingPoint.OP_MEDIUM, u"central", [u"up", u"down"])

reader.load(calibration, BT.JetFlavor.FLAV_B,   "comb");
reader.load(calibration, BT.JetFlavor.FLAV_C,   "comb");
reader.load(calibration, BT.JetFlavor.FLAV_UDSG,"incl");
###############################################
#Load data
print "loading data"


ana_obj = analysis_setup()
print process_yields(ana_obj.df, ana_obj.df_da, scales=ana_obj.scales)

ana_obj.apply_pre_cuts() 
ana_obj.apply_flat_jet_correction() 

scales = ana_obj.scales
df = ana_obj.df
df_da = ana_obj.df_da
df_ww = ana_obj.df_ww
df_ggww = ana_obj.df_ggww
rf = ana_obj.rfs


#temp = df[rf["features_fTT"]]
#temp = temp.replace([np.inf,-np.inf], 0)
#df[rf["features_fTT"]] = temp
print process_yields(rf_ana(df), rf_ana(df_da), scales=scales)
"""
df = load_origMC(columns=columns_jets_unc)#?Need to think about this
df_da = load_origDATA(columns=columns_jets_unc)
df_ww = rp.read_root(data_path+"/ww_complete.root", columns=columns_jets_unc)
rf = load_randomForest()
"""



###############################################
def calc_bjet_weight(df, reader, up_down=u"central", flv_interest= -1):
  """
  Go through each jet for every event and apply the appropriate bjet scaling factor to the event.
  """
  b_weights = np.array( [1.]*df.shape[0] )

  for jet_ in range(6):
    if (jet_ + 1 <= 6):
      jet = str(jet_ + 1)
 
      csv_vec = list(df["jet"+jet+"_csv"].values)
      flv_vec = list(df["jet"+jet+"_flv"].values)
      eta_vec = list(df["jet"+jet+"_eta"].values)
      pt_vec  = list(df["jet"+jet+"_pt"].values)

      weights = np.array(BT.bjetScale(up_down,  csv_vec, flv_vec, eta_vec, pt_vec, reader, flv_interest=flv_interest))
      b_weights *= weights


  return b_weights
###############################################



#Set up data with correct ht, numb_jets values
#jet_scale_shift_flat(df, jet_pt_shift=1.0, pt_cut=30, rf=rf)
#jet_scale_shift_flat(df_da, jet_pt_shift=1.0, pt_cut=30, rf=rf)

#Preselection
#df_da = pre_cuts(df_da, diff_charge= False)
#df = pre_cuts(df, diff_charge= False)

#Cross section stuff
print "Original weights: ", df.weight.sum()

nominal =  calc_bjet_weight(df[(df.numb_BJet == 0)  ], reader)
df["weight"] = df.weight.values * nominal
print "Nominal weights: ", df.weight.sum()

#?Make this stuff per jet bin a function 
def x_per_jet_bin( dic, df, df_da, df_ww, df_ggww, calc_func, scales, args= None): 
  if args == None:
    dic["tot"] = calc_func(df, df_da, scales=scales) 
    print dic["tot"], scales["WW"]
    dic["j0"] =  calc_func(df[df.numb_jets == 0], df_da[df_da.numb_jets == 0], scales=scales )
    dic["j1"] =  calc_func(df[df.numb_jets == 1], df_da[df_da.numb_jets == 1], scales=scales )
    dic["j2"] =  calc_func(df[df.numb_jets == 2], df_da[df_da.numb_jets == 2], scales=scales )
  else:
    dic["tot"] = calc_func(df, df_da, df_ww, df_ggww, scales=scales, fiducial=args["fiducial"], pseudo=args["pseudo"]["tot"]) 
    print dic["tot"]
    dic["j0"] =  calc_func(df[df.numb_jets == 0], df_da[df_da.numb_jets == 0], df_ww[df_ww.numb_jets == 0], df_ggww[df_ggww.numb_jets == 0], scales=scales, fiducial=args["fiducial"], pseudo=args["pseudo"]["j0"] )
    dic["j1"] =  calc_func(df[df.numb_jets == 1], df_da[df_da.numb_jets == 1], df_ww[df_ww.numb_jets == 1], df_ggww[df_ggww.numb_jets == 1], scales=scales, fiducial=args["fiducial"], pseudo=args["pseudo"]["j1"] )
    dic["j2"] =  calc_func(df[df.numb_jets == 2], df_da[df_da.numb_jets == 2], df_ww[df_ww.numb_jets == 2], df_ggww[df_ggww.numb_jets == 2], scales=scales, fiducial=args["fiducial"], pseudo=args["pseudo"]["j2"] )



def pretty_print_dic(dic):
  return "\t".join([ k + ": "+str(dic[k]) for k in dic])



def compute_tagging_results(df, df_da, df_ww, df_ggww, flavor="", reader=reader, nominal=nominal):

  pseudo = {}
  x_per_jet_bin(pseudo, rf_ana(df), rf_ana(df_da), df_ww, df_ggww, pseudo_data_yield_sum, scales)

  orig = {"tot":0, "j0":0, "j1":0, "j2":0}
  x_per_jet_bin(orig, df, df_da, df_ww, df_ggww, cross_calc, scales, args={"fiducial":True, "pseudo":pseudo})
  print orig

  print "Nominal WW\n",process_yields(rf_ana(df), rf_ana(df_da), scales=scales)
  print "Nominal DY\n",process_yields(rf_ana_DY(df), rf_ana_DY(df_da), scales=scales)
  print "Nominal TT\n",process_yields(rf_ana_TT(df), rf_ana_TT(df_da), scales=scales)
  #########################
  #Fits
  nominal_fit_result = fit.comprehensive_fit(df, df_da, "metMod", scales)
  #########################

  ###################################################
  ##        B/C    TAGGING 
  up =  calc_bjet_weight(df[(df.numb_BJet == 0)], reader, u"up", flv_interest=5)
  df["weight"] = df.weight.values / nominal * up
  print "\nB/C Up weights: ", df.weight.sum()

  upper = {"tot":0, "j0":0, "j1":0, "j2":0}
  x_per_jet_bin(upper, df, df_da, df_ww, df_ggww, cross_calc, scales, args={"fiducial":True, "pseudo":pseudo})
  #########################
  #Fits
  bc_up_fit_result = fit.comprehensive_fit(df, df_da, "metMod", scales)
  print "Up WW\n",process_yields(rf_ana(df), rf_ana(df_da), scales=scales)
  print "Up DY\n",process_yields(rf_ana_DY(df), rf_ana_DY(df_da), scales=scales)
  print "Up TT\n",process_yields(rf_ana_TT(df), rf_ana_TT(df_da), scales=scales)
  #########################


  down =  calc_bjet_weight(df[(df.numb_BJet == 0) ], reader, u"down", flv_interest=5)
  df["weight"] = df.weight.values / up * down
  print "\nB/C Up weights: ", df.weight.sum()

  downer = {"tot":0, "j0":0, "j1":0, "j2":0}
  x_per_jet_bin(downer, df, df_da, df_ww, df_ggww, cross_calc, scales, args={"fiducial":True, "pseudo":pseudo})
  #########################
  #Fits
  bc_down_fit_result = fit.comprehensive_fit(df, df_da, "metMod", scales)
  print "Down WW\n",process_yields(rf_ana(df), rf_ana(df_da), scales=scales)
  print "Down DY\n",process_yields(rf_ana_DY(df), rf_ana_DY(df_da), scales=scales)
  print "Down TT\n",process_yields(rf_ana_TT(df), rf_ana_TT(df_da), scales=scales)
  #########################

  print orig
  print "BTAG B/C up", upper
  print "BTAB B/C down", downer
  bc_cs_result = ["Orig:\n"+pretty_print_dic(orig)+"\n", "Up:\n"+pretty_print_dic(upper)+"\n", "Down:\n"+pretty_print_dic(downer)+"\n"]  
  bc_cs_unc = abs(orig["tot"] - upper["tot"]) / (2*orig["tot"]) * 100
  ###################################################



  ###################################################
  ##        LIGHT    TAGGING 
  up =  calc_bjet_weight(df[(df.numb_BJet == 0)  ], reader, u"up", flv_interest=0)
  df["weight"] = df.weight.values / down * up
  print "\nLight Up weights: ", df.weight.sum()

  upper = {"tot":0, "j0":0, "j1":0, "j2":0}
  x_per_jet_bin(upper, df, df_da, df_ww, df_ggww, cross_calc, scales, args={"fiducial":True, "pseudo":pseudo})
  #########################
  #Fits
  l_up_fit_result = fit.comprehensive_fit(df, df_da, "metMod", scales)
  #########################


  down =  calc_bjet_weight(df[(df.numb_BJet == 0) ], reader, u"down", flv_interest=0)
  df["weight"] = df.weight.values / up * down
  print "\nLight Down weights: ", df.weight.sum(), "UP/DOWN avg.: ", up.mean(), down.mean()
  #########################
  #Fits
  l_down_fit_result = fit.comprehensive_fit(df, df_da, "metMod", scales)
  #########################

  downer = {"tot":0, "j0":0, "j1":0, "j2":0}
  x_per_jet_bin(downer, df, df_da, df_ww, df_ggww, cross_calc, scales, args={"fiducial":True, "pseudo":pseudo})
  print "BTAG L up", upper
  print "BTAB L down", downer
  l_cs_result = ["Orign:\n"+pretty_print_dic(orig), "Up:\n"+pretty_print_dic(upper), "Down:\n"+pretty_print_dic(downer)]  
  l_cs_unc = abs(orig["tot"] - upper["tot"]) / (2*orig["tot"]) * 100
  ###################################################



  print "///////////////////////////////////////"
  print "///////////////////////////////////////"
  results = {}
  for process in df.process.unique():
    results[process] = {}
    nominal= 0
    up_b   = 0
    down_b = 0
    up_l   = 0
    down_l = 0
    for process_decay in df[df.process == process].process_decay.unique():
      nominal+=  calc_bjet_weight(df[(df.numb_BJet == 0) & (df.process_decay == process_decay) ], reader).sum()
      up_b   +=  calc_bjet_weight(df[(df.numb_BJet == 0) & (df.process_decay == process_decay) ], reader, u"up",   flv_interest= 5).sum()
      down_b +=  calc_bjet_weight(df[(df.numb_BJet == 0) & (df.process_decay == process_decay) ], reader, u"down", flv_interest= 5).sum()
      if "ttbar" in process_decay:
        print process_decay, "\nNominal",  calc_bjet_weight(df[(df.numb_BJet == 0) & (df.process_decay == process_decay) ], reader), "\nLight up",  calc_bjet_weight(df[(df.numb_BJet == 0) & (df.process_decay == process_decay) ], reader, u"up", flv_interest= 0)
      up_l   +=  calc_bjet_weight(df[(df.numb_BJet == 0) & (df.process_decay == process_decay) ], reader, u"up",   flv_interest= 0).sum()
      down_l +=  calc_bjet_weight(df[(df.numb_BJet == 0) & (df.process_decay == process_decay) ], reader, u"down", flv_interest= 0).sum()
    results[process]["up_l"] = up_l / nominal
    results[process]["up_b"] = up_b / nominal
    results[process]["down_l"] = down_l / nominal
    results[process]["down_b"] = down_b / nominal


  for process in results:
    print process, results[process]


  #########################
  print "Fit results"
  fit_processes = ["WW", "Top", "DY"] 
  print fit_processes
  print "nominal: ", nominal_fit_result.x
  print "b/c up: ", bc_up_fit_result.x
  print "b/c down: ", bc_down_fit_result.x
  print "l up: ", l_up_fit_result.x
  print "l down: ", l_down_fit_result.x
  print "b/c Diffs(%)"
  for it, ele in enumerate(nominal_fit_result.x):
    print  fit_processes[it], (abs(ele - bc_up_fit_result.x[it]) + abs(ele - bc_down_fit_result.x[it]))/2. * 100
  print "l Diffs(%)"
  for it, ele in enumerate(nominal_fit_result.x):
    print  fit_processes[it], (abs(ele - l_up_fit_result.x[it]) + abs(ele - l_down_fit_result.x[it]))/2. * 100
  #########################




  date = datetime.date.today() 
  f = open("results/mar/btag_"+ str(date.month) + "_" + str(date.day)+ "_" + flavor + ".txt", "w")
  f.write("\nBC unc(%): " + str(bc_cs_unc) +"\n")
  f.write("B/C Cross-Section\n")
  for i in bc_cs_result:
    f.write(i)
  f.write("\nL unc(%): " + str(l_cs_unc) +"\n")
  f.write("L Cross-Section\n")
  f.write( "\t".join([i for i in l_cs_result]))
  f.write("\nProcess Results\n")
  for i in [process + "\n" + "\t".join([k + " " +str(results[process][k]) for k in results[process]]) for process in results]: 
    f.write("\n"+i)
  #f.write(str([process + "\n" + "\n".join([k + str(results[process][k]) for k in results[process]]) for process in results]))
  #########################
  f.write("\n\n\nFit results\n")
  f.write(str(fit_processes) + "\n")
  f.write("nominal: "+ str(nominal_fit_result.x) + "\n")
  f.write("b/c up: "+ str(bc_up_fit_result.x) + "\n")
  f.write("b/c down: "+ str(bc_down_fit_result.x) + "\n")
  f.write("l up: "+ str(l_up_fit_result.x) + "\n")
  f.write("l down: "+ str(l_down_fit_result.x) + "\n")
  f.write("b/c Diffs(%)" + "\n")
  for it, ele in enumerate(nominal_fit_result.x):
    f.write(str(fit_processes[it]) + str((abs(ele - bc_up_fit_result.x[it]) + abs(ele - bc_down_fit_result.x[it]))/2. * 100.))
  f.write("\nl Diffs(%)" + "\n")
  for it, ele in enumerate(nominal_fit_result.x):
    f.write(str(fit_processes[it]) + " " + str((abs(ele - l_up_fit_result.x[it]) + abs(ele - l_down_fit_result.x[it]))/2. * 100. ) + "\n")
  #########################
  f.close()




if __name__ == "__main__":
  for flavor in ["", "diff", "same"]:
    print "Computing", flavor
    if flavor == "":
      compute_tagging_results(df, df_da, df_ww, df_ggww, flavor=flavor)
    elif flavor == "diff":
      compute_tagging_results(df[(df.lep1_type != df.lep2_type )],\
                              df_da[(df_da.lep1_type != df_da.lep2_type )],\
                              df_ww[(df_ww.lep1_type != df_ww.lep2_type )], df_ggww[(df_ggww.lep1_type != df_ggww.lep2_type )], flavor=flavor, nominal=nominal[df.lep1_type != df.lep2_type])
    elif flavor == "same":
      compute_tagging_results(  df[(df.lep1_type == df.lep2_type )],\
                                df_da[(df_da.lep1_type == df_da.lep2_type )],\
                                df_ww[(df_ww.lep1_type == df_ww.lep2_type )], df_ggww[(df_ggww.lep1_type != df_ggww.lep2_type )], flavor=flavor, nominal=nominal[df.lep1_type == df.lep2_type])



