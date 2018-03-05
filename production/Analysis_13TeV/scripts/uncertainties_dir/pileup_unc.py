import os, sys, time
sys.path.append(os.getcwd() + "/../../")
sys.path.append(os.getcwd() + "/../../tools/")
sys.path.append(os.getcwd() + "/../../tools/fit/")

from prep_ana_II import * 
from cross_section_calc import calc_cross_stuff, cross_calc, stat_unc_calc, normalization_unc_calc 
import fit
import pile_up
import warnings
np.random.seed(108)


###############################################################
###############################################################
def pileup_unc( ana_obj, flavor='both', mc_flag=1, w_flag=1):

  warnings.filterwarnings('ignore')
  scales = ana_obj.scales

  df = ana_obj.df 
  df_da = ana_obj.df_da 
  df_ww = ana_obj.df_ww 
  df_ggww = ana_obj.df_ggww 
  

  if flavor == "diff":
    df = df[df.lep1_type != df.lep2_type]
    df_da = df_da[df_da.lep1_type != df_da.lep2_type]
    df_ww = df_ww[df_ww.lep1_type != df_ww.lep2_type]
    df_ggww = df_ggww[df_ggww.lep1_type != df_ggww.lep2_type]
  if flavor == "same":
    df = df[df.lep1_type == df.lep2_type]
    df_da = df_da[df_da.lep1_type == df_da.lep2_type]
    df_ww = df_ww[df_ww.lep1_type == df_ww.lep2_type]
    df_ggww = df_ggww[df_ggww.lep1_type == df_ggww.lep2_type]


  orig = cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=True)
  difs = []
  fit_diffs = []
  #########################
  #Fits
  nominal_fit_result = fit.comprehensive_fit(df, df_da, "metMod", scales)
  #########################

  #print "Flags mc:", mc_flag, "weight:", w_flag
  if mc_flag == 0: 
    df      = df[df.metFilter_flag == 0]
    df_da   = df_da[df_da.metFilter_flag == 0]
    df_ww   = df_ww[df_ww.metFilter_flag == 0]
    df_ggww = df_ggww[df_ggww.metFilter_flag == 0]
    

  #print "orig", orig
  for i in [ pile_up.Min_Bias.HIGH, pile_up.Min_Bias.MID, pile_up.Min_Bias.LOW]:

    warnings.filterwarnings('ignore')
    gw = np.ones(df_ww.shape[0])
    gw *= df_ww.gen_weight.values
    l = list(df_ww.gen_pu.values)
    w = np.array(pile_up.pileUpFunction(l, i, w_flag))
    df_ww["w_alt"] = w * gw
    warnings.filterwarnings('ignore')

    gw = np.ones(df_ggww.shape[0])
    gw *= df_ggww.gen_weight.values
    l = list(df_ggww.gen_pu.values)
    w = np.array(pile_up.pileUpFunction(l, i, w_flag))
    df_ggww["w_alt"] = w * gw
    warnings.filterwarnings('ignore')


    gw = np.ones(df.shape[0])
    gw *= df.gen_weight.values
    l = list(df.gen_pu.values)
    w = np.array(pile_up.pileUpFunction(l, i, w_flag))
    df["w_alt"] = w * gw
    warnings.filterwarnings('ignore')


    kwargs = calc_cross_stuff(df, df_da, df_ww, df_ggww, scales, weights="w_alt")#, flavor=flavor)
    #print i, cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=True, **kwargs)

    if i != pile_up.Min_Bias.MID: 
      difs.append( cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=True, **kwargs))
      #########################
      #Fits
      #df["weight"] = df.w_alt
      fit_diffs.append(fit.comprehensive_fit(df, df_da, "metMod", scales))
      #########################
    else: 
      orig = cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=True, **kwargs)
      #########################
      #Fits
      #df["weight"] = df.w_alt
      nominal_fit_result = fit.comprehensive_fit(df, df_da, "metMod", scales)
      #########################

  print "orig Xs:", orig
  print "different Xs:", difs
  #########################
  print "\n\n\nFit results"
  fit_processes = ["WW", "Top", "DY"] 
  print fit_processes
  print "nominal: ", nominal_fit_result.x
  print "up: ", fit_diffs[0].x
  print "down: ", fit_diffs[1].x
  print "Diffs(%)"
  for it, ele in enumerate(nominal_fit_result.x):
    print  fit_processes[it], (abs(ele - fit_diffs[0].x[it]) + abs(ele - fit_diffs[1].x[it]))/2. * 100
  #########################



  f = open("results/jan/pileup_" + flavor +".txt", "w")
  f.write("Orig Xs: " + str(orig) + "\n")
  f.write("Up/Down Xs: " + str(difs) + "\n")
  f.write("Unc(%): " + str(sum([abs(i - orig) for i in difs]) / (2. * orig) * 100.))
  f.write("Percent difference(%): " + str((difs[0] - difs[1]) / (2. * (sum(difs)/2.)) * 100))
  #########################
  f.write("Fit results\n")
  f.write(str(fit_processes) + "\n")
  f.write("nominal: "+ str(nominal_fit_result.x) + "\n")
  f.write("up: "+ str(fit_diffs[0].x) + "\n")
  f.write("down: "+ str(fit_diffs[1].x) + "\n")
  f.write("Diffs(%)" + "\n")
  for it, ele in enumerate(nominal_fit_result.x):
    f.write(str(fit_processes[it]) + str((abs(ele - fit_diffs[0].x[it]) + abs(ele - fit_diffs[1].x[it]))/2. * 100.))
  #########################
  return sum([abs(i - orig) for i in difs]) / (2 * orig) * 100.

###############################################################
###############################################################



def pile_up_unc_study(df=load_presel_w_fDY_fTT_MC(), mc_flag=0, w_flag=0 ):
  difs = []

  print "mc", mc_flag, "w", w_flag
  if mc_flag == 0:
    df = df[df.metFilter_flag == 0]
    

  #print "orig", orig
  results_list_1 = []
  results_list_2 = []
  for i in [ pile_up.Min_Bias.HIGH, pile_up.Min_Bias.MID, pile_up.Min_Bias.LOW]:

    orig_w = np.ones(df.shape[0])
    l = list(df.gen_pu.values)
    w = np.array(pile_up.pileUpFunction(l, i, w_flag))
    df["w_alt"] = w

    print i, "Pre selection", df.w_alt.values.sum(), "Post Selection", pre_cuts(df[(df.pred_fDY_WW > .9) & (df.pred_fTT_WW > .6)], diff_charge= False)["w_alt"].values.sum() 

    results_list_1.append(df.w_alt.values.sum())
    results_list_2.append(pre_cuts(df[(df.pred_fDY_WW > .9) & (df.pred_fTT_WW > .6)], diff_charge= False)["w_alt"].values.sum())


  print "Preselection"
  print pile_up.Min_Bias.HIGH, results_list_1[0] / results_list_1[1] 
  print pile_up.Min_Bias.LOW, results_list_1[2] / results_list_1[1] 


  print "RF selection"
  print pile_up.Min_Bias.HIGH, results_list_2[0] / results_list_2[1] 
  print pile_up.Min_Bias.LOW,  results_list_2[2] / results_list_2[1] 

  print "Diffs"
  print pile_up.Min_Bias.HIGH, results_list_1[0] / results_list_1[1] -  results_list_2[0] / results_list_2[1] 
  print pile_up.Min_Bias.LOW, results_list_1[2] / results_list_1[1]  -  results_list_2[2] / results_list_2[1] 





if __name__ == "__main__":
  #pile_up_unc_lep(mc_flag=0, w_flag=1)
  print "Pile up unc"
  ana_obj = analysis_setup()
  ana_obj.apply_pre_cuts() 
  #ana_obj.apply_flat_jet_correction() 
  for flavor in ["both"]:#, "same", "diff"]:
    print "\n\n\n\nFlavor:", flavor
    print flavor, pileup_unc(ana_obj, mc_flag=0, w_flag=1, flavor=flavor), "%"
