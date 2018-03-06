import os, sys, time
#Clean this shit up
sys.path.append(os.getcwd() + "/../")
sys.path.append(os.getcwd() + "/../../")
sys.path.append(os.getcwd() + "/../tools/")
sys.path.append(os.getcwd() + "/../../tools/")
sys.path.append(os.getcwd() + "/../../tools/fit/")
sys.path.append(os.getcwd() + "/../scripts/uncertainties_dir/")

from prep_ana_II import *
from cross_section_calc import calc_cross_stuff, cross_calc, stat_unc_calc, normalization_unc_calc, pseudo_data_yield_sum
from jet_scale import jec_setup, jet_scale_shift_flat
from sklearn.externals import joblib
import fit
import warnings
warnings.filterwarnings("ignore")

#ToDo:
#?Effect on yields messes up with working with same and different flavors

def score_df(df, rfs):
    
    features_fDY = rfs["features_fDY"]
    features_fTT = rfs["features_fTT"]
    clf_fDY = rfs["clf_fDY"]
    clf_fTT = rfs["clf_fTT"]
    
    #Remove undesirable measurements TT
    temp = df[features_fTT]
    temp = temp.replace([np.inf, -np.inf], 0)

    #Predict TT
    pred_fTT = clf_fTT.predict_proba(np.float32(temp.values))
    df["pred_fTT_WW"] = pred_fTT[:,0]

    #Remove undesirable measurements DY
    temp = df[features_fDY]
    temp = temp.replace([np.inf, -np.inf], 0)
    pred_fDY = clf_fDY.predict_proba(np.float32(temp.values))
    df["pred_fDY_WW"] = pred_fDY[:,0]

def yield_string(df, string):
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
        string += process + "\t" +  str(j0) + "\t" + \
                                    str(j1) + "\t" + \
                                    str(j2) + "\n"
    return string, raw_numbers


def compute_ule(flavor, ana_obj, scales=scales):

  scales_ = scales
  print scales_["WW"]
  scales = ana_obj.scales
  df     = ana_obj.df
  df_da  = ana_obj.df_da
  df_ww  = ana_obj.df_ww


  #WW up and down data frames
  df_ww_up   = rp.read_root(data_path+"/ww_up_complete.root", columns=columns)
  df_ww_down = rp.read_root(data_path+"/ww_down_complete.root", columns=columns)


  #jet_scale_shift_flat(df_ww_up,   jet_pt_shift=1.0, pt_cut=30, rf=ana_obj.rfs)
  #jet_scale_shift_flat(df_ww_down, jet_pt_shift=1.0, pt_cut=30, rf=ana_obj.rfs)

  if flavor == "diff":
    df    = df[   df.lep1_type    != df.lep2_type]
    df_da = df_da[df_da.lep1_type != df_da.lep2_type]
    df_ww = df_ww[df_ww.lep1_type != df_ww.lep2_type]
    df_ww_up   = df_ww_up[df_ww_up.lep1_type     != df_ww_up.lep2_type]
    df_ww_down = df_ww_down[df_ww_down.lep1_type != df_ww_down.lep2_type]
  elif flavor == "same":
    df    = df[   df.lep1_type    == df.lep2_type]
    df_da = df_da[df_da.lep1_type == df_da.lep2_type]
    df_ww = df_ww[df_ww.lep1_type == df_ww.lep2_type]
    df_ww_up   = df_ww_up[  df_ww_up.lep1_type   == df_ww_up.lep2_type]
    df_ww_down = df_ww_down[df_ww_down.lep1_type == df_ww_down.lep2_type]



  #Unpack Random forest 
  random_forests = ana_obj.rfs



  df_ww_up["process"] = "WW"
  df_ww_down["process"] = "WW"
  df_ww_up["process_decay"] = "WW"
  df_ww_down["process_decay"] = "WW"

  #df_ww_up["gen_weight"] = df_ww_up.gen_weight.values/ np.abs(df_ww_up.gen_weight.values)
  #df_ww_up["weight"] = df_ww_up.weight.values * df_ww_up.gen_weight.values

  #df_ww_down["gen_weight"] = df_ww_down.gen_weight.values/ np.abs(df_ww_down.gen_weight.values)
  #df_ww_down["weight"]     = df_ww_down.weight.values * df_ww_down.gen_weight.values

  score_df(df_ww_up, random_forests)
  score_df(df_ww_down, random_forests)



  #produce pseudo data
  rf_selection_cuts = "pred_fDY_WW > .9 & pred_fTT_WW > .6"
  pseudo = pseudo_data_yield_sum(df.query(rf_selection_cuts), df_da.query(rf_selection_cuts), scales=scales)
  cross_orig   =  cross_calc(df, df_da, df_ww, None, scales,  fiducial=True, pseudo=pseudo)
  #cross_orig   =  cross_calc(df, df_da, df_ww, ana_obj.df_ggww, scales,  fiducial=True, pseudo=pseudo)
  cross_sections = {}
  yields = {}
  yield_eff = {}

  ww_selection = (df.process_decay == "WW") 
  ww_yield_orig = df[ww_selection & (df.lep1_Charge != df.lep2_Charge)].query(rf_selection_cuts).weight.sum() * scales["WW"]

  print "Xs", cross_orig
  print "Scaled WW yields: Post cuts:", ww_yield_orig, df_ww[df_ww.process_decay == "WW"].weight.sum() * scales_["WW"],  "Eff", ww_yield_orig / (df_ww[df_ww.process_decay == "WW"].weight.sum() * scales_["WW"]) 
  print "Raw WW yields: Post cuts:", ww_yield_orig/ scales["WW"], "Pre cuts:", df_ww[df_ww.process_decay == "WW"].weight.sum(), "RAW eff:", ( ww_yield_orig/ scales["WW"]) /  df_ww[df_ww.process_decay == "WW"].weight.sum() 
  ww_orig_eff = ww_yield_orig / (df_ww[df_ww.process_decay == "WW"].weight.sum() * scales_["WW"]) 
  #print process_yields(df, df_da, scales=scales)

  nominal_ww_total_yields = df_ww[df_ww.process_decay == "WW"].weight.sum() * scales_["WW"] 
  #########################
  #Fits
  nominal_fit_result = fit.comprehensive_fit(df, ana_obj.df_da, "metMod", scales)
  #########################
  fit_results = []

  for it, ww_variant in enumerate([df_ww_down, df_ww_up]):
      df_ = df[df.process_decay != "WW"]
      df_ = pre_cuts(pd.concat([df_, ww_variant[ww_variant.metFilter_flag == 0]]), diff_charge= False)
      
      ww_var_selection = (df_.process_decay == "WW") 

      ww_type = ""
      if it == 0:
          ww_type = "down"
          scales["WW"] = 35.9e3 * 118.7 * (3*.108)**2 / 1987956.
          cross_sections[ww_type] = cross_calc(df_, df_da, ww_variant, None, scales,  fiducial=True, pseudo=pseudo, scale_ww_tot= scales["WW"] )
          down_ww_total_yields = ww_variant.weight.sum() * scales["WW"] 
      else:
          ww_type = "up"
          scales["WW"] = 35.9e3 * 118.7 * (3*.108)**2 / 1866452.
          cross_sections[ww_type] = cross_calc(df_, df_da, ww_variant, None, scales,  fiducial=True, pseudo=pseudo, scale_ww_tot= scales["WW"] )
          up_ww_total_yields = ww_variant.weight.sum() * scales["WW"] 
      
      yields[ww_type] = df_[ww_var_selection & (df_.lep1_Charge != df_.lep2_Charge)].query(rf_selection_cuts).weight.sum() * scales["WW"]
      yield_eff[ww_type] = df_[ww_var_selection & (df_.lep1_Charge != df_.lep2_Charge)].query(rf_selection_cuts).weight.sum() / ww_variant.weight.sum()  
      #########################
      #Fits
      print ww_type, "fit imminent"
      fit_results.append(fit.comprehensive_fit(df_, df_da, "metMod", scales))
      #########################
      
      if it == 1:
        cross_sections[ww_type] = cross_sections[ww_type] 

      print "\nXs", cross_sections[ww_type]
      print "Yield type", ww_type, "Post cuts:", df_[ww_var_selection].query(rf_selection_cuts).weight.sum()*scales["WW"] , "yield_eff", yield_eff[ww_type] , yields[ww_type] 
      print "Raw WW yields", ww_type, "Post cuts:", df_[ww_var_selection].query(rf_selection_cuts).weight.sum(), "Pre cuts", ww_variant.weight.sum()


  print "Orig Xs", cross_orig
  print "Xs:", cross_sections

  print "\n\n"
  print "Results"
  print "Uncertainty on cross section", sum([abs(cross_sections[c] - cross_orig) for c in cross_sections]) / (2 * cross_orig) * 100, "%"
  avg = sum([ cross_sections[c] for c in cross_sections]) / 2.
  print "Uncertainty on cross section 2 ", sum([abs(cross_sections[c] - avg) for c in cross_sections]) / (2 * avg) * 100, "%"
  print "Effect on WW eff: ", sum([abs(float(yield_eff[y]) - ww_orig_eff ) for y in yield_eff]) / (2 * ww_orig_eff) * 100, "%"

  #########################
  print "Fit results"
  fit_processes = ["WW", "Top", "DY"] 
  print fit_processes
  print "nominal: ", nominal_fit_result.x
  print "down: ", fit_results[0].x
  print "up: ", fit_results[1].x
  print "Diffs(%)"
  for it, ele in enumerate(nominal_fit_result.x):
    if fit_processes[it].lower() == "ww":
      print  fit_processes[it], (abs(ele - fit_results[0].x[it] * nominal_ww_total_yields / down_ww_total_yields)  +\
                                 abs(ele - fit_results[1].x[it] * nominal_ww_total_yields / up_ww_total_yields))/2. * 100
    else:
      print  fit_processes[it], (abs(ele - fit_results[0].x[it]) +\
                                 abs(ele - fit_results[1].x[it]))/2. * 100
  #########################
  print "UP totals", nominal_ww_total_yields , up_ww_total_yields, nominal_ww_total_yields / up_ww_total_yields
  print "Down totals", nominal_ww_total_yields , down_ww_total_yields, nominal_ww_total_yields / down_ww_total_yields
  #save to file
  if flavor != "":
    flavor = "_" + flavor
  if flavor == "_both":
    flavor = "_tot"
  f = open("results/jan//ule_DATE" + flavor +".txt", "w")
  f.write( "Orig Xs "+ str(cross_orig) + "\t")
  f.write( "Xs: " + str(cross_sections) + "\n")
  f.write( "Uncertainty on cross section "+ str(sum([abs(cross_sections[c] - cross_orig) for c in cross_sections]) / (2 * cross_orig) * 100) + " %\n")
  f.write( "Uncertainty on cross section 2 "+ str(sum([abs(cross_sections[c] - avg) for c in cross_sections]) / (2 * avg) * 100) + " %\n")
  f.write( "Effect on WW eff:  "+ str(sum([abs(float(yield_eff[y]) - ww_orig_eff) for y in yield_eff]) / (2 * ww_orig_eff) * 100)+ " %")


if __name__ == "__main__":
  ana_obj = analysis_setup(unc="jet")
  ana_obj.apply_pre_cuts() 
  #ana_obj.apply_flat_jet_correction() 
  for flavor in ["both"]:#, "same", "diff"]:
    print flavor
    compute_ule(flavor, ana_obj)




