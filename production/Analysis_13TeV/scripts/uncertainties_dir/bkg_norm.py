#Thoth Gunter
import os, sys, time
sys.path.append(os.getcwd() + "/../../")
sys.path.append(os.getcwd() + "/../../tools/")

from prep_ana_II import * 
from cross_section_calc import calc_cross_stuff, cross_calc, stat_unc_calc, normalization_unc_calc, pseudo_data_yield_sum  
import pile_up

np.random.seed(108)


def background_uncertainty( flavor = 'both'):
  ana_obj = analysis_setup()
  ana_obj.apply_pre_cuts() 
  ana_obj.apply_flat_jet_correction() 

  scales = ana_obj.scales
  df     = ana_obj.df
  df_da  = ana_obj.df_da
  df_ww  = ana_obj.df_ww
  df_ggww  = ana_obj.df_ggww

  if flavor == "diff":
    df =    df[   df.lep1_type    != df.lep2_type]
    df_da = df_da[df_da.lep1_type != df_da.lep2_type]
    df_ww = df_ww[df_ww.lep1_type != df_ww.lep2_type]
  if flavor == "same":
    df =    df[   df.lep1_type    == df.lep2_type]
    df_da = df_da[df_da.lep1_type == df_da.lep2_type]
    df_ww = df_ww[df_ww.lep1_type == df_ww.lep2_type]


  pseudo = {}
  pseudo["tot"] = pseudo_data_yield_sum(rf_ana(df), rf_ana(df_da), scales=scales)

  normalization_unc = normalization_unc_calc(df, df_da, df_ww, df_ggww, scales, unc_mc_process, fiducial=True, flavor=flavor)
  normalization_percent = normalization_unc / cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=True, flavor=flavor, pseudo=pseudo['tot']) * 100
  print "Systematics ", flavor, "unc: ", normalization_unc, "unc (%)", normalization_percent, "Xs:", cross_calc(df, df_da, df_ww, df_ggww, scales, fiducial=True, pseudo=pseudo['tot'])  
  return normalization_unc, normalization_percent 


if __name__ == "__main__":
  for flavor in ["both"]:#, "same", "diff"]:
    f = open("results/jan/bkg_normalization_" + flavor + ".txt", "w")
    unc = background_uncertainty(flavor= flavor)
    f.write("Unc " + str(unc[0]) )
    f.write("Unc(%) " + str(unc[1]) )
    f.close()
