#Thoth Gunter
import os, sys, time
sys.path.append(os.getcwd() + "/../../")
sys.path.append(os.getcwd() + "/../../tools/")

from prep_ana_II import * 
from cross_section_calc import calc_cross_stuff, cross_calc, stat_unc_calc, normalization_unc_calc, pseudo_data_yield_sum   
import pile_up
np.random.seed(108)



def stat_uncertainty( flavor = ''):
  ana_obj = analysis_setup()
  ana_obj.apply_pre_cuts() 
  #ana_obj.apply_flat_jet_correction() 

  scales   = ana_obj.scales
  df       = ana_obj.df
  df_da    = ana_obj.df_da
  df_ww    = ana_obj.df_ww
  df_ggww  = ana_obj.df_ggww

  if flavor == "diff":
    df    =   df[   df.lep1_type        != df.lep2_type]
    df_da =   df_da[df_da.lep1_type     != df_da.lep2_type]
    df_ww =   df_ww[df_ww.lep1_type     != df_ww.lep2_type]
    df_ggww = df_ggww[df_ggww.lep1_type != df_ggww.lep2_type]
  if flavor == "same":
    df      = df[   df.lep1_type        == df.lep2_type]
    df_da   = df_da[df_da.lep1_type     == df_da.lep2_type]
    df_ww   = df_ww[df_ww.lep1_type     == df_ww.lep2_type]
    df_ggww = df_ggww[df_ggww.lep1_type == df_ggww.lep2_type]


  pseudo = {}
  pseudo["tot"] = pseudo_data_yield_sum(rf_ana(df), rf_ana(df_da), scales=scales)




  signal_cut = "pred_fDY_WW > .9 & pred_fTT_WW > .6" 
  Xs = cross_calc(df, df_da, df_ww, df_ggww, scales, flavor="both", fiducial=True, pseudo=pseudo['tot'])
  stat = stat_unc_calc(df.query(signal_cut), df_da.query(signal_cut), df_ww, df_ggww, scales, unc_mc_process, flavor="both", fiducial=True)
  print "Xs:", Xs, "stat(%):", stat * 100, "stat:", stat*Xs
  return Xs, stat * 100, stat*Xs


if __name__ == "__main__":
  for flavor in ["", "diff", "same"]:
    unc = stat_uncertainty(flavor)
    f = open("results/stat_" + flavor + ".txt", "w")
    f.write("Xs " + str(unc[0]) )
    f.write("Unc " + str(unc[2]) )
    f.write("Unc(%) " + str(unc[1]) )
    f.close()

