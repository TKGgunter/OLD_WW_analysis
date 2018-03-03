#Thoth Gunter
import os, sys

sys.path.append(os.getcwd() + "/../")

from prep_ana_II import *


if __name__ == "__main__":
  pre_unc = ""
  if len(sys.argv) > 1:
     pre_unc = sys.argv[1]
  print pre_unc

  if pre_unc == "all":
    for i in ["", "jet", "lep", "lhe"]:
      write_preselMC(i)
      write_preselDATA(i)
  else:
    write_preselMC(pre_unc)
    write_preselDATA(pre_unc)







#def pre_cuts( df ): #Charge cut removed
#  """
#    Initial selection apply before preforming analysis. Initial selection includes a same sign Z peak cut, 3 or more lepton cut, b jet cut and invariant mass cut of 30.  
#    pre_cuts( df, diff_charge)
#    df: data frame
#    diff_charge: bool (Apply opposite sign charge cut)
#  """
#  dif_lep = df.lep_Type > 0
#  sam_lep = df.lep_Type < 0
#  z_mass = (df.mll < 81) | (df.mll > 101 ) #altered from (df.mll < 76) | (df.mll > 106 ) for uncertainty calculations
#  nBjet = df.numb_BJet == 0
#  extra_lep = df.numbExtraLep == 0
#  quality_cuts = df.mll > 28 #altered from df.mll > 30 for uncertainty calculations
#
#  lep_pt = df.lep2_pt > 20
#
#  s_df = df[ sam_lep & z_mass & nBjet & extra_lep & quality_cuts & lep_pt]
#  d_df = df[ dif_lep & nBjet & extra_lep & quality_cuts & lep_pt]
#
#
#  return pd.concat([ s_df, d_df ])
#
##MonteCarlo
#df = load_origMC( )
#pre_cuts(df).to_hdf("../data/preselMC.hdf", 'table', complevel=3)
#
##Data
#df_da = load_origDATA()
#pre_cuts(df_da).to_hdf("../data/preselDATA.hdf", 'table', complevel=3)
