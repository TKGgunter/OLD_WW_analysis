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
    for i in ["jet", "lep", "", "lhe"]:
      print "Currently running over", i
      write_preselRFMC(i)
      write_preselRFDATA(i)
  else:
    write_preselRFMC(pre_unc)
    write_preselRFDATA(pre_unc)


#df = load_preselMC()
#df_da = load_preselDATA()
#
#
#from sklearn.externals import joblib
#
#random_forests = joblib.load("../data/RF/sep_05_fDY_fTT.jbl")
#features_fDY = random_forests["features_fDY"]
#clf_fDY = random_forests["clf_fDY"]
#features_fTT = random_forests["features_fTT"]
#clf_fTT = random_forests["clf_fTT"]
#
#
##Predict MC
#pred_fTT = clf_fTT.predict_proba(np.float32(df[features_fTT].values))
#df["pred_fTT_WW"] = pred_fTT[:,0]
#
#
#temp = df[features_fDY] 
#temp = temp.replace([np.inf, -np.inf], 0)
#pred_fDY = clf_fDY.predict_proba(np.float32(temp.values))
#df["pred_fDY_WW"] = pred_fDY[:,0]
#
##Predict Data
#pred_fTT = clf_fTT.predict_proba(np.float32(df_da[features_fTT].values))
#df_da["pred_fTT_WW"] = pred_fTT[:,0]
#
#temp = df_da[features_fDY] 
#temp = temp.replace([np.inf, -np.inf], 0)
#pred_fDY = clf_fDY.predict_proba(np.float32(temp.values))
#df_da["pred_fDY_WW"] = pred_fDY[:,0]
#
#
##dump data
#df.to_hdf("../data/preselMC_rf.hdf", "table", complevel=3)
#df_da.to_hdf("../data/preselDATA_rf.hdf", "table", complevel=3)





