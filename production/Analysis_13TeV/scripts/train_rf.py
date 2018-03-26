import os, sys

sys.path.append(os.getcwd() + "/../")

import time
from prep_ana_II import *
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
from validation_plots import set_month

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






def train_rfs():
  #Top quark 
  df_tt_train = pre_cuts(load_tt_trainset(), diff_charge=False)
  features_fTT = ['recoil', 'HT', 'jet1_csv', 'qT', 'numb_jets', 'metMod', 'dPhiLLJet', 'dPhiLLMET', 'dPhiMETJet']
  labels_fTT = np.empty( df_tt_train.shape[0] )

  labels_fTT[(df_tt_train.process == "WW").as_matrix()] = 1
  labels_fTT[(df_tt_train.process == "Top").as_matrix()] = 2


  clf_fTT = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
              max_depth=20, max_features='auto', max_leaf_nodes=None,
              min_impurity_split=1e-07, min_samples_leaf=1,
              min_samples_split=50, min_weight_fraction_leaf=0.0,
              n_estimators=50, n_jobs=-1, oob_score=False, random_state=None,
              verbose=0, warm_start=False)
  print "TT"
  print "nan", np.any(np.isnan(df_tt_train[features_fTT].values)), "inf", np.any(np.isinf(df_tt_train[features_fTT].values))
  clf_fTT.fit(np.float32(df_tt_train[features_fTT].values), np.float32(labels_fTT))



  #DY 
  df_dy_train = pre_cuts(load_dy_trainset(), diff_charge=False)
  features_fDY =['lep_Type', 'metMod', 'METProj', 'qT', 'mllMET', 'dPhiLL', 'mll', 'dPhiLLMET', 'lep2_pt', 'recoil',]

  labels_fDY = np.empty( df_dy_train.shape[0] )

  labels_fDY[(df_dy_train.process == "WW").as_matrix()] = 1
  labels_fDY[(df_dy_train.process == "DY").as_matrix()] = 2


  clf_fDY = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
              max_depth=20, max_features='auto', max_leaf_nodes=None,
              min_impurity_split=1e-07, min_samples_leaf=1,
              min_samples_split=50, min_weight_fraction_leaf=0.0,
              n_estimators=50, n_jobs=-1, oob_score=False, random_state=None,
              verbose=0, warm_start=False)
  print "DY"
  print "nan", np.any(np.isnan(df_dy_train[features_fDY].values)), "inf", np.any(np.isinf(df_dy_train[features_fDY].values))
  print np.where(np.isinf(df_dy_train[features_fDY].values))
  temp = df_dy_train[features_fDY] 
  temp = temp.replace([np.inf, -np.inf], 0)

  #clf_fDY.fit(np.float32(df_dy_train[features_fDY].values), np.float32(labels_fDY))
  clf_fDY.fit(np.float32(temp.values), np.float32(labels_fDY))


  month = set_month(time.strftime("%d/%m/%Y").split("/")[1]) # time.strftime returns date string

  rfs   = {"clf_fDY":clf_fDY, "features_fDY":features_fDY,\
                "clf_fTT":clf_fTT, "features_fTT":features_fTT}
  joblib.dump( rfs, "../data/RF/"+month+"_fDY_fTT.jbl", compress=3)
  return rfs




def score_testsets(pre_unc, rfs):
  scales, df_test = load_testset(pre_unc)
  score_df(df_test, rfs) 

  if pre_unc != "":
    pre_unc += "_"

  df_test.to_hdf(post_data_dir+pre_unc+"test.hdf", 'table', complevel=3)
  print "Completed scoring", df_test[df_test.process == "WW"].shape, df_test[(df_test.process == "WW") & (df_test.pred_fDY_WW > .9) & (df_test.pred_fTT_WW > .6)].shape


def score_datasets(pre_unc, rfs):

  write_preselRFDATA(pre_unc)


if __name__ == "__main__":
  pre_unc = "all"
  if len(sys.argv) > 1:
     pre_unc = sys.argv[1]
  print pre_unc

  print "Train random forest"
  rfs = train_rfs()
  rfs  = analysis_setup().rfs
  print "Score testset data"
  if pre_unc == "all":
    for unc in ["", "jet", "lep", "lhe"]:
      score_testsets(unc, rfs)
      score_datasets(unc, rfs)




