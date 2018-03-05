#Thoth Gunter
import os, sys, time
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier

sys.path.append(os.getcwd() + "/../")

from prep_ana_II import *

path = "../data/"


def train_test_split( unc="" ):
  df_train_test = load_preselMC(unc)
  print df_train_test.shape[0]

  train_WW = df_train_test[df_train_test.process=="WW"].sample(frac=0.4)
  train_DY = df_train_test[df_train_test.process=="DY"].sample(n=train_WW.shape[0]*2)
  train_TT = df_train_test[df_train_test.process=="Top"].sample(n=train_WW.shape[0]*2)

  training_indices = []
  for i in train_DY.index:
    training_indices.append(i)

  for i in train_TT.index:
    training_indices.append(i)

  for i in train_WW.index:
    training_indices.append(i)


  train_fTT = pd.concat([train_WW, train_TT])
  train_fDY = pd.concat([train_WW, train_DY])

  print "Total number of Drell Yan events in Drell Yan training sample:", train_DY.shape[0]
  print "Total number of Top events in Top training sample:", train_TT.shape[0]
  print "Total number of WW events in WW sample:", train_WW.shape[0]

  
  test = df_train_test.drop(train_WW.index)
  test = test.drop(train_TT.index)           #Dropping Top training sample and WW training sample
  test = test.drop(train_DY.index)           #Dropping Drell Yan training sample and WW training sample


  train_fDY.to_hdf(path + unc +"train_dy.hdf", "table", complevel=3)
  train_fTT.to_hdf(path + unc +"train_top.hdf", "table", complevel=3)
  test.to_hdf(path + unc + "test.hdf", "table", complevel=3)


  print"\n"
  print "Total number of Drell Yan in testing sample:", test[test.process == "DY"].shape[0]
  print "Total number of Top events in testing sample:", test[test.process == "Top"].shape[0]
  print "Total number of WW events in testing sample:", test[test.process == "WW"].shape[0]
  print "Total number of events in testing sample:", test.shape[0]
  print "Total number of events in testing sample:", df_train_test.drop(training_indices).shape[0], len(training_indices)


  scales_test = {}
  for key in scales.keys():
    if key in df_train_test.process_decay.unique():
      scales_test[key] = scales[key] * ( float(df_train_test[ df_train_test.process_decay == key].shape[0]) / float(test[ test.process_decay == key].shape[0]) ) 
      if key == "WW": print scales[key] , scales_test[key], float(df_train_test[ df_train_test.process_decay == key].shape[0]) , float(test[ test.process_decay == key].shape[0]) 
    else:
      scales_test[key] = scales[key]

  import pickle
  pickle.dump(scales_test,      open(path + "scales_test.pkl", "w"))
  pickle.dump(training_indices, open(path + "training_indices.pkl", "w"))




def test_split_alt_samples():
  training_indices = pickle.load(open(path + "/training_indices.pkl", "r"))

  for unc in  ["", "jet", "lep", "lhe"]:
    print "Loading the:", unc, "data set"
    test = load_preselMC(unc)
    test = test.drop(training_indices)
    print "test shape", test.shape[0], "Number of removed events:", len(training_indices), "\n"
    print "number of WW shape", test[test.process == "WW"].shape[0]
    print "numb of DY shape", test[test.process == "DY"].shape[0]
    print "numb of Top shape", test[test.process == "Top"].shape[0]
    if unc != "":
      unc += "_"
    test.to_hdf(path + unc + "test.hdf", "table", complevel=3)
   

if __name__ == "__main__":

  print "/////////////////////////"
  print "Creating train test split"
  print "/////////////////////////"
  train_test_split()
  print "\n/////////////////////////"
  print "Creating alt samples"
  print "/////////////////////////"
  test_split_alt_samples()











