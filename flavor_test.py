# Thoth Gunter
from prep_ana import *

from sklearn.ensemble import RandomForestClassifier
import pandas as pd

drop_mll = np.where(np.isinf(df[[u'mll']].values))[0]
            
df =df.drop(df.index[drop_mll])

df_s = df[ df.lep_Type < 0 ]
print "WW", df_s[df_s.WW ==1].shape
print df_s[df_s.TT ==1].shape
print df_s[df_s.DY ==1].shape
test_s = df_s[(df_s.WW ==1) | (df_s.DY == 1) | (df_s.TT == 1) | (df_s.ZZ == 1) | (df_s.WZ == 1) ]
print test_s.shape

df_d = df[ df.lep_Type > 0 ]
print "WW",df_d[df_d.WW ==1].shape
print df_d[df_d.TT ==1].shape
print df_d[df_d.DY ==1].shape
test_d = df_d[(df_d.WW ==1) | (df_d.DY == 1) | (df_d.TT == 1) | (df_d.ZZ == 1) | (df_d.WZ == 1) ]
print test_d.shape

train_s = test_s[(test_s.WW ==1) | (test_s.DY == 1) | (test_s.TT == 1)]
train_d = test_d[(test_d.WW ==1) | (test_d.DY == 1) | (test_d.TT == 1)] 

train = pd.concat([ train_s, train_d ])
train = pd.concat( [train[train.WW ==1].sample(n=50000), train[train.DY == 1].sample(n=70000), train[train.TT == 1].sample(n=33000) ] )

train_d = train[ train.lep_Type > 0]
train_s = train[ train.lep_Type < 0]

labels_d = np.empty( train_d.shape[0] )
labels_s = np.empty( train_s.shape[0] )
labels = np.empty( train.shape[0] )

labels_d[(train_d.WW == 1).as_matrix()] = 1
labels_d[(train_d.DY == 1).as_matrix()] = 2
labels_d[(train_d.TT == 1).as_matrix()] = 3

labels_s[(train_s.WW == 1).as_matrix()] = 1
labels_s[(train_s.DY == 1).as_matrix()] = 2
labels_s[(train_s.TT == 1).as_matrix()] = 3

labels[(train.WW == 1).as_matrix()] = 1
labels[(train.DY == 1).as_matrix()] = 2
labels[(train.TT == 1).as_matrix()] = 3


weights_d = np.empty(train_d.shape[0] )
weights_s = np.empty(train_s.shape[0] )
weights = np.empty(train.shape[0] )

weights_d[(train_d.WW == 1).as_matrix()] = 1.
weights_d[(train_d.DY == 1).as_matrix()] = 1.#.01
weights_d[(train_d.TT == 1).as_matrix()] = 1.#.01

weights_s[(train_s.WW == 1).as_matrix()] = 1.
weights_s[(train_s.DY == 1).as_matrix()] = 1.#.01
weights_s[(train_s.TT == 1).as_matrix()] = 1.#.01

weights[(train.WW == 1).as_matrix()] = 1.
weights[(train.DY == 1).as_matrix()] = 1.#.01
weights[(train.TT == 1).as_matrix()] = 1.#.01

train_d["Truth"] = labels_d
train_s["Truth"] = labels_s


features = ['nBJet', 'numb_jets', 'lep1_pt', 'lep2_pt', 'METProj', 'qT', 'mll', 'metMod', 'dPhiLLMET', 'HT', 'lep_Type'] +\
            ['jetPt'+str(i) for i in range(1,7)]


# some if statement 
clf = RandomForestClassifier(n_estimators=100, n_jobs=-1, min_samples_split=10, max_depth=15, max_features='sqrt')
clf = clf.fit( np.float32(train[features].values) , np.float32(labels), sample_weight=weights)

pred = clf.predict_proba(np.float32(df[features].values))
df["pred_WW"] = pred[:,0]
df["pred_DY"] = pred[:,1]
df["pred_TT"] = pred[:,2]


df.to_pickle("all_labels")

'''
for process in scales.keys():
    if process in df.keys():
        print process, df[ (df[process] == 1) & (df.pred > .97912) & (df.lep_Type > 0)].shape[0] *  scales[process],\
        df[ (df[process] == 1) & (df.pred > .97912) & (df.lep_Type < 0)].shape[0] *  scales[process]
'''


