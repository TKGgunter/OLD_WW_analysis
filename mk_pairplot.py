#Thoth Gunter
from prep_ana import *

import seaborn as sns

df_WW = rp.read_root("data/WW_tot_complete.root", "trees_vec")
df_DY = rp.read_root("data/DY_ll_complete.root", "trees_vec")
df_TT = rp.read_root("data/TT_ll_complete.root", "trees_vec")

df = pd.concat( [df_WW, df_TT, df_DY] )
df = df.reset_index()

train = df.sample( frac= 0.33 )
test_df = df.drop( train.index )

features = ['nBJet', 'numb_jets', 'lep1_pt', 'lep2_pt', 'METProj_sin', 'qT', 'mll', 'metMod', 'dPhiLLMET', 'dPhiLLJet', 'HT', 'lep_Type', 'mllMET', 'recoil' ] +\
            ['jetPt'+str(i) for i in range(1,3)] + [ 'jet'+str(i)+'_csv' for i in range(1,3) ]


g = sns.pairplot(train[ (train.metMod < 400) & (train.mll > 50) & (train.mll < 500) ][features + ['process']], hue='process' )

g.savefig("pairPlots.png")
