#Thoth Gunter
#pileup test example
import os, sys
  
sys.path.append(os.getcwd() + "/../")

from prep_ana_II import * 

import pile_up

df = load_preselMC()

pile_up.pileUpFunction( list(df.gen_pu.values), pile_up.Min_Bias.HIGH, 2)
