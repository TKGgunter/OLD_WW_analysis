#Thoth Gunter
import numpy as np
import root_pandas as rp
import pdf_qcd_tools 

#a = rp.read_root("/home/gunter/WW_analysis/data_alpha/ww_complete.root", columns=["lhe_weight_string"])
a = rp.read_root("/home/gunter/WW_analysis/data_alpha/dyjetstoll_m-50_complete.root", columns=["lhe_weight_string"])
temp_list = list(a.lhe_weight_string.values) 

#print pdf_qcd_tools.get_pdf_weights(temp_list, 2)
for i in range(0, 9):
  b =  np.array(pdf_qcd_tools.get_pdf_weights(temp_list, i))
  print "I", i, b.max(), b[b > 120].shape, b.shape[0]
  print b.sum(), a.shape[0], b.sum() / a.shape[0]
  b[b > 120] = b[b < 120].mean()
  print b.sum(), a.shape[0], b.sum() / a.shape[0]
#pdf_qcd_tools.calc_qcd_unc(temp_list)
#b = np.array(pdf_qcd_tools.calc_qcd_unc_min_max(temp_list, 1))
#b[b > 20] = b[b < 20].mean()
#print b.mean()



