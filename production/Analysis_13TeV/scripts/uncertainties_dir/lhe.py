#Thoth Gunter

import os, sys, pickle, copy
import datetime
sys.path.append(os.getcwd() + "/../../")
sys.path.append(os.getcwd() + "/../../tools/fit/")
sys.path.append(os.getcwd() + "/../../tools/PdfQcd_unc_tools")
from prep_ana_II import *
from pdf_qcd_tools import get_pdf_weights, calc_unc, calc_qcd_unc
from jet_scale import jec_setup, jet_scale_shift_flat, pseudo_data_yield_sum 
from cross_section_calc import calc_cross_stuff, cross_calc, stat_unc_calc, normalization_unc_calc
from tools import x_per_jet_bin, yield_string
import fit
import warnings
warnings.filterwarnings("ignore")


#ToDo:
#?ADD phi and eta jet elements to the default column 
#?Calc pdf variations 
  #?Calc cross sections with new weights and take the rmsd for the PDF
#?Calc qcd uncertainty
  #?



def calc_lhe_weight(df, element):
  weights = np.array(get_pdf_weights(list(df.lhe_weight_string.values), element))
  return weights


def calc_rmsd(some_dic, nominal, jet_bin = 'tot'):
  result = 0
  it = 0.
  for pdf in some_dic:
    print pdf, some_dic[pdf][jet_bin]
    it += 1.
    result += (nominal - some_dic[pdf][jet_bin]) **2.

  result = (result / it)**.5
  return result



def pdf_unc_calc(flavor, jer_obj):
  df = jer_obj.df
  df_ww = jer_obj.df_ww
  df_ggww = jer_obj.df_ggww

  scales = jer_obj.scales

  query_str = None

  if flavor == "diff":
    query_str = "lep1_type != lep2_type"

  if flavor == "same":
    query_str = "lep1_type == lep2_type"


  if query_str != None:
    df = df.query(query_string)
    jer_obj.df_da = jer_obj.df_da.query(query_string)
    df_ww = df_ww.query(query_string)
    df_ggww = df_ggww.query(query_string)

  static_weights      = df.weight.values.copy() 
  static_weights_ww   = df_ww.weight.values.copy()
  static_weights_ggww = df_ggww.weight.values.copy() 
  pseudo = {}
  x_per_jet_bin(pseudo, rf_ana(df), rf_ana(jer_obj.df_da), jer_obj.df_ww, jer_obj.df_ggww, pseudo_data_yield_sum, scales=scales)
  print pseudo

  nominal_dic = {}
  x_per_jet_bin(nominal_dic, df, jer_obj.df_da, jer_obj.df_ww, jer_obj.df_ggww, cross_calc, scales=scales, args={"fiducial":True, "pseudo":pseudo})
  #########################
  #Fits
  nominal_fit_result = fit.comprehensive_fit(df, ana_obj.df_da, "metMod", scales)
  #########################
  print "nominal dic", nominal_dic
  
  diffs = []
  lib_dics = {}
  fit_results = {}
  for process in df.process.unique():
    print "process", process
    post_dic = {}

    ##### 
    df.weight.values[df.process == process] = df.weight.values[df.process == process] * (1.0 + df.pdf_weight.values[df.process == process])
    if process == "WW":
      print "WW process"
      df_ww.weight.values[df_ww.process == process] = df_ww.weight.values[df_ww.process == process] * (1.0 + df_ww.pdf_weight.values[df_ww.process == process])
      df_ggww.weight.values[df_ggww.process == process] = df_ggww.weight.values[df_ggww.process == process] * (1.0 + df_ggww.pdf_weight.values[df_ggww.process == process])
    
    x_per_jet_bin(post_dic, df, jer_obj.df_da, df_ww, df_ggww, cross_calc, scales=scales, args={"fiducial":True, "pseudo":pseudo})
    diffs.append((nominal_dic["tot"] - post_dic["tot"])**2.)
    lib_dics[process] = post_dic["tot"]
    print process, post_dic

    #########################
    #Fits
    fit_results[process] = fit.comprehensive_fit(df, ana_obj.df_da, "metMod", scales)
    #########################



    df["weight"]      = static_weights 
    df_ww["weight"]   = static_weights_ww
    df_ggww["weight"] = static_weights_ggww 
    

  print "Nominal dic", nominal_dic, "\n"
  print "Lib dics", lib_dics, "\n"
  print "Diffs", diffs, "\n"
  print "quad_sum: ", sum(diffs)**.5, "\n"
  print sum(diffs)**.5 / nominal_dic["tot"] * 100.
  #########################
  print "Fit results"
  fit_processes = ["WW", "Top", "DY"] 
  print fit_processes
  print "nominal: ", nominal_fit_result.x
  quad_sum = 0
  for process in fit_results:
    if process == "WW" : continue 
    print process, fit_results[process].x
    quad_sum += (nominal_fit_result.x[0] - fit_results[process].x[0])**2
  print "Diffs(%)"
  print quad_sum**.5 * 100. 
  #########################
  f = open("results/jan/pdf" + flavor + ".txt", "w") 
  f.write("Nominal: "+str(nominal_dic))
  f.write("Process effects:\n"+str(lib_dics))
  f.write("UNC: "+str(sum(diffs)**.5) + "\tUNC(%): " + str(sum(diffs)**.5/nominal_dic["tot"] * 100.))

def pdf_try2(flavor, jer_obj):
  df = jer_obj.df
  df_ww = jer_obj.df_ww
  df_ggww = jer_obj.df_ggww

  if flavor == "diff":
    query_str = "lep1_type != lep2_type"

  if flavor == "same":
    query_str = "lep1_type == lep2_type"


  if query_str != None:
    df = df.query(query_string)
    jer_obj.df_da = jer_obj.df_da.query(query_string)
    df_ww = df_ww.query(query_string)
    df_ggww = df_ggww.query(query_string)

  jer_obj.df = df
  jer_obj.df_ww = df_ww
  jer_obj.df_ggww = df_ggww


  #Create copy of weights for future use
  static_weights_df = copy.copy(df.weight)
  static_weights_ww = copy.copy(jer_obj.df_ww.weight)
  static_weights_ggww = copy.copy(jer_obj.df_ggww.weight)


  pseudo = {}
  x_per_jet_bin(pseudo, df, jer_obj.df_da, jer_obj.df_ww, jer_obj.df_ggww, pseudo_data_yield_sum)
  print pseudo

  nominal_dic = {}
  x_per_jet_bin(nominal_dic, df, jer_obj.df_da, jer_obj.df_ww, jer_obj.df_ggww, cross_calc, args={"scales": scales, "fiducial":True, "pseudo":pseudo})

  ##
  impacts = []
  diffs = []
  for process in df.process.unique():
    print "process", process
    post_dic = {}
    mean_square_diff = 0.

    #Apply Weights
    for element in range(10, 110):
      temp_weights = np.array(get_pdf_weights(list(df[df.process == process].lhe_weight_string), element))
      temp_weights[temp_weights > 20] = temp_weights[temp_weights < 20].mean()
      df.weight.values[df.process == process] = df.weight.values[df.process == process] * temp_weights

      if process == "WW":
        temp_weights_ww = np.array(get_pdf_weights(list(df_ww[df_ww.process == process].lhe_weight_string), element))
        temp_weights_ww[temp_weights_ww > 20] = temp_weights_ww[temp_weights_ww < 20].mean()
        jer_obj.df_ww["weight"] = jer_obj.df_ww.weight * temp_weights_ww

        temp_weights_ggww = np.array(get_pdf_weights(list(df_ww[df_ggww.process == process].lhe_weight_string), element))
        temp_weights_ggww[temp_weights_ggww > 20] = temp_weights_ggww[temp_weights_ggww < 20].mean()
        jer_obj.df_ggww["weight"] = jer_obj.df_ggww.weight * temp_weights_ggww

      #Calc cross section
      x_per_jet_bin(post_dic, df, jer_obj.df_da, jer_obj.df_ww, jer_obj.df_ggww, cross_calc, args={"scales": scales, "fiducial":True, "pseudo":pseudo})

      df["weight"] = static_weights_df
      jer_obj.df_ww["weight"] = static_weights_ww
      jer_obj.df_ggww["weight"] = static_weights_ggww

      #print process, element, post_dic["tot"], (nominal_dic["tot"] - post_dic["tot"]) 
      diffs.append(nominal_dic["tot"] - post_dic["tot"])


      mean_square_diff += (nominal_dic["tot"] - post_dic["tot"])**2 / 100.

    print "nominal", nominal_dic["tot"], "post", post_dic["tot"], "diffs", np.array(diffs).std(),\
          "mean sqr diff", mean_square_diff, "impact:", (mean_square_diff)**.5 / nominal_dic["tot"]
    impacts.append((mean_square_diff)**.5 / nominal_dic["tot"])


  sumquad_impacts = sum([impact**2 for impact in impacts])**.5
  print "FIN", sumquad_impacts


  if flavor != "":
    flavor = "_" + flavor

  f = open("results/pdf" + flavor + ".txt", "w") 
  f.write(str(sumquad_impacts))
  plt.hist(diffs, bins=50)
  #plt.show()




def qcd_calc_print( flavor, jer_obj ):
  df = jer_obj.df
  df_ww = jer_obj.df_ww
  df_ggww = jer_obj.df_ggww
  df_da = jer_obj.df_da
  scales = jer_obj.scales


  if flavor == "diff":
    df =    df[   df.lep1_type    != df.lep2_type]
    df_ww = df_ww[df_ww.lep1_type != df_ww.lep2_type]
    df_ggww = df_ggww[df_ggww.lep1_type != df_ggww.lep2_type]
    jer_obj.df_da = jer_obj.df_da[jer_obj.df_da.lep1_type != jer_obj.df_da.lep2_type]

    jer_obj.df = df
    jer_obj.df_ww = df_ww
    jer_obj.df_ggww = df_ggww
    
  if flavor == "same":
    df =    df[   df.lep1_type    == df.lep2_type]
    df_ww = df_ww[df_ww.lep1_type == df_ww.lep2_type]
    df_ggww = df_ggww[df_ggww.lep1_type == df_ggww.lep2_type]
    jer_obj.df_da = jer_obj.df_da[jer_obj.df_da.lep1_type == jer_obj.df_da.lep2_type]

    jer_obj.df = df
    jer_obj.df_ww = df_ww
    jer_obj.df_ggww = df_ggww

  #Create copy of weights for future use
  static_weights_df = copy.copy(df.weight)
  static_weights_ww = copy.copy(jer_obj.df_ww.weight)
  static_weights_ggww = copy.copy(jer_obj.df_ggww.weight)


  pseudo = {}
  x_per_jet_bin(pseudo, rf_ana(df), rf_ana(df_da), df_ww, df_ggww, pseudo_data_yield_sum, scales=scales)
  print pseudo

  nominal_dic = {}
  x_per_jet_bin(nominal_dic, df, jer_obj.df_da, df_ww, df_ggww, cross_calc, scales=scales, args={"fiducial":True, "pseudo":pseudo})

  #########################
  #Fits
  nominal_fit_result = fit.comprehensive_fit(df, ana_obj.df_da, "metMod", scales)
  #########################

  ##
  impacts = []
  results = ""
  fit_results = {}
  for process in df.process.unique():
    print "process", process
    post_dic = {}
    diff_cross_section = []

    

    _weights = df[df.process == process].qcd_weight1.values

    for qcd_weight in ["qcd_weight"+str(i) for i in range(1,9)]:
      #print "Looking at ", qcd_weight
      if "6" in qcd_weight or "8" in qcd_weight: continue
      temp_weights = df[df.process == process][qcd_weight].values
      _weights[ np.abs(1 - _weights) < np.abs(1 - temp_weights) ] = temp_weights[ np.abs(1 - _weights) < np.abs(1 - temp_weights) ] 

    df.weight.values[(df.process == process)] *= _weights

    if process == "WG": continue
    if process == "WW": continue
    if process == "WW":
      ######pp->WW->2l2n
      _weights = df_ww[df_ww.process == process].qcd_weight.values

      for qcd_weight in ["qcd_weight"+str(i) for i in range(1,9)]:
        temp_weights = df_ww[df_ww.process == process][qcd_weight].values
        _weights[ np.abs(1 - _weights) < np.abs(1 - temp_weights) ] = temp_weights[ np.abs(1 - _weights) < np.abs(1 - temp_weights) ] 
        df_ww.weight.values[(df_ww.process == process)] *= _weights

      ######GG->WW->2l2n
      _weights = df_ggww[df_ggww.process == process].qcd_weight.values

      for qcd_weight in ["qcd_weight"+str(i) for i in range(1,9)]:
        temp_weights = df_ggww[df_ggww.process == process][qcd_weight].values
        _weights[ np.abs(1 - _weights) < np.abs(1 - temp_weights) ] = temp_weights[ np.abs(1 - _weights) < np.abs(1 - temp_weights) ] 
        df_ggww.weight.values[(df_ggww.process == process)] *= _weights

    plt.hist(_weights, bins=50, )
    plt.yscale("log")
    plt.title("Max diff")
    #plt.show()


    x_per_jet_bin(post_dic, df, df_da, df_ww, df_ggww, cross_calc, scales=scales, args={"fiducial":True, "pseudo":pseudo})
    impacts.append((nominal_dic["tot"] - post_dic["tot"]) / nominal_dic["tot"])
    print process, "Mean weights: ", _weights.mean(), "Xs", nominal_dic["tot"], "Xs'", post_dic["tot"], "%Diff", (nominal_dic["tot"] - post_dic["tot"]) / nominal_dic["tot"] * 100
    results +=  process +  " Mean weights: " + str( _weights.mean() ) + " Xs "+ str(nominal_dic["tot"]) + " Xs " + str(post_dic["tot"]) + " %Diff " +  str( (nominal_dic["tot"] - post_dic["tot"]) / nominal_dic["tot"] * 100 ) + "\n"

    #########################
    #Fits
    fit_results[process] = fit.comprehensive_fit(df, ana_obj.df_da, "metMod", scales)
    #########################

    df["weight"] = static_weights_df
    df_ww["weight"] = static_weights_ww  
    df_ggww["weight"] = static_weights_ggww



  sumquad_impacts = sum([impact**2 for impact in impacts])**.5
  print "FIN", sumquad_impacts, impacts
  #########################
  print "Fit results"
  fit_processes = ["WW", "Top", "DY"] 
  print fit_processes
  print "nominal: ", nominal_fit_result.x
  quad_sum = 0
  for process in fit_results:
    #if process == "WW" : continue 
    print process, fit_results[process].x
    quad_sum += (nominal_fit_result.x[0] - fit_results[process].x[0])**2
  print "Diffs(%)"
  print quad_sum**.5 * 100. 
  #########################

  if flavor != "":
    flavor = "_" + flavor

  f = open("results/jan/qcd" + flavor + ".txt", "w") 
  f.write( "UNC: " + str(sumquad_impacts) + " UNC(%) " + str(sumquad_impacts * 100) + "\n")
  f.write(results)



if __name__ == "__main__":
  for flavor in [""]:#, "same", "diff"]:
    ana_obj = analysis_setup("lhe")
    ana_obj.apply_pre_cuts()
    ana_obj.apply_flat_jet_correction() 
    print "PDF", flavor
    pdf_unc_calc(flavor, ana_obj)
    print "\n\n\n", "QCD"
    #qcd_calc_print(flavor, ana_obj)















