#Thoth Gunter
"""
ToDo:
+Use analysis frame work instead of placing reselection cuts by hand
"""
import os, sys, time

sys.path.append(os.getcwd() + "/../")
from prep_ana_II import * 
from physics_selections import pre_cuts

plt.ioff()

def set_month( month):
  dic_month = {"01":"jan", "02": "feb", "03": "mar", "04": "apr", "05": "may", "06": "june", "07": "jul", "08": "aug", "09": "sep", "10": "oct", "11": "nov", "12": "dec"}
  return dic_month[month]

def main(scales=scales):
  if len(sys.argv) < 2:
    print("\x1B[31mError you did not give an argument. cut_level: 'full', 'pre', ...; output_type: 'plot', 'table', 'plot_table' .\x1B[0m")
    exit()

  cut_level = sys.argv[1].lower()


  print("You've selected ", cut_level, "\nWill begin making valedation plots for this cut level.")
  if cut_level == "full":
    df    = load_origMC()
    df_da = load_origDATA()
    q     = "lep2_pt > 25 & mll > 30 & metFilter_flag == 0" 
    df    = df.query(q)
    df_da = df_da.query(q)

  elif cut_level == "pre":
    df = pre_cuts(load_preselMC(), False)
    df_da = pre_cuts(load_preselDATA(), False)
    df = df[df.lep2_pt > 25]
    df_da = df_da[df_da.lep2_pt > 25]

  elif cut_level == "rf":
    ana_obj = analysis_setup("lep")
    scales = ana_obj.scales
    ana_obj.apply_pre_cuts()

    df = ana_obj.df 
    df_da = ana_obj.df_da 

  else:
    print("Please: full or pre or rf")
    exit()

  #Set-up recalc numb of jets, HT for jet pt cut of 30 pt
  #jet_scale_shift_flat(df)
  #jet_scale_shift_flat(df_da)

  def rfDY_ana( df, diff_charge=True ):
    return df[(df.pred_fDY_WW < .6) & (df.pred_fTT_WW > .6)]

  def rfTT_ana( df, diff_charge=True ):
    return df[(df.pred_fDY_WW > .60) & (df.pred_fTT_WW < .60)]

  def rfWW_ana( df, diff_charge=True ):
    return df[(df.pred_fDY_WW > .9) & (df.pred_fTT_WW > .6) ] 

  def rfWW_tt( df, diff_charge=True ):
    return df[(df.pred_fTT_WW > .6) ] 

  def rfWW_dy( df, diff_charge=True ):
    return df[(df.pred_fDY_WW > .9) ] 

  def tot_ana(df, diff_charge=True):
    return df

#
#
  control_regions = None
  #control_regions = {"WW": WW_ana, "TT": TT_ana, "DY": DY_ana, "Z_tt": Z_tt_ana, "same": same_ana, "diff": diff_ana, "full": full_ana, "same_sign": same_sign_ana}
  #make plots 
  #Add number of jets table region
  if "plot" in sys.argv[2].lower():
    month = set_month(time.strftime("%d/%m/%Y").split("/")[1]) # time.strftime returns date string
    day = time.strftime("%d/%m/%Y").split("/")[0]  # time.strftime returns date string
    print("Month: "+ month)

    if cut_level == "rf":

      control_regions = {"rfDY": rfDY_ana, "rfTT": rfTT_ana, "rfWW": rfWW_ana, "tot": tot_ana, "rfWW_tt": rfWW_tt, "rfWW_dy": rfWW_dy}
    make_control_plots(df[df.mll > 30], df_da[df_da.mll > 30], month+"_"+day, cut_level, energy_dir= "13TeV", control_regions= control_regions, scales=scales) 


  if "table" in sys.argv[2].lower():
    #make tables
    if control_regions == None:
      control_regions = {"WW": WW_ana, "TT": TT_ana, "DY": DY_ana, "Z_tt": Z_tt_ana, "same": same_ana, "diff": diff_ana, "full": full_ana}
    if cut_level == "rf":
      control_regions = {"rfDY": rfDY_ana, "rfTT": rfTT_ana, "rfWW": rfWW_ana, "tot": tot_ana, "rfWW_tt": rfWW_tt, "rfWW_dy": rfWW_dy}
 

    for region in control_regions: 
      latex_string = process_yields(control_regions[region](df), control_regions[region](df_da), scales=scales).to_latex(columns=["Process", "Diff Flavor", "Same Flavor"], index=False)
      latex_string = latex_string.split("\n")
      latex_string.insert(0, "\\begin{table}[ht]\n\t\centering\n\t\\topcaption{"+region+" control region for 13\TeV.}" )

      re_dic = {"toprule": "", "midrule": "\t\t\hline", "bottomrule": "\n"}
      for it, ls in enumerate(latex_string):
        for rm_str in re_dic:
          #Remove crap we don't want
          if rm_str in ls:
            latex_string[it] = re_dic[rm_str]
        #Add style to total
        if "Total" in ls:
          latex_string[it] = "\hline\n" + ls


      latex_string = "\n".join(latex_string) 
      latex_string = latex_string + "\n\end{table}"

      month = set_month(time.strftime("%d/%m/%Y").split("/")[1])

      f = open("../tables/"+cut_level+"_"+month+"_"+region+"_table.tex", "w")
      f.write(latex_string)
      f.close()
  

if __name__ == "__main__":
  main()

