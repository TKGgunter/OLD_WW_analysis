# Thoth Gunter

import os, sys, pickle, copy
import datetime
sys.path.append(os.getcwd() + "/../../")
sys.path.append(os.getcwd() + "/../../tools/")
sys.path.append(os.getcwd() + "/../../tools/fit/")

from prep_ana_II import *
from jet_scale import jec_setup, jet_scale_shift_flat, pseudo_data_yield_sum 
from cross_section_calc import calc_cross_stuff, cross_calc, stat_unc_calc, normalization_unc_calc
import fit


def create_ana_obj(flavor, dataset):
    ana_obj = analysis_setup(dataset)
    ana_obj.apply_pre_cuts() 
    ana_obj.apply_flat_jet_correction() 

    
    if flavor == "same":
        ana_obj.df = ana_obj.df[ana_obj.df.lep1_type == ana_obj.df.lep2_type]
        ana_obj.df_da = ana_obj.df_da[ana_obj.df_da.lep1_type == ana_obj.df_da.lep2_type]
        ana_obj.df_ww = ana_obj.df_ww[ana_obj.df_ww.lep1_type == ana_obj.df_ww.lep2_type]
    elif flavor == "diff":
        ana_obj.df = ana_obj.df[ana_obj.df.lep1_type != ana_obj.df.lep2_type]
        ana_obj.df_da = ana_obj.df_da[ana_obj.df_da.lep1_type != ana_obj.df_da.lep2_type]
        ana_obj.df_ww = ana_obj.df_ww[ana_obj.df_ww.lep1_type != ana_obj.df_ww.lep2_type]

    return ana_obj


def UncFrameWork( compute_unc= None, flavor=None, dataset=None ):
    if compute_unc == None:
        return "Error"

    #NOTE
    #Checking to make sure user pushed a function
    if type(compute_unc) != type(UncFrameWork):
        return "Error did not supply a function"
    if type(flavor) == type("string"):
        if not (flavor.lower() == "same" or flavor.lower() == "diff" or flavor.lower() == "both"):
            return "Error did not give a proper lepton flavor"

    if type(dataset) != type("ASDF"):
        return "Error did not give a proper dataset"



    ana_obj = create_ana_obj(flavor, dataset)
    


    pseudo = {}
    pseudo["tot"] = pseudo_data_yield_sum(rf_ana(ana_obj.df), rf_ana(ana_obj.df_da), scales=scales)
    pseudo["j0"] = pseudo_data_yield_sum(rf_ana(ana_obj.df), rf_ana(ana_obj.df_da), scales=scales, query="numb_jets == 0")
    pseudo["j1"] = pseudo_data_yield_sum(rf_ana(ana_obj.df), rf_ana(ana_obj.df_da), scales=scales, query="numb_jets == 1")
    pseudo["j2"] = pseudo_data_yield_sum(rf_ana(ana_obj.df), rf_ana(ana_obj.df_da), scales=scales, query="numb_jets == 2")

    cross_results = {}
    fit_results = {}
    yields = ""
    for nominal_up_down in ["nominal", "up", "down"]:
        compute_unc(ana_obj, nominal_up_down)

        cross_results[nominal_up_down] = {"tot":0, "j0":0, "j1":0, "j2":0}
        cross_results[nominal_up_down]["tot"]=  cross_calc(ana_obj.df, ana_obj.df_da, ana_obj.df_ww, ana_obj.df_ggww, scales,  fiducial=True, pseudo=pseudo["tot"])  
        cross_results[nominal_up_down]["j0"] =  cross_calc(ana_obj.df, ana_obj.df_da, ana_obj.df_ww, ana_obj.df_ggww, scales,  fiducial=True, pseudo=pseudo["j0"], query="numb_jets == 0")  
        cross_results[nominal_up_down]["j1"] =  cross_calc(ana_obj.df, ana_obj.df_da, ana_obj.df_ww, ana_obj.df_ggww, scales,  fiducial=True, pseudo=pseudo["j1"], query="numb_jets == 1")  
        cross_results[nominal_up_down]["j2"] =  cross_calc(ana_obj.df, ana_obj.df_da, ana_obj.df_ww, ana_obj.df_ggww, scales,  fiducial=True, pseudo=pseudo["j2"], query="numb_jets == 2")  
        #########################
        #Fits
        fit_results[nominal_up_down] = fit.comprehensive_fit(ana_obj.df, ana_obj.df_da, "pred_fTT_WW;pred_fDY_WW", scales)
        yields += "\n" + nominal_up_down + "\n" + str(process_yields(rf_ana(ana_obj.df), rf_ana(ana_obj.df_da), scales=scales))
        #########################
        ana_obj = create_ana_obj(flavor, dataset)
  

    #Compute the uncertainty with fit and direct cross section results 
    cross_unc_up   =  np.abs(cross_results["up"]["tot"] - cross_results["nominal"]["tot"] ) /  cross_results["nominal"]["tot"]
    cross_unc_down =  np.abs(cross_results["down"]["tot"] - cross_results["nominal"]["tot"] ) /  cross_results["nominal"]["tot"]

    cross_unc_tot  =  np.abs(cross_unc_down - cross_unc_up) / 2.0


    fit_unc_up     =  np.abs(fit_results["up"].x[0]   - fit_results["nominal"].x[0]) / fit_results["nominal"].x[0]
    fit_unc_down   =  np.abs(fit_results["down"].x[0] - fit_results["nominal"].x[0]) / fit_results["nominal"].x[0]

    fit_unc_tot    =  (fit_unc_up - fit_unc_down) / 2.0
    
    print  "Cross unc: ", cross_unc_down, cross_unc_up, cross_unc_tot
    print  "Fit unc: ", fit_unc_down, fit_unc_up, fit_unc_tot

    
    f = open('TEST_ASDF.txt', 'w')
    f.write(yields)




if __name__ == "__main__":

    print "Muon id unc"
    from muon_idiso import muon_computeFW_unc
    print UncFrameWork( compute_unc=muon_computeFW_unc, flavor="both", dataset="lep" )
    




 
