//Thoth Gunter
#include "TGraph.h"
#include "TFile.h"



#include "TROOT.h"
#include "TObject.h"                                                          
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h" 
#include "TGraphErrors.h"    
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"  




#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 


namespace py = pybind11;

struct LeptonStruct{
  float Eta;
  float Pt;

};


py::list muonEff( py::list pt_list, py::list eta_list ){

    TGraphAsymmErrors *_muSF_BCDEF_ID_DATA[4];
    TGraphAsymmErrors *_muSF_GH_ID_DATA[4];

    //Load Tight muon id efficiencies
    std::string idFileName = "/home/gunter/WW_analysis/production/Analysis_13TeV/data/leptonEfficiencies/MuonID_EfficienciesAndSF_BCDEF.root";
    TFile* f_muRecoSF2012_ID = new TFile(idFileName.c_str(), "OPEN");

    std::string filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
    _muSF_BCDEF_ID_DATA[0] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
    _muSF_BCDEF_ID_DATA[1] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
    _muSF_BCDEF_ID_DATA[2] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
    _muSF_BCDEF_ID_DATA[3] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin3_DATA").c_str());


    std::string idFileName_ = "/home/gunter/WW_analysis/production/Analysis_13TeV/data/leptonEfficiencies/MuonID_EfficienciesAndSF_GH.root";
    TFile* f_muRecoSF2012_ID_ = new TFile(idFileName_.c_str(), "OPEN");

    filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
    _muSF_GH_ID_DATA[0] = (TGraphAsymmErrors*)f_muRecoSF2012_ID_->Get((filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
    _muSF_GH_ID_DATA[1] = (TGraphAsymmErrors*)f_muRecoSF2012_ID_->Get((filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
    _muSF_GH_ID_DATA[2] = (TGraphAsymmErrors*)f_muRecoSF2012_ID_->Get((filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
    _muSF_GH_ID_DATA[3] = (TGraphAsymmErrors*)f_muRecoSF2012_ID_->Get((filePath + "pt_PLOT_abseta_bin3_DATA").c_str());


    py::list list_;
    for (int it=0; it < pt_list.size(); it++)
    {
        LeptonStruct muon = {eta_list[it].cast<float>(), pt_list[it].cast<float>()};

        float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
        int etaBin = 0; 
        for (int i = 0; i < 4; ++i) {
            if (fabs(muon.Eta) > binningEta[i] && fabs(muon.Eta) <= binningEta[i+1]) {
                etaBin = i;
                break;
            }
        }

        float binningPt[] = {20., 25, 30, 40, 50, 60, 200};
        int ptBin = 0;
        for (int i = 0; i < 6; ++i) {
            if (fabs(muon.Pt) > binningPt[i] && fabs(muon.Pt) <= binningPt[i+1]) {
                ptBin = i;
                break;
            }
        }


        float weight = 1;
        float high = 0;
        float low = 0;

        auto error = [](float a, float b, float a_, float b_){
          return pow( pow(1./b * a_, 2) + pow( a/b * b_, 2) , .5);
        };

        if (rand() > .75) { //Replace with actual ratio lf B-F/ (B-f + GH) 
          float w_data = _muSF_BCDEF_ID_DATA[etaBin]->Eval(muon.Pt);
          weight   *= w_data ;

          high = _muSF_BCDEF_ID_DATA[etaBin]->GetErrorYhigh(ptBin);
          low  = _muSF_BCDEF_ID_DATA[etaBin]->GetErrorYlow(ptBin);
        }
        else{
          float w_data = _muSF_GH_ID_DATA[etaBin]->Eval(muon.Pt);
          weight  *=  w_data;

          high =  _muSF_GH_ID_DATA[etaBin]->GetErrorYhigh(ptBin);
          low  = _muSF_GH_ID_DATA[etaBin]->GetErrorYlow(ptBin);
        }

        list_.append(py::make_tuple(weight, high, low));

    }

    return list_;
}




py::list electronEff( py::list pt_list, py::list eta_list ){

    std::string el_idFileName = "/home/gunter/WW_analysis/production/Analysis_13TeV/data/leptonEfficiencies/egamma_tightSF.root";
    TFile* f_elSF2012_ID = new TFile( el_idFileName.c_str(), "OPEN");
    TH2D* _elSF2012 = (TH2D*)f_elSF2012_ID->Get("EGamma_EffData2D");


    py::list list_;
    for (int it=0; it < pt_list.size(); it++)
    {
        LeptonStruct electron = {eta_list[it].cast<float>(), pt_list[it].cast<float>()};
        float binningEta[] =  {-2.5, -2.0, -1.56, -1.4442, -1.0, 0, 1.0, 1.4442, 1.56, 2.0, 2.5}; //{0., 0.8, 1.442, 1.556, 2., 2.5};
        int etaBin = 0; 
        for (int i = 0; i < 10; ++i) { 
            if (fabs(electron.Eta) > binningEta[i] && fabs(electron.Eta) <= binningEta[i+1]) {
                etaBin = i+1;
                break;
            }
        }

        float binningPt[] = { 10, 20, 30, 40, 50, 2000 };//{10., 15., 20., 30, 40, 50, 200};
        int ptBin = 0;
        for (int i = 0; i < 5; ++i) { 
            if (fabs(electron.Pt) > binningPt[i] && fabs(electron.Pt) <= binningPt[i+1]) {
                ptBin = i+1;
                break;
            }
        }

        float weight = 1;
        float high = 0;
        float low = 0;
        if (electron.Pt < 2000 ){
          weight   *= _elSF2012->GetBinContent(etaBin, ptBin);
          high = _elSF2012->GetBinError(etaBin, ptBin);
          low = high;
        }
        else{ 
          weight *= _elSF2012->GetBinContent(etaBin, ptBin);
          high    = _elSF2012->GetBinError(etaBin, ptBin);
          low = high;
        }
        list_.append(py::make_tuple(weight, high, low));
    }

    return list_;
}



void test(){
  printf("This is a test.");
}


PYBIND11_MODULE(lepton_eff, m)
{

    m.doc() = "pybind11 pileup plugin";
    m.def("muonEff", &muonEff, "A function that calculates the muon data efficiences for 13TeV.");
    m.def("electronEff", &electronEff, "A function that calculates the electron data efficiences for 13TeV.");
    m.def("test", &test, "This is a test to test the functionality of pybind11");

}

