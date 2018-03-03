//Thoth Gunrer
#include "TGraph.h"
#include "TFile.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 

enum Min_Bias{
  LOW,
  MID,
  HIGH
};


namespace py = pybind11;

float tempFunction( float gen_pu ){
  float temp[] = {0.0, 0.3147786463115499, 0.8796341236710382, 1.061972001528117, 0.8802282960472673, 1.058114604115144, 1.0826811014683957, 0.7376762489819555, 0.4640722020021043, 0.6911186023764881, 0.832269092316236, 0.913473399623423, 1.0161972002864539, 1.0710041569048485, 1.1243837419706226, 1.1499659410244831, 1.1587551709463875, 1.1554041424116381, 1.144042056182669, 1.1094289579481478, 1.0668274383185186, 1.0404346109901768, 1.029891628406652, 1.0339130061432469, 1.037733889430177, 1.0382417236700447, 1.0511368608367437, 1.068552128905686, 1.0842553859643556, 1.0998779588806775, 1.1193274579437247, 1.109747024568451, 1.1034763199861946, 1.0650194923185534, 1.0104285944506537, 0.9367447964708603, 0.8487355413398124, 0.7435939748906176, 0.6368959509755441, 0.5267835991314975, 0.42666684010567874, 0.32858860190287276, 0.24249304884695433, 0.17484210250062987, 0.12216327717436092, 0.08399232324547708, 0.05571696626397032, 0.03493211395714428, 0.022131963049736067, 0.013622917258736898};
  if( gen_pu > 50 ) return temp[49];

  for( int i =0; i < 50; i++){

    if ((gen_pu <= i+1) && (gen_pu > i)){
      return temp[i];
    }

  }

  return -999;
}

py::list pileUpFunction( py::list list,  Min_Bias type, int metFilter_flag=1 ){

  std::string puFileName = "";
  if (metFilter_flag == 1){
    printf("metFilter on\n");
    //Load File
    puFileName = "/home/gunter/WW_analysis/production/Analysis_13TeV/tools/pileup_data/pileup_sf_2016_full_69216_bins75_metFilter_feb2018.root";
    if( type == Min_Bias::LOW){
      puFileName =  "/home/gunter/WW_analysis/production/Analysis_13TeV/tools/pileup_data/pileup_sf_2016_full_66013_bins75_metFilter_feb2018.root";
    }
    if( type == Min_Bias::HIGH){
      puFileName = "/home/gunter/WW_analysis/production/Analysis_13TeV/tools/pileup_data/pileup_sf_2016_full_72386_bins75_metFilter_feb2018.root";
    }
  }

  else{
    printf("metFilter off\n");
    //Load File
    puFileName = "/home/gunter/WW_analysis/production/Analysis_13TeV/tools/pileup_data/pileup_sf_2016_full_69216_bins75_in.root";
    if( type == Min_Bias::LOW){
      puFileName =  "/home/gunter/WW_analysis/production/Analysis_13TeV/tools/pileup_data/pileup_sf_2016_full_66013_bins75_in.root";
    }
    if( type == Min_Bias::HIGH){
      puFileName = "/home/gunter/WW_analysis/production/Analysis_13TeV/tools/pileup_data/pileup_sf_2016_full_72386_bins75_in.root";
    }
  }

  std::cout << "File name:" << puFileName << std::endl;
  TFile* puFile = new TFile(puFileName.c_str(), "OPEN");
  //Get tgraph
  TGraph* _puReweight = (TGraph*)puFile->Get("pileup_sf");


  py::list list_;

  //Loop over list
  for( auto i : list ){
    //append new list  list value to using eval method
    float gen_pu = i.cast<float>();
    float new_weight = _puReweight->Eval( gen_pu );
    //if( fabs(new_weight - tempFunction(gen_pu)) / new_weight > .05 ) printf( "warning(input, orig, temp): %f %f %f\n", gen_pu, new_weight, tempFunction(gen_pu));

    //new_weight = tempFunction(gen_pu);
    
    list_.append(new_weight);
  }
  //Close File
  puFile->Close();


  return list_;
}



void test(){
  printf("This is a test.");
}
 
PYBIND11_MODULE(pile_up, m)
{
    using namespace py;

    m.doc() = "pybind11 pileup plugin";
    m.def("pileUpFunction", &pileUpFunction, "A function that recalculateds pile up weights.");
    m.def("test", &test, "This is a test to test the functionality of pybind11");

    py::enum_<Min_Bias>(m, "Min_Bias")
      .value("LOW", Min_Bias::LOW)
      .value("MID", Min_Bias::MID)
      .value("HIGH", Min_Bias::HIGH);

    //return m.ptr();
}









