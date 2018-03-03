#include <cmath>
#include <TRandom.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <Riostream.h>
#include <stdint.h>

int test_tree() {
  // create a TTree
  TTree *tree = new TTree("tree","treelibrated tree");
  float_t arr[110];
  std::string s;
  // create a branch with energy
  tree->Branch("arr",&arr, "&arr[110]/F");
  tree->Branch("s",&s);
  // fill some events with random numbers
  Int_t nevent=10;
  for (Int_t iev=0;iev<nevent;iev++) {
    float_t ea;
    s = "";

    for(int i=0; i < 5; i++){
      ea = fabs( 1 + gRandom->Gaus(0.,.1));  // identical to a.t but a gaussian
      arr[i] = ea;
      uint16_t temp = (uint16_t)round(ea*500);
      s += std::string((char*)&temp, sizeof(uint16_t));
      printf("%f %i ", ea, temp);
    }
    std::cout << s << std::endl;
    
    printf("\n");
    //Decodes string 
    for (int i=0; i<s.length()/2; ++i)
    {
      std::cout << ((uint8_t)s[i*2] | (s[i*2+1] << 8)) << " from " <<
       (uint16_t) s[i*2] << " and " << (uint8_t)  s[i*2 + 1] <<"\n";
    } 
    tree->Fill();  // fill the tree with the current event
  }
  TFile f1("file1.root","recreate");
  tree->Write();
  f1.Close();
  return 0;
}
