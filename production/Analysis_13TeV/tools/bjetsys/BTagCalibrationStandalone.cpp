//https://github.com/cms-sw/cmssw/tree/CMSSW_8_0_X/CondTools/BTau/test
#include "BTagCalibrationStandalone.h"
#include <iostream>
#include <exception>
#include <algorithm>
#include <sstream>


BTagEntry::Parameters::Parameters(
  OperatingPoint op,
  std::string measurement_type,
  std::string sys_type,
  JetFlavor jf,
  float eta_min,
  float eta_max,
  float pt_min,
  float pt_max,
  float discr_min,
  float discr_max
):
  operatingPoint(op),
  measurementType(measurement_type),
  sysType(sys_type),
  jetFlavor(jf),
  etaMin(eta_min),
  etaMax(eta_max),
  ptMin(pt_min),
  ptMax(pt_max),
  discrMin(discr_min),
  discrMax(discr_max)
{
  std::transform(measurementType.begin(), measurementType.end(),
                 measurementType.begin(), ::tolower);
  std::transform(sysType.begin(), sysType.end(),
                 sysType.begin(), ::tolower);
}

BTagEntry::BTagEntry(const std::string &csvLine)
{
  // make tokens
  std::stringstream buff(csvLine);
  std::vector<std::string> vec;
  std::string token;
  while (std::getline(buff, token, ","[0])) {
    token = BTagEntry::trimStr(token);
    if (token.empty()) {
      continue;
    }
    vec.push_back(token);
  }
  if (vec.size() != 11) {
std::cerr << "ERROR in BTagCalibration: "
          << "Invalid csv line; num tokens != 11: "
          << csvLine;
throw std::exception();
  }

  // clean string values
  char chars[] = " \"\n";
  for (unsigned int i = 0; i < strlen(chars); ++i) {
    vec[1].erase(remove(vec[1].begin(),vec[1].end(),chars[i]),vec[1].end());
    vec[2].erase(remove(vec[2].begin(),vec[2].end(),chars[i]),vec[2].end());
    vec[10].erase(remove(vec[10].begin(),vec[10].end(),chars[i]),vec[10].end());
  }

  // make formula
  formula = vec[10];
  TF1 f1("", formula.c_str());  // compile formula to check validity
  if (f1.IsZombie()) {
std::cerr << "ERROR in BTagCalibration: "
          << "Invalid csv line; formula does not compile: "
          << csvLine;
throw std::exception();
  }

  // make parameters
  unsigned op = stoi(vec[0]);
  if (op > 3) {
std::cerr << "ERROR in BTagCalibration: "
          << "Invalid csv line; OperatingPoint > 3: "
          << csvLine;
throw std::exception();
  }
  unsigned jf = stoi(vec[3]);
  if (jf > 2) {
std::cerr << "ERROR in BTagCalibration: "
          << "Invalid csv line; JetFlavor > 2: "
          << csvLine;
throw std::exception();
  }
  params = BTagEntry::Parameters(
    BTagEntry::OperatingPoint(op),
    vec[1],
    vec[2],
    BTagEntry::JetFlavor(jf),
    stof(vec[4]),
    stof(vec[5]),
    stof(vec[6]),
    stof(vec[7]),
    stof(vec[8]),
    stof(vec[9])
  );
}

BTagEntry::BTagEntry(const std::string &func, BTagEntry::Parameters p):
  formula(func),
  params(p)
{
  TF1 f1("", formula.c_str());  // compile formula to check validity
  if (f1.IsZombie()) {
std::cerr << "ERROR in BTagCalibration: "
          << "Invalid func string; formula does not compile: "
          << func;
throw std::exception();
  }
}

BTagEntry::BTagEntry(const TF1* func, BTagEntry::Parameters p):
  formula(std::string(func->GetExpFormula("p").Data())),
  params(p)
{
  if (func->IsZombie()) {
std::cerr << "ERROR in BTagCalibration: "
          << "Invalid TF1 function; function is zombie: "
          << func->GetName();
throw std::exception();
  }
}

// Creates chained step functions like this:
// "<prevous_bin> : x<bin_high_bound ? bin_value : <next_bin>"
// e.g. "x<0 ? 1 : x<1 ? 2 : x<2 ? 3 : 4"
std::string th1ToFormulaLin(const TH1* hist) {
  int nbins = hist->GetNbinsX();
  TAxis const* axis = hist->GetXaxis();
  std::stringstream buff;
  buff << "x<" << axis->GetBinLowEdge(1) << " ? 0. : ";  // default value
  for (int i=1; i<nbins+1; ++i) {
    char tmp_buff[50];
    sprintf(tmp_buff,
            "x<%g ? %g : ",  // %g is the smaller one of %e or %f
            axis->GetBinUpEdge(i),
            hist->GetBinContent(i));
    buff << tmp_buff;
  }
  buff << 0.;  // default value
  return buff.str();
}

// Creates step functions making a binary search tree:
// "x<mid_bin_bound ? (<left side tree>) : (<right side tree>)"
// e.g. "x<2 ? (x<1 ? (x<0 ? 0:0.1) : (1)) : (x<4 ? (x<3 ? 2:3) : (0))"
std::string th1ToFormulaBinTree(const TH1* hist, int start=0, int end=-1) {
  if (end == -1) {                      // initialize
    start = 0.;
    end = hist->GetNbinsX()+1;
    TH1* h2 = (TH1*) hist->Clone();
    h2->SetBinContent(start, 0);  // kill underflow
    h2->SetBinContent(end, 0);    // kill overflow
    std::string res = th1ToFormulaBinTree(h2, start, end);
    delete h2;
    return res;
  }
  if (start == end) {                   // leave is reached
    char tmp_buff[20];
    sprintf(tmp_buff, "%g", hist->GetBinContent(start));
    return std::string(tmp_buff);
  }
  if (start == end - 1) {               // no parenthesis for neighbors
    char tmp_buff[70];
    sprintf(tmp_buff,
            "x<%g ? %g:%g",
            hist->GetXaxis()->GetBinUpEdge(start),
            hist->GetBinContent(start),
            hist->GetBinContent(end));
    return std::string(tmp_buff);
  }

  // top-down recursion
  std::stringstream buff;
  int mid = (end-start)/2 + start;
  char tmp_buff[25];
  sprintf(tmp_buff,
          "x<%g ? (",
          hist->GetXaxis()->GetBinUpEdge(mid));
  buff << tmp_buff
       << th1ToFormulaBinTree(hist, start, mid)
       << ") : ("
       << th1ToFormulaBinTree(hist, mid+1, end)
       << ")";
  return buff.str();
}

BTagEntry::BTagEntry(const TH1* hist, BTagEntry::Parameters p):
  params(p)
{
  int nbins = hist->GetNbinsX();
  TAxis const* axis = hist->GetXaxis();

  // overwrite bounds with histo values
  if (params.operatingPoint == BTagEntry::OP_RESHAPING) {
    params.discrMin = axis->GetBinLowEdge(1);
    params.discrMax = axis->GetBinUpEdge(nbins);
  } else {
    params.ptMin = axis->GetBinLowEdge(1);
    params.ptMax = axis->GetBinUpEdge(nbins);
  }

  // balanced full binary tree height = ceil(log(2*n_leaves)/log(2))
  // breakes even around 10, but lower values are more propable in pt-spectrum
  if (nbins < 15) {
    formula = th1ToFormulaLin(hist);
  } else {
    formula = th1ToFormulaBinTree(hist);
  }

  // compile formula to check validity
  TF1 f1("", formula.c_str());
  if (f1.IsZombie()) {
std::cerr << "ERROR in BTagCalibration: "
          << "Invalid histogram; formula does not compile (>150 bins?): "
          << hist->GetName();
throw std::exception();
  }
}

std::string BTagEntry::makeCSVHeader()
{
  return "OperatingPoint, "
         "measurementType, "
         "sysType, "
         "jetFlavor, "
         "etaMin, "
         "etaMax, "
         "ptMin, "
         "ptMax, "
         "discrMin, "
         "discrMax, "
         "formula \n";
}

std::string BTagEntry::makeCSVLine() const
{
  std::stringstream buff;
  buff << params.operatingPoint
       << ", " << params.measurementType
       << ", " << params.sysType
       << ", " << params.jetFlavor
       << ", " << params.etaMin
       << ", " << params.etaMax
       << ", " << params.ptMin
       << ", " << params.ptMax
       << ", " << params.discrMin
       << ", " << params.discrMax
       << ", \"" << formula
       << "\" \n";
  return buff.str();
}

std::string BTagEntry::trimStr(std::string str) {
  size_t s = str.find_first_not_of(" \n\r\t");
  size_t e = str.find_last_not_of (" \n\r\t");

  if((std::string::npos == s) || (std::string::npos == e))
    return "";
  else
    return str.substr(s, e-s+1);
}


#include <fstream>
#include <sstream>



BTagCalibration::BTagCalibration(const std::string &taggr):
  tagger_(taggr)
{}

BTagCalibration::BTagCalibration(const std::string &taggr,
                                 const std::string &filename):
  tagger_(taggr)
{
  std::ifstream ifs(filename);
  if (!ifs.good()) {
std::cerr << "ERROR in BTagCalibration: "
          << "input file not available: "
          << filename;
throw std::exception();
  }
  readCSV(ifs);
  ifs.close();
}

void BTagCalibration::addEntry(const BTagEntry &entry)
{
  data_[token(entry.params)].push_back(entry);
}

const std::vector<BTagEntry>& BTagCalibration::getEntries(
  const BTagEntry::Parameters &par) const
{
  std::string tok = token(par);
  if (!data_.count(tok)) {
std::cerr << "ERROR in BTagCalibration: "
          << "(OperatingPoint, measurementType, sysType) not available: "
          << tok;
throw std::exception();
  }
  return data_.at(tok);
}

void BTagCalibration::readCSV(const std::string &s)
{
  std::stringstream buff(s);
  readCSV(buff);
}

void BTagCalibration::readCSV(std::istream &s)
{
  std::string line;

  // firstline might be the header
  getline(s,line);
  if (line.find("OperatingPoint") == std::string::npos) {
    addEntry(BTagEntry(line));
  }

  while (getline(s,line)) {
    line = BTagEntry::trimStr(line);
    if (line.empty()) {  // skip empty lines
      continue;
    }
    addEntry(BTagEntry(line));
  }
}

void BTagCalibration::makeCSV(std::ostream &s) const
{
  s << tagger_ << ";" << BTagEntry::makeCSVHeader();
  for (std::map<std::string, std::vector<BTagEntry> >::const_iterator i
           = data_.cbegin(); i != data_.cend(); ++i) {
    const std::vector<BTagEntry> &vec = i->second;
    for (std::vector<BTagEntry>::const_iterator j
             = vec.cbegin(); j != vec.cend(); ++j) {
      s << j->makeCSVLine();
    }
  }
}

std::string BTagCalibration::makeCSV() const
{
  std::stringstream buff;
  makeCSV(buff);
  return buff.str();
}

std::string BTagCalibration::token(const BTagEntry::Parameters &par)
{
  std::stringstream buff;
  buff << par.operatingPoint << ", "
       << par.measurementType << ", "
       << par.sysType;
  return buff.str();
}




class BTagCalibrationReader::BTagCalibrationReaderImpl
{
  friend class BTagCalibrationReader;

public:
  struct TmpEntry {
    float etaMin;
    float etaMax;
    float ptMin;
    float ptMax;
    float discrMin;
    float discrMax;
    TF1 func;
  };

private:
  BTagCalibrationReaderImpl(BTagEntry::OperatingPoint op,
                            const std::string & sysType,
                            const std::vector<std::string> & otherSysTypes={});

  BTagCalibrationReaderImpl(BTagEntry::OperatingPoint op,
                            const std::string & sysType,
                            const pybind11::list & otherSysTypes);

  void load(const BTagCalibration & c,
            BTagEntry::JetFlavor jf,
            std::string measurementType);

  double eval(BTagEntry::JetFlavor jf,
              float eta,
              float pt,
              float discr) const;

  double eval_auto_bounds(const std::string & sys,
                          BTagEntry::JetFlavor jf,
                          float eta,
                          float pt,
                          float discr) const;

  std::pair<float, float> min_max_pt(BTagEntry::JetFlavor jf,
                                     float eta,
                                     float discr) const;

  BTagEntry::OperatingPoint op_;
  std::string sysType_;
  std::vector<std::vector<TmpEntry> > tmpData_;  // first index: jetFlavor
  std::vector<bool> useAbsEta_;                  // first index: jetFlavor
  std::map<std::string, std::shared_ptr<BTagCalibrationReaderImpl>> otherSysTypeReaders_;
};


BTagCalibrationReader::BTagCalibrationReaderImpl::BTagCalibrationReaderImpl(
                                             BTagEntry::OperatingPoint op,
                                             const std::string & sysType,
                                             const std::vector<std::string> & otherSysTypes):
  op_(op),
  sysType_(sysType),
  tmpData_(3),
  useAbsEta_(3, true)
{
  for (const std::string & ost : otherSysTypes) {
    if (otherSysTypeReaders_.count(ost)) {
std::cerr << "ERROR in BTagCalibration: "
            << "Every otherSysType should only be given once. Duplicate: "
            << ost;
throw std::exception();
    }
    otherSysTypeReaders_[ost] = std::auto_ptr<BTagCalibrationReaderImpl>(
        new BTagCalibrationReaderImpl(op, ost)
    );
  }
}





//Custom Thoth Gunter
BTagCalibrationReader::BTagCalibrationReaderImpl::BTagCalibrationReaderImpl(
                                             BTagEntry::OperatingPoint op,
                                             const std::string & sysType,
                                             const pybind11::list & otherSysTypes): 
  op_(op),
  sysType_(sysType),
  tmpData_(3),
  useAbsEta_(3, true)
{
  for (auto ost : otherSysTypes) {
    if (otherSysTypeReaders_.count(ost.cast<std::string>())) {
std::cerr << "ERROR in BTagCalibration: "
            << "Every otherSysType should only be given once. Duplicate: "
            << ost.cast<std::string>();
throw std::exception();
    }
    otherSysTypeReaders_[ost.cast<std::string>()] = std::auto_ptr<BTagCalibrationReaderImpl>(
        new BTagCalibrationReaderImpl(op, ost.cast<std::string>())
    );
  }
}






void BTagCalibrationReader::BTagCalibrationReaderImpl::load(
                                             const BTagCalibration & c,
                                             BTagEntry::JetFlavor jf,
                                             std::string measurementType)
{
  if (tmpData_[jf].size()) {
std::cerr << "ERROR in BTagCalibration: "
          << "Data for this jet-flavor is already loaded: "
          << jf;
throw std::exception();
  }

  BTagEntry::Parameters params(op_, measurementType, sysType_);
  const std::vector<BTagEntry> &entries = c.getEntries(params);

  for (const auto &be : entries) {
    if (be.params.jetFlavor != jf) {
      continue;
    }

    TmpEntry te;
    te.etaMin = be.params.etaMin;
    te.etaMax = be.params.etaMax;
    te.ptMin = be.params.ptMin;
    te.ptMax = be.params.ptMax;
    te.discrMin = be.params.discrMin;
    te.discrMax = be.params.discrMax;

    if (op_ == BTagEntry::OP_RESHAPING) {
      te.func = TF1("", be.formula.c_str(),
                    be.params.discrMin, be.params.discrMax);
    } else {
      te.func = TF1("", be.formula.c_str(),
                    be.params.ptMin, be.params.ptMax);
    }

    tmpData_[be.params.jetFlavor].push_back(te);
    if (te.etaMin < 0) {
      useAbsEta_[be.params.jetFlavor] = false;
    }
  }

  for (auto & p : otherSysTypeReaders_) {
    p.second->load(c, jf, measurementType);
  }
}

double BTagCalibrationReader::BTagCalibrationReaderImpl::eval(
                                             BTagEntry::JetFlavor jf,
                                             float eta,
                                             float pt,
                                             float discr) const
{
  bool use_discr = (op_ == BTagEntry::OP_RESHAPING);
  if (useAbsEta_[jf] && eta < 0) {
    eta = -eta;
  }

  // search linearly through eta, pt and discr ranges and eval
  // future: find some clever data structure based on intervals
  const auto &entries = tmpData_.at(jf);
  for (unsigned i=0; i<entries.size(); ++i) {
    const auto &e = entries.at(i);
    if (
      e.etaMin <= eta && eta < e.etaMax                   // find eta
      && e.ptMin < pt && pt <= e.ptMax                    // check pt
    ){
      if (use_discr) {                                    // discr. reshaping?
        if (e.discrMin <= discr && discr < e.discrMax) {  // check discr
          return e.func.Eval(discr);
        }
      } else {
        return e.func.Eval(pt);
      }
    }
  }

  return 0.;  // default value
}

double BTagCalibrationReader::BTagCalibrationReaderImpl::eval_auto_bounds(
                                             const std::string & sys,
                                             BTagEntry::JetFlavor jf,
                                             float eta,
                                             float pt,
                                             float discr) const
{
  auto sf_bounds = min_max_pt(jf, eta, discr);
  float pt_for_eval = pt;
  bool is_out_of_bounds = false;

  if (pt < sf_bounds.first) {
    pt_for_eval = sf_bounds.first + .0001;
    is_out_of_bounds = true;
  } else if (pt > sf_bounds.second) {
    pt_for_eval = sf_bounds.second - .0001;
    is_out_of_bounds = true;
  }

  // get central SF (and maybe return)
  double sf = eval(jf, eta, pt_for_eval, discr);
  if (sys == sysType_) {
    return sf;
  }

  // get sys SF (and maybe return)
  if (!otherSysTypeReaders_.count(sys)) {
std::cerr << "ERROR in BTagCalibration: "
        << "sysType not available (maybe not loaded?): "
        << sys;
throw std::exception();
  }
  double sf_err = otherSysTypeReaders_.at(sys)->eval(jf, eta, pt_for_eval, discr);
  if (!is_out_of_bounds) {
    return sf_err;
  }

  // double uncertainty on out-of-bounds and return
  sf_err = sf + 2*(sf_err - sf);
  return sf_err;
}

std::pair<float, float> BTagCalibrationReader::BTagCalibrationReaderImpl::min_max_pt(
                                               BTagEntry::JetFlavor jf,
                                               float eta,
                                               float discr) const
{
  bool use_discr = (op_ == BTagEntry::OP_RESHAPING);
  if (useAbsEta_[jf] && eta < 0) {
    eta = -eta;
  }

  const auto &entries = tmpData_.at(jf);
  float min_pt = -1., max_pt = -1.;
  for (const auto & e: entries) {
    if (
      e.etaMin <= eta && eta < e.etaMax                   // find eta
    ){
      if (min_pt < 0.) {                                  // init
        min_pt = e.ptMin;
        max_pt = e.ptMax;
        continue;
      }

      if (use_discr) {                                    // discr. reshaping?
        if (e.discrMin <= discr && discr < e.discrMax) {  // check discr
          min_pt = min_pt < e.ptMin ? min_pt : e.ptMin;
          max_pt = max_pt > e.ptMax ? max_pt : e.ptMax;
        }
      } else {
        min_pt = min_pt < e.ptMin ? min_pt : e.ptMin;
        max_pt = max_pt > e.ptMax ? max_pt : e.ptMax;
      }
    }
  }

  return std::make_pair(min_pt, max_pt);
}


BTagCalibrationReader::BTagCalibrationReader(BTagEntry::OperatingPoint op,
                                             const std::string & sysType,
                                             const std::vector<std::string> & otherSysTypes):
  pimpl(new BTagCalibrationReaderImpl(op, sysType, otherSysTypes)) {}


//Custom cause we will not wrap vectors
BTagCalibrationReader::BTagCalibrationReader(BTagEntry::OperatingPoint op,
                                             const std::string & sysType,
                                             const pybind11::list & otherSysTypes): 
  pimpl(new BTagCalibrationReaderImpl(op, sysType, otherSysTypes)) {}
  
  /*{
    std::vector<std::string> str_otherSysTypes;
    for(auto other_sys_type: otherSysTypes){
      str_otherSysTypes.push_back( other_sys_type.cast<std::string>() );
    }

  pimpl = std::make_shared<BTagCalibrationReaderImpl>(op, sysType, str_otherSysTypes);
  }*/



void BTagCalibrationReader::load(const BTagCalibration & c,
                                 BTagEntry::JetFlavor jf,
                                 const std::string & measurementType)
{
  pimpl->load(c, jf, measurementType);
}

double BTagCalibrationReader::eval(BTagEntry::JetFlavor jf,
                                   float eta,
                                   float pt,
                                   float discr) const
{
  return pimpl->eval(jf, eta, pt, discr);
}

double BTagCalibrationReader::eval_auto_bounds(const std::string & sys,
                                               BTagEntry::JetFlavor jf,
                                               float eta,
                                               float pt,
                                               float discr) const
{
  return pimpl->eval_auto_bounds(sys, jf, eta, pt, discr);
}

std::pair<float, float> BTagCalibrationReader::min_max_pt(BTagEntry::JetFlavor jf,
                                                          float eta,
                                                          float discr) const
{
  return pimpl->min_max_pt(jf, eta, discr);
}

//#include <pybind11/pybind11.h>
//#include <pybind11/stl.h> 
namespace py = pybind11;

py::list bjetScale( std::string& sys_type, py::list csvlist, py::list flvlist, py::list etalist, py::list ptlist, BTagCalibrationReader& bt_obj, float cut_csv = 0.8484, int flv_interest = -1){

  py::list list_;

  if(!(ptlist.size() == etalist.size() && ptlist.size() == csvlist.size() && csvlist.size() == flvlist.size())){
    printf("!!!BJET SCALE ERROR!!!\npt list and eta list are not the same size!");
    return list_;
  }


  float btag_eff[] = {
  0.5542576310070667,
  0.5606102430319312,
  0.5914605031773699,
  0.6180488071650029,
  0.6344665209139178,
  0.6479497556938635,
  0.6608028574597213,
  0.6644499981441845,
  0.6721201396126378,
  0.6394896941827578,
  0.6479497556938635,
  0.6608028574597213,
  0.6644499981441845,
  0.6721201396126378,
  0.6394896941827578,
  0.616956077630235,
  0.5740658362989324,
  0.5450061652281134,
  0.4839797639123103,
  0.437125748502994,
  0.4032258064516129,
  0.25};

  float cmistag_eff[] = {
  0.098510927868065,
  0.098510927868065,
  0.09679576483700195,
  0.09573161369568556,
  0.09555001039717197,
  0.09905397885364496,
  0.11177688201329121,
  0.10397376543209877,
  0.12468683754095201,
  0.11891679748822606,
  0.11082024432809773,
  0.11702127659574468,
  0.11559139784946236,
  0.1115702479338843,
  0.1,
  0.19047619047619047,
  0.11764705882352941};

  float lmistag_eff[] = {
  0.010224535206605491,
  0.010224535206605491,
  0.0090752110514198,
  0.00806412198817262,
  0.008831463864830328,
  0.008478437431905073,
  0.008668121599606553,
  0.009386300787066582,
  0.011345999464811346,
  0.010502471169686986,
  0.013101476893758932,
  0.018541930046354824,
  0.01507537688442211,
  0.030280649926144758,
  0.015306122448979591,
  0.030054644808743168,
  0.012345679012345678};

   
  float ptmin[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
  float ptmax[] = {30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800, 1000};


  float return_weight;
  float mc_eff_;
  float sf;
  BTagEntry::JetFlavor jet_flavor;
  for(int i =0; i < ptlist.size(); i++){

    sf = 0.0;
    return_weight = 1.0;
    float flv = flvlist[i].cast<int>();
    if(flv == 5) jet_flavor = BTagEntry::JetFlavor::FLAV_B; 
    else if(flv == 4) jet_flavor = BTagEntry::JetFlavor::FLAV_C; 
    else if(flv == 0) jet_flavor = BTagEntry::JetFlavor::FLAV_UDSG; 
    else {
      list_.append(1);
      continue;  
    }


    mc_eff_ = -99;
    for(int it_pt=0; it_pt < 16; it_pt++){

      if( ptlist[i].cast<float>() > ptmin[it_pt] && ptlist[i].cast<float>() <= ptmax[it_pt] )
      {
        if(jet_flavor == BTagEntry::JetFlavor::FLAV_B ) mc_eff_ = btag_eff[it_pt]; 
        if(jet_flavor == BTagEntry::JetFlavor::FLAV_C ) mc_eff_ = cmistag_eff[it_pt]; 
        if(jet_flavor == BTagEntry::JetFlavor::FLAV_UDSG ) mc_eff_ = lmistag_eff[it_pt]; 
      }
    }

    if( ptlist[i].cast<float>() < ptmax[0] ){
      if(jet_flavor == BTagEntry::JetFlavor::FLAV_B ) mc_eff_ = btag_eff[0]; 
      if(jet_flavor == BTagEntry::JetFlavor::FLAV_C ) mc_eff_ = cmistag_eff[0]; 
      if(jet_flavor == BTagEntry::JetFlavor::FLAV_UDSG ) mc_eff_ = lmistag_eff[0]; 
    }
    if( ptlist[i].cast<float>() > ptmax[15] ){ 
      if(jet_flavor == BTagEntry::JetFlavor::FLAV_B ) mc_eff_ = btag_eff[15]; 
      if(jet_flavor == BTagEntry::JetFlavor::FLAV_C ) mc_eff_ = cmistag_eff[15]; 
      if(jet_flavor == BTagEntry::JetFlavor::FLAV_UDSG ) mc_eff_ = lmistag_eff[15]; 
    }

    bool mistagged = true;
    float csv = csvlist[i].cast<float>();
    if( csv > cut_csv ){
      mistagged = false;
      if(jet_flavor != BTagEntry::JetFlavor::FLAV_B) mistagged = true;
    }
    if( csv < cut_csv ) mistagged = true;



      //Central value defaults    
      sf = bt_obj.eval_auto_bounds("central", jet_flavor, etalist[i].cast<float>(), ptlist[i].cast<float>(), 0.); 

    //Default
    if (flv_interest == -1){
      sf = bt_obj.eval_auto_bounds(sys_type, jet_flavor, etalist[i].cast<float>(), ptlist[i].cast<float>(), 0.); 
    }

    // B/C selection
    if(flv >= 4 && flv_interest >= 4){ 
      sf = bt_obj.eval_auto_bounds(sys_type, jet_flavor, etalist[i].cast<float>(), ptlist[i].cast<float>(), 0.); 
    }

    // Light 
    if(flv < 4 && flv_interest == 0){
      sf = bt_obj.eval_auto_bounds(sys_type, jet_flavor, etalist[i].cast<float>(), ptlist[i].cast<float>(), 0.);
      //if (i % 111 == 0 ) std::cout<< "light sf " << sf << std::endl;
    }


    //
    if(mistagged == true){
      if(sf != 0) return_weight = (1 - sf * mc_eff_) / (1 - mc_eff_) ; 
    }
    else{
      if(sf != 0) return_weight = sf;
    }
    if(sf == 0) return_weight = 1;


    //std::cout <<  fabs(return_weight - 1) << " flavor " << flv << " flavor of interst " << flv_interest << " flavor " << jet_flavor << std::endl;
    //std::cout << " sf " << sf << " mc eff " << mc_eff_ << "\n\n" << std::endl;
    

    if(flv < 4 && fabs(return_weight - 1) > 0.15) std::cout << "Return weight: " << return_weight << " sf " << sf << " mc_eff " << mc_eff_ << std::endl;
    list_.append(return_weight);

  }

  return list_;
}




PYBIND11_MODULE(BTagCalibrationStandalone, m)
{

    m.doc() = "BJet scale factors";

 

    py::class_<BTagEntry>(m, "BTagEntry") 
      .def(py::init<>())
      ;

 
    py::enum_<BTagEntry::OperatingPoint>(m, "OperatingPoint")
      .value("OP_LOOSE" , BTagEntry::OperatingPoint::OP_LOOSE)
      .value("OP_MEDIUM", BTagEntry::OperatingPoint::OP_MEDIUM)
      .value("OP_TIGHT" , BTagEntry::OperatingPoint::OP_TIGHT)
      .value("OP_RESHAPING", BTagEntry::OperatingPoint::OP_RESHAPING)
      ;


    py::enum_<BTagEntry::JetFlavor>(m, "JetFlavor")
      .value("FLAV_B" , BTagEntry::JetFlavor::FLAV_B)
      .value("FLAV_C" , BTagEntry::JetFlavor::FLAV_C)
      .value("FLAV_UDSG" , BTagEntry::JetFlavor::FLAV_UDSG)
      ;


    py::class_<BTagCalibrationReader>(m, "BTagCalibrationReader") 
      .def(py::init<>()) 
      .def(py::init<BTagEntry::OperatingPoint&, 
                    const std::string&,
                    const std::vector<std::string>& >(), 
                    py::arg("op") = BTagEntry::OperatingPoint(),
                    py::arg("sysType") = "central",
                    py::arg("otherSysTypes") = py::list() )
      .def("eval_auto_bounds", &BTagCalibrationReader::eval_auto_bounds)
      .def("load", &BTagCalibrationReader::load)
      ;


    //Complete calibration 
    py::class_<BTagCalibration>(m, "BTagCalibration") 
      .def(py::init<const std::string&>(),
                    py::arg("tagger"))
      .def(py::init<const std::string&,
                    const std::string&>(),
                    py::arg("tagger"),
                    py::arg("filename"))
      .def("readCSV", static_cast< void (BTagCalibration::*)(const std::string&)>(&BTagCalibration::readCSV))
      ;

//py::list bjetScale( std::string& sys_type, py::list csvlist, py::list flvlist, py::list etalist, py::list ptlist, BTagCalibrationReader& bt_obj, float cut_csv = 0.8484){
    m.def("bjetScale", bjetScale,
                        py::arg("sys_type"), 
                        py::arg("csvlist"), 
                        py::arg("flvlist"), 
                        py::arg("etalist"), 
                        py::arg("ptlist"), 
                        py::arg("bt_obj"), 
                        py::arg("cut_csv") = 0.8484,
                        py::arg("flv_interest") = -1);
}
















