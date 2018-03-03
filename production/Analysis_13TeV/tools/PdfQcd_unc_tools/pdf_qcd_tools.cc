//Thoth Gunter
#include <iostream>   // std::cout
#include <algorithm>
#include <string> 
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include <stdint.h>

namespace py = pybind11;


//Should move to py::array_t<string>
//http://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html
py::list get_gen_settings( py::list list_string, int element){
  if(element >= 110) element = 109;
  if(element < 0) element = 0;
   
  std::string delimiter = ";";
  py::list ret_list;
  for(auto handle_string : list_string ){


    std::string string = handle_string.cast<std::string>();
    float weight=1.;
    int iter = 0;
    //From stack over flow
    //https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
    size_t pos = 0;
    std::string token;
    while ((pos = string.find(delimiter)) != std::string::npos) {
        token = string.substr(0, pos);
        //std::cout << token << std::endl;
        weight = std::stof(token);

        string.erase(0, pos + delimiter.length());
        if(iter == element) break;
    }
    ret_list.append(weight);  
  }
  return ret_list; 
}







float rmsv( py::list list_string){
  float result = 0.;

  for(int element =9; element < 110; element++){
    float sum_inner = 0;
    int numb_events = 0;

    for(auto handle_string : list_string){
      float temp_float = 1.;
      std::string string = handle_string.cast<std::string>();

      if (string.length() != 0){
        numb_events += 1;
        uint16_t temp_int = (uint16_t) ((uint8_t)string[2 * element]) | ((uint8_t)string[2 * element + 1] << 8);  
        temp_float *= (float)temp_int / 500.;
      }

      sum_inner += temp_float;
    }

   result = pow(float(numb_events) - sum_inner, 2.); 
  }


  return pow(result / 100., .5);
}











py::list get_qcd_weights( py::list list_string, int element){
  py::list list;
  if(element >= 9) element = 8;
  if(element < 0 ) element = 0;
/*Replace above
 ((uint8_t)s[i*2] | (s[i*2+1] << 8)) 
*/
  for(auto handle_string : list_string){
    std::string string = handle_string.cast<std::string>();
    float temp_float = 1.;
    if (string.length() != 0){
      uint16_t temp_int = (uint16_t) ((uint8_t)string[2 * element]) | ((uint8_t)string[2 * element + 1] << 8);  
      temp_float *= (float)temp_int / 500.;
      //if (temp_float > 2.) std::cout << temp_float << " " << temp_int << " " << (uint8_t)string[2 * element] << " " << (uint8_t)string[2 * element + 1] << std::endl;
    }
    list.append(temp_float); 
  }
  return list;
}
py::list get_pdf_weights( py::list list_string, int element){
  py::list list;
  if(element >= 110) element = 109;
  if(element < 0 ) element = 109;
/*Replace above
 ((uint8_t)s[i*2] | (s[i*2+1] << 8)) 
*/
  for(auto handle_string : list_string){
    std::string string = handle_string.cast<std::string>();
    float temp_float = 1.;
    if (string.length() != 0){
      uint16_t temp_int = (uint16_t) ((uint8_t)string[2 * element]) | ((uint8_t)string[2 * element + 1] << 8);  
      temp_float *= (float)temp_int / 500.;
      //if (temp_float > 2.) std::cout << temp_float << " " << temp_int << " " << (uint8_t)string[2 * element] << " " << (uint8_t)string[2 * element + 1] << std::endl;
    }
    list.append(temp_float); 
  }
  return list;
}

py::list calc_unc(py::list list_string){
  py::list list;

  for(auto handle_string : list_string){
    std::string string = handle_string.cast<std::string>();
    uint16_t temp_int = (uint16_t) ((uint8_t)string[0]) | ((uint8_t)string[1] << 8);  
    float nominal = (float)temp_int / 500.;

    float weight = 0;
    for (int element=9; element < string.length()/2; element++ ){
      if(string.length() < 100)
        break;
      float temp_float = 0;

      if (string.length() != 0){
        temp_int = (uint16_t) ((uint8_t)string[2 * element]) | ((uint8_t)string[2 * element + 1] << 8);  
        temp_float = (float)temp_int / 500.;
        weight += pow(nominal - temp_float, 2);
        //if ((nominal - temp_float) / nominal > 1.5) std::cout << nominal << " " << temp_float << " " << std::endl;
      }
    }
    weight = pow(weight / (99. * nominal*nominal), .5); 
    if (weight > 2 ){  
      temp_int = (uint16_t) ((uint8_t)string[20]) | ((uint8_t)string[21] << 8);
      //std::cout << weight << " " << nominal << " " << temp_int << " " << (float)temp_int/500 << std::endl;
    }
    list.append(weight); 
  }
  return list;
}


py::list calc_qcd_unc(py::list list_string){
  py::list list;

  int iterator = 0;
  for(auto handle_string : list_string){
    std::string string = handle_string.cast<std::string>();
    uint16_t temp_int = (uint16_t) ((uint8_t)string[0]) | ((uint8_t)string[1] << 8);  
    float nominal = (float)temp_int / 500.;

    float weight = nominal;
    //std::cout << "event number: " << ++iterator << std::endl;
    for (int element=1; element < 10; element++ ){
      if(string.length() < 100)
        break;
      if(element == 5) continue;
      if(element == 7) continue;
      float temp_float = nominal;

      if (string.length() != 0){
        temp_int = (uint16_t) ((uint8_t)string[2 * element]) | ((uint8_t)string[2 * element + 1] << 8);  
        temp_float = (float)temp_int / 500.;
        //std::cout << nominal << " temp_weight " << temp_float << " " << weight << " " << fabs(nominal - temp_float) << std::endl;
        weight = (fabs(nominal - temp_float) > fabs(nominal - weight)) ? temp_float : weight ;
      }
    }
    list.append(weight); 
  }
  return list;
}

py::list calc_qcd_unc_min_max(py::list list_string, int min_max=0){
  py::list list;

  int iterator = 0;
  for(auto handle_string : list_string){
    std::string string = handle_string.cast<std::string>();
    uint16_t temp_int = (uint16_t) ((uint8_t)string[0]) | ((uint8_t)string[1] << 8);  
    float nominal = (float)temp_int / 500.;

    float weight = nominal;
    //std::cout << "event number: " << ++iterator << std::endl;
    std::vector<float> vector_weights;
    for (int element=1; element < 10; element++ ){
      if(string.length() < 100)
        break;
      if(element == 5) continue;
      if(element == 7) continue;
      float temp_float = nominal;

      if (string.length() != 0){
        temp_int = (uint16_t) ((uint8_t)string[2 * element]) | ((uint8_t)string[2 * element + 1] << 8);  
        temp_float = (float)temp_int / 500.;
        
        vector_weights.push_back(temp_float);
      }
    }
    if(min_max == 0){
      list.append(*std::min_element(vector_weights.begin(), vector_weights.end())); 
    }
    if(min_max == 1){
      list.append(*std::max_element(vector_weights.begin(), vector_weights.end())); 
    }
  }
  return list;
}



PYBIND11_MODULE(pdf_qcd_tools, m)
{
    using namespace py;

    m.doc() = "pybind11 list of string-floats to list of floats";
    m.def("get_gen_settings", &get_gen_settings, "A function returns a list of gen settings.");
    m.def("get_pdf_weights", &get_pdf_weights, "A function returns a list of pdf weights.");
    m.def("get_qcd_weights", &get_qcd_weights, "A function returns a list of qcd weights.");
    m.def("calc_unc", &calc_unc, "A function returns a list of pdf unc weights.");
    m.def("calc_qcd_unc", &calc_qcd_unc, "A function returns a list of qcd unc weights.");
    m.def("calc_qcd_unc_min_max", &calc_qcd_unc_min_max, "A function returns a min max weights.");

}




