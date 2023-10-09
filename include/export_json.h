#ifndef EXPORT_TXT_TO_JSON
#define EXPORT_TXT_TO_JSON
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

void json_export (std::ifstream &is, std::ofstream &os) {
  json J;
  std::vector<std::string> variable_names;
  size_t num_var = 0;
  std::string line;
  std::getline(is, line);
  std::stringstream header(line);
  while (!header.eof()){
    variable_names.push_back("");
    header >> variable_names[num_var];
    J[variable_names[num_var]] = std::vector<double>();
    num_var++;
  }

  double num;
  size_t count;
  while (!is.eof()){
    for (size_t i = 0; i < num_var; ++i){
      if (is >> num)
        J[variable_names[i]].push_back(num);
    }
  }
  os << std::setw(4) << J;
  return;
}

void txt2json(std::string const & ifile, std::string const & ofile) {
      std::ifstream is;
      std::ofstream os;
      is.open(ifile);
      os.open(ofile);
      json_export(is, os);
      is.close();
      os.close();
    
}

#endif