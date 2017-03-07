#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include "progressBar.h"

// [[Rcpp::plugins(cpp11)]]s

// Macro genotypes to encode them as a char so
// comp opt will generate a fast jumping table
#define A1 1
#define A2 2
#define A3 3
#define A4 4
#define B15 5
#define B26 6
#define B37 7
#define C8 8
#define D19 9
#define D110 10
#define D111 11
#define D112 12
#define D113 13
#define D214 14
#define D215 15
#define D216 16
#define D217 17
#define D218 18

static std::map<std::string, char> segregation = {
  {"A.1", A1}, {"A.2", A2}, {"A.3", A3}, {"A.4", A4},
  {"B1.5", B15}, {"B2.6", B26}, {"B3.7", B37},
  {"C.8", C8}, {"D1.9", D19}, {"D1.10", D110},
  {"D1.11", D111}, {"D1.12", D112}, {"D1.13", D113},
  {"D2.14", D214}, {"D2.15", D215}, {"D2.16", D216},
  {"D2.17", D217}, {"D2.18", D218}
};

static std::map<std::pair<char, std::string>, unsigned short> codify_geno = {
  {{A1, "ac"}, 1}, {{A1, "ad"}, 2}, {{A1, "bc"}, 3}, {{A1, "bd"}, 4}, {{A1, "-"}, 0},
  {{A2, "a"}, 1},  {{A2, "ac"}, 2}, {{A2, "ba"}, 3}, {{A2, "bc"}, 4}, {{A2, "-"}, 0},
  {{A3, "ac"}, 1}, {{A3, "a"}, 2},  {{A3, "bc"}, 3}, {{A3, "b"}, 4}, {{A3, "-"}, 0},
  {{A4, "ab"}, 1}, {{A4, "a"}, 2},  {{A4, "b"}, 3},  {{A4, "o"}, 4}, {{A4, "-"}, 0},
  {{B15, "a"}, 1}, {{B15, "ab"}, 2}, {{B15, "b"}, 3}, {{B15, "-"}, 0},
  {{B26, "a"}, 1}, {{B26, "ab"}, 2}, {{B26, "b"}, 3}, {{B26, "-"}, 0},
  {{B37, "a"}, 1}, {{B37, "ab"}, 2}, {{B37, "b"}, 3}, {{B37, "-"}, 0},
  {{C8, "a"}, 1}, {{C8, "o"}, 2}, {{C8, "-"}, 0},
  {{D19, "ac"}, 1}, {{D19, "bc"}, 2}, {{D19, "-"}, 0},
  {{D110, "a"}, 1}, {{D110, "ab"}, 2}, {{D110, "-"}, 0},
  {{D111, "a"}, 1}, {{D111, "b"}, 2}, {{D111, "-"}, 0},
  {{D112, "ab"}, 1}, {{D112, "a"}, 2}, {{D112, "-"}, 0},
  {{D113, "a"}, 1}, {{D113, "o"}, 2}, {{D113, "-"}, 0},
  {{D214, "ac"}, 1}, {{D214, "bc"}, 2}, {{D214, "-"}, 0},
  {{D215, "a"}, 1}, {{D215, "ab"}, 2}, {{D215, "-"}, 0},
  {{D216, "a"}, 1}, {{D216, "b"}, 2}, {{D216, "-"}, 0},
  {{D217, "ab"}, 1}, {{D217, "a"}, 2}, {{D217, "-"}, 0},
  {{D218, "a"}, 1}, {{D218, "o"}, 2}, {{D218, "-"}, 0}
};

unsigned short codify_segregation(char segregation)
{

  switch (segregation)
  {
    case A1:
    return 1;
    case A2:
    return 1;
    case A3:
    return 1;
    case A4:
    return 1;
    case B15:
    return 2;
    case B26:
    return 3;
    case B37:
    return 4;
    case C8:
    return 5;
    case D19:
    return 6;
    case D110:
    return 6;
    case D111:
    return 6;
    case D112:
    return 6;
    case D113:
    return 6;
    case D214:
    return 7;
    case D215:
    return 7;
    case D216:
    return 7;
    case D217:
    return 7;
    case D218:
    return 7;
    default:
    return 0;
  }
  return 0; //never called
}

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
SEXP READ_OUTCROSS(SEXP file)
{
  std::vector<std::string> segr_type;
  std::vector<unsigned short> segr_type_num;
  std::vector<unsigned short> geno;
  std::vector<std::string> markers;
  // Get nymber of markers, individuals and phenotypes
  unsigned long n_mar = 0;
  unsigned long n_ind = 0;
  unsigned long n_phe = 0;

  std::string inf = as<std::string>(file);

  std::ifstream ifs(inf.c_str(), std::ios::in);

  std::string x;

  ifs >> x; n_ind = std::stoul(x);
  ifs >> x; n_mar = std::stoul(x);
  ifs >> x; n_phe = std::stoul(x);
  Rcpp::Rcout << "Reading data... ";
  progressBar prog_bar(40, n_mar);

  bool first = true;

  // Start to read markers

  unsigned long lcount = 1; // Lines read
  unsigned long mcount = 0; // Already counted markers

  std::string line;
  while(std::getline(ifs, line))
  {
    lcount++;
    if(first)
    {
      first = false; continue;
    }
    if(line.size() == 0) continue;
    if(line.at(0) == '#') continue;
    std::stringstream ss(line);
    char segregation_chr;
    if(line.at(0) == '*')
    {
      mcount++;
      prog_bar.update(mcount);
      if(mcount > n_mar) break; // We ignore phenotype information
      ss.ignore(1);

      std::string marker;
      ss >> marker;

      std::string segregation_str;
      ss >> segregation_str;

      segregation_chr = segregation.at(segregation_str);
      unsigned short segregation_us = codify_segregation(segregation_chr);

      std::string gt;
      ss >> gt;

      std::stringstream ss2(gt);

      std::string genotype;
      while(std::getline(ss2, genotype, ','))
      {
        geno.push_back(codify_geno.at(std::pair<char, std::string>(segregation_chr,genotype)));
      }

      segr_type.push_back(segregation_str);
      segr_type_num.push_back(segregation_us);
      markers.push_back(marker);
    }
    else
    {
      std::string gt;
      ss >> gt;

      std::stringstream ss2(gt);

      std::string genotype;
      while(std::getline(ss2, genotype, ','))
      {
        geno.push_back(codify_geno.at(std::pair<char, std::string>(segregation_chr,genotype)));
      }
    }
  }

    return Rcpp::List::create(Rcpp::Named("geno") = wrap(geno),
                            Rcpp::Named("marker") = wrap(markers),
                            Rcpp::Named("n.ind") = wrap(n_ind),
                            Rcpp::Named("n.mar") = wrap(n_mar),
                            Rcpp::Named("n.phe") = wrap(n_phe),
                            Rcpp::Named("segr.type") = wrap(segr_type),
                            Rcpp::Named("segr.type.num") = wrap(segr_type_num));

}

