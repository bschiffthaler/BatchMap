#include <Rcpp.h>
#include "flip_phase.h"
// [[Rcpp::export]]
Rcpp::IntegerVector flip_phases(Rcpp::IntegerVector seq_type,
                                Rcpp::IntegerVector phases) {
  Rcpp::IntegerVector new_phases(phases.size());
  for(unsigned int i = 0; i < (seq_type.size() - 1); i++)
  {
    int s1 = seq_type( i );
    int s2 = seq_type( i + 1 );
    if(s1 < s2){
      int s_tmp = s1;
      s1 = s2; s2 = s_tmp;
    }
    int p_new;

    switch(s1){

    case 4: //B3
      switch(s2){
      case 4:
        p_new = flip_B3_B3(phases(i)); break;
      case 6:
        p_new = flip_B3_D1(phases(i)); break;
      case 7:
        p_new = flip_B3_D2(phases(i)); break;
      } break; // END switch(s2)

    case 6: //D1
      switch(s2){
      case 6:
        p_new = flip_D1_D1(phases(i)); break;
      } break; // END switch(s2)

    case 7: //D2
      switch(s2){
      case 7:
        p_new = flip_D2_D2(phases(i)); break;
      } break; // END switch(s2)

    } // END switch(s1)

    new_phases(i) = p_new;
  } // END  for(unsigned int i = 0; i < (seq_type.size() - 1); i++)

  return new_phases;
}
