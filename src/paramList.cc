/* C++11 */

/**
 * @file   paramList.cc
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * This file is part of the "DosExtraction" project.
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 */

#include "paramList.h"

using namespace constants;

ParamList::ParamList(const RowVectorXr & list)
{
  assert( list.size() == PARAMS_NO );
  
  assert( list( 0) >  0   );
  assert( list( 1) >  0.0 );
  assert( list( 2) >  0.0 );
  assert( list( 7) >= 0.0 );
  assert( list( 9) >= 0.0 );
  assert( list(12) >= 0.0 );
  assert( list(15) >= 0.0 );
  assert( list(22) >  0   );
  assert( list(23) >  0   );
  
  simulationNo_ = list( 0)          ;
  t_semic_      = list( 1)          ;
  t_ins_        = list( 2)          ;
  eps_semic_    = list( 3)          ;
  eps_ins_      = list( 4)          ;
  Wf_           = list( 5) * Q      ;
  Ea_           = list( 6) * Q      ;
  N0_           = list( 7)          ;
  sigma_        = list( 8) * K_B * T;
  N0_2_         = list( 9)          ;
  sigma_2_      = list(10) * K_B * T;
  shift_2_      = list(11) * ( - Q) ;
  N0_3_         = list(12)          ;
  sigma_3_      = list(13) * K_B * T;
  shift_3_      = list(14) * ( - Q) ;
  N0_4_         = list(15)          ;
  sigma_4_      = list(16) * K_B * T;
  shift_4_      = list(17) * ( - Q) ;
  N0_exp_       = list(18)          ;
  lambda_exp_   = list(19) * K_B * T;
  A_semic_      = list(20)          ;
  C_sb_         = list(21)          ;
  nNodes_       = list(22)          ;
  nSteps_       = list(23)          ;
  V_min_        = list(24)          ;
  V_max_        = list(25)          ;
}
