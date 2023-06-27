/*
*                    y=1_____________
*                    / |           /|
*                  /   | [3]     /  |
*                /_____|_______/    |
*                |     |    [5|     |
*                | [0] |      | [1] |
*                |     |Y     |     |
*    x=y=z=0     |     |__x___|_____| x=1
*                |   z/       |    /
*                |  /    [2]  |  /
*                |/___________|/
*               z=1
!               !  Neumann     du/dn = 0
!               !  Dirichlet   u = 1
!               !  no_slip     u = 0
*/

#pragma once
#include <tuple>
#include <vector>

using NumDirType = std::pair<std::vector<std::vector<double> >,
                             std::vector<std::vector<double> > >;

class BoundaryConditionSelector {
 public:
 
   std::vector<std::vector<double> > num;
   std::vector<std::vector<double> > dir;


  BoundaryConditionSelector() {
    dir.resize(6, std::vector<double>(3, 0.0));
    num.resize(6, std::vector<double>(3, 0.0));
   }

   auto CavityFlow() {
#include "cavity.hpp";
   }

  auto ChannelFlow(){
#include "channel_flow.hpp"
  }

  auto FlowPassing(){
#include "flow_passing.hpp"
  }

  NumDirType
  operator()() {
    return std::make_pair(num, dir);
  }
};