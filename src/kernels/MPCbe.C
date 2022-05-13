//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MPCbe.h"
// #include "AddConstraintAction.h"
#include "LinearNodalConstraint.h"
#include "MooseMesh.h"
#include "FEProblem.h"

#include <sstream>
#include <stdexcept>
#include "libmesh/libmesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fe.h"

registerMooseAction("moose_behApp", MPCbe, "add_constraint");

InputParameters
MPCbe::validParams()
{
  InputParameters params = Action::validParams();
  params.addParam<std::vector<unsigned int>>("primary_node_ids", "The primary node IDs.");
  params.addParam<std::vector<unsigned int>>("secondary_node_ids",
                                           "The list of secondary node ids");

  return params;
}
MPCbe::MPCbe(const InputParameters & params)
  : Action(params),
primary_node_ids(getParam<std::vector<unsigned int>>("primary_node_ids")),
secondary_node_ids(getParam<std::vector<unsigned int>>("secondary_node_ids"))

{
}
void
MPCbe::act()

{
for (unsigned cur_num = 0; cur_num<primary_node_ids.size() ; cur_num++)
{

  std::vector<Real> weights_in (1);
  std::vector<unsigned int> primary_node_ids_in (1);
  std::vector<unsigned int> secondary_node_ids_in (1);

  weights_in.at(0) = 1;
  primary_node_ids_in.at(0) = primary_node_ids[cur_num];
  secondary_node_ids_in.at(0) = secondary_node_ids[cur_num];

  InputParameters params = _factory.getValidParams("LinearNodalConstraint");
  params.set<NonlinearVariableName>("variable") = "porepressure";
  params.set<std::vector<Real>>("weights")= weights_in;
  params.set<std::vector<unsigned int>>("primary") = primary_node_ids_in;
  params.set<std::vector<unsigned int>>("secondary_node_ids") = secondary_node_ids_in;
  params.set<Real>("penalty") = 1e10;
  _problem->addConstraint("LinearNodalConstraint", "MPCbe" + Moose::stringify(cur_num), params);

}
}
