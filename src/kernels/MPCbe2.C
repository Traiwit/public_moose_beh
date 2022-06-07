//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MPCbe2.h"
// #include "AddConstraintAction.h"
#include "LinearNodalConstraint.h"
#include "MooseMesh.h"
#include "FEProblem.h"


#include <sstream>
#include <stdexcept>
#include <iostream>
using namespace std;
#include "libmesh/libmesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fe.h"

#include "gtest/gtest.h"

// MOOSE includes
#include "DelimitedFileReader.h"
#include "MooseException.h"
#include "libmesh/point.h"
#include "DelimitedFileReader.h"
#include "MooseUtils.h"
#include "MooseError.h"
#include "pcrecpp.h"

registerMooseAction("MOOSE_BEHApp", MPCbe2, "add_constraint");

InputParameters
MPCbe2::validParams()
{
  InputParameters params = Action::validParams();
  params.addRequiredParam<FileName>("csv_file",
                                    "The name of the CSV file to read. Currently, with "
                                    "the exception of the header row, only numeric "
                                    "values are supported.");
params.addParam<bool>("header",
                              "When true it is assumed that the first row contains column headers, these "
                              "headers are used as the VectorPostprocessor vector names. If false the "
                              "file is assumed to contain only numbers and the vectors are named "
                             "automatically based on the column number (e.g., 'column_0000', "
                              "'column_0001'). If not supplied the reader attempts to auto detect the "
                              "headers.");

  return params;
}

MPCbe2::MPCbe2(const InputParameters & params)
  : Action(params)
{
}

void
MPCbe2::act()
{
  MooseUtils::DelimitedFileReader reader(getParam<FileName>("csv_file"), &_communicator);
  if (isParamValid("header"))
  reader.setHeaderFlag(getParam<bool>("header")
                               ? MooseUtils::DelimitedFileReader::HeaderFlag::ON
                               : MooseUtils::DelimitedFileReader::HeaderFlag::OFF);

  reader.read();
  const std::vector<std::vector<double>> & data = reader.getData();

  for (unsigned cur_num = 0; cur_num<10000000 ; cur_num++)
  {

    std::vector<Real> weights_in (1);
    std::vector<unsigned int> primary_node_ids_in (1);
    std::vector<unsigned int> secondary_node_ids_in(1);

    weights_in.at(0) = 1;
    primary_node_ids_in.at(0) = data[1][cur_num];
    secondary_node_ids_in.at(0) =  data[0][cur_num];

          if (primary_node_ids_in.at(0) == 0){
              break;
          }

    InputParameters params = _factory.getValidParams("LinearNodalConstraint");
    params.set<NonlinearVariableName>("variable") = "porepressure";
    params.set<std::vector<Real>>("weights")= weights_in;
    params.set<std::vector<unsigned int>>("primary") = primary_node_ids_in;
    params.set<std::vector<unsigned int>>("secondary_node_ids") = secondary_node_ids_in;
    params.set<Real>("penalty") = 1e10;
    _problem->addConstraint("LinearNodalConstraint", "MPCbe2" + Moose::stringify(cur_num), params);

  }
}
