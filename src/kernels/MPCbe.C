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

// STL includes
#include <fstream>

// MOOSE includes
#include "CSVReader.h"
#include "MooseUtils.h"
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <math.h>
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
#include "GeneralVectorPostprocessor.h"

registerMooseAction("MOOSE_BEHApp", MPCbe, "add_constraint");

InputParameters
MPCbe::validParams()
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

MPCbe::MPCbe(const InputParameters & params)
  : Action(params)
{
}
void
MPCbe::act()

{
  MooseUtils::DelimitedFileReader csv_reader(getParam<FileName>("csv_file"), &_communicator);
  if (isParamValid("header"))
  csv_reader.setHeaderFlag(getParam<bool>("header")
                               ? MooseUtils::DelimitedFileReader::HeaderFlag::ON
                               : MooseUtils::DelimitedFileReader::HeaderFlag::OFF);

  csv_reader.read();
    const std::vector<std::string> & names = csv_reader.getNames();
    const std::vector<std::vector<double>> & data = csv_reader.getData();


    int row = data.size();


    _console << "number of column   " <<  row << std::endl;

for (unsigned cur_num = 0; cur_num<10000000 ; cur_num++)
{

  std::vector<Real> weights_in (1);
  std::vector<unsigned int> primary_node_ids_in (1);
  std::vector<unsigned int> secondary_node_ids_in_1 (1);
  std::vector
  <unsigned int> secondary_node_ids_in_2 (1);

  weights_in.at(0) = 1;
  primary_node_ids_in.at(0) = data[0][cur_num];
  secondary_node_ids_in_1.at(0) =  data[1][cur_num];
  std::vector<unsigned int> secondary_node_ids_in_final(secondary_node_ids_in_1);

  for (unsigned n = 2; n<row ; n++)
  {

  secondary_node_ids_in_2.at(0) =  data[n][cur_num];

  secondary_node_ids_in_final.insert(secondary_node_ids_in_final.end(), secondary_node_ids_in_2.begin(), secondary_node_ids_in_2.end());
  }

  sort( secondary_node_ids_in_final.begin(), secondary_node_ids_in_final.end() );
  secondary_node_ids_in_final.erase( unique( secondary_node_ids_in_final.begin(), secondary_node_ids_in_final.end() ), secondary_node_ids_in_final.end() );

  if (primary_node_ids_in.at(0) == 0){
      break;
  }

  InputParameters params = _factory.getValidParams("LinearNodalConstraint");
  params.set<NonlinearVariableName>("variable") = "porepressure";
  params.set<std::vector<Real>>("weights")= weights_in;
  params.set<std::vector<unsigned int>>("primary") = primary_node_ids_in;
  params.set<std::vector<unsigned int>>("secondary_node_ids") = secondary_node_ids_in_final;
  params.set<Real>("penalty") = 1e10;
  _problem->addConstraint("LinearNodalConstraint", "MPCbe" + Moose::stringify(cur_num), params);

}


_console << "MPC done!" << std::endl;
}
