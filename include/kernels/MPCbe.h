//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Action.h"


class MPCbe : public Action
{
public:
  static InputParameters validParams();

  MPCbe(const InputParameters & params);

  virtual void act() override;


    // Holds the primary node ids
  std::vector<unsigned int> primary_node_ids;
  // Holds the list of secondary node ids
  std::vector<unsigned int> secondary_node_ids;

  // std::vector<Real> weights;



};
