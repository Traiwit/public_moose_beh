//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowPorosityBase.h"
#include "ElementPropertyReadFile.h"
#include "Material.h"

/**
 * Material to provide a constant value of porosity. This can be specified
 * by either a constant value in the input file, or taken from an aux variable.
 * Note: this material assumes that the porosity remains constant throughout a
 * simulation, so the coupled aux variable porosity must also remain constant.
 */
class PorousFlowPorosityConstCSV : public PorousFlowPorosityBase
{
public:
  static InputParameters validParams();

  PorousFlowPorosityConstCSV(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Constant porosity
  // const VariableValue & _input_porosity;

  const ElementPropertyReadFile * const _read_prop_user_object;

  MaterialProperty<Real> & _mat_prop;
  
  const MaterialProperty<Real> & _damage;

  const MaterialProperty<Real> & _excav_poro;

};
