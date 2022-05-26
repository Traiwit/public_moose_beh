//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowPermeabilityBase.h"
#include "ElementPropertyReadFile.h"
#include "RankTwoTensor.h"
#include "Material.h"

/**
 * Material designed to provide a constant permeability tensor
 */
class PorousFlowPermeabilityConstCSV : public PorousFlowPermeabilityBase
{
public:
  static InputParameters validParams();

  PorousFlowPermeabilityConstCSV(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  /// Constant value of permeability tensor
  const RealTensorValue _input_permeability;

  const ElementPropertyReadFile * const _read_prop_user_object;

  MaterialProperty<RealVectorValue> & _mat_prop;

  const MaterialProperty<Real> & _damage;
  const MaterialProperty<Real> & _excav_perm;

  MaterialProperty<RealVectorValue> & _aw;

  MaterialProperty<RealVectorValue> & _kmax;

  MaterialProperty<RealVectorValue> & _RMD_max;

  MaterialProperty<RealVectorValue> & _perm;






};
