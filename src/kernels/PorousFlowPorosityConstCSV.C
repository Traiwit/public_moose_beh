//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPorosityConstCSV.h"

registerMooseObject("MOOSE_BEHApp", PorousFlowPorosityConstCSV);

InputParameters
PorousFlowPorosityConstCSV::validParams()
{
  InputParameters params = PorousFlowPorosityBase::validParams();
  params.addParam<UserObjectName>("read_prop_user_object",
                                  "The ElementReadPropertyFile "
                                  "GeneralUserObject to read element "
                                  "specific property values from file");
  params.addClassDescription("This Material calculates the porosity assuming it is constant");
  params.addParam<MaterialPropertyName>(
      "excav_poro", 0, "The coupled variable which provides the force");
  params.addParam<MaterialPropertyName>("damage", 0, "damage");
  return params;
}

PorousFlowPorosityConstCSV::PorousFlowPorosityConstCSV(const InputParameters & parameters)
  : PorousFlowPorosityBase(parameters),
    _read_prop_user_object(isParamValid("read_prop_user_object")
                               ? &getUserObject<ElementPropertyReadFile>("read_prop_user_object")
                               : nullptr),
    _mat_prop(declareProperty<Real>("container")),
    _damage(getMaterialProperty<Real>("damage")),
    _excav_poro(getMaterialProperty<Real>("excav_poro"))
{
}

void
PorousFlowPorosityConstCSV::initQpStatefulProperties()
{

  if (_current_elem->id() > _mesh.nElem())
    _console << _current_elem->true_centroid();

  if (_damage[_qp] < -100)
  {
    _porosity[_qp] = 1;
  }
  else
  {
    _porosity[_qp] = _read_prop_user_object->getData(_current_elem, 0) + _excav_poro[_qp];
  }
}

void
PorousFlowPorosityConstCSV::computeQpProperties()
{
  initQpStatefulProperties();

  // The derivatives are zero for all time
  _dporosity_dvar[_qp].assign(_num_var, 0.0);
  _dporosity_dgradvar[_qp].assign(_num_var, RealGradient());
}
