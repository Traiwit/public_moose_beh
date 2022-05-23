//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPermeabilityConstCSV.h"
#include <iostream>
#include <algorithm>

registerMooseObject("MOOSE_BEHApp", PorousFlowPermeabilityConstCSV);

InputParameters
PorousFlowPermeabilityConstCSV::validParams()
{
  InputParameters params = PorousFlowPermeabilityBase::validParams();
  params.addParam<RealTensorValue>(
      "permeability",
      "The permeability tensor (usually in m^2), which is assumed constant for this material");
  params.addParam<UserObjectName>("read_prop_user_object",
                                      "The ElementReadPropertyFile "
                                      "GeneralUserObject to read element "
                                      "specific property values from file");
  params.addClassDescription(
      "This Material calculates the permeability tensor assuming it is constant");
  params.addParam<MaterialPropertyName>("damage",0,"damage");
  params.addParam<MaterialPropertyName>("excav_perm", 0 ,"The coupled variable which provides the force");
  return params;
}

PorousFlowPermeabilityConstCSV::PorousFlowPermeabilityConstCSV(const InputParameters & parameters)
  : PorousFlowPermeabilityBase(parameters), _damage(getMaterialProperty<Real>("damage")),
  _read_prop_user_object(isParamValid("read_prop_user_object")
                             ? &getUserObject<ElementPropertyReadFile>("read_prop_user_object")
                             : nullptr),
  _mat_prop(declareProperty<RealVectorValue>("Euler_angles")),
  _aw(declareProperty<RealVectorValue>("aw")),
  _kmax(declareProperty<RealVectorValue>("kmax")),
  _RMD_max(declareProperty<RealVectorValue>("RMD_MAX")),
  _perm(declareProperty<RealVectorValue>("perm_tensor")),
  _excav_perm(getMaterialProperty<Real>("excav_perm"))
{

}

void
PorousFlowPermeabilityConstCSV::computeQpProperties()
{

/// RMD_max is given
//  read values from .CSV file
//   _mat_prop[_qp](0) = _read_prop_user_object->getData(_current_elem, 0);
//   _mat_prop[_qp](1) = _read_prop_user_object->getData(_current_elem, 1);
//   _mat_prop[_qp](2) = _read_prop_user_object->getData(_current_elem, 2);
//   _mat_prop[_qp](3) = _read_prop_user_object->getData(_current_elem, 3);
//   _mat_prop[_qp](4) = _read_prop_user_object->getData(_current_elem, 4);
//   _mat_prop[_qp](5) = _read_prop_user_object->getData(_current_elem, 5);
//   _RMD_max[_qp](0) = _read_prop_user_object->getData(_current_elem, 6);
//   _kmax[_qp](0)    = _read_prop_user_object->getData(_current_elem, 7);
// // calculate _aw from input data
//   _aw[_qp](0) = log(_kmax[_qp](0)/_mat_prop[_qp](0))/_RMD_max[_qp](0);
//   _aw[_qp](1) = log(_kmax[_qp](0)/_mat_prop[_qp](1))/_RMD_max[_qp](0);
//   _aw[_qp](2) = log(_kmax[_qp](0)/_mat_prop[_qp](2))/_RMD_max[_qp](0);
//
//   if (_mat_prop[_qp](3) > 0){
//   _aw[_qp](3) = log(_kmax[_qp](0)/_mat_prop[_qp](3))/_RMD_max[_qp](0);
//   }
//   else {
//     _aw[_qp](3) = 0;
//   }
//   if (_mat_prop[_qp](4) > 0){
//   _aw[_qp](4) = log(_kmax[_qp](0)/_mat_prop[_qp](4))/_RMD_max[_qp](0);
//   }
//   else {
//     _aw[_qp](3) = 0;
//   }
//   if (_mat_prop[_qp](5) > 0){
//   _aw[_qp](5) = log(_kmax[_qp](0)/_mat_prop[_qp](5))/_RMD_max[_qp](0);
//   }
//   else {
//     _aw[_qp](5) = 0;
//   }
//
// // calculate perm tensor
//     _perm[_qp](0) = std::min(exp(_damage[_qp]*_aw[_qp](0))*_mat_prop[_qp](0),_kmax[_qp](0));
//     _perm[_qp](1) = std::min(exp(_damage[_qp]*_aw[_qp](1))*_mat_prop[_qp](1),_kmax[_qp](0));
//     _perm[_qp](2) = std::min(exp(_damage[_qp]*_aw[_qp](2))*_mat_prop[_qp](2),_kmax[_qp](0));
//     _perm[_qp](3) = std::min(exp(_damage[_qp]*_aw[_qp](3))*_mat_prop[_qp](3),_kmax[_qp](0));
//     _perm[_qp](4) = std::min(exp(_damage[_qp]*_aw[_qp](4))*_mat_prop[_qp](4),_kmax[_qp](0));
//     _perm[_qp](5) = std::min(exp(_damage[_qp]*_aw[_qp](5))*_mat_prop[_qp](5),_kmax[_qp](0));



////// Aw is given
//  read values from .CSV file
  _mat_prop[_qp](0) = _read_prop_user_object->getData(_current_elem, 0);
  _mat_prop[_qp](1) = _read_prop_user_object->getData(_current_elem, 1);
  _mat_prop[_qp](2) = _read_prop_user_object->getData(_current_elem, 2);

  //_mat_prop[_qp](3) = _read_prop_user_object->getData(_current_elem, 3);
  //_mat_prop[_qp](4) = _read_prop_user_object->getData(_current_elem, 4);
  //_mat_prop[_qp](5) = _read_prop_user_object->getData(_current_elem, 5);
  //
  //
  _aw[_qp](0) = _read_prop_user_object->getData(_current_elem, 6);
  _kmax[_qp](0)    = _read_prop_user_object->getData(_current_elem, 7);
  // // calculate perm tensor
  if (_damage[_qp] > -50){
    _perm[_qp](0) = std::min(exp(_damage[_qp]*_aw[_qp](0))*_mat_prop[_qp](0)*1e-7,_kmax[_qp](0)*1e-7);
    _perm[_qp](1) = std::min(exp(_damage[_qp]*_aw[_qp](0))*_mat_prop[_qp](1)*1e-7,_kmax[_qp](0)*1e-7);
    _perm[_qp](2) = std::min(exp(_damage[_qp]*_aw[_qp](0))*_mat_prop[_qp](2)*1e-7,_kmax[_qp](0)*1e-7);
  }

    else if (_damage[_qp] > -110) {
       _perm[_qp](0) = pow(10,_damage[_qp]+100)*1e-7;
       _perm[_qp](1) = pow(10,_damage[_qp]+100)*1e-7;
       _perm[_qp](2) = pow(10,_damage[_qp]+100)*1e-7;
     }


    else {
       _perm[_qp](0) = pow(10,_damage[_qp]+101)*1e-7;
       _perm[_qp](1) = pow(10,_damage[_qp]+101)*1e-7;
       _perm[_qp](2) = pow(10,_damage[_qp]+101)*1e-7;
     }




RealTensorValue permeability(  _perm[_qp](0),
                               0,
                               0,
                               0,
                               _perm[_qp](1),
                               0,
                               0,
                               0,
                               _perm[_qp](2));

_permeability_qp[_qp]  = permeability;

  _dpermeability_qp_dvar[_qp].assign(_num_var, RealTensorValue());
  _dpermeability_qp_dgradvar[_qp].resize(LIBMESH_DIM);
  for (unsigned i = 0; i < LIBMESH_DIM; ++i)
    _dpermeability_qp_dgradvar[_qp][i].assign(_num_var, RealTensorValue());
}
