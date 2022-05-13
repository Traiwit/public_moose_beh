//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GenericConstantArrayBEH.h"
#include "libmesh/quadrature.h"

registerMooseObject("moose_behApp", GenericConstantArrayBEH);


InputParameters
GenericConstantArrayBEH::validParams()
{

  InputParameters params = Material::validParams();
  params.addParam<std::string>("prop_name",
                                       "The name of the property this material will have");
  // params.addParam<RealEigenVector>("prop_value",
  //                                          "The values associated with the named property");
  // params.declareControllable("prop_value");
  params.addClassDescription(
      "A material evaluating one material property in type of RealEigenVector");
  // params.set<MooseEnum>("constant_on") = "SUBDOMAIN";
  params.addParam<UserObjectName>("read_prop_user_object",
                                      "The ElementReadPropertyFile "
                                      "GeneralUserObject to read element "
                                      "specific property values from file");
  return params;
}
GenericConstantArrayBEH::GenericConstantArrayBEH(const InputParameters & parameters)
  : Material(parameters),
    _prop_name(getParam<std::string>("prop_name")),
    // _prop_value(getParam<RealEigenVector>("prop_value")),
    _property(declareProperty<Real>(_prop_name)),
    _mat_prop(declareProperty<Real>("prop_value")),
    _read_prop_user_object(isParamValid("read_prop_user_object")
                               ? &getUserObject<ElementPropertyReadFile>("read_prop_user_object")
                               : nullptr)
{
}

void
GenericConstantArrayBEH::initQpStatefulProperties()
{
  computeQpProperties();
}

void
GenericConstantArrayBEH::computeQpProperties()
{
   _mat_prop[_qp] = _read_prop_user_object->getData(_current_elem, 0);
  _property[_qp] = _mat_prop[_qp];
}
