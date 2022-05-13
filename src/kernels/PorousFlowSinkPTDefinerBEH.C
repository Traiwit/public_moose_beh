//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowSinkPTDefinerBEH.h"

InputParameters
PorousFlowSinkPTDefinerBEH::validParams()
{
  InputParameters params = PorousFlowSinkBEH::validParams();
  params.addCoupledVar("PT_shift",
                       0.0,
                       "Whenever the sink is an explicit function of porepressure "
                       "(such as a PiecewiseLinear function) the argument of the "
                       "function is set to P - PT_shift instead of simply P.  "
                       "Similarly for temperature.  PT_shift does not enter into "
                       "any use_* calculations.");
  return params;
}

PorousFlowSinkPTDefinerBEH::PorousFlowSinkPTDefinerBEH(const InputParameters & parameters)
  : PorousFlowSinkBEH(parameters),
    _pp(_involves_fluid ? &getMaterialProperty<std::vector<Real>>("PorousFlow_porepressure_nodal")
                        : nullptr),
    _dpp_dvar(_involves_fluid ? &getMaterialProperty<std::vector<std::vector<Real>>>(
                                    "dPorousFlow_porepressure_nodal_dvar")
                              : nullptr),
    _temp(!_involves_fluid ? &getMaterialProperty<Real>("PorousFlow_temperature_nodal") : nullptr),
    _dtemp_dvar(!_involves_fluid
                    ? &getMaterialProperty<std::vector<Real>>("dPorousFlow_temperature_nodal_dvar")
                    : nullptr),
    _pt_shift(coupledDofValues("PT_shift"))
{
  if (_involves_fluid && (_pp == nullptr || _dpp_dvar == nullptr))
    mooseError("PorousFlowSinkBEH: There is no porepressure Material");
  if (!_involves_fluid && (_temp == nullptr || _dtemp_dvar == nullptr))
    mooseError("PorousFlowSinkBEH: There is no temperature Material");
}

Real
PorousFlowSinkPTDefinerBEH::ptVar() const
{
  if (_involves_fluid)
    return (*_pp)[_i][_ph] - _pt_shift[_i];
  return (*_temp)[_i] - _pt_shift[_i];
}

Real
PorousFlowSinkPTDefinerBEH::dptVar(unsigned pvar) const
{
  if (_involves_fluid)
    return (*_dpp_dvar)[_i][_ph][pvar];
  return (*_dtemp_dvar)[_i][pvar];
}
