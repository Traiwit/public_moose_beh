//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPiecewiseLinearSinkBEH.h"

registerMooseObject("moose_behApp", PorousFlowPiecewiseLinearSinkBEH);

InputParameters
PorousFlowPiecewiseLinearSinkBEH::validParams()
{
  InputParameters params = PorousFlowSinkPTDefinerBEH::validParams();
  params.addRequiredParam<std::vector<Real>>(
      "pt_vals",
      "Tuple of pressure values (for the fluid_phase specified).  Must be monotonically "
      "increasing.  For heat fluxes that don't involve fluids, these are temperature "
      "values");
  params.addRequiredParam<std::vector<Real>>(
      "multipliers", "Tuple of multiplying values.  The flux values are multiplied by these.");
  params.addClassDescription("Applies a flux sink to a boundary. The base flux defined by "
                             "PorousFlowSink is multiplied by a piecewise linear function.");
  return params;
}

PorousFlowPiecewiseLinearSinkBEH::PorousFlowPiecewiseLinearSinkBEH(const InputParameters & parameters)
  : PorousFlowSinkPTDefinerBEH(parameters),
    _sink_func(getParam<std::vector<Real>>("pt_vals"), getParam<std::vector<Real>>("multipliers"))
{
}

Real
PorousFlowPiecewiseLinearSinkBEH::multiplier() const
{
  return PorousFlowSinkBEH::multiplier() * _sink_func.sample(ptVar());
}

Real
PorousFlowPiecewiseLinearSinkBEH::dmultiplier_dvar(unsigned int pvar) const
{
  return PorousFlowSinkBEH::dmultiplier_dvar(pvar) * _sink_func.sample(ptVar()) +
         PorousFlowSinkBEH::multiplier() * _sink_func.sampleDerivative(ptVar()) * dptVar(pvar);
}
