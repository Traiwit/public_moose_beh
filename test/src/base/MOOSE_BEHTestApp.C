//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "MOOSE_BEHTestApp.h"
#include "MOOSE_BEHApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
MOOSE_BEHTestApp::validParams()
{
  InputParameters params = MOOSE_BEHApp::validParams();
  return params;
}

MOOSE_BEHTestApp::MOOSE_BEHTestApp(InputParameters parameters) : MooseApp(parameters)
{
  MOOSE_BEHTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

MOOSE_BEHTestApp::~MOOSE_BEHTestApp() {}

void
MOOSE_BEHTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  MOOSE_BEHApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"MOOSE_BEHTestApp"});
    Registry::registerActionsTo(af, {"MOOSE_BEHTestApp"});
  }
}

void
MOOSE_BEHTestApp::registerApps()
{
  registerApp(MOOSE_BEHApp);
  registerApp(MOOSE_BEHTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
MOOSE_BEHTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  MOOSE_BEHTestApp::registerAll(f, af, s);
}
extern "C" void
MOOSE_BEHTestApp__registerApps()
{
  MOOSE_BEHTestApp::registerApps();
}
