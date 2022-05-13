#include "MOOSE_BEHApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
MOOSE_BEHApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

MOOSE_BEHApp::MOOSE_BEHApp(InputParameters parameters) : MooseApp(parameters)
{
  MOOSE_BEHApp::registerAll(_factory, _action_factory, _syntax);
}

MOOSE_BEHApp::~MOOSE_BEHApp() {}

void
MOOSE_BEHApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"MOOSE_BEHApp"});
  Registry::registerActionsTo(af, {"MOOSE_BEHApp"});

  /* register custom execute flags, action syntax, etc. here */

  auto & syntax = s;  // for resiterSyntax macros
  auto & factory = f; // for resiterSyntax macros

  // registerExecFlag(EXEC_JUST_GO);
  registerSyntax("MPCbe", "MPCbe");
  registerSyntax("MPCbe2", "MPCbe2");

  /* register custom execute flags, action syntax, etc. here */
}

void
MOOSE_BEHApp::registerApps()
{
  registerApp(MOOSE_BEHApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
MOOSE_BEHApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  MOOSE_BEHApp::registerAll(f, af, s);
}
extern "C" void
MOOSE_BEHApp__registerApps()
{
  MOOSE_BEHApp::registerApps();
}
