/******************************************************************************
*******************************************************************************
***   Interpolation... classes                                              ***
***                                                                         ***
***  ... can interpolate nodal or Gauss-pt values to points                 ***
***                                                                         ***
***  IMPORTANT NOTE:                                                        ***
***  These classes use a reference of the ptWeights arrays (containing      ***
***  node/element numbers and weight factors) that have been supplied to    ***
***  them initially.                                                        ***
***  DON'T delete / free those arrays as long as you use the interpolation  ***
***  object!                                                                ***
***                                                                         ***
*******************************************************************************
******************************************************************************/

#ifndef INTERPOLATION_DEFINED
#define INTERPOLATION_DEFINED

#include <stdio.h>
#include <vector>

/* Python integration */
#include <Python.h>
#include <structmember.h>   // needed for e.g. T_INT

#include "fieldValueMap.h"

class InterpolationBase
{
public:
  // m_gaussptsPerPoint needs to be declared in InterpolationBase because it
  // will be accessed through the generic pointer "interpolation"
  // (of type InterpolationBase*).
  // It holds the nb of Gauss points according to what's been supplied by
  // the constructor. The constructor is typically called by
  // Communicator_storeRemapWeights.
  int m_gaussptsPerPoint;

public:
  InterpolationBase(int gaussptsPerPoint=0);
  virtual ~InterpolationBase() {};

  // void error(const char* fmt, ...);
  virtual int interpolateToPyList(FieldValueMapBase* field,
                                  PyObject* result) = 0;
  virtual int interpolateToPyListPieceWiseConst
    (FieldValueMapBase* field, PyObject* result) = 0;
  // virtual int getNbOfGaussptsPerPoint();
};

class InterpolationFromNodal : public InterpolationBase
{
protected:
  int m_nodesPerPoint;
 public:    // for debugging - output, otherwise protected would be fine...
  int m_nbPoints;
  char* m_ptWeights;
  int m_stridesPts;
  int m_stridesNodes;
public:
  InterpolationFromNodal(int nodesPerPoint, int nbPoints,
			 char* ptWeights);
  virtual ~InterpolationFromNodal() {};
  virtual int interpolateToPyList(FieldValueMapBase* field, PyObject* result);
  virtual int interpolateToPyListPieceWiseConst
    (FieldValueMapBase* field, PyObject* result);
};

class InterpolationFromGaussPts : public InterpolationBase
{
protected:
  int m_nbPoints;
  char* m_ptWeights;
  int m_stridesPts;

public:
  InterpolationFromGaussPts(int gaussptsPerPoint, int nbPoints,
			    char* ptWeights);
  virtual ~InterpolationFromGaussPts() {};
  virtual int getNbOfGaussptsPerPoint();

  /* get the element number that has been stored in the ptWeights array
   * for a specific pointIdx */
  virtual int getElem(int pointIdx);
  virtual int interpolateToPyList(FieldValueMapBase* field, PyObject* result);
  virtual int interpolateToPyListPieceWiseConst
    (FieldValueMapBase* field, PyObject* result);
};

class InterpolationFromElem : public InterpolationBase
{
protected:
  int m_nbPoints;
  char* m_ptWeights;
  int m_stridesPts;
public:
  InterpolationFromElem(int gaussptsPerPoint, int nbPoints,
			char* ptWeights);
  virtual ~InterpolationFromElem() {};

  /* get the element number that has been stored in the ptWeights array
   * for a specific pointIdx */
  virtual int getElem(int pointIdx);
  virtual int interpolateToPyList(FieldValueMapBase* field, PyObject* result);
  virtual int interpolateToPyListPieceWiseConst
    (FieldValueMapBase* field, PyObject* result);
};


#endif
