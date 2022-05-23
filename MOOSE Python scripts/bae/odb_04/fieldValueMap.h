/******************************************************************************
*******************************************************************************
***   FieldValueMap... classes                                              ***
***                                                                         ***
***   dictionaries to hold field values in dictionaries (maps)              ***

fieldValueMap-classes  ... maps node/elem-label -> scalar/vector/vec-of-vectors
 - FieldValueMapVector:         nodal, whole elem; vectors
 - FieldValueMapVectorToVec:    integr. pt       ; vectors
 - FieldValueMapVecComp:        nodal, whole elem; scalars
 - FieldValueMapVecCompToVec:   integr. pt       ; scalars

*******************************************************************************
******************************************************************************/

#ifndef FIELDVALUEMAP_DEFINED
#define FIELDVALUEMAP_DEFINED

#include <stdio.h>

#include <map>
#include <vector>

/* Python integration */
#include <Python.h>
#include <structmember.h>   // needed for e.g. T_INT

/******************************************************************************
 * class FieldValueMapBase
 *
 * Base class for all following classes to define the common interface.
 * Not to be instantiated directly.
 */

class FieldValueMapBase
{
public:
  int m_nbComponents;

  FieldValueMapBase(int nbComponents=-1);
  virtual ~FieldValueMapBase();  // essential to make the destructor virtual
  // void error(const char* fmt, ...);
  virtual std::size_t size() = 0;
  virtual int store(int key, float *val);
  virtual void storeNotDefined(int key, float * val);
  virtual int getNbComp();
  virtual float getScalar(int key);
  virtual const float* getVector(int key);
  virtual const std::vector<float>& getVecOfScalars(int key);
  virtual const std::vector<const float*>& getVecOfVectors(int key);

  // for values at nodes
  virtual void initResult();
  virtual void addWeighted(int key, float weight);
  virtual PyObject* getPyResult();

  // for values at integration points
  virtual PyObject* getPyResult(int key, float* weights);
  virtual PyObject* getItemAsPyObj(int key, int intPtIdx);

  // for values at whole elements
  virtual PyObject* getPyResult(int key);
};


/******************************************************************************
 * class FieldValueMapVector
 *
 * For storing vector and tensor values at nodes and whole elements
 * (anything but integration points).
 * For values at integration points use class FieldValueMapVectorToVec.
 */

class FieldValueMapVector : public FieldValueMapBase,
  public std::map<int, const float*>
{
protected:
  std::vector<float>* m_notDefinedValuePtr;
  std::vector<float> m_resVecFloat;

public:
  FieldValueMapVector(int nbComponents);
  virtual ~FieldValueMapVector();
  std::size_t size();
  int store(int key, float *val);
  void storeNotDefined(int key, float *val);
  int getNbComp();
  const float* getVector(int key);

  // for values at nodes
  void initResult();
  void addWeighted(int key, float weight);
  PyObject* getPyResult();

  // for values at whole elements
  // and for nodal fields without interpolation
  PyObject* getPyResult(int key);
};


/******************************************************************************
 * class FieldValueMapVectorToVec
 *
 * For storing vector and tensor values at integration points.
 */

class FieldValueMapVectorToVec : public FieldValueMapBase,
  public std::map<int, std::vector<const float*> >
{
protected:
  std::vector<const float*> m_notDefinedValue;
  std::vector<float>* m_notDefinedValuePtr;  // could be protected
  std::vector<float> m_resVecFloat;

public:
  FieldValueMapVectorToVec(int nbComponents);
  virtual ~FieldValueMapVectorToVec();
  std::size_t size();
  int store(int key, float *val);
  void storeNotDefined(int key, float * val);
  int getNbComp();
  const std::vector<const float*>& getVecOfVectors(int key);

  PyObject* getPyResult(int key, float* weights);
  PyObject* getItemAsPyObj(int key, int intPtIdx);
};


/******************************************************************************
 * class FieldValueMapVecComp
 *
 * For storing scalars or vector and tensor components, i.e. scalar values
 * at nodes and whole elements (anything but integration points).
 * For values at integration points use class FieldValueMapVecCompToVec.
 */

class FieldValueMapVecComp : public FieldValueMapBase,
  public std::map<int, float>
{
protected:
  float m_notDefinedValue;
  float m_resFloat;

public:
  FieldValueMapVecComp();
  std::size_t size();
  int store(int key, float *val);
  void storeNotDefined(int key, float * val);
  float getScalar(int key);

  // for values at nodes
  void initResult();
  void addWeighted(int key, float weight);
  PyObject* getPyResult();

  // for values at whole elements
  // and for nodal fields without interpolation
  PyObject* getPyResult(int key);
};


/******************************************************************************
 * class FieldValueMapVecCompToVec
 *
 * For storing scalars or vector and tensor components, i.e. scalar values
 * at integration points.
 */

class FieldValueMapVecCompToVec : public FieldValueMapBase,
  public std::map<int, std::vector<float> >
{
protected:
  std::vector<float> m_notDefinedValue;

public:
  std::size_t size();
  int store(int key, float *val);
  void storeNotDefined(int key, float * val);
  const std::vector<float>& getVecOfScalars(int key);
  PyObject* getPyResult(int key, float* weights);
  PyObject* getItemAsPyObj(int key, int intPtIdx);
};


/******************************************************************************
 * class FieldValueMapInt
 *
 * For storing integers (scalars) at whole elements. Nodal values can be stored
 * as well but methods for interpolating nodal values are not implemented.
 * Not suitable for integration point values.
 *
 * Note that some functions have the same name as in FieldValueMapBase but they
 * are not overloaded because they have different parameters.
 * E.g. store, storeNotDefined.
 * Method getScalar is missing alltogether because a different return type
 * would not separate it from the base class method and leads to a conflict.
 *
 * Why then inherits FieldValueMapInt from FieldValueMapBase at all?
 * Because InterpolationFromElem.interpolateToPyList requires a pointer to a
 * FieldValueMapBase object as argument.
 */

class FieldValueMapInt : public FieldValueMapBase,
  public std::map<int, int>
{
protected:
  int m_notDefinedValue;

public:
  FieldValueMapInt();
  std::size_t size();
  int store(int key, int val);
  void storeNotDefined(int key, int val);
  int getInt(int key);

  // for values at whole elements
  PyObject* getPyResult(int key);
};

#endif
