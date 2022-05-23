/******************************************************************************
*******************************************************************************
***   FieldValueMap... classes                                            ***
***                                                                         ***
***   dictionaries to hold field data from the odb                          ***

fieldValueMap-classes  ... maps node/elem-label -> scalar/vector/vec-of-vectors
 - FieldValueMapVector:         nodal, whole elem; vectors
 - FieldValueMapVectorToVec:    integr. pt       ; vectors
 - FieldValueMapVecComp:        nodal, whole elem; scalars
 - FieldValueMapVecCompToVec:   integr. pt       ; scalars

*******************************************************************************
******************************************************************************/

#include "fieldValueMap.h"
#include <cstdlib>
#include <cstdarg>

/* Python integration */
#include <Python.h>
#include <structmember.h>   // needed for e.g. T_INT

/******************************************************************************
 * class FieldValueMapBase
 *
 * Base class for all following classes to define the common interface.
 * Not to be instantiated directly.
 */
FieldValueMapBase::FieldValueMapBase(int nbComponents) :
  m_nbComponents(nbComponents) {}

FieldValueMapBase::~FieldValueMapBase() {}

int FieldValueMapBase::getNbComp() {
  return 0;
}

int FieldValueMapBase::store(int key, float *val) {
  fprintf(stderr, "FieldValueMapBase: int store() not implemented.\n");
  exit(255);
}

void FieldValueMapBase::storeNotDefined(int key, float * val) {
  fprintf(stderr, "FieldValueMapBase: storeNotDefined() not implemented.\n");
  exit(255);
}


float FieldValueMapBase::getScalar(int key) {
  fprintf(stderr, "FieldValueMapBase: float getScalar() not implemented.\n");
  exit(255);
}

const float* FieldValueMapBase::getVector(int key) {
  fprintf(stderr, "FieldValueMapBase: float* getVector() not implemented.\n");
  exit(255);
}

const std::vector<float>& FieldValueMapBase::getVecOfScalars
(int key) {
  fprintf(stderr, "FieldValueMapBase: std::vector<float>& getVecOfScalars()"
          " not implemented.\n");
  exit(255);
}

const std::vector<const float*>& FieldValueMapBase::getVecOfVectors
(int key) {
  fprintf(stderr, "FieldValueMapBase: std::vector<float*>& getVecOfVectors()"
          " not implemented.\n");
  exit(255);
}

// for nodal fields
void FieldValueMapBase::initResult() {
  fprintf(stderr, "FieldValueMapBase: initResult() not implemented.\n");
  exit(255);
}

void FieldValueMapBase::addWeighted(int key, float weight) {
  fprintf(stderr, "FieldValueMapBase: addWeighted(float weight)"
          " not implemented.\n");
  exit(255);
}

PyObject* FieldValueMapBase::getPyResult() {
  fprintf(stderr, "FieldValueMapBase: getPyResult() not implemented.\n");
  exit(255);
}

// for integration pt fields
PyObject* FieldValueMapBase::getPyResult(int key, float* weights) {
  fprintf(stderr, "FieldValueMapBase: getPyResult(int key, float* weights)"
          " not implemented.\n");
  exit(255);
}

PyObject* FieldValueMapBase::getItemAsPyObj(int key, int intPtIdx) {
  fprintf(stderr, "FieldValueMapBase: getItemAsPyObj(int key, int intPtIdx)"
          " not implemented.\n");
  exit(255);
}

// for whole element fields
PyObject* FieldValueMapBase::getPyResult(int key) {
  fprintf(stderr, "FieldValueMapBase: getPyResult(int key)"
          " not implemented.\n");
  exit(255);
}



/******************************************************************************
 * class FieldValueMapVector
 *
 * For storing vector and tensor values at nodes and whole elements
 * (anything but integration points).
 * For values at integration points use class FieldValueMapVectorToVec.
 */
FieldValueMapVector::FieldValueMapVector(int nbComponents) :
  FieldValueMapBase(nbComponents),
  m_notDefinedValuePtr(NULL),
  m_resVecFloat(nbComponents) {}

FieldValueMapVector::~FieldValueMapVector() {
  delete m_notDefinedValuePtr;
}

std::size_t FieldValueMapVector::size() {
  return std::map<int, const float*>::size();
}

int FieldValueMapVector::store(int key, float *val) {
  // The following would be similar. Except we could not detect overwrite.
  // (*this)[key] = val;

  // declaration of local var rv: return value of map::insert
  // which is (iterator, flag)
  std::pair< std::map<int, const float*>::iterator, bool > rv;
  rv = insert( std::make_pair(key, val) );
  if (!rv.second) {
    PyErr_Format(PyExc_RuntimeError,
                 "FieldValueMapVector::store() called second time for the same"
                 " key (node/element label) %d", key);
    return 1;
  }
  return 0;
}

void FieldValueMapVector::storeNotDefined(int key, float *val) {
  /* currently key (for not defined) is not used in any way */

  // store value-vector
  delete m_notDefinedValuePtr;
  m_notDefinedValuePtr = new std::vector<float>(val, val+m_nbComponents);
}

int FieldValueMapVector::getNbComp() {
  return m_nbComponents;
}

const float* FieldValueMapVector::getVector(int key) {
  // old: return (*this)[key];
  std::map<int, const float*>::const_iterator it = find( key );
  if (it == end()) {
    // return the vector *m_notDefinedValuePtr as float*
    return (const float*) &((*m_notDefinedValuePtr)[0]);
  }
  else {
    return it->second;
  }
}

// for nodal fields
void FieldValueMapVector::initResult() {
  for (int i=0; i<m_nbComponents; ++i) {
    m_resVecFloat[i] = 0.0;
  }
}

void FieldValueMapVector::addWeighted
(int key, float weight) {

  // find requested entry
  std::map<int, const float*>::const_iterator it = find( key );
  if (it == end()) {
    // add m_notDefinedValuePtr to resVecFloat
    for (int i=0; i<m_nbComponents; ++i) {
      m_resVecFloat[i] += weight * (*m_notDefinedValuePtr)[i];
    }
  }
  else {
    // add requested value to resVecFloat
    for (int i=0; i<m_nbComponents; ++i) {
      m_resVecFloat[i] += weight * (it->second)[i];
    }
  }
}

PyObject* FieldValueMapVector::getPyResult() {
  PyObject* resVal;
  resVal = PyList_New(m_nbComponents);
  if (resVal == NULL) {
    PyErr_SetString(PyExc_RuntimeError,
                        "FieldValueMapVector::getPyResult:"
                        " Could not create result vector.");
    return NULL;
  }
  for (int i=0; i<m_nbComponents; ++i) {
    PyList_SET_ITEM(resVal, i, PyFloat_FromDouble((double) m_resVecFloat[i]));
  }
  return resVal;
}

// for whole element fields
// and for nodal fields without interpolation
PyObject* FieldValueMapVector::getPyResult(int key) {
  PyObject* resVal;
  resVal = PyList_New(m_nbComponents);
  if (resVal == NULL) {
    PyErr_SetString(PyExc_RuntimeError,
                        "FieldValueMapVector::getPyResult:"
                        " Could not create result vector.");
    return NULL;
  }

  // find requested entry
  std::map<int, const float*>::const_iterator it = find( key );
  if (it == end()) {
    // store the vector *m_notDefinedValuePtr in resVal
    for (int i=0; i<m_nbComponents; ++i) {
      PyList_SET_ITEM(resVal, i, PyFloat_FromDouble((double) (*m_notDefinedValuePtr)[i]));
    }
  }
  else {
    // store element value in resVal
    for (int i=0; i<m_nbComponents; ++i) {
      PyList_SET_ITEM(resVal, i, PyFloat_FromDouble((double) (it->second)[i]));
    }
  }
  return resVal;
}


/******************************************************************************
 * class FieldValueMapVectorToVec
 *
 * For storing vector and tensor values at integration points.
 *
 * Typical use case:
 *   FieldValueMapVector fieldValMap(nbComp);
 *   for (int i=0; i<nbVals; ++i) {
 *      fieldValMap.store(labels[i], values + i*nbComp);
 *   }
 *   fieldValMap.storeNotDefined(notDefinedKey, notDefinedValue);
 *   value = fieldValMap.getVecOfVectors(key);
 *
 * Note that storeNotDefined() must be called *after* something has been
 * stored with store() for storeNotDefined() being able to determine the
 * number of integration points per key.
 */
FieldValueMapVectorToVec::FieldValueMapVectorToVec(int nbComponents) :
  FieldValueMapBase(nbComponents),
  m_notDefinedValuePtr(NULL),
  m_resVecFloat(nbComponents) {}

FieldValueMapVectorToVec::~FieldValueMapVectorToVec() {
  delete m_notDefinedValuePtr;
}

std::size_t FieldValueMapVectorToVec::size() {
  return std::map<int, std::vector<const float*> >::size();
}

int FieldValueMapVectorToVec::store(int key, float *val) {
  (*this)[key].push_back(val);
  return 0;
}

void FieldValueMapVectorToVec::storeNotDefined(int key, float * val) {
  /* Store the value that is to be returned for keys that are not in the
     FieldValueMap.

     Must be called *after* store(), see class description.

     Note: Currently key (for not defined) is not used in any way
   */

  // store value-vector as a copy of what's in val
  // This is the primary store to hold the notDefinedValue!
  delete m_notDefinedValuePtr;
  m_notDefinedValuePtr = new std::vector<float>(val, val+m_nbComponents);

  // create a vector of float-pointers into the array m_notDefinedValuePtr,
  // for each integration point one component
  // This is created for the method getVecOfVectors() to return a
  // tuple (vector) of values for each integration point.
  int nbIntegrationPts = begin()->second.size();
  m_notDefinedValue.clear();
  for (int i=0; i<nbIntegrationPts; ++i) {
    m_notDefinedValue.push_back(&((*m_notDefinedValuePtr)[0]));
  }
}

int FieldValueMapVectorToVec::getNbComp() {
  return m_nbComponents;
}

const std::vector<const float*>& FieldValueMapVectorToVec::getVecOfVectors
(int key) {
  // old: return (*this)[key];
  std::map<int, std::vector<const float*> >::const_iterator it = find( key );
  if (it == end()) {
    return m_notDefinedValue;
  }
  else {
    return it->second;
  }
}

PyObject* FieldValueMapVectorToVec::getPyResult(int key, float* weights) {

  // result is a list (representing a vector or tensor)
  PyObject* resVal;
  resVal = PyList_New(m_nbComponents);
  if (resVal == NULL) {
    PyErr_SetString(PyExc_RuntimeError,
                    "FieldValueMapVectorToVec::getPyResult:"
                    " Could not create result vector.");
    return NULL;
  }

  // find requested entry
  std::map<int, std::vector<const float*> >::const_iterator it = find( key );
  if (it == end()) {
    // not found ...
    // store m_notDefinedValue in resVal
    for (int i=0; i<m_nbComponents; ++i) {
      PyList_SET_ITEM(resVal, i, PyFloat_FromDouble((double) (*m_notDefinedValuePtr)[i]));
    }
  }
  else {

    std::vector<float> resVecFloat(m_nbComponents);
    for (int i=0; i<m_nbComponents; ++i) {
      resVecFloat[i] = 0.0;
    }
    
    // add requested value to resVecFloat
    int nbGaussPts = (int) (it->second).size();
    for (int i=0; i<nbGaussPts; ++i) {
      for (int j=0; j<m_nbComponents; ++j) {
        resVecFloat[j] += (it->second)[i][j] * weights[i];
      }
    }

    // convert resVecFloat to Python list resVal
    for (int i=0; i<m_nbComponents; ++i) {
      PyList_SET_ITEM(resVal, i, PyFloat_FromDouble((double) resVecFloat[i]));
    }

  }
  return resVal;
}

PyObject* FieldValueMapVectorToVec::getItemAsPyObj(int key, int intPtIdx) {

  // result is a list (representing a vector or tensor)
  PyObject* resVal;
  resVal = PyList_New(m_nbComponents);
  if (resVal == NULL) {
    PyErr_SetString(PyExc_RuntimeError,
                    "FieldValueMapVectorToVec::getPyResult:"
                    " Could not create result vector.");
    return NULL;
  }

  // find requested entry
  std::map<int, std::vector<const float*> >::const_iterator it = find( key );
  if (it == end()) {
    // not found ...
    // store m_notDefinedValue in resVal
    for (int i=0; i<m_nbComponents; ++i) {
      PyList_SET_ITEM(resVal, i, PyFloat_FromDouble((double) (*m_notDefinedValuePtr)[i]));
    }
  }
  else {

    // retrieve the value at the requested integration point
    const float* val = (it->second)[intPtIdx];
    
    // store requested value in Python list resVal
    for (int i=0; i<m_nbComponents; ++i) {
      PyList_SET_ITEM(resVal, i, PyFloat_FromDouble((double) val[i]));
    }
  }
  return resVal;
}


/******************************************************************************
 * class FieldValueMapVecComp
 *
 * For storing scalars or vector and tensor components, i.e. scalar values
 * at nodes and whole elements (anything but integration points).
 * For values at integration points use class FieldValueMapVecCompToVec.
 */
FieldValueMapVecComp::FieldValueMapVecComp() :
  m_notDefinedValue(0.0) {}

std::size_t FieldValueMapVecComp::size() {
  return std::map<int, float>::size();
}

int FieldValueMapVecComp::store(int key, float *val) {
  // The following would be similar. Except we could not detect overwrite.
  // (*this)[key] = *val;

  // declaration of local var rv: return value of map::insert
  // which is (iterator, flag)
  std::pair< std::map<int, float>::iterator, bool > rv;
  rv = insert( std::make_pair(key, *val) );
  if (!rv.second) {
    PyErr_Format(PyExc_RuntimeError,
                 "FieldValueMapVecComp::store() called second time for the same key"
                 " (node/element label) %d", key);
    return 1;
  }
  return 0;
}

void FieldValueMapVecComp::storeNotDefined(int key, float * val) {
  /* currently key (for not defined) is not used in any way */
  m_notDefinedValue = *val;
}

float FieldValueMapVecComp::getScalar(int key) {
  // old: return (*this)[key];
  std::map<int, float>::const_iterator it = find( key );
  if (it == end()) {
    return m_notDefinedValue;
  }
  else {
    return it->second;
  }
}

// for nodal fields
void FieldValueMapVecComp::initResult() {
  m_resFloat = 0.0;
}

void FieldValueMapVecComp::addWeighted
(int key, float weight) {

  // find requested entry
  std::map<int, float>::const_iterator it = find( key );
  if (it == end()) {
    // add m_notDefinedValue to resFloat
    m_resFloat += weight * m_notDefinedValue;
  }
  else {
    // add requested value to resFloat
    m_resFloat += weight * (it->second);
  }
}

PyObject* FieldValueMapVecComp::getPyResult() {
  PyObject* resVal;
  resVal = PyFloat_FromDouble((double) m_resFloat);
  return resVal;
}

// for whole element fields
// and for nodal fields without interpolation
PyObject* FieldValueMapVecComp::getPyResult(int key) {

  // find requested entry
  std::map<int, float>::const_iterator it = find( key );
  if (it == end()) {
    // return m_notDefinedValue
    return PyFloat_FromDouble((double) m_notDefinedValue);
  }
  else {
    // return element value
    return PyFloat_FromDouble((double) it->second);
  }
}


/******************************************************************************
 * class FieldValueMapVecCompToVec
 *
 * For storing scalars or vector and tensor components, i.e. scalar values
 * at integration points.
 */
std::size_t FieldValueMapVecCompToVec::size() {
  return std::map<int, std::vector<float> >::size();
}
int FieldValueMapVecCompToVec::store(int key, float *val) {
  (*this)[key].push_back(*val);
  return 0;
}
void FieldValueMapVecCompToVec::storeNotDefined(int key, float * val) {
  /* currently key (for not defined) is not used in any way */

  // create a vector, for each integration point one component
  int nbIntegrationPts = begin()->second.size();
  if (!nbIntegrationPts) nbIntegrationPts = 1;  // store at least one value!
  m_notDefinedValue.clear();
  for (int i=0; i<nbIntegrationPts; ++i) {
    m_notDefinedValue.push_back(*val);
  }
}
const std::vector<float>& FieldValueMapVecCompToVec::getVecOfScalars
(int key) {
  // old: return (*this)[key];
  std::map<int, std::vector<float> >::const_iterator it = find( key );
  if (it == end()) {
    return m_notDefinedValue;
  }
  else {
    return it->second;
  }
}

PyObject* FieldValueMapVecCompToVec::getPyResult(int key, float* weights) {
  PyObject* resVal;

  // find requested entry
  std::map<int, std::vector<float> >::const_iterator it = find( key );
  if (it == end()) {
    // store m_notDefinedValue in resVal
    resVal = PyFloat_FromDouble((double) m_notDefinedValue[0]);
  }
  else {
    int nbGaussPts = (int) (it->second).size();
    float resFloat = 0.0;

    // add weighted requested values to resFloat
    for (int i=0; i<nbGaussPts; ++i) {
      resFloat += (it->second)[i] * weights[i];
    }
    resVal = PyFloat_FromDouble((double) resFloat);
  }
  return resVal;
}

PyObject* FieldValueMapVecCompToVec::getItemAsPyObj(int key, int intPtIdx) {
  PyObject* resVal;

  // find requested entry
  std::map<int, std::vector<float> >::const_iterator it = find( key );
  if (it == end()) {
    // store m_notDefinedValue in resVal
    resVal = PyFloat_FromDouble((double) m_notDefinedValue[0]);
  }
  else {
    // retrieve the value at the requested integration point and convert to PyObject
    resVal = PyFloat_FromDouble((double) (it->second)[intPtIdx]);
  }
  return resVal;
}


/******************************************************************************
 * class FieldValueMapInt
 *
 * For storing integers (scalars) at nodes and whole elements (anything but
 * integration points).
 */
FieldValueMapInt::FieldValueMapInt() :
  m_notDefinedValue(0) {}

std::size_t FieldValueMapInt::size() {
  return std::map<int, int>::size();
}

int FieldValueMapInt::store(int key, int val) {
  // The following would be similar. Except we could not detect overwrite.
  // (*this)[key] = *val;

  // declaration of local var rv: return value of map::insert
  // which is (iterator, flag)
  std::pair< std::map<int, int>::iterator, bool > rv;
  rv = insert( std::make_pair(key, val) );
  if (!rv.second) {
    PyErr_Format(PyExc_RuntimeError,
                 "FieldValueMapInt::store() called second time for the same key"
                 " (node/element label) %d", key);
    return 1;
  }
  return 0;
}

void FieldValueMapInt::storeNotDefined(int key, int val) {
  /* currently key (for not defined) is not used in any way */
  m_notDefinedValue = val;
}

int FieldValueMapInt::getInt(int key) {
  // old: return (*this)[key];
  std::map<int, int>::const_iterator it = find( key );
  if (it == end()) {
    return m_notDefinedValue;
  }
  else {
    return it->second;
  }
}

// for "interpolating" whole element fields
PyObject* FieldValueMapInt::getPyResult(int key) {

  // find requested entry
  std::map<int, int>::const_iterator it = find( key );
  if (it == end()) {
    // return m_notDefinedValue
    return PyInt_FromLong((long) m_notDefinedValue);
  }
  else {
    // return element value
    return PyInt_FromLong((long) it->second);
  }
}
