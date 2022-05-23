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

#include "interpolation.h"
#include <cstdarg>
#include <vector>
#include <algorithm>

/* Debugging */
// #define DEBUG
#undef DEBUG

#ifdef DEBUG
#define debug_print(fmt, ...) fprintf(stderr, fmt, ##__VA_ARGS__)
#else
#define debug_print(fmt, ...) do {} while (0)
#endif

/******************************************************************************
 * class InterpolationBase
 *
 *
 */
InterpolationBase::InterpolationBase(int gaussptsPerPoint) :
  m_gaussptsPerPoint(gaussptsPerPoint) {}
// void InterpolationBase::error(const char* fmt, ...) {
//   va_list args;
//   va_start(args, fmt);
//   vfprintf(stderr, fmt, args);
//   va_end(args);
//   throw odb_Exception(odb_TEXT_MESSAGE);
// }
// int InterpolationBase::getNbOfGaussptsPerPoint() {
//   error("Interpolation::getNbOfGaussptsPerPoint() called for other than"
//         " InterpolationFromGaussPts.\n"
//         "This is an internal program error. Ask the programmer.\n");
// }


/******************************************************************************
 * class InterpolationFromNodal
 *
 *
 */
InterpolationFromNodal::InterpolationFromNodal(int nodesPerPoint, int nbPoints,
                                               char* ptWeights)
  : m_nodesPerPoint(nodesPerPoint),
    m_nbPoints(nbPoints),
    m_ptWeights(ptWeights),
    m_stridesNodes(sizeof(int)+sizeof(float))
{
  m_stridesPts = (sizeof(int)+sizeof(float))*nodesPerPoint;
  debug_print("InterpolationFromNodal-constructor:"
              " nodesPerPoint=%d, nbPoints=%d, m_stridesPts=%d\n",
              nodesPerPoint, nbPoints, m_stridesPts);
}

int InterpolationFromNodal::interpolateToPyList
(FieldValueMapBase* field, PyObject* result) {

  int errorFlag = 0;
  int node;
  float weight;
  char* ptData = m_ptWeights;
  PyObject* resVal;
  debug_print("InterpolationFromNodal.interpolateToPyList: nbPoints=%d,"
              " nbComp=%d\n", m_nbPoints, field->getNbComp());

  for (int ptIdx=0; ptIdx<m_nbPoints; ++ptIdx) {
    field->initResult();
    for (int i=0; i<m_nodesPerPoint; ++i) {
      node = *( (int*) ptData);
      weight = *( (float*) (ptData+sizeof(int)));
      field->addWeighted(node, weight);
      ptData += m_stridesNodes;
    }
    resVal = field->getPyResult();
    errorFlag = PyList_Append(result, resVal);
    Py_XDECREF(resVal);
    if (errorFlag) break;
    debug_print("InterpolationFromNodal.interpolateToPyList(vec): stored"
                " value for ptIdx %d in result list.\n", ptIdx);
  }
  debug_print("InterpolationFromNodal.interpolateToPyList: ending,"
              " errorFlag=%d\n", errorFlag);
  return errorFlag;
}

int InterpolationFromNodal::interpolateToPyListPieceWiseConst
(FieldValueMapBase* field, PyObject* result) {

  int errorFlag = 0;
  int node;
  float weight;
  float maxweight;
  char* ptData = m_ptWeights;
  PyObject* resVal;
  debug_print("InterpolationFromNodal.interpolateToPyListPieceWiseConst:"
              " nbPoints=%d, nbComp=%d\n", m_nbPoints, field->getNbComp());

  for (int ptIdx=0; ptIdx<m_nbPoints; ++ptIdx) {

    // initialize search for max weight
    node = *( (int*) ptData);
    maxweight = *( (float*) (ptData+sizeof(int)));
    ptData += m_stridesNodes;

    // search node with max weight, start with second
    for (int i=1; i<m_nodesPerPoint; ++i) {
      weight = *( (float*) (ptData+sizeof(int)));
      if (weight>maxweight) {
        maxweight = weight;
        node = *( (int*) ptData);
      }
      ptData += m_stridesNodes;
    }

    // get value from that node and store it in list
    resVal = field->getPyResult(node);
    errorFlag = PyList_Append(result, resVal);
    Py_XDECREF(resVal);
    if (errorFlag) break;
    debug_print("InterpolationFromNodal.interpolateToPyListPieceWiseConst(vec):"
                " stored value for ptIdx %d in result list.\n", ptIdx);
  }
  debug_print("InterpolationFromNodal.interpolateToPyListPieceWiseConst:"
              " ending, errorFlag=%d\n", errorFlag);
  return errorFlag;
}



/******************************************************************************
 * class InterpolationFromGaussPts
 *
 *
 */
InterpolationFromGaussPts::InterpolationFromGaussPts
(int gaussptsPerPoint, int nbPoints, char* ptWeights)
  : InterpolationBase(gaussptsPerPoint),
    m_nbPoints(nbPoints),
    m_ptWeights(ptWeights)
{
  m_stridesPts = (sizeof(int)+sizeof(float)*gaussptsPerPoint);
}

int InterpolationFromGaussPts::getNbOfGaussptsPerPoint() {
  return m_gaussptsPerPoint;
}

/* get the element number that has been stored in the ptWeights array
 * for a specific pointIdx */
int InterpolationFromGaussPts::getElem(int pointIdx) {
  return *((int*) (m_ptWeights+(pointIdx*m_stridesPts)));
};

int InterpolationFromGaussPts::interpolateToPyList
(FieldValueMapBase* field, PyObject* result) {

  int errorFlag = 0;
  int elem;
  float* weights;
  char* ptData = m_ptWeights;
  PyObject* resVal;

  for (int ptIdx=0; ptIdx<m_nbPoints; ++ptIdx) {
    elem = *((int*) ptData);
    ptData += sizeof(int);
    weights = (float*) ptData;
    ptData += m_gaussptsPerPoint*sizeof(float);
    resVal = field->getPyResult(elem, weights);
    errorFlag = PyList_Append(result, resVal);
    Py_DECREF(resVal);
    if (errorFlag) break;
  }
  debug_print("InterpolationFromGaussPts_interpolateToPyList:"
              " Interpolated results for %d points.\n", m_nbPoints);

  return errorFlag;
}

int InterpolationFromGaussPts::interpolateToPyListPieceWiseConst
(FieldValueMapBase* field, PyObject* result) {

  int errorFlag = 0;
  int elem;
  float* weights;
  char* ptData = m_ptWeights;
  int maxIdx;
  PyObject* resVal;

  for (int ptIdx=0; ptIdx<m_nbPoints; ++ptIdx) {
    elem = *((int*) ptData);
    ptData += sizeof(int);
    weights = (float*) ptData;
    ptData += m_gaussptsPerPoint*sizeof(float);
    maxIdx = std::max_element(weights, (float*) ptData) - weights;
    resVal = field->getItemAsPyObj(elem, maxIdx);
    errorFlag = PyList_Append(result, resVal);
    Py_DECREF(resVal);
    if (errorFlag) break;
  }
  debug_print("InterpolationFromGaussPts_interpolateToPyListPieceWiseConst:"
              " Interpolated results for %d points.\n", m_nbPoints);

  return errorFlag;
}

/******************************************************************************
 * class InterpolationFromElem
 *
 *
 */
InterpolationFromElem::InterpolationFromElem(int gaussptsPerPoint, int nbPoints,
                                             char* ptWeights)
  : m_nbPoints(nbPoints),
    m_ptWeights(ptWeights)
{
  m_stridesPts = (sizeof(int)+sizeof(float)*gaussptsPerPoint);
}

/* get the element number that has been stored in the ptWeights array
 * for a specific pointIdx */
int InterpolationFromElem::getElem(int pointIdx) {
  return *((int*) (m_ptWeights+(pointIdx*m_stridesPts)));
}

int InterpolationFromElem::interpolateToPyList
(FieldValueMapBase* field, PyObject* result) {

  int errorFlag = 0;
  int elem;
  char* ptData = m_ptWeights;
  PyObject* resVal;

  for (int ptIdx=0; ptIdx<m_nbPoints && !errorFlag; ++ptIdx) {

    /* store element value in result list */
    elem = *((int*) ptData);
    resVal = field->getPyResult(elem);
    errorFlag = PyList_Append(result, resVal);
    Py_DECREF(resVal);
    ptData += m_stridesPts;
  }
  return errorFlag;
}

int InterpolationFromElem::interpolateToPyListPieceWiseConst
(FieldValueMapBase* field, PyObject* result) {
  return interpolateToPyList(field, result);
}
