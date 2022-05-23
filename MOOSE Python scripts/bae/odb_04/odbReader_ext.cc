/*****************************************************************************
 *  Help recieving and processing data from the odbAccessServer through a    *
 *  pipe.                                                                    *
 *  This is Python / C extension module functionality that is independent    *
 *  from the Abaqus API.                                                     *
 *****************************************************************************/

/*****************************************************************************/
static const char
module_version[] = "1.08";
/* version history:

 1.07: first version as split-off from the old Abq-API odbReader_ext.cc
 1.08: GP added interpolationType argument for
       Communicator_updateStructPtFieldValues for piecewise constant
       interpolation, e.g. for the status field.
 *****************************************************************************/

/* structure

interpolation-classes

 subclasses:
 - InterpolationFromNodal:
 - InterpolationFromGaussPts:
 - InterpolationFromElem:

 common methods:
 - interpolateToPyList(fieldValueMap* field, result-array)
   . different method depending on subclass (nodal, integration pt, whole elem)
   . constructor takes appropriate point,weights-array
   . currently constant gaussptsPerPoint! in whole elem and integration pt

 */

#include <stdio.h>
#include <string.h>
#include <map>
#include <vector>
#include <string>

// for basic I/O using integer file handles
// needed to read from and write to pipes
#include <unistd.h>

// #include <stdlib.h>
// #include <assert.h>

/* Python integration */
#include <Python.h>
#include <structmember.h>   // needed for e.g. T_INT

/* utility modules */
#include "fieldValueMap.h"
#include "interpolation.h"

/* Debugging */
// #define DEBUG
#undef DEBUG

// #define DEBUG_MEMORYLEAK
#undef DEBUG_MEMORYLEAK

#ifdef DEBUG
#define debug_print(fmt, ...) fprintf(stderr, fmt, ##__VA_ARGS__)
#else
#define debug_print(fmt, ...) do {} while (0)
#endif

#ifdef DEBUG_MEMORYLEAK
/* prints /proc/self/statm
 * from man 5 proc:
              Provides information about memory usage, measured in pages.  The
              columns are:

                  size       (1) total program size
                             (same as VmSize in /proc/[pid]/status)
                  resident   (2) resident set size
                             (same as VmRSS in /proc/[pid]/status)
                  share      (3) shared pages (i.e., backed by a file)
                  text       (4) text (code)
                  lib        (5) library (unused in Linux 2.6)
                  data       (6) data + stack
                  dt         (7) dirty pages (unused in Linux 2.6)
 */
void debug_print_mem(char* where) {
  FILE* f;
  char status[1000];
  int i;

  f=fopen("/proc/self/statm", "r");
  fgets(status, 1000, f);
  fclose(f);

  for (i=0; i<1000; ++i) {
    if (status[i] == '\n') {
      status[i] = 0;
    }
  }
  fprintf(stderr, "%s: current /proc/self/statm: <%s>\n", where, status);
}
#else
#define debug_print_mem(where) do {} while (0)
#endif



/******************************************************************************
*******************************************************************************
***  Communicator class                                                     ***
*******************************************************************************
******************************************************************************/

/* data structure for the Communicator object
 */
typedef struct {
  PyObject_HEAD

  int dataChannelIn;
  int dataChannelOut;


  int nbPoints;
  /* in order to not need to copy the actual data in the python-string objects
   * the OdbReader object keeps a reference of the actual python strings */
  PyObject* pyStrNodalweights;
  PyObject* pyStrGaussPtsweights;
  InterpolationFromNodal* interpolationFromNodal;
  InterpolationFromGaussPts* interpolationFromGaussPts;
  InterpolationFromElem* interpolationFromElem;

} Communicator;

static PyMemberDef Communicator_members[] = {
  {const_cast<char *>("nbPoints"),
   T_INT, offsetof(Communicator, nbPoints), 0, const_cast<char *>
   ("Number of output points. Will be set by storeRemapWeights.")},
  {NULL}  /* Sentinel */
};


/******************************************************************************
*******************************************************************************
***  functions for the communication to odbAccessServer through the pipe    ***
*******************************************************************************
******************************************************************************/

/* Read single integer from the pipe.
 *
 * return zero on success.
 */
int readOneInt(const int fd, int *result) {
  int bytesRead;

  bytesRead = read(fd, result, sizeof(int));
  if (bytesRead<0) {
    PyErr_SetString(PyExc_IOError, "readOneInt: read from pipe failed.");
    return errno;
  }
  if (bytesRead != ((int) sizeof(int))) {
    PyErr_SetString(PyExc_ValueError, "readOneInt: value read from pipe is not"
                    " an int.");
    return 1;
  }
  return 0;
}


/******************************************************************************
*******************************************************************************
***  functions to get model data (nodes, elements, ...) from the odb        ***
*******************************************************************************
******************************************************************************/

/* Communicator_updateNodeCoordsFromPipe
 *
 * get nodes from the odb and store in the supplied dict
 * On the other side of the pipe it's OdbReader.getNodeCoordsFromOdb
 * writing the data to the pipe."
 */
static PyObject*
Communicator_updateNodeCoordsFromPipe(Communicator* self, PyObject *args) {
  debug_print("This is Communicator_updateNodeCoordsFromPipe...\n");

  /* get and check argument */
  PyObject* nodeCoords;  /* should be a PyDictObject */
  if (!PyArg_ParseTuple(args, "O:updateNodeCoordsFromPipe", &nodeCoords)) {
    return NULL;
  }
  if (!PyDict_Check(nodeCoords)) {
    PyErr_SetString(PyExc_ValueError,
        	    "Communicator_updateNodeCoordsFromPipe:"
        	    " nodeCoords argument must be a dict" );
    return NULL;
  }
  Py_XINCREF(nodeCoords);

  Py_BEGIN_ALLOW_THREADS

  /* determine dimensionality of the model and nb of nodes */
  int dim;
  int size;
  debug_print("Communicator_updateNodeCoordsFromPipe: now reading dim from"
              " datachannel fd=%d...\n", self->dataChannelIn);
  if (readOneInt(self->dataChannelIn, &dim)) {
    PyErr_SetString(PyExc_RuntimeError,
         "Communicator_updateNodeCoordsFromPipe: Couldn't determine"
         " dimensionality.");
    return NULL;
  }
  debug_print("Communicator_updateNodeCoordsFromPipe: got dim %d, now reading"
              " size...\n", dim);
  if (readOneInt(self->dataChannelIn, &size)) {
    PyErr_SetString(PyExc_RuntimeError,
         "Communicator_updateNodeCoordsFromPipe: Couldn't determine nb of"
         " nodes.");
    return NULL;
  }
  debug_print("Communicator_updateNodeCoordsFromPipe: dim %d, size %d\n",
              dim, size);

  int rc;
  int nodeLabel;
  PyObject* pyObjLabel;
  float coordinates[dim];
  PyObject* pyObjCoords;

  Py_BLOCK_THREADS
  for (int i=0; i<size; ++i) {

    if (readOneInt(self->dataChannelIn, &nodeLabel)) {
      PyErr_SetString(PyExc_RuntimeError,
           "Communicator_updateNodeCoordsFromPipe: Couldn't read node label.");
      return NULL;
    }
    rc = read(self->dataChannelIn, coordinates, dim*sizeof(float));
    if (rc != dim*( (int) sizeof(float) )) {
      PyErr_SetString(PyExc_RuntimeError,
           "Communicator_updateElNodesFromPipe: Couldn't read element"
           " coordinates.");
      return NULL;
    }
    
    pyObjCoords = PyList_New(dim);
    for (Py_ssize_t i=0; i<dim; ++i) {
      PyList_SET_ITEM(pyObjCoords, i,  // note comment on PyList_SET_ITEM below
        	      PyFloat_FromDouble((double) coordinates[i]));
    }
    pyObjLabel = PyInt_FromLong((long) nodeLabel);
    if (0 != PyDict_SetItem(nodeCoords, pyObjLabel, pyObjCoords)) {
      PyErr_SetString(PyExc_RuntimeError,
           "Communicator_updateElNodesFromPipe: Couldn't update node"
           " coordinates.");
      Py_DECREF(pyObjLabel);
      Py_DECREF(pyObjCoords);
      return NULL;
    }
    Py_DECREF(pyObjLabel);
    Py_DECREF(pyObjCoords);
  }
  Py_UNBLOCK_THREADS
  debug_print("Communicator_updateNodeCoordsFromPipe, ending.\n");

  Py_END_ALLOW_THREADS
  Py_XDECREF(nodeCoords);
  Py_RETURN_NONE;
}


/* Communicator_updateElNodesFromPipe
 *
 * get elements from the odb and store in the supplied dicts
 * arguments are in that order: elNodes, elType, typeEl
 */
static PyObject*
Communicator_updateElNodesFromPipe(Communicator* self, PyObject *args) {
  debug_print("This is Communicator_updateElNodesFromPipe...\n");

  /* get and check arguments */
  PyObject* elNodes;  /* should be a PyDictObject */
  PyObject* elType;  /* should be a PyDictObject */
  PyObject* typeEl;  /* should be a PyDictObject */
  if (!PyArg_ParseTuple(args, "OOO:updateElNodesFromPipe",
			&elNodes, &elType, &typeEl)) {
    return NULL;
  }
  if (!PyDict_Check(elNodes)) {
    PyErr_SetString(PyExc_ValueError,
       "Communicator_updateElNodesFromPipe: elNodes argument must be a dict");
  }
  if (!PyDict_Check(elType)) {
    PyErr_SetString(PyExc_ValueError,
       "Communicator_updateElNodesFromPipe: elType argument must be a dict");
  }
  if (!PyDict_Check(typeEl)) {
    PyErr_SetString(PyExc_ValueError,
       "Communicator_updateElNodesFromPipe: typeEl argument must be a dict");
  }

  /* get a private reference so they can't be deleted in the meantime */
  Py_XINCREF(elNodes);
  Py_XINCREF(elType);
  Py_XINCREF(typeEl);

  int errorFlag = 0;

  /* determine nb of elements */
  int size;
  if (readOneInt(self->dataChannelIn, &size)) {
    PyErr_SetString
      (PyExc_RuntimeError,
       "Communicator_updateElNodesFromPipe: Couldn't determine nb of"
       " elements.");

    /* release private references */
    Py_XDECREF(elNodes);
    Py_XDECREF(elType);
    Py_XDECREF(typeEl);
    errorFlag = 1;
    return NULL;
 }

  // data buffer for one element
  int elemLabel;
  int nbNodes;
  int nodeLabel;
  int sizeType;
  const int maxSizeType = 1000;
  char elemType[maxSizeType];
  int rc;

  PyObject* pyObjLabel;
  PyObject* pyObjNodes;
  PyObject* pyObjType;
  PyObject* pyObjElset;

  for (int i=0; i<size; ++i) {

    /* create element label (element number) */
    if (readOneInt(self->dataChannelIn, &elemLabel)) {
      PyErr_SetString
        (PyExc_RuntimeError,
         "Communicator_updateElNodesFromPipe: Couldn't determine element"
         " label.");
      errorFlag = 1;
      break;
    }
    pyObjLabel = PyInt_FromLong((long) elemLabel);

    /* create list of node numbers */
    if (readOneInt(self->dataChannelIn, &nbNodes)) {
      errorFlag = 1;
      break;
    }
    pyObjNodes = PyList_New(nbNodes);
    if (pyObjNodes == NULL) {
      PyErr_SetString(PyExc_RuntimeError,
         "Communicator_updateElNodesFromPipe: Could not create nodes list.");
      errorFlag = 1;
      Py_DECREF(pyObjLabel);
      break;
    }
    for (Py_ssize_t j=0; j<nbNodes; ++j) {
      if (readOneInt(self->dataChannelIn, &nodeLabel)) {
        errorFlag = 1;
        Py_DECREF(pyObjLabel);
        break;
      }
      PyList_SET_ITEM(pyObjNodes, j,  // note comment on PyList_SET_ITEM below
      		      PyInt_FromLong((long) nodeLabel));
    }
    if (errorFlag) break;

    /* store nodes in elNodes dict */
    if (0 != PyDict_SetItem(elNodes, pyObjLabel, pyObjNodes)) {
      errorFlag = 1;
      Py_DECREF(pyObjLabel);
      Py_DECREF(pyObjNodes);
      break;
    }
    Py_DECREF(pyObjNodes);

    /* read element type from pipe */
    if (readOneInt(self->dataChannelIn, &sizeType) || sizeType>maxSizeType-1) {
      errorFlag = 1;
      Py_DECREF(pyObjLabel);
      break;
    }
    rc = read(self->dataChannelIn, elemType, sizeType);
    if (rc != sizeType) {
      PyErr_SetString(PyExc_RuntimeError,
         "Communicator_updateElNodesFromPipe: Couldn't read element type.");
      errorFlag = 1;
      Py_DECREF(pyObjLabel);
      break;
    }
    elemType[sizeType] = 0;  /* string end marker */

    /* store element type in elType dict */
    pyObjType = PyString_FromString(elemType);
    if (pyObjType == NULL) {
      PyErr_SetString(PyExc_RuntimeError,
         "Communicator_updateElNodesFromPipe: Could not create type string.");
      errorFlag = 1;
      Py_DECREF(pyObjLabel);
      break;
    }
    if (0 != PyDict_SetItem(elType, pyObjLabel, pyObjType)) {
      errorFlag = 1;
      Py_DECREF(pyObjLabel);
      Py_DECREF(pyObjType);
      break;
    }

    /* store element in typeEl */
    pyObjElset = PyDict_GetItemString(typeEl, elemType);
    if (pyObjElset != NULL) {
      if (!PySet_Check(pyObjElset)) {
	PyErr_SetString(PyExc_RuntimeError,
	   "Communicator_updateElNodesFromPipe: Found non-set in typeEl");
	errorFlag = 1;
	Py_DECREF(pyObjLabel);
	Py_DECREF(pyObjType);
	break;
      }
      // PyDict_GetItemString returns borrowed reference, get private reference
      // so that we can give it back at the end as in the case that it's a
      // newly created set object
      Py_INCREF(pyObjElset);
    }
    else {
      pyObjElset = PySet_New(NULL);
      if (pyObjElset == NULL) {
	PyErr_SetString(PyExc_RuntimeError,
         "Communicator_updateElNodesFromPipe: Could not create set for new type.");
	errorFlag = 1;
	Py_DECREF(pyObjLabel);
	Py_DECREF(pyObjType);
	break;
      }
      if (0 != PyDict_SetItemString(typeEl, elemType, pyObjElset)) {
	errorFlag = 1;
	Py_DECREF(pyObjLabel);
	Py_DECREF(pyObjType);
	Py_DECREF(pyObjElset);
	break;
      }
    }
    if (0 != PySet_Add(pyObjElset, pyObjLabel)) {
      errorFlag = 1;
      Py_DECREF(pyObjLabel);
      Py_DECREF(pyObjType);
      Py_DECREF(pyObjElset);
      break;
    }

    Py_DECREF(pyObjLabel);
    Py_DECREF(pyObjType);
    Py_DECREF(pyObjElset);
  }
  debug_print("Communicator_updateElNodesFromPipe, ending, errorFlag=%d.\n",
	      errorFlag);


  /* release private references */
  Py_XDECREF(elNodes);
  Py_XDECREF(elType);
  Py_XDECREF(typeEl);

  if (errorFlag) {
    return NULL;
  }
  Py_RETURN_NONE;
}
/* comment on PyList_SET_ITEM vs. PyList_SetItem
 * from http://stackoverflow.com/questions/10305327/pylist-setitem-vs-pylist-setitem, see also: http://stackoverflow.com/questions/17635782/difference-between-pyint-fromlong-and-py-buildvalue
 *
 * PyList_SET_ITEM is an unsafe macro that basically sticks an object into the list's internal pointer array without any bound checks. If anything non-NULL is in the ith position of the list, a reference leak will occur. PyList_SET_ITEM steals the reference to the object you put in the list. PyList_SetItem also steals the reference, but it checks bounds and decrefs anything which may be in the ith position. The rule-of-thumb is use PyList_SET_ITEM to initialize lists you've just created and PyList_SetItem otherwise.
 */


/******************************************************************************
*******************************************************************************
***  functions for preparation of interpolation                             ***
*******************************************************************************
******************************************************************************/

/* Communicator_storeRemapWeights
 *
 * store the weighting factors for nodal and Gauss-pt interpolation
 * create InterpolationFromNodal object self->interpolationFromNodal
 * and InterpolationFromGaussPts object self->interpolationFromGaussPts
 *
 * Call this function to specify the points at which you want to get field
 * value output by subsequent calls to getFieldDataForPoints.
 *
 * Arguments:
 * nbPoints: number of points for which field output is being queried.
      Determines the length of the following two arrays.
 * nodalWeights: A string that will be interpreted as a binary array of node
 *    numbers and associated weight factors for each of the output points.
 *    See bae.odb_03.odbReader.OdbReaderWorker.initOutputPoints() and
 *    bae.mesh_01.InterpolMeshToPoints.writeRemapParam for further details.
 * gaussPointsWeights: A string that will be interpreted as a binary array of
 *    element numbers and associated weight factors for each of the output
 *    points.
 */
static PyObject*
Communicator_storeRemapWeights(Communicator* self, PyObject *args) {

  debug_print("This is Communicator_storeRemapWeights...\n");

  /* clear OdbReader object from possible old data that is supposed to be
   * initialized by this function */
  Py_XDECREF(self->pyStrNodalweights);
  Py_XDECREF(self->pyStrGaussPtsweights);
  self->pyStrNodalweights = NULL;
  self->pyStrGaussPtsweights = NULL;
  if (self->interpolationFromNodal) {
    delete self->interpolationFromNodal;
    self->interpolationFromNodal = NULL;
  }
  if (self->interpolationFromGaussPts) {
    delete self->interpolationFromGaussPts;
    self->interpolationFromGaussPts = NULL;
  }
  if (self->interpolationFromElem) {
    delete self->interpolationFromElem;
    self->interpolationFromElem = NULL;
  }

  /* get and check arguments */
  int notDefinedKey = -1;
  if (!PyArg_ParseTuple(args, "iOO|i:storeRemapWeights",
			&(self->nbPoints),
			&(self->pyStrNodalweights),
			&(self->pyStrGaussPtsweights),
			&notDefinedKey)) {
    return NULL;
  }
  Py_XINCREF(self->pyStrNodalweights);
  Py_XINCREF(self->pyStrGaussPtsweights);

  /* get data pointer and length of nodal weights string */
  Py_ssize_t length;
  char* dataPointer;
  int nodesPerPoint;
  if (-1 == PyString_AsStringAndSize
      (self->pyStrNodalweights, &dataPointer, &length)) {
    PyErr_SetString(PyExc_ValueError,
       "Communicator_storeRemapWeights: nodalweights argument must be a string.");
    return NULL;
  }
  nodesPerPoint = length / (sizeof(int)+sizeof(float)) / self->nbPoints;
  debug_print("Communicator_storeRemapWeights: pyStrNodalweights has %d bytes"
	      " for %d points, so we have %d nodes per element.\n",
	      (int) length, self->nbPoints, nodesPerPoint);

  /* create interpolationFromNodal */
  self->interpolationFromNodal = new InterpolationFromNodal
    (nodesPerPoint, self->nbPoints, dataPointer);
  if (!self->interpolationFromNodal) {
    PyErr_SetString(PyExc_ValueError,
		    "Communicator_storeRemapWeights: could not create"
		    " interpolationFromNodal.");
    return NULL;
  }
  debug_print("Communicator_storeRemapWeights:"
              " Created InterpolationFromNodal for %d points.\n",
	      self->interpolationFromNodal->m_nbPoints);

  /* get data pointer and length of Gauss pts weights string */
  int gaussptsPerPoint;
  if (-1 == PyString_AsStringAndSize
      (self->pyStrGaussPtsweights, &dataPointer, &length)) {
    PyErr_SetString(PyExc_ValueError,
       "Communicator_storeRemapWeights: gaussPtsweights argument must be a"
       " string.");
    return NULL;
  }
  gaussptsPerPoint = (length/self->nbPoints - sizeof(int)) / sizeof(float);
  debug_print("Communicator_storeRemapWeights: pyStrGaussPtsweights has %d bytes"
	      " for %d points, so we have %d Gauss pts per element\n",
	      (int) length, self->nbPoints, gaussptsPerPoint);

  /* create interpolationFromGaussPts */
  self->interpolationFromGaussPts = new InterpolationFromGaussPts
    (gaussptsPerPoint, self->nbPoints, dataPointer);
  if (!self->interpolationFromGaussPts) {
    PyErr_SetString(PyExc_ValueError,
		    "Communicator_storeRemapWeights: could not create"
		    " interpolationFromGaussPts.");
    return NULL;
  }

  /* create interpolationFromElem */
  self->interpolationFromElem = new InterpolationFromElem
    (gaussptsPerPoint, self->nbPoints, dataPointer);
  if (!self->interpolationFromElem) {
    PyErr_SetString(PyExc_ValueError,
		    "Communicator_storeRemapWeights: could not create"
		    " interpolationFromElem.");
    return NULL;
  }

  debug_print("Communicator_storeRemapWeights, ending.\n");
  Py_RETURN_NONE;
}


/******************************************************************************
*******************************************************************************
***  functions for receiving field values from the pipe                     ***
*******************************************************************************
******************************************************************************/

/* service functions for Communicator_updateFieldValues
 */

/* updateFieldValues_storeIntPtSectPt:
 *
 * This function stores the intPtSectPt array in fieldvalues -the python list
 * object (of type bae.field_01.Field) to be updated by updateFieldValues.
 * It is called only if intPtSectPt is periodic and consistent with all
 * nb of Gauss-pts * nb of section pts.
 */
int updateFieldValues_storeIntPtSectPt
(PyObject* fieldvalues,
 int intPtSectPt[], int nbIntPtSectPt)
{
  int errorFlag;
  PyObject* resVal;

  //... create the list object
  PyObject* pyIntPtSectPt;
  pyIntPtSectPt = PyList_New(nbIntPtSectPt);
  if (!pyIntPtSectPt) {
    PyErr_SetString(PyExc_RuntimeError,
		    "Communicator_updateFieldValues:"
		    " Could not create intPtSectPt list attribute.");
    return 1;
  }

  //... attach the list object intPtSectPt to the result object fieldvalues
  errorFlag = PyObject_SetAttrString
    (fieldvalues, "intPtSectPt", pyIntPtSectPt);
  if (errorFlag) {
    PyErr_SetString(PyExc_RuntimeError,
		    "Communicator_updateFieldValues: Could not assign"
		    " intPtSectPt attribute to result field.");
    Py_XDECREF(pyIntPtSectPt);
    return 2;
  }

  //... fill the list objects with (Gauss pt, section pt) tuples
  for (int i=0; i<nbIntPtSectPt; ++i) {
    resVal = Py_BuildValue("(ii)",
			   intPtSectPt[i*2 + 0], intPtSectPt[i*2 + 1]);
    if (!resVal) {
      PyErr_Format(PyExc_RuntimeError,
		   "Communicator_updateFieldValues: Could not create"
		   " %dth tuple item for intPtSectPt.", i);
      Py_XDECREF(pyIntPtSectPt);
      return 3;
    }
    PyList_SET_ITEM(pyIntPtSectPt, i, resVal);
  }
  return 0;
}


/* updateFieldValues_storeItem:
 *
 * For the i-th value from the odbAccessServer find its place in the result
 * dictionary fieldvalues and store the corresponding value item.
 */
int updateFieldValues_storeItem
(PyObject* fieldvalues,
 PyObject* item,
 int i, int labels[],
 int intPtSectPt[],
 int nbIntPtSectPt,
 bool periodicIntPtSectPt,
 std::map<int,int>* elemToNbVals
 ) {

  // create the label (dict-key)
  PyObject* label = PyInt_FromLong((long) labels[i]);
  if (!label) {
    return 1;
  }

  // find dict entry for label (this is a list, for sure)
  // resVal is a borrowed reference, don't DECREF resVal!
  PyObject* resVal = PyDict_GetItem(fieldvalues, label);
  Py_XDECREF(label);
  if (!resVal) {
    PyErr_Format
      (PyExc_RuntimeError,
       "Communicator_updateFieldValues: For the %dth field value element"
       " %d there is no entry in the result dict. Ask the programmer.",
       i+1, labels[i]);
    return 1;
  }

  // determine Gauss pt number
  int intPt = intPtSectPt[(i%nbIntPtSectPt)*2 + 0];
  if (intPt < 1 || intPt > PyList_GET_SIZE(resVal)) {
    PyErr_Format
      (PyExc_ValueError,
       "Communicator_updateFieldValues: For the %dth field value"
       " expected Gauss pt nb between 1 and %d but got %d instead.",
       i+1, (int) PyList_GET_SIZE(resVal), intPt);
    return 1;
  }

  // determine section pt number
  int sectPt = intPtSectPt[(i%nbIntPtSectPt)*2 + 1];

  // determine index in result list
  int resIdx;
  if (periodicIntPtSectPt) {
    resIdx = i%nbIntPtSectPt;
  }
  else {   // not periodicIntPtSectPt
    if (sectPt > 0) {
      resIdx = ((*elemToNbVals)[labels[i]])++;
    }
    else {
      // we have no section points -> index in result = index of Gauss pt
      resIdx = intPt-1;
    }
  }

  // and finally store the item
  int errorFlag = PyList_SetItem(resVal, resIdx, item);
  debug_print("Communicator_updateFieldValues:"
	      " position elemIP, elem %d, intPt %d, sectPt %d"
	      " in fieldvalues dict, error:%d.\n",
	      labels[i], intPt, sectPt, errorFlag);
  
  return errorFlag;
}


/* Communicator_updateFieldValues
 *
 * Receive field values from the pipe and update the field values object.
 *
 * If possible create a intPtSectPt attribute:
 * list of (Gauss pt number, section point number) tuples.
 *
 * Arguments (function accepts keyword arguments):
 * position: "node", "element" or "elemIP"
 * fieldvalues: a bae.field_01.Field / dict object to take the results.
 */
static PyObject*
Communicator_updateFieldValues
(Communicator* self, PyObject *args, PyObject *keywds) {

  size_t rc;
  int* intPtSectPt = 0;
  int* labels = 0;
  float* values = 0;
  debug_print("This is Communicator_updateFieldValues...\n");

  /* get and check arguments */
  const char* position = 0;
  PyObject* fieldvalues = 0;
  static const char *kwlist[] = {
    "position", "fieldvalues", NULL};
  if (!PyArg_ParseTupleAndKeywords
      (args, keywds,
       "sO:updateFieldValues", const_cast<char **>(kwlist),
       &position, &fieldvalues)) {
    return NULL;
  }
  if (!PyDict_Check(fieldvalues)) {
    PyErr_SetString(PyExc_ValueError,
		    "Communicator_updateFieldValues: Expected dict"
		    " object as second argument fieldvalues.");
    return NULL;
  }
  Py_XINCREF(fieldvalues);

  /* Read data from the pipe:
   * (I): integration point / section point numbers
   *
   * - nbIntPtSectPt: number of integration point / section point tuples in the
   *   following array. An integer. E.g. 1 for nodal values and 4 for SDV1 on
   *   C3D10M elements.
   * - (integration point, section point) - tuples: 2*nbIntPtSectPt integers.
   */
  int nbIntPtSectPt;
  if (readOneInt(self->dataChannelIn, &nbIntPtSectPt)) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Communicator_updateFieldValues: Couldn't read"
                    " nbIntPtSectPt.");
    Py_XDECREF(fieldvalues);
    return NULL;
  }

  intPtSectPt = new int[2*nbIntPtSectPt];
  if (!intPtSectPt) {
    PyErr_SetString
      (PyExc_MemoryError,
       "Communicator_updateFieldValues: Couldn't allocate"
       " memory for list of integration- and section point numbers.");
    Py_XDECREF(fieldvalues);
    return NULL;
  }
  rc = read(self->dataChannelIn, intPtSectPt,
            2*nbIntPtSectPt*sizeof(int));
  if (rc != 2*nbIntPtSectPt*sizeof(int)) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Communicator_updateFieldValues: Couldn't read"
                    " list of integration- and section point numbers.");
    Py_XDECREF(fieldvalues);
    delete[] intPtSectPt;
    return NULL;
  }
  debug_print("Communicator_updateFieldValues:"
              " Read %d items for intPtSectPt from pipe\n", nbIntPtSectPt);
  ///#############################################################  DEBUG
  // for (int i=0; i<nbIntPtSectPt; ++i) {
  //   debug_print("Communicator_updateFieldValues:"
  //               " Read intPtSectPt from pipe: intpt %d, sectpt %d\n",
  //               intPtSectPt[i*2+0], intPtSectPt[i*2+1]);
  // }
  ///#############################################################  DEBUG

  /* Read data from the pipe:
   * (II): field values
   *
   * Creates:
   * - nbComp: nb of components of each field value: 1 for scalars, 3 for
   *   vectors, 6 for tensors
   * - nbVals: nb of field values (for nodal: nb of nodes, for Gauss-pt: nb of
   *   elems*4 (C3D10M), in general nb of elems * nb of Gauss-pts per elem * nb
   *   of section pts per Gauss-pt)
   * - labels: integer array of node or element labels
   * - values: float array of field values.
   */
  int nbComp;
  int nbVals;
  if (readOneInt(self->dataChannelIn, &nbComp)) {
    PyErr_SetString(PyExc_RuntimeError,
		    "Communicator_updateFieldValues: Couldn't read"
		    " number of components of each field value.");
    Py_XDECREF(fieldvalues);
    delete[] intPtSectPt;
    return NULL;
  }
  if (readOneInt(self->dataChannelIn, &nbVals)) {
    PyErr_SetString(PyExc_RuntimeError,
		    "Communicator_updateFieldValues: Couldn't read"
		    " number of field values.");
    Py_XDECREF(fieldvalues);
    delete[] intPtSectPt;
    return NULL;
  }
  debug_print("Communicator_updateFieldValues:"
              " Number of components: %d. Number of field values: %d.\n",
              nbComp, nbVals);

  labels = new int[nbVals];
  if (!labels) {
    PyErr_SetString(PyExc_MemoryError,
                    "Communicator_updateFieldValues: Couldn't allocate"
                    " memory for node- or element- labels.");
    Py_XDECREF(fieldvalues);
    delete[] intPtSectPt;
    return NULL;
  }
  rc = read(self->dataChannelIn, labels, nbVals*sizeof(int));
  if (rc != nbVals*sizeof(int)) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Communicator_updateFieldValues: Couldn't read"
                    " node- or element- labels.");
    Py_XDECREF(fieldvalues);
    delete[] intPtSectPt;
    delete[] labels;
    return NULL;
  }
  values = new float[nbVals*nbComp];
  if (!values) {
    PyErr_SetString
      (PyExc_MemoryError,
       "Communicator_updateFieldValues: Couldn't allocate"
       " memory for field values.");
    Py_XDECREF(fieldvalues);
    delete[] intPtSectPt;
    delete[] labels;
    return NULL;
  }
  rc = read(self->dataChannelIn, values, nbVals*nbComp*sizeof(float));
  if (rc != nbVals*nbComp*sizeof(float)) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Communicator_updateFieldValues: Couldn't read"
                    " field values from the pipe.");
    Py_XDECREF(fieldvalues);
    delete[] intPtSectPt;
    delete[] labels;
    delete[] values;
    return NULL;
  }
  debug_print("Communicator_updateFieldValues:"
              " Finished reading from the pipe: %d labels: %d...,\n"
              "    ...   %d value-floats, %d components per value: %g....\n",
              nbVals, labels[0], nbVals*nbComp, nbComp, values[0]);
  ///#############################################################  DEBUG
  // for (int i=0; i<nbVals; ++i) {
  //   debug_print("label %d, value %g, intpt %d, sectpt %d\n",
  //               labels[i], values[i],
  //               intPtSectPt[(i%nbIntPtSectPt)*2+0],
  //               intPtSectPt[(i%nbIntPtSectPt)*2+1]);
  // }
  ///#############################################################  DEBUG

  int errorFlag = 0;

  /* select field-value-map type and initialize values
   * type and procedure depending on on position argument and nbComp */
  if (strcmp(position, "node")==0 || strcmp(position, "element")==0) {

    /* check Gauss-pt and section pt number */
    if (nbIntPtSectPt!=1) {
      PyErr_Format
        (PyExc_NotImplementedError,
         "Communicator_updateFieldValues: For nodal and whole"
         " element values expect nbIntPtSectPt==1 but got %d.",
         nbIntPtSectPt);
      errorFlag = 1;
    }
    if (!errorFlag && intPtSectPt[0]!=0) {
      PyErr_Format
        (PyExc_NotImplementedError,
         "Communicator_updateFieldValues: For nodal and whole"
         " element values expect Gauss pt nb 0 but got %d instead.",
         intPtSectPt[0]);
      errorFlag = 1;
    }
    if (!errorFlag && intPtSectPt[1]>0) {
      PyErr_Format
        (PyExc_NotImplementedError,
         "Communicator_updateFieldValues: For nodal or whole"
         " element field value got section point nb %d. Not"
         " implemented!", intPtSectPt[1]);
      errorFlag = 1;
    }
    if (errorFlag) {
      Py_XDECREF(fieldvalues);
      delete[] intPtSectPt;
      delete[] labels;
      delete[] values;
      return NULL;
    }
    debug_print("Communicator_updateFieldValues:"
                " position node or element, checked intPtSectPt: ok.\n");
    
    /* store values in dict fieldvalues */
    if (nbComp>1) {  // vector/tensor values
      PyObject* label;
      PyObject* resVal;
      PyObject* item;
      for (int i=0; i<nbVals && !errorFlag; ++i) {
      
        // create the label (dict-key)
        label = PyInt_FromLong((long) labels[i]);
        if (!label) {
          errorFlag = 1;
          break;
        }
      
        // create the vector/tensor (list) - value
        resVal = PyList_New(nbComp);
        if (!resVal) {
          PyErr_SetString(PyExc_RuntimeError,
                          "Communicator_updateFieldValues:"
                          " Could not create result vector.");
          Py_XDECREF(label);
          errorFlag = 1;
          break;
        }
        for (int j=0; j<nbComp && !errorFlag; ++j) {
          item = PyFloat_FromDouble((double) values[i*nbComp+j]);
          if (!item) {
            errorFlag = 1;
            break;  // exit only the inner for j... loop
          }
          // PyList_SET_ITEM is the version *without* error checking!
          // PyList_SET_ITEM steals a reference => don't DECREF item!
          PyList_SET_ITEM(resVal, j, item);
        }
        if (errorFlag) {
          Py_XDECREF(label);
          Py_XDECREF(resVal);
          break;
        }
      
        // store value in dict
        errorFlag = PyDict_SetItem(fieldvalues, label, resVal);
        Py_XDECREF(label);
        Py_XDECREF(resVal);
      }
    }
    else {  // nbComp==1, scalar values
      PyObject* label = 0;
      PyObject* resVal = 0;
      for (int i=0; i<nbVals && !errorFlag; ++i) {
    
        // create the label (dict-key)
        label = PyInt_FromLong((long) labels[i]);
        if (!label) {
          errorFlag = 1;
          break;
        }
    
        // create scalar value
        resVal = PyFloat_FromDouble((double) values[i]);
        if (!resVal) {
          Py_XDECREF(label);
          errorFlag = 1;
          break;
        }
      
        // store value in dict
        errorFlag = PyDict_SetItem(fieldvalues, label, resVal);
        debug_print("Communicator_updateFieldValues:"
                    "after setitem, errorFlag=%d\n", errorFlag);
        Py_XDECREF(label);
        Py_XDECREF(resVal);
      }
    }
    if (errorFlag) {
      Py_XDECREF(fieldvalues);
      delete[] intPtSectPt;
      delete[] labels;
      delete[] values;
      return NULL;
    }
    debug_print("Communicator_updateFieldValues:"
                " position node or element, stored values in fieldvalues"
                " (nbComp=%d), error:%d.\n", nbComp, errorFlag);
  }

  else if (strcmp(position, "elemIP")==0) {
    PyObject* label = 0;
    PyObject* resVal = 0;

    // create elemToNbVals:
    // number of values per element (nb of Gauss pts * section pts)
    // elemToNbVals is on the stack because it's got millions of items
    std::map<int,int>* elemToNbVals = 0;
    elemToNbVals = new std::map<int,int>;
    if (!elemToNbVals) {
      Py_XDECREF(fieldvalues);
      delete[] intPtSectPt;
      delete[] labels;
      delete[] values;
      return NULL;
    }
    for (int i=0; i<nbVals; ++i) {
      (*elemToNbVals)[labels[i]] += 1;
    }
    debug_print("Communicator_updateFieldValues:"
                " position elemIP, the %d values refer to %d elements.\n",
                nbVals, (int) elemToNbVals->size());

    // check if all elements have the same number of values (same nb of
    // Gauss-pts * nb of sections pts) and if this is number is identical to
    // nbIntPtSectPt (nb of items in intPtSectPt array)
    bool periodicIntPtSectPt = true;
    for (std::map<int,int>::iterator it=elemToNbVals->begin();
	 it != elemToNbVals->end() && periodicIntPtSectPt;
         ++it) {
      periodicIntPtSectPt &= (it->second == nbIntPtSectPt);
      // debug_print("Communicator_updateFieldValues: Checking periodicity:"
      // 		  " elemToNbVals[el=%d]=%d, nbIntPtSectPt %d.\n",
      // 		  it->first, it->second, nbIntPtSectPt);
    }
    debug_print("Communicator_updateFieldValues: Checked consistency of"
                " elemToNbVals with length of intPtSectPt. It is %s.\n",
                periodicIntPtSectPt
		? "consistent and intPtSectPt is periodic"
		: "not consistent or intPtSectPt is not periodic");

    if (periodicIntPtSectPt) {
      // If intPtSectPt is periodic and consistent with all
      // nb of Gauss-pts * nb of section pts
      // then store intPtSectPt array in result object fieldvalues
      errorFlag = updateFieldValues_storeIntPtSectPt
	(fieldvalues, intPtSectPt, nbIntPtSectPt);
      if (errorFlag) {
	Py_XDECREF(fieldvalues);
	delete[] intPtSectPt;
	delete[] labels;
	delete[] values;
	return NULL;
      }
      debug_print("Communicator_updateFieldValues: stored intPtSectPt (1).\n");
    }
    else {
      // create condensedIntPtSectPt:
      // condensedIntPtSectPt on the heap because we stop it eventually
      bool consistentIntPtSectPt = true;
      std::vector<int> condensedIntPtSectPt;
      std::map<int,int>* condensedElemToNbVals = 0;
      condensedElemToNbVals = new std::map<int,int>;
      if (!condensedElemToNbVals) {
	Py_XDECREF(fieldvalues);
	delete[] intPtSectPt;
	delete[] labels;
	delete[] values;
	return NULL;
      }
      for (int i=0; i<nbIntPtSectPt; ++i) {
	// get current number of values for this element as idxForElem
	int idxForElem = (*condensedElemToNbVals)[labels[i]];
	(*condensedElemToNbVals)[labels[i]] += 1;

	if (idxForElem<((int) condensedIntPtSectPt.size())/2) {
	  // we already have an item in condensedIntPtSectPt for this idxForElem
	  // check that intPtSectPt is still consistent with
	  // condensedIntPtSectPt
	  if (intPtSectPt[i*2+0] != condensedIntPtSectPt[idxForElem*2+0]
	      or intPtSectPt[i*2+1] != condensedIntPtSectPt[idxForElem*2+1]) {
	    debug_print("Communicator_updateFieldValues:"
			" For the %d. field value (of %d to be checked)"
			" --the %d. val for elem %d-- found inconsistency:"
			" intPtSectPt: %d, %d  / condensedIntPtSectPt: %d, %d\n",
			i+1, nbIntPtSectPt,
			idxForElem+1, labels[i],
			intPtSectPt[i*2+0], intPtSectPt[i*2+1],
			condensedIntPtSectPt[idxForElem*2+0],
			condensedIntPtSectPt[idxForElem*2+1]);
	    
	    consistentIntPtSectPt = false;
	    break;
	  }
	}
	else {
	  // no corresponding entry in condensedIntPtSectPt yet
	  // -> store a new item
	  condensedIntPtSectPt.push_back(intPtSectPt[i*2+0]);
	  condensedIntPtSectPt.push_back(intPtSectPt[i*2+1]);
	}
      } // end for i=0 to nbIntPtSectPt
      debug_print("Communicator_updateFieldValues: Checked consistency of"
		  " intPtSectPt for all elements. It is %sconsistent.\n",
		  consistentIntPtSectPt ? "" : "not ");

      if (consistentIntPtSectPt) {
	// If intPtSectPt is consistent for all elements (but not periodic in
	// the sense checked earlier) ...
	// then store intPtSectPt array in result object fieldvalues
	errorFlag = updateFieldValues_storeIntPtSectPt
	  (fieldvalues,
	   condensedIntPtSectPt.data(), condensedIntPtSectPt.size()/2);
	if (errorFlag) {
	  Py_XDECREF(fieldvalues);
	  delete condensedElemToNbVals;
	  delete[] intPtSectPt;
	  delete[] labels;
	  delete[] values;
	  return NULL;
	}
	debug_print("Communicator_updateFieldValues:"
		    " stored intPtSectPt (2).\n");
      }
      delete condensedElemToNbVals;
    }  // end if (periodicIntPtSectPt) ... else

    // create list entries in fieldvalues
    // after that we have the dictionary filled with all keys but only empty
    // lists as values so far.
    for (std::map<int,int>::iterator it=elemToNbVals->begin();
         it != elemToNbVals->end() && !errorFlag;
         ++it) {

      // create the label (dict-key)
      label = PyInt_FromLong((long) it->first);
      if (!label) {
        errorFlag = 1;
        break;
      }

      // create the list for IP values
      // it->second is the number of values per element
      resVal = PyList_New(it->second);
      if (!resVal) {
        PyErr_SetString(PyExc_RuntimeError,
                        "Communicator_updateFieldValues:"
                        " Could not create result vector.");
        Py_XDECREF(label);
        errorFlag = 1;
        break;
      }

      // store empty list in dict
      errorFlag = PyDict_SetItem(fieldvalues, label, resVal);
      Py_XDECREF(label);
      Py_XDECREF(resVal);
    }
    if (errorFlag) {
      Py_XDECREF(fieldvalues);
      delete elemToNbVals;
      delete[] intPtSectPt;
      delete[] labels;
      delete[] values;
      return NULL;
    }
    debug_print("Communicator_updateFieldValues:"
                " created dict entries.\n");

    if (!periodicIntPtSectPt) {
      elemToNbVals->clear();
    }

    /* store values in dict fieldvalues */
    if (nbComp>1) {  // vector/tensor values
      PyObject* vec;
      PyObject* item;
      for (int i=0; i<nbVals && !errorFlag; ++i) {

        // create the vector/tensor (list) - value
        vec = PyList_New(nbComp);
        if (!vec) {
          PyErr_SetString(PyExc_RuntimeError,
                          "Communicator_updateFieldValues:"
                          " Could not create result vector.");
          errorFlag = 1;
          break;
        }
        for (int j=0; j<nbComp && !errorFlag; ++j) {
          item = PyFloat_FromDouble((double) values[i*nbComp+j]);
          if (!item) {
            errorFlag = 1;
            break;  // exit only the inner for j... loop
          }
          // PyList_SET_ITEM is the version *without* error checking!
          // PyList_SET_ITEM steals a reference => don't DECREF item!
          PyList_SET_ITEM(vec, j, item);
        }
        if (errorFlag) {
          Py_XDECREF(vec);
          break;
        }
    
        // store value in list of integration pt values
        // updateFieldValues_storeItem steals a reference => don't DECREF vec!
        errorFlag = updateFieldValues_storeItem
	  (fieldvalues, vec, i, labels, intPtSectPt, nbIntPtSectPt,
	   periodicIntPtSectPt, elemToNbVals);
      }  // end for i... loop
    }
    
    else {  // nbComp==1, scalar values
      PyObject* item;
      for (int i=0; i<nbVals && !errorFlag; ++i) {

        // create the scalar value to be stored
        item = PyFloat_FromDouble((double) values[i]);
        if (!item) {
          errorFlag = 1;
          break;
        }
    
        // store value in list of integration pt values
        // PyList_SetItem is the version with error checking!
        // PyList_SetItem steals a reference => don't DECREF item!
        errorFlag = updateFieldValues_storeItem
	  (fieldvalues, item, i, labels, intPtSectPt, nbIntPtSectPt,
	   periodicIntPtSectPt, elemToNbVals);
        debug_print("Communicator_updateFieldValues:"
                    " scalar value: %g\n", values[i]);
      }  // end for i... loop
    }
    debug_print("Communicator_updateFieldValues:"
                " position elemIP, stored values in fieldvalues dict"
                " (nbComp=%d), error:%d.\n",
                nbComp, errorFlag);
    delete elemToNbVals;
    elemToNbVals = 0;
  }
  else {
    PyErr_Format(PyExc_ValueError,
                 "Communicator_updateFieldValues: Unrecognized"
                 " position argument <%s>.", position);
    errorFlag = 1;
  }

  debug_print("Communicator_updateFieldValues, ending, error=%d.\n",
              errorFlag);
  Py_XDECREF(fieldvalues);
  delete[] intPtSectPt;
  delete[] labels;
  delete[] values;

  if (errorFlag) {
    return NULL;
  }
  else {
    Py_RETURN_NONE;
  }
}

/* Communicator_updateFieldValuesConstInt
 *
 * Receive integer scalars for elements from the pipe and update the given
 * field values object (a dict).
 *
 *
 * Arguments (function accepts keyword arguments):
 * fieldvalues: a bae.field_01.Field / dict object to take the results.
 */
static PyObject*
Communicator_updateFieldValuesConstInt
(Communicator* self, PyObject *args, PyObject *keywds) {

  size_t rc;
  int* labels = 0;
  int* values = 0;
  debug_print("This is Communicator_updateFieldValuesConstInt...\n");

  /* get and check arguments */
  PyObject* fieldvalues = 0;
  static const char *kwlist[] = {"fieldvalues", NULL};
  if (!PyArg_ParseTupleAndKeywords
      (args, keywds,
       "O:updateFieldValuesConstInt", const_cast<char **>(kwlist),
       &fieldvalues)) {
    return NULL;
  }
  if (!PyDict_Check(fieldvalues)) {
    PyErr_SetString(PyExc_ValueError,
		    "Communicator_updateFieldValuesConstInt: Expected dict"
		    " object as first argument fieldvalues.");
    return NULL;
  }
  Py_XINCREF(fieldvalues);

  /* Read data from the pipe:
   *
   * - nb of field values (=nb of elements)
   * - integer array of element labels
   * - integer array of field values.
   */
  int nbVals;
  if (readOneInt(self->dataChannelIn, &nbVals)) {
    PyErr_SetString(PyExc_RuntimeError,
		    "Communicator_updateFieldValuesConstInt: Couldn't read"
		    " number of field values.");
    Py_XDECREF(fieldvalues);
    return NULL;
  }
  debug_print("Communicator_updateFieldValuesConstInt:"
              " Number of field values: %d.\n", nbVals);

  labels = new int[nbVals];
  if (!labels) {
    PyErr_SetString(PyExc_MemoryError,
                    "Communicator_updateFieldValuesConstInt: Couldn't allocate"
                    " memory for element labels.");
    Py_XDECREF(fieldvalues);
    return NULL;
  }
  rc = read(self->dataChannelIn, labels, nbVals*sizeof(int));
  if (rc != nbVals*sizeof(int)) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Communicator_updateFieldValuesConstInt: Couldn't read"
                    " element labels.");
    Py_XDECREF(fieldvalues);
    delete[] labels;
    return NULL;
  }
  values = new int[nbVals];
  if (!values) {
    PyErr_SetString
      (PyExc_MemoryError,
       "Communicator_updateFieldValuesConstInt: Couldn't allocate"
       " memory for field values.");
    Py_XDECREF(fieldvalues);
    delete[] labels;
    return NULL;
  }
  rc = read(self->dataChannelIn, values, nbVals*sizeof(int));
  if (rc != nbVals*sizeof(int)) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Communicator_updateFieldValuesConstInt: Couldn't read"
                    " field values from the pipe.");
    Py_XDECREF(fieldvalues);
    delete[] labels;
    delete[] values;
    return NULL;
  }
  debug_print("Communicator_updateFieldValuesConstInt:"
              " Finished reading from the pipe: %d labels: %d...,\n"
              "    ... and int values: %d....\n",
              nbVals, labels[0], values[0]);
  ///#############################################################  DEBUG
  // for (int i=0; i<nbVals; ++i) {
  //   debug_print("label %d, value %d\n", labels[i], values[i]);
  // }
  ///#############################################################  DEBUG

  int errorFlag = 0;

  /* store values in dict fieldvalues */
  PyObject* label = 0;
  PyObject* resVal = 0;
  for (int i=0; i<nbVals && !errorFlag; ++i) {
  
    // create the label (dict-key)
    label = PyInt_FromLong((long) labels[i]);
    if (!label) {
      errorFlag = 1;
      break;
    }
  
    // create scalar value
    resVal = PyInt_FromLong((long) values[i]);
    if (!resVal) {
      Py_XDECREF(label);
      errorFlag = 1;
      break;
    }
  
    // store value in dict
    errorFlag = PyDict_SetItem(fieldvalues, label, resVal);
    debug_print("Communicator_updateFieldValuesConstInt:"
                "after setitem, errorFlag=%d\n", errorFlag);
    Py_XDECREF(label);
    Py_XDECREF(resVal);
  }
  if (errorFlag) {
    Py_XDECREF(fieldvalues);
    delete[] labels;
    delete[] values;
    return NULL;
  }

  debug_print("Communicator_updateFieldValuesConstInt, ending, error=%d.\n",
              errorFlag);
  Py_XDECREF(fieldvalues);
  delete[] labels;
  delete[] values;
  Py_RETURN_NONE;
}

/* Communicator_updateStructPtFieldValues
 *
 * Receive field values from the pipe interpolate them to the grid points and
 * updated the field values object.
 *
 *
 * Arguments (function accepts keyword arguments):
 * position: "node", "element" or "elemIP"
 * fieldvalues: a bae.field_01.Field / list object to take the results.
 * notDefinedValue (optional): This value will be returned for points that are
 *    outside the mesh and that are therefore not defined. It must be an array
 *    of length nbComp. Note: For scalar fields it must be an array of length
 *    one!
 * notDefinedKey (optional, default=-1): Use this as key in the dictionary
 *    of field values for the not-defined-value. This must be guaranteed to
 *    never coincide with a real key, i.e. node or element label.
 * interpolationType (optional, default=1): type of interpolation, an int
 *    . 1 a.k.a. "default": depending on the element type, for C3D10M: linear
 *        for integration point values, quadratic or nodal values
 *    . 2 a.k.a. "const": piecewise constant: take the value of the closest
 *        node or integration point.
 */
static PyObject*
Communicator_updateStructPtFieldValues
(Communicator* self, PyObject *args, PyObject *keywds) {

  FieldValueMapBase* fieldValMap = 0;
  size_t rc;
  debug_print("This is Communicator_updateStructPtFieldValues...\n");

  /* get and check arguments */
  const char* position;
  PyObject* fieldvalues;
  float* notDefinedValue = NULL;
  Py_ssize_t notDefinedValueLength = 0;
  int notDefinedKey = -1;
  int interpolationType = 1;
  static const char *kwlist[] = {
    "position", "fieldvalues", "notDefinedValue", "notDefinedKey",
    "interpolationType", NULL};
  if (!PyArg_ParseTupleAndKeywords
      (args, keywds,
       "sO|w#ii:updateStructPtFieldValues", const_cast<char **>(kwlist),
       &position, &fieldvalues, &notDefinedValue, &notDefinedValueLength,
       &notDefinedKey, &interpolationType)) {
    return NULL;
  }
  if (!PyList_Check(fieldvalues)) {
    PyErr_SetString(PyExc_ValueError,
		    "Communicator_updateStructPtFieldValues: Expected list"
		    " object as second argument.");
    return NULL;
  }
  Py_XINCREF(fieldvalues);
  notDefinedValueLength /= sizeof( (*notDefinedValue) );
  debug_print("Communicator_updateStructPtFieldValues; position=%s\n",
              position);

  /* select interpolation type based on position argument
   * Note: The interpolationType argument refers to something else: namely
   * wether we are using interpolation functions of an order as high as
   * possible or if we want piecewise constant interpolation e.g. for the
   * status field.
   * This question is not dealt with here but way further down where the
   * actual interpolation is performed by either the interpolation-method
   * interpolateToPyList or interpolateToPyListPieceWiseConst. */
  InterpolationBase* interpolation;
  if (strcmp(position, "node")==0) {
    debug_print("Communicator_updateStructPtFieldValues, chose"
                " interpolationFromNodal.\n");
    interpolation = self->interpolationFromNodal;
  }
  else if (strcmp(position, "element")==0) {
    debug_print("Communicator_updateStructPtFieldValues, chose"
                " interpolationFromElem.\n");
    interpolation = self->interpolationFromElem;
  }
  else if (strcmp(position, "elemIP")==0) {
    debug_print("Communicator_updateStructPtFieldValues, chose"
                " interpolationFromGaussPts.\n");
    interpolation = self->interpolationFromGaussPts;
  }
  else {
    PyErr_Format(PyExc_ValueError,
                 "Communicator_updateStructPtFieldValues: Unrecognized"
                 " position argument <%s>.", position);
    Py_XDECREF(fieldvalues);
    return NULL;
  }

  /* clear the result-list / make sure it's empty at the beginning */
  Py_ssize_t size = PyList_GET_SIZE(fieldvalues);
  if (PyList_SetSlice(fieldvalues, 0, size, NULL)) {
    PyErr_SetString(PyExc_RuntimeError,
		    "Communicator_updateStructPtFieldValues: Couldn't clear"
		    " the fieldvalues list.");
    Py_XDECREF(fieldvalues);
    return NULL;
  }

  /* Read data from the pipe:
   * (I): integration point / section point numbers
   *
   * - nbIntPtSectPt: number of integration point / section point tuples in the
   *   following array. An integer. E.g. 1 for nodal values and 4 for SDV1 on
   *   C3D10M elements.
   * - (integration point, section point) - tuples: 2*nbIntPtSectPt integers.
   */
  int nbIntPtSectPt;
  if (readOneInt(self->dataChannelIn, &nbIntPtSectPt)) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Communicator_updateStructPtFieldValues: Couldn't read"
                    " nbIntPtSectPt.");
    Py_XDECREF(fieldvalues);
    return NULL;
  }

  int* intPtSectPt = new int[2*nbIntPtSectPt];
  if (!intPtSectPt) {
    PyErr_SetString
      (PyExc_MemoryError,
       "Communicator_updateStructPtFieldValues: Couldn't allocate"
       " memory for list of integration- and section point numbers.");
    Py_XDECREF(fieldvalues);
    return NULL;
  }
  rc = read(self->dataChannelIn, intPtSectPt, 2*nbIntPtSectPt*sizeof(int));
  if (rc != 2*nbIntPtSectPt*sizeof(int)) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Communicator_updateStructPtFieldValues: Couldn't read"
                    " list of integration- and section point numbers.");
    Py_XDECREF(fieldvalues);
    delete[] intPtSectPt;
    return NULL;
  }
  // ///#############################################################  DEBUG
  // for (int i=0; i<nbIntPtSectPt; ++i) {
  //   debug_print("Communicator_updateStructPtFieldValues:"
  //               " Read intPtSectPt from pipe: intpt %d, sectpt %d\n",
  //               intPtSectPt[i*2+0], intPtSectPt[i*2+1]);
  // }
  // ///#############################################################  DEBUG

  /* Read data from the pipe:
   * (II): field values
   *
   * - nb of components of each field value: 1 for scalars, 3 for vectors, 6 for
   *   tensors
   * - nb of field values
   * - integer array of node or element labels
   * - float array of field values.
   */
  int nbComp;
  int nbVals;
  if (readOneInt(self->dataChannelIn, &nbComp)) {
    PyErr_SetString(PyExc_RuntimeError,
		    "Communicator_updateStructPtFieldValues: Couldn't read"
		    " number of components of each field value.");
    Py_XDECREF(fieldvalues);
    delete[] intPtSectPt;
    return NULL;
  }
  if (readOneInt(self->dataChannelIn, &nbVals)) {
    PyErr_SetString(PyExc_RuntimeError,
		    "Communicator_updateStructPtFieldValues: Couldn't read"
		    " number of field values.");
    Py_XDECREF(fieldvalues);
    delete[] intPtSectPt;
    return NULL;
  }
  if (nbComp!=notDefinedValueLength) {
    PyErr_Format(PyExc_RuntimeError,
                 "Communicator_updateStructPtFieldValues: expected %d"
                 " components in notDefinedValue but got %ld.",
                 nbComp, notDefinedValueLength);
    Py_XDECREF(fieldvalues);
    delete[] intPtSectPt;
    return NULL;
  }
  debug_print("Communicator_updateStructPtFieldValues:"
              " Number of components: %d. Number of field values: %d.\n",
              nbComp, nbVals);

  int* labels = new int[nbVals];
  if (!labels) {
    PyErr_SetString(PyExc_MemoryError,
                    "Communicator_updateStructPtFieldValues: Couldn't allocate"
                    " memory for node- or element- labels.");
    Py_XDECREF(fieldvalues);
    delete[] intPtSectPt;
    return NULL;
  }
  rc = read(self->dataChannelIn, labels, nbVals*sizeof(int));
  if (rc != nbVals*sizeof(int)) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Communicator_updateStructPtFieldValues: Couldn't read"
                    " node- or element- labels.");
    Py_XDECREF(fieldvalues);
    delete[] intPtSectPt;
    delete[] labels;
    return NULL;
  }
  float* values = new float[nbVals*nbComp];
  if (!values) {
    PyErr_SetString
      (PyExc_MemoryError,
       "Communicator_updateStructPtFieldValues: Couldn't allocate"
       " memory for field values.");
    Py_XDECREF(fieldvalues);
    delete[] intPtSectPt;
    delete[] labels;
    return NULL;
  }
  rc = read(self->dataChannelIn, values, nbVals*nbComp*sizeof(float));
  if (rc != nbVals*nbComp*sizeof(float)) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Communicator_updateStructPtFieldValues: Couldn't read"
                    " field values from the pipe.");
    Py_XDECREF(fieldvalues);
    delete[] intPtSectPt;
    delete[] labels;
    delete[] values;
    return NULL;
  }
  debug_print("Communicator_updateStructPtFieldValues:"
              " Got %d labels: %d...,\n"
              "    ...   %d value-floats, %d components per value: %g....\n",
              nbVals, labels[0], nbVals*nbComp, nbComp, values[0]);
  // ///#############################################################  DEBUG
  // for (int i=0; i<nbVals; ++i) {
  //   debug_print("label %d, value %g\n", labels[i], values[i]);
  // }
  // ///#############################################################  DEBUG

  /* select field-value-map type and initialize values
   * type and procedure depending on on position argument and nbComp */
  int errorFlag = 0;
  if (strcmp(position, "node")==0 || strcmp(position, "element")==0) {

    /* check Gauss-pt and section pt number */
    if (nbIntPtSectPt!=1) {
      PyErr_Format
        (PyExc_NotImplementedError,
         "Communicator_updateStructPtFieldValues: For nodal and whole"
         " element values expect nbIntPtSectPt==1 but got %d.",
         nbIntPtSectPt);
      errorFlag = 1;
    }
    if (!errorFlag && intPtSectPt[0]!=0) {
      PyErr_Format
        (PyExc_NotImplementedError,
         "Communicator_updateStructPtFieldValues: For nodal and whole"
         " element values expect Gauss pt nb 0 but got %d instead.",
         intPtSectPt[0]);
      errorFlag = 1;
    }
    if (!errorFlag && intPtSectPt[1]>0) {
      PyErr_Format
        (PyExc_NotImplementedError,
         "Communicator_updateStructPtFieldValues: For nodal or whole"
         " element field value got section point nb %d. Not"
         " implemented!", intPtSectPt[1]);
      errorFlag = 1;
    }
    if (errorFlag) {
      Py_XDECREF(fieldvalues);
      delete[] intPtSectPt;
      delete[] labels;
    delete[] values;
    return NULL;
    }
    debug_print("Communicator_updateStructPtFieldValues:"
                " position node or element, checked intPtSectPt: ok.\n");

    /* store values in map */
    if (nbComp>1) {
      fieldValMap = new FieldValueMapVector(nbComp);
      for (int i=0; i<nbVals && !errorFlag; ++i) {
        errorFlag = fieldValMap->store(labels[i], values + i*nbComp);
      }
      fieldValMap->storeNotDefined(notDefinedKey, notDefinedValue);
    }
    else {  // nbComp==1
      fieldValMap = new FieldValueMapVecComp();
      for (int i=0; i<nbVals && !errorFlag; ++i) {
        errorFlag = fieldValMap->store(labels[i], values+i);
      }
      fieldValMap->storeNotDefined(notDefinedKey, notDefinedValue);
    }
    debug_print("Communicator_updateStructPtFieldValues:"
                " position node or element, stored values in fieldValMap"
                " (nbComp=%d), error:%d.\n",
                fieldValMap->getNbComp(), errorFlag);

  }
  else if (strcmp(position, "elemIP")==0) {
    int intPt, sectPt;
    int interpol_nbIntPt = interpolation->m_gaussptsPerPoint;
    if (interpol_nbIntPt<1 && !errorFlag) {
      PyErr_Format(PyExc_ValueError,
                   "Communicator_updateStructPtFieldValues: The number of"
                   " Gauss points must be greater than zero but it's %d.",
                   interpol_nbIntPt);
      errorFlag = 1;
    }
    if (nbIntPtSectPt != interpol_nbIntPt && !errorFlag) {
      PyErr_Format(PyExc_ValueError,
                   "Communicator_updateStructPtFieldValues: We have a"
                   " mismatch between the number of Gauss points according"
                   " to the 'remap'-weights (%d) and the number of Gauss"
                   " points and section points received from the odb (%d)."
                   " The reason might be that section points are not being"
                   " treated properly by this framework, yet.",
                   interpol_nbIntPt, nbIntPtSectPt);
      errorFlag = 1;
    }

    /* check Gauss-pt and section pt number */
    for (int i=0; i<nbIntPtSectPt && !errorFlag; ++i) {
      intPt = intPtSectPt[2*i+0];
      sectPt = intPtSectPt[2*i+1];
      if (intPt != (i%interpol_nbIntPt)+1) {
        PyErr_Format(PyExc_ValueError,
                     "Communicator_updateStructPtFieldValues: For the %dth"
                     " field value expected Gauss pt nb %d but got %d instead.",
                     i+1, (i%interpol_nbIntPt)+1, intPt);
        errorFlag = 1;
        break;
      }
      if (sectPt > 0) {
        PyErr_Format(PyExc_NotImplementedError,
                     "Communicator_updateStructPtFieldValues: For the %dth"
                     " field value got section point nb %d. Not implemented!",
                     i+1, sectPt);
        errorFlag = 1;
        break;
      }
    }
    debug_print("Communicator_updateStructPtFieldValues:"
                " position elemIP, checked intPtSectPt: ok.\n");

    if (!errorFlag) {
      /* store values in map */
      if (nbComp>1) {
        fieldValMap = new FieldValueMapVectorToVec(nbComp);
        for (int i=0; i<nbVals && !errorFlag; ++i) {
          errorFlag = fieldValMap->store(labels[i], values + i*nbComp);
        }
        fieldValMap->storeNotDefined(notDefinedKey, notDefinedValue);
      }
      else {  // nbComp==1
        fieldValMap = new FieldValueMapVecCompToVec();
        for (int i=0; i<nbVals && !errorFlag; ++i) {
          errorFlag = fieldValMap->store(labels[i], values + i);
        }
        fieldValMap->storeNotDefined(notDefinedKey, notDefinedValue);
      }
      // ///########################## DEBUG
      // debug_print("Communicator_updateStructPtFieldValues: notDefinedValue"
      //             " (%ld values):\n   ... [", notDefinedValueLength);
      // for (int i=0; i<notDefinedValueLength; ++i) {
      //   debug_print("%g, ", notDefinedValue[i]);
      // }
      // debug_print("]\n");
      // ///########################## DEBUG
      debug_print("Communicator_updateStructPtFieldValues:"
                  " position elemIP, stored values in fieldValMap"
                  " (nbComp=%d), error:%d.\n",
                  fieldValMap->getNbComp(), errorFlag);
    }
  }
  else {
    PyErr_Format(PyExc_ValueError,
                 "Communicator_updateStructPtFieldValues: Unrecognized"
                 " position argument <%s>.", position);
    errorFlag = 1;
  }
  if (errorFlag) {
    Py_XDECREF(fieldvalues);
    delete fieldValMap;
    delete[] intPtSectPt;
    delete[] labels;
    delete[] values;
    return NULL;
  }
  ///#############################################################  DEBUG
  // 
  // for (std::map<int, float>::iterator it=((FieldValueMapVecComp*) fieldValMap)->begin();
  //      it!=((FieldValueMapVecComp*) fieldValMap)->end(); ++it) {
  //   debug_print("in map: label %d, value %g\n", it->first, it->second);
  // }
  ///#############################################################  DEBUG

  debug_print("Communicator_updateStructPtFieldValues:"
              " starting interpolation.\n");
  if (interpolationType==1) {       // =="default"
    errorFlag = interpolation->interpolateToPyList(fieldValMap, fieldvalues);
  }
  else if (interpolationType==2) {  // =="const"
    errorFlag = interpolation->interpolateToPyListPieceWiseConst
      (fieldValMap, fieldvalues);
  }
  else {
    PyErr_Format(PyExc_ValueError,
                 "Communicator_updateStructPtFieldValues: Unrecognized"
                 " interpolationType argument <%d>.", interpolationType);
    errorFlag = 1;
  }
  
  debug_print("Communicator_updateStructPtFieldValues:"
              " finished interpolating: errorFlag: %d.\n", errorFlag);
  if (errorFlag) {
    Py_XDECREF(fieldvalues);
    delete fieldValMap;
    delete[] intPtSectPt;
    delete[] labels;
    delete[] values;
    return NULL;
  }    

  debug_print("Communicator_updateStructPtFieldValues, ending.\n");
  Py_XDECREF(fieldvalues);
  delete fieldValMap;
  delete[] intPtSectPt;
  delete[] labels;
  delete[] values;
  Py_RETURN_NONE;
}


/* Communicator_updateStructPtFieldValuesConstInt
 *
 * Receive integer scalars for elements from the pipe and assign them to
 * grid points. Updates a given list object as field values result.
 *
 * Note: This is not the same as Communicator_updateStructPtFieldValues()
 * with interpolationType="const". This function is intended for
 * OdbReaderInterpolateToPoints.getMatFieldFromOdb().
 *
 *
 * Arguments (function accepts keyword arguments):
 * result: a bae.field_01.Field / list object to take the results.
 * notDefinedValue (optional, default=-1): An integer value that will be
 *    returned for points that are outside the mesh and that are therefore
 *    not defined.
 * notDefinedKey (optional, default=-1): Use this as key in the dictionary
 *    of field values for the not-defined-value. This must be guaranteed to
 *    never coincide with a real key, i.e. node or element label.
 */
static PyObject*
Communicator_updateStructPtFieldValuesConstInt
(Communicator* self, PyObject *args, PyObject *keywds) {

  size_t rc;
  debug_print("This is Communicator_updateStructPtFieldValuesConstInt...\n");
  
  /* get and check arguments */
  PyObject* result;
  int notDefinedValue = -1;
  int notDefinedKey = -1;
  static const char *kwlist[] = {
    "result", "notDefinedValue", "notDefinedKey", NULL};
  if (!PyArg_ParseTupleAndKeywords
      (args, keywds,
       "O|ii:updateStructPtFieldValuesConstInt", const_cast<char **>(kwlist),
       &result, &notDefinedValue, &notDefinedKey)) {
    return NULL;
  }
  if (!PyList_Check(result)) {
    PyErr_SetString(PyExc_ValueError,
		    "Communicator_updateStructPtFieldValuesConstInt:"
		    " Expected list object as first argument.");
    return NULL;
  }
  Py_XINCREF(result);

  /* clear the result-list / make sure it's empty at the beginning */
  Py_ssize_t size = PyList_GET_SIZE(result);
  if (PyList_SetSlice(result, 0, size, NULL)) {
    PyErr_SetString(PyExc_RuntimeError,
		    "Communicator_updateStructPtFieldValuesConstInt:"
		    " Couldn't clear the result list.");
    Py_XDECREF(result);
    return NULL;
  }

  /* Read data from the pipe:
   *
   * - nb of values
   * - integer array of node or element labels
   * - integer array of values.
   */
  int nbVals;
  if (readOneInt(self->dataChannelIn, &nbVals)) {
    PyErr_SetString(PyExc_RuntimeError,
		    "Communicator_updateStructPtFieldValuesConstInt:"
		    " Couldn't read number of field values.");
    Py_XDECREF(result);
    return NULL;
  }
  debug_print("Communicator_updateStructPtFieldValuesConstInt:"
              " Number of values to read: %d.\n", nbVals);

  int* labels = new int[nbVals];
  if (!labels) {
    PyErr_SetString
      (PyExc_MemoryError,
       "Communicator_updateStructPtFieldValuesConstInt: Couldn't allocate"
       " memory for node- or element- labels.");
    Py_XDECREF(result);
    return NULL;
  }
  rc = read(self->dataChannelIn, labels, nbVals*sizeof(int));
  if (rc != nbVals*sizeof(int)) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Communicator_updateStructPtFieldValuesConstInt:"
                    " Couldn't read node- or element- labels.");
    Py_XDECREF(result);
    delete[] labels;
    return NULL;
  }
  int* values = new int[nbVals];
  if (!values) {
    PyErr_SetString
      (PyExc_MemoryError,
       "Communicator_updateStructPtFieldValuesConstInt: Couldn't allocate"
       " memory for field values.");
    Py_XDECREF(result);
    delete[] labels;
    return NULL;
  }
  rc = read(self->dataChannelIn, values, nbVals*sizeof(int));
  if (rc != nbVals*sizeof(int)) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Communicator_updateStructPtFieldValuesConstInt:"
                    " Couldn't read field values from the pipe.");
    Py_XDECREF(result);
    delete[] labels;
    delete[] values;
    return NULL;
  }
  debug_print("Communicator_updateStructPtFieldValuesConstInt:"
              " Got %d labels: %d..., and values: %d....\n",
              nbVals, labels[0], values[0]);
  // ///#############################################################  DEBUG
  // for (int i=0; i<nbVals; ++i) {
  //   debug_print("label %d, value %d\n", labels[i], values[i]);
  // }
  // ///#############################################################  DEBUG

  /* select field-value-map type and initialize values
   * type and procedure depending on on position argument and nbComp */
  int errorFlag = 0;

  /* store values in map */
  FieldValueMapInt fieldValMap;
  for (int i=0; i<nbVals && !errorFlag; ++i) {
    errorFlag = fieldValMap.store(labels[i], values[i]);
  }
  fieldValMap.storeNotDefined(notDefinedKey, notDefinedValue);
  debug_print("Communicator_updateStructPtFieldValuesConstInt:"
              " stored values in fieldValMap, error:%d.\n", errorFlag);

  if (errorFlag) {
    Py_XDECREF(result);
    delete[] labels;
    delete[] values;
    return NULL;
  }
  ///#############################################################  DEBUG
  // 
  // for (std::map<int, int>::iterator it=fieldValMap.begin();
  //      it!=fieldValMap.end(); ++it) {
  //   debug_print("in map: label %d, value %d\n", it->first, it->second);
  // }
  ///#############################################################  DEBUG

  debug_print("Communicator_updateStructPtFieldValuesConstInt:"
              " starting interpolation.\n");
  errorFlag = self->interpolationFromElem->interpolateToPyList(&fieldValMap, result);
  debug_print("Communicator_updateStructPtFieldValuesConstInt:"
              " finished interpolating: errorFlag: %d.\n", errorFlag);
  if (errorFlag) {
    Py_XDECREF(result);
    delete[] labels;
    delete[] values;
    return NULL;
  }    

  debug_print("Communicator_updateStructPtFieldValuesConstInt, ending.\n");
  Py_XDECREF(result);
  delete[] labels;
  delete[] values;
  Py_RETURN_NONE;
}


/******************************************************************************
*******************************************************************************
***  function for initialize, open, close and such...                       ***
*******************************************************************************
******************************************************************************/


/* Communicator destructor
 *
 * Finalizes the odb-API if there is no other Communicator object present
 * anymore. And frees self.
 */
static void
Communicator_dealloc(Communicator* self)
{
  debug_print("In dealloc (1).\n");

  // free interpolation objects
  if (self->interpolationFromNodal) {
    delete self->interpolationFromNodal;
    self->interpolationFromNodal = NULL;
  }
  if (self->interpolationFromGaussPts) {
    delete self->interpolationFromGaussPts;
    self->interpolationFromGaussPts = NULL;
  }
  if (self->interpolationFromElem) {
    delete self->interpolationFromElem;
    self->interpolationFromElem = NULL;
  }

  // decrease python reference of raw data for weights (PyString)
  Py_XDECREF(self->pyStrNodalweights);
  Py_XDECREF(self->pyStrGaussPtsweights);

  // free self-object
  self->ob_type->tp_free((PyObject*)self);
  debug_print("In dealloc (6).\n");
}


/* Communicator constructor: __new__()
 *
 * Initializes the odb-API if that's not been done before.
 * Increases nbActiveOdbs.
 */
static PyObject *
Communicator_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  Communicator *self;

  self = (Communicator *)type->tp_alloc(type, 0);
  if (self != NULL) {

    // data channels not specified yet, will be done by open()
    self->dataChannelIn = 0;
    self->dataChannelOut = 0;
    
  }
  return (PyObject *)self;
}


/* Communicator constructor: __init__()
 *
 * Opens the odb.
 *
 * @param odb: file name of the odb.
 * @param dataChannelIn: file descriptor (int) of pipe for data to the worker
 *                       (this process)
 * @param dataChannelOut: file descriptor (int) of pipe for data output from
 *          the worker to the proxy (or its caller) in the parent process
 */
static int
Communicator_init(Communicator *self, PyObject *args, PyObject *kwds)
{
  static char const *kwlist[] = {"dataChannelIn", "dataChannelOut", NULL};
  if (! PyArg_ParseTupleAndKeywords
      ( args, kwds, "ii", const_cast<char **>(kwlist),
        &(self->dataChannelIn), &(self->dataChannelOut)))
    return -1;
  return 0;
}


/******************************************************************************
*******************************************************************************
***  Python API stuff                                                       ***
*******************************************************************************
******************************************************************************/


/* Communicator methods
 */
static PyMethodDef Communicator_methods[] = {
  {"updateNodeCoordsFromPipe",
   (PyCFunction)Communicator_updateNodeCoordsFromPipe, METH_VARARGS,
   "Read node definitions from the pipe and store in the supplied dict."
   " On the other side of the pipe it's OdbReader.getNodeCoordsFromOdb"
   " writing the data to the pipe."
  },
  {"updateElNodesFromPipe",
   (PyCFunction)Communicator_updateElNodesFromPipe, METH_VARARGS,
   "Read element definitions from the odb and store in the supplied dicts"
   " elNodes, elType and typeEl."
  },
  {"storeRemapWeights",
   (PyCFunction)Communicator_storeRemapWeights, METH_VARARGS,
   "Store the weighting factors for nodal and Gauss-pt interpolation,"
   " create InterpolationFromNodal object self->interpolationFromNodal"
   " and InterpolationFromGaussPts object self->interpolationFromGaussPts.\n"
   "Call this function to specify the points at which you want to get field"
   " value output by subsequent calls to getFieldDataForPoints.\n\n"
   "Arguments:\n"
   "nbPoints: number of points for which field output is being queried."
   " Determines the length of the following two arrays.\n"
   "nodalWeights: A string that will be interpreted as a binary array of node"
   " numbers and associated weight factors for each of the output points."
   " See bae.odb_03.OdbReaderInterpolateToPoints.initOutputPoints() and"
   " bae.mesh_01.InterpolMeshToPoints.writeRemapParam for further details.\n"
   "gaussPointsWeights: A string that will be interpreted as a binary array of"
   " element numbers and associated weight factors for each of the output"
   " points."
  },
  {"updateFieldValues",
   (PyCFunction)Communicator_updateFieldValues,
   METH_VARARGS | METH_KEYWORDS,
   "Receive field values from the pipe and update the supplied field values"
   " L{bae.field_01.Field}-dict-object.\n\n"
   "Arguments:\n"
   "position: 'node', 'element' or 'elemIP'\n"
   "fieldvalues: a L{bae.field_01.Field} / dict object to take the results."
  },
  {"updateFieldValuesConstInt",
   (PyCFunction)Communicator_updateFieldValuesConstInt,
   METH_VARARGS | METH_KEYWORDS,
   "Receive integer values for elements from the pipe and update the supplied"
   " field values L{bae.field_01.Field}-dict-object.\n\n"
   "Arguments:\n"
   "fieldvalues: a L{bae.field_01.Field} / dict object to take the results."
  },
  {"updateStructPtFieldValues",
   (PyCFunction)Communicator_updateStructPtFieldValues,
   METH_VARARGS | METH_KEYWORDS,
   "Receive field values from the pipe, interpolate them to the grid points and"
   " update the supplied field values L{bae.field_01.Field}-list-object.\n\n"
   "Arguments:\n"
   "position: 'node', 'element' or 'elemIP'\n"
   "fieldvalues: a L{bae.field_01.Field} / list object to take the results.\n"
   "notDefinedValue (optional): This value will be returned for points that are"
   " outside the mesh and that are therefore not defined. It must be an array"
   " of length nbComp. Note: For scalar fields it must be an array of length"
   " one!\n"
   "notDefinedKey (optional, default=-1): Use this as key in the dictionary"
   " of field values for the not-defined-value. This must be guaranteed to"
   " never coincide with a real key, i.e. node or element label."
  },
  {"updateStructPtFieldValuesConstInt",
   (PyCFunction)Communicator_updateStructPtFieldValuesConstInt,
   METH_VARARGS | METH_KEYWORDS,
   "Receive integer values for elements from the pipe, assign them to the"
   " grid points and updated the supplied field values"
   " L{bae.field_01.Field}-list-object.\n\n"
   "Arguments:\n"
   "result: a bae.field_01.Field / list object to take the results.\n"
   "notDefinedValue (optional, default=-1): An integer value that will be"
   " returned for points that are outside the mesh and that are therefore"
   " not defined.\n"
   "notDefinedKey (optional, default=-1): Use this as key in the dictionary"
   " of field values for the not-defined-value. This must be guaranteed to"
   " never coincide with a real key, i.e. node or element label."
  },
  {NULL}  /* Sentinel */
};


/* Communicator type description
 */
static PyTypeObject CommunicatorType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "odbReader_ext.Communicator", /*tp_name*/
    sizeof(Communicator),         /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)Communicator_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Communicator class"          /* tp_doc */
    "\n"
    "\nCommunicator(dataChannelIn, dataChannelOut)"
    "\n"
    "\nStores the data channels (pipe connections) for I/O."
    "\n"
    "\n@param dataChannelIn: file descriptor (int) of pipe for data input from"
    "\n         the worker in the child process running the Abaqus odb-API."
    "\n@param dataChannelOut: file descriptor (int) of pipe for data output"
    "\n         from the proxy (or its caller) to the worker.",
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    Communicator_methods,      /* tp_methods */
    Communicator_members,      /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Communicator_init,  /* tp_init */
    0,                         /* tp_alloc */
    Communicator_new,                 /* tp_new */
};


/*   Python module stuff
 */
static PyMethodDef ModuleMethods[] = {
  /* ... */
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

/* Note: PYEXTMODULEINITFUNC will be defined in setup.py
 * ... it's something like initodbReader_ext_abq613
 * the same applies to PYEXTMODULENAME, it's for example
 * "odbReader_ext_abq613" (including the quotes)
 */
PyMODINIT_FUNC
initodbReader_ext(void)
{
  debug_print("Initializing odbReader_ext module...\n");

  // Initialize the module
  PyObject* m;
  m = Py_InitModule3("odbReader_ext", ModuleMethods,
  		     "Extension module for the (binary-) Communicator class.");
  if (m == NULL)
    return;

  // add version string to module 
  if (PyModule_AddStringConstant(m, "__version__", module_version) != 0)
    return;

  // add Communicator class to module
  if (PyType_Ready(&CommunicatorType) < 0)
    return;
  Py_INCREF(&CommunicatorType);
  PyModule_AddObject(m, "Communicator", (PyObject *)&CommunicatorType);

}
