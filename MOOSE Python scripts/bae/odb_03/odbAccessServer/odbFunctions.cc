/*****************************************************************************
 * server functions to access the odb


Server functions implementing possible requests to the odbAccessServer stick
to the following conventions:
* The functions take no arguments.
* They return an int return code:
  . 0 on success,
  . 1 if they raise a recoverable exception,
  . 2 in case of an unrecoverable error
* When adding a new server function log an entry in funcMap by extending the
  function initFuncMap() at the very end of this file odbFunctions.cc.
* Arguments for the function as well as return values are being passed
  through the pipe. Each server function takes care of this communication
  on it's own. Various variants of pipeReadXXX() and pipeWrite() functions
  are provided for this purpose (and defined in communications.cc).
  Typically a server function first reads all arguments from the pipe at
  the very beginning of the function. Plausibility checks are highly
  recommended to follow immediately.
* Raising a recoverable exception will cause this exception to be raised
  in the calling python script. The odbAccessServer stays in a fully
  functional state and can serve further requests. The corresponding code
  in the server function looks something like this:
   >>>   if (! fieldRepos.isMember(fieldName.c_str()) ) {
   >>>     pipeWrite("");  // pipe dummy return values according to the
   >>>     pipeWrite("");  // function interface. Here: position, dataType
   >>>     pipeException
   >>>       ("ValueError",
   >>>        "getFieldProps: Field %s does not exist in step %d frame %d.",
   >>>        fieldName.c_str() , currentStep, currentFrame);
   >>>     return 1;
   >>>   }
  Note the dummy return values to avoid obfuscating the pipe. Then call
  pipeException() and return 1.
* To signal an unrecoverable error simply call error("Errormessage.")
  This will stop the odbAccessServer.
*****************************************************************************/

#include <map>
#include <set>
#include <algorithm>

#include <odb_MaterialTypes.h>
#include <odb_SectionTypes.h>

/******************************************************************************
*******************************************************************************
*** global variables                                                        ***
*******************************************************************************
******************************************************************************/

odb_Odb* ptOdb = NULL;

odb_Instance* ptInstance = NULL;
int currentStep = -1;
int currentFrame = -1;
odb_Step* ptOdbStep = NULL;
odb_Frame* ptOdbFrame = NULL;

/* filterNsetName, filterElsetName: keys for filter sets in the
 * instance.nodeSets() and instance.elementSets() repositories. */
odb_String* filterNsetName = NULL;
odb_String* filterElsetName = NULL;

/******************************************************************************
 * type GetInvariantFuncPtr
 * and lookup table invariantIntToFunc: lookup table for functions
 *
 * Used by server function getFieldData.
 * Must be consistent with OdbReaderBase._odbInvariantNameToNb defined in
 * odbReaderBase.py.
 */
typedef float (odb_FieldValue:: *GetInvariantFuncPtr)() const;
/* array size is arbitrary, we'll not have more than 20 */
const static GetInvariantFuncPtr invariantIntToFunc[] = {
  NULL,
  &odb_FieldValue::magnitude,
  &odb_FieldValue::mises,
  &odb_FieldValue::tresca,
  &odb_FieldValue::press,
  &odb_FieldValue::inv3,
  &odb_FieldValue::maxPrincipal,
  &odb_FieldValue::midPrincipal,
  &odb_FieldValue::minPrincipal,
  &odb_FieldValue::maxInPlanePrincipal,
  &odb_FieldValue::minInPlanePrincipal,
  &odb_FieldValue::outOfPlanePrincipal,
};

/******************************************************************************
*******************************************************************************
*** server functions                                                        ***
*******************************************************************************
******************************************************************************/

/* The following server functions all are of type void func().
 *
 * Input arguments are listed in the correct order in the header after "Input:"
 * Those input arguments are read from the ctrl-pipe (stdin, fd 0).
 *
 * Return values are listed in the correct order in the header after "Output:"
 * Those output items are written to the ctrl-pipe (stdout, fd 1).
 */


/* Test function
 *
 * Call this test with the following command from within a method of
 * OdbReaderInterpolateToPoints:
 * f,s = self._callmethod("test", 3, 3.14, "teststring", resultTypeCodes="fs")
 *
 * input:
 *   int, float, string
 * output:
 *   float, string
 */
int test() {
  debug_print("Hi this is the test function.\n");
  int valInt = pipeReadInt();
  float valFloat = pipeReadFloat();
  std::string valString = pipeReadString();
  debug_print("test function got args: %d, %f, %s\n",
              valInt, valFloat, valString.c_str());

  valFloat = valInt*0.01;
  valString = "ResulTsTr";
  debug_print("test function writing: %f, %s\n", valFloat, valString.c_str());
  pipeWrite(valFloat);
  pipeWrite(valString);
  debug_print("test function fini.\n");
  return 0;
}


/* init function
 *
 * To initialize the odb functionality of the odbAccessServer.
 *
 * Why don't we do it on odbAccessServer-startup? It's more a matter of taste:
 * We've got all odb related stuff in this odbFunctions.cc.
 * And why not in openOdb? Good point, again a matter of taste: Maybe
 * eventually we have a reason to separate init from opening the odb, dunno.
 * This function also serves as a means to test communication to the
 * odbAccessServer right after it has been started.
 */
int init() {
  filterNsetName = new odb_String();
  filterElsetName = new odb_String();
  return 0;
}


/* open the odb
 *
 * input:
 * - odbPath: path of the odb
 */
int closeOdb();  // declare to be used in openOdb
int openOdb() {
  std::string odbPath = pipeReadString();

  // Open odb read-only:
  static const odb_String name = "";
  static const odb_String odb_path = odbPath.c_str();
  debug_print("Opening odb %s ...\n", odbPath.c_str());
  closeOdb();
  ptOdb = &openOdb(name, odb_path, true);
  debug_print("Opening odb done. ptOdb: %p\n", (void*) ptOdb);
  return 0;
}


/* close the odb
 */
int closeOdb() {

  delete ptInstance; ptInstance = NULL;

  currentFrame = -1;
  delete ptOdbFrame; ptOdbFrame = NULL;
  currentStep = -1;
  delete ptOdbStep; ptOdbStep = NULL;

  if (ptOdb) {
    debug_print("In closeOdb. ptOdb: %p\n", (void*) ptOdb);
    ptOdb->close();
    ptOdb = NULL;
  }
  return 0;
}

/******************************************************************************
*** read model data (mesh, elsets...)                                       ***
******************************************************************************/

/* setPtInstance
 *
 * initialize the global attribute ptInstance
 */
void setPtInstance() {
  debug_print("This is setPtInstance...\n");

  /* get first part instance of the root assembly */
  odb_InstanceRepositoryIT instIter( ptOdb->rootAssembly().instances() );
  instIter.first();
  // the following does not work, because iter.currentValue() yields const obj
  // need to access the repository via key to get the non-const instance obj
  // ptInstance = &(instIter.currentValue());
  delete ptInstance;
  ptInstance = new odb_Instance
    (ptOdb->rootAssembly().instances()[instIter.currentKey()]);

  debug_print("Found %s as part instance\n", instIter.currentKey().CStr());
}


/* getNodeCoordsFromOdb
 *
 * get nodes from the odb and send the raw data through the pipe
 */
int getNodeCoordsFromOdb() {
  debug_print("This is getNodeCoordsFromOdb...\n");

  if (ptInstance == NULL) {
    setPtInstance();
  }

  /* determine dimensionality of the model */
  int dim;
  switch (ptInstance->embeddedSpace()) {
  case (odb_Enum::THREE_D):
    dim = 3;
    break;
  case odb_Enum::TWO_D_PLANAR:
  case odb_Enum::AXISYMMETRIC:
    dim = 2;
    break;
  case odb_Enum::UNKNOWN_DIMENSION:
  default:
    error("updateNodeCoordsFromOdb: Couldn't determine dimensionality.");
    return 2;  // can't reach this point
  }

  odb_SequenceNode& nodes( ptInstance->nodes() );

  int size = nodes.size();

  // write nb of nodes to be transfered through the pipe
  debug_print("getNodeCoordsFromOdb: dim=%d, size=%d, sending to pipe\n",
              dim, size);
  pipeWrite(dim);
  pipeWrite(size);
  debug_print("getNodeCoordsFromOdb: sent dim and size to pipe\n");

  // now transfer node data through the pipe
  int label;
  for (int i=0; i<size; ++i) {
    label = nodes[i].label();
    pipeWrite(label);
    pipeWrite(nodes[i].coordinates(), dim);
  }
  debug_print("getNodeCoordsFromOdb: sent node data to pipe, ending.\n");
  return 0;
}


/* getElementsFromOdb
 *
 * get elements from the odb and send the raw data through the pipe.
 */
inline void storeInt(char*& dataPtr, const int data) {
  *((int*) dataPtr) = data;
  dataPtr += sizeof(int);
}

int getElementsFromOdb() {
  debug_print("This is getElementsFromOdb...\n");

  if (ptInstance == NULL) {
    setPtInstance();
  }

  odb_SequenceElement& elements( ptInstance->elements() );

  // write nb of elements to be transferred through the pipe
  int size = elements.size();
  pipeWrite(size);

  // data buffer for one element
  const int maxSizeElemData = 1000;
  char elemData[maxSizeElemData];
  char *dataPtr;
  int sizeElType;

  for (int i=0; i<size; ++i) {

    dataPtr = elemData;

    /* store element label in elemData */
    storeInt(dataPtr, elements[i].label());

    /* list of node numbers */
    const odb_SequenceInt& nodes = elements[i].connectivity();

    /* ... store in elemData */
    storeInt(dataPtr, nodes.size());
    for (int j=0; j<nodes.size(); ++j) {
      storeInt(dataPtr, nodes[j]);
    }

    /* store element type in elemData */
    sizeElType = elements[i].type().length();
    storeInt(dataPtr, sizeElType);
    strncpy(dataPtr, elements[i].type().cStr(), sizeElType);

    /* write elemData to the pipe */
    pipeWrite(elemData, (dataPtr-elemData)+sizeElType);

  }

  debug_print("getElementsFromOdb, ending.\n");

  return 0;
}


/* getElsetsFromOdb
 *
 * get elsets from the odb and send the raw data through the pipe.
 *
 * Note: elsets are being read from the first part instance and from the
 * root assembly as well. If the same elset name occurs twice the elset
 * in the root assembly is ignored.
 */
void pipeWriteOdbElemSet(const odb_Set& odbSet) {
  //
  pipeWrite((std::string) odbSet.name().cStr());

  const odb_SequenceString& instanceNames = odbSet.instanceNames();
  int nbInstances = instanceNames.size();
  const odb_SequenceElement& odbSeqElems =
    (nbInstances ?
     odbSet.elements(instanceNames[0]) :
     odbSet.elements() );
  int nbElems = odbSeqElems.size();
  int* elems = new int[nbElems];
  for (int i=0; i<nbElems; ++i) {
    elems[i] = odbSeqElems[i].label();
  }
  pipeWrite(nbElems);  // nb of elements
  pipeWrite(elems, nbElems);  // array of element labels
  delete[] elems;
}

int getElsetsFromOdb() {
  debug_print("This is getElsetsFromOdb...\n");

  if (ptInstance == NULL) {
    setPtInstance();
  }

  std::set<std::string> elsetsInInstance;
  {
    const odb_SetRepository& odbSetRepository = ptInstance->elementSets();
    debug_print("getElsetsFromOdb: Found %d elsets in the first part instance\n",
                odbSetRepository.size());
    odb_SetRepositoryIT elsetIter(odbSetRepository);
    for (elsetIter.first(); !elsetIter.isDone(); elsetIter.next()) {
      const odb_Set& odbSet = elsetIter.currentValue();
      elsetsInInstance.insert( odbSet.name().cStr() );
      pipeWriteOdbElemSet(odbSet);
    }
    debug_print("getElsetsFromOdb: Sent %d elsets from the first part"
                " instance\n", elsetsInInstance.size());
  }

  {
    const odb_SetRepository& odbSetRepository
      = ptOdb->rootAssembly().elementSets();
    debug_print("getElsetsFromOdb: Found %d elsets in the root assembly\n",
                odbSetRepository.size());
    odb_SetRepositoryIT elsetIter(odbSetRepository);
    for (elsetIter.first(); !elsetIter.isDone(); elsetIter.next()) {
      const odb_Set& odbSet = elsetIter.currentValue();
      if (elsetsInInstance.count( odbSet.name().cStr() )) {
        // ignore this set if we already found it earlier
        // (std::set::count(key) returns number of occurences of key in set.
        //  might be 0 or 1)
        continue;
      }
      pipeWriteOdbElemSet(odbSet);
    }
    debug_print("getElsetsFromOdb: Sent elsets from the root assembly.\n");
  }

  pipeWrite(0);  // send end marker to the pipe

  debug_print("getElsetsFromOdb, ending.\n");
  return 0;
}


/******************************************************************************
***  functions to specify the odb frame for subsequent queries              ***
******************************************************************************/

/* setFrame
 *
 * set (possibly change) current frame
 * input arguments: stepNb, frameNb (2x int)
 *
 * This needs to be called before calling any function that gets field data
 * from the odb like getFieldProps and getFieldData.
 */
int setFrame() {
  debug_print("This is setFrame...\n");

  /* get and check arguments */
  int stepNb = pipeReadInt();
  int frameNb = pipeReadInt();
  int framesCnt;

  /* update step */
  if (stepNb != currentStep) {
    debug_print("setFrame: changing from step %d to step %d...\n",
        	currentStep, stepNb);
    odb_StepRepositoryIT stepIter( ptOdb->steps() );
    stepIter.first();
    for (int i=1; i<stepNb && !stepIter.isDone(); ++i) {
      stepIter.next();
    }
    if (stepIter.isDone()) {
      currentStep = -1;
      delete ptOdbFrame; ptOdbStep = NULL;
      currentFrame = -1;
      delete ptOdbFrame; ptOdbFrame = NULL;
      pipeException
        ("ValueError",
         "setFrame: There is no step %d in the odb, last step is"
         " %d.", stepNb, ptOdb->steps().size());
      return 1;
    }
    else {
      currentStep = stepNb;
      delete ptOdbStep;
      ptOdbStep = new odb_Step(stepIter.currentValue());
      currentFrame = -1;
      delete ptOdbFrame; ptOdbFrame = NULL;
    }
  }

  /* update frame */
  if (ptOdbStep) {
    framesCnt = ptOdbStep->frames().size();
    if (framesCnt < 1) {
      currentFrame = -1;
      delete ptOdbFrame; ptOdbFrame = NULL;
      pipeException
        ("ValueError",
         "setFrame: There is no frame in step %d.", stepNb);
      return 1;
    }

    /* negative frame numbers count from the end */
    while (frameNb < 0) {
      frameNb += framesCnt;
    }
  }

  if (ptOdbStep && framesCnt > 0 && frameNb != currentFrame) {
    debug_print("setFrame: changing from frame %d to frame %d...\n",
        	currentFrame, frameNb);

    /* check that frequested frame actually exists */
    if (frameNb >= framesCnt) {
      currentFrame = -1;
      delete ptOdbFrame; ptOdbFrame = NULL;
      pipeException
        ("ValueError",
         "setFrame: There is no frame %d in step %d, last frame is"
         " %d.", frameNb, currentStep, framesCnt-1);
      return 1;
    }
    else if (frameNb != currentFrame) {
      currentFrame = frameNb;
      delete ptOdbFrame;
      ptOdbFrame = new odb_Frame(ptOdbStep->frames(frameNb));
    }
  }
  debug_print("setFrame, ending.\n");
  return 0;
}


/******************************************************************************
***  functions to set a filter nset and a filter elset                      ***
******************************************************************************/

/* createFilterNset
 *
 * create a new filter nset from node numbers
 * input arguments (to be passed first):
 *  1. nb of nodes (int)
 *  2. node numbers: int array
 *
 * output argument (to be read after passing the input args):
 *  - nset name
 */
int createFilterNset() {
  debug_print("This is createFilterNset...\n");

  /* read number of nodes from the pipe */
  int N = pipeReadInt();
  debug_print("createFilterNset: Will read %d node labels from the pipe.\n", N);

  /* no nodes means: disable filter */
  if (N==0) {
    debug_print("createFilterNset: disabling filter.\n");
    *(filterNsetName) = "";
    pipeWrite((int) 0);  /* write filterNsetName: no filter, no name */
    debug_print("createFilterNset ending.\n");
    return 0;
  }

  if (ptInstance == NULL) {
    setPtInstance();
  }

  /* create filter set name */
  int i = 0;
  do {
    i++;
    *(filterNsetName) = "filterNset_";
    filterNsetName->append(odb_String(i));
  } while (ptInstance->nodeSets().isMember(*(filterNsetName)));

  /* read arguments from pipe: node labels */
  odb_SequenceInt listOfNodes(N);
  for (i=0; i<N; ++i) {
    listOfNodes.append(pipeReadInt());
  }
  debug_print("createFilterNset: got %d node labels from the pipe.\n", N);

  /* create filter nset */
  ptInstance->NodeSet(*(filterNsetName), listOfNodes);
  debug_print("createFilterNset: created node set %s with %d nodes.\n",
              filterNsetName->cStr(), N);

  /* write nset name to pipe */
  pipeWrite((int) filterNsetName->length());
  pipeWrite(filterNsetName->cStr(), filterNsetName->length());
  debug_print("createFilterNset ending.\n");

  return 0;
}

/* setFilterNset
 *
 * activate filter nset already defined in the odb
 * input argument:
 *  - nset name
 *
 * output argument (to be read after passing the input args):
 *  - nb of nodes (int)
 */
int setFilterNset() {
  debug_print("This is setFilterNset...\n");

  /* get and check argument */
  std::string nsetName = pipeReadString();

  /* empty string means: disable filter */
  if (0 == nsetName.size()) {
    debug_print("setFilterNset: disabling filter.\n");
    *(filterNsetName) = "";
    pipeWrite((int) 0);  /* no filter, nb of nodes=0 */
    debug_print("setFilterNset ending.\n");
    return 0;
  }

  /* store nset name */
  *(filterNsetName) = nsetName.c_str();

  /* check if this nset exists... */
  if (ptInstance == NULL) {
    setPtInstance();
  }
  if (!ptInstance->nodeSets().isMember(*(filterNsetName))) {
    debug_print("setFilterNset: nset not found, disabling filter.\n");
    *(filterNsetName) = "";
    pipeWrite((int) -1);  /* nset not found, no filter */
    debug_print("setFilterNset ending.\n");
    return 0;
  }

  /* determine number of nodes */
  int nbItems = ptInstance->nodeSets().get(*(filterNsetName)).nodes().size();
  pipeWrite(nbItems);  /* nb of nodes */
  debug_print("setFilterNset ending.\n");
  return 0;
}

/* createFilterElset
 *
 * create a new filter elset from element numbers
 * input arguments (to be passed first):
 *  1. nb of elements (int)
 *  2. element numbers (labels): int array
 *
 * output argument (to be read after passing the input args):
 *  - elset name
 */
int createFilterElset() {
  debug_print("This is createFilterElset...\n");

  /* read number of elements from the pipe */
  int N = pipeReadInt();
  debug_print("createFilterElset:"
              " Will read %d element labels from the pipe.\n", N);

  /* no elements means: disable filter */
  if (N==0) {
    debug_print("createFilterElset: disabling filter.\n");
    *(filterElsetName) = "";
    pipeWrite((int) 0);  /* write filterElsetName: no filter, no name */
    debug_print("createFilterElset ending.\n");
    return 0;
  }

  if (ptInstance == NULL) {
    setPtInstance();
  }

  /* create filter set name */
  int i = 0;
  do {
    i++;
    *(filterElsetName) = "filterElset_";
    filterElsetName->append(odb_String(i));
  } while (ptInstance->elementSets().isMember(*(filterElsetName)));

  /* read arguments from pipe: element labels from pipe */
  odb_SequenceInt listOfElems(N);
  for (i=0; i<N; ++i) {
    listOfElems.append(pipeReadInt());
  }
  debug_print("createFilterElset: got %d element labels from the pipe.\n", N);

  /* create filter elset */
  ptInstance->ElementSet(*(filterElsetName), listOfElems);
  debug_print("createFilterElset: created element set %s with %d elements.\n",
              filterElsetName->cStr(), N);

  /* write elset name to pipe */
  pipeWrite((int) filterElsetName->length());
  pipeWrite(filterElsetName->cStr(), filterElsetName->length());
  debug_print("createFilterElset ending.\n");

  return 0;
}

/* setFilterElset
 *
 * activate filter elset already defined in the odb
 * input argument:
 *  - elset name
 *
 * output argument (to be read after passing the input args):
 *  - nb of elements (int)
 */
int setFilterElset() {
  debug_print("This is setFilterElset...\n");

  /* get and check argument */
  std::string elsetName = pipeReadString();

  /* empty string means: disable filter */
  if (0 == elsetName.size()) {
    debug_print("setFilterElset: disabling filter.\n");
    *(filterElsetName) = "";
    pipeWrite((int) 0);  /* no filter, nb of elements=0 */
    debug_print("setFilterElset ending.\n");
    return 0;
  }

  /* store elset name */
  *(filterElsetName) = elsetName.c_str();

  /* check if this elset exists... */
  if (ptInstance == NULL) {
    setPtInstance();
  }
  if (!ptInstance->elementSets().isMember(*(filterElsetName))) {
    debug_print("setFilterElset: elset not found, disabling filter.\n");
    *(filterElsetName) = "";
    pipeWrite((int) -1);  /* elset not found, no filter */
    debug_print("setFilterElset ending.\n");
    return 0;
  }

  /* determine number of elements */
  int nbItems = ptInstance->elementSets()[*(filterElsetName)].elements().size();
  pipeWrite(nbItems);  /* nb of elements */
  debug_print("setFilterElset ending.\n");
  return 0;
}

/******************************************************************************
***  functions for retrieving field data from the odb                       ***
******************************************************************************/

/* getFieldProps
 *
 * get field properties from fieldName
 * input argument: fieldName as in the odb (string) E.g. "S", "U", ...
 *
 * output arguments: position, dataType
 * where position can be "node", "elemIP" or "element"
 * and dataType can be "scalar", "vector" or "tensor"
 */
int getFieldProps() {
  debug_print("This is getFieldProps...\n");

  /* get and check arguments */
  std::string fieldName = pipeReadString();

  if (!ptOdbFrame) {
    pipeWrite("");  /* position */
    pipeWrite("");  /* dataType */
    error("getFieldProps: Frame not set.\n");
  }

  /* get odb_FieldOutput object from fieldName */
  odb_FieldOutputRepository& fieldRepos = ptOdbFrame->fieldOutputs();
  if (! fieldRepos.isMember(fieldName.c_str()) ) {
    pipeWrite("");  /* position */
    pipeWrite("");  /* dataType */
    pipeException
      ("ValueError",
       "getFieldProps: Field %s does not exist in step %d frame %d.",
       fieldName.c_str() , currentStep, currentFrame);
    return 1;
  }
  const odb_FieldOutput& odbField = fieldRepos[fieldName.c_str()];

  /* get position, precision, type from the odb_FieldValue object */
  odb_Enum::odb_ResultPositionEnum odbFldPos;
  odb_Enum::odb_PrecisionEnum odbFldPrec;
  odb_Enum::odb_DataTypeEnum odbFldType;

  const odb_FieldValue& fieldValue = odbField.values(0);
  odbFldPos = fieldValue.position();
  odbFldPrec = fieldValue.precision();
  odbFldType = fieldValue.type();

  /* get position from odb */

  /* position: node, elemIP, element */
  switch (odbFldPos) {
  case odb_Enum::NODAL:
    pipeWrite("node");
    break;
  case odb_Enum::INTEGRATION_POINT:
    pipeWrite("elemIP");
    break;
  case odb_Enum::CENTROID:
  case odb_Enum::WHOLE_ELEMENT:
    pipeWrite("element");
    break;
  default:
    pipeWrite("");  /* position */
    pipeWrite("");  /* dataType */
    error("getFieldProps: position %d for field %s not"
          " implemented. Search for odb_ResultPositionEnum in"
          " /usr/abaqus/6.13-2/code/include/odb_Enum.h",
          odbFldPos, fieldName.c_str());
  }

  /* check precision */
  /* precision */
  // /*** pointer to struct member   ***/
  // typedef void*(odb_FieldValue:: *FieldDataPtr);
  // FieldDataPtr* fieldDataPtr
  // *fieldDataPtr = &odb_FieldValue::data;
  // *fieldDataPtr = &odb_FieldValue::dataDouble;
  // data = (double) odbField.*fieldDataPtr();
  // switch (fieldValue.precision()) {
  // case odb_Enum::SINGLE_PRECISION:

  // case odb_Enum::DOUBLE_PRECISION:
  // }
  if (odbFldPrec != odb_Enum::SINGLE_PRECISION) {
    pipeWrite("");  /* dataType */
    error("getFieldProps: double precision not implemented.");
  }

  /* get dataType from odb */
  switch (odbFldType) {
  case odb_Enum::SCALAR:
    pipeWrite("scalar");
    break;
  case odb_Enum::VECTOR:
    pipeWrite("vector");
    break;
  case odb_Enum::TENSOR_3D_FULL:
  case odb_Enum::TENSOR_3D_PLANAR:
  case odb_Enum::TENSOR_3D_SURFACE:
    pipeWrite("tensor");
    break;
  default:
    pipeWrite("");  /* dataType */
    error("getFieldProps: dataType %d for field %s not"
          " implemented. Search for odb_DataTypeEnum in"
          " /usr/abaqus/6.13-2/code/include/odb_Enum.h",
          odbFldType, fieldName.c_str());
  }

  return 0;
}


/* getFieldData
 *
 * this is to read arbitrary fields from the odb
 * (Without interpolation, no remap functionality. This may be done on the
 * receiving side of the pipe.)
 *
 * If createFilterNset or setFilterNset and/or createFilterElset or
 * setFilterElset have been called before only those nodes and elements will
 * be considered that belong to those sets.
 *
 * input arguments: fieldName (string), component (int), invariant (int)
 *
 * component is a zero based component index. if -1 it's not a component
 * but either the value as is (float*) or an invariant
 * invariant is an int with the following meaning:
 * 0 = no invariant (so either the value as is (float*) or a component)
 * 1 = magnitude
 * 2 = mises
 * 3 = tresca
 * 4 = press
 * 5 = inv3
 * 6 = maxPrincipal
 * 7 = midPrincipal
 * 8 = minPrincipal
 * 9 = maxInPlanePrincipal
 * 10 = minInPlanePrincipal
 * 11 = outOfPlanePrincipal
 *
 *
 * Output:
 *
 * Whether the requested field is defined at the nodes or integration points
 * and whether it's scalar or vector has to be determined by means of the
 * method getFieldProps beforehand.
 * The first output of getFieldData is a sequence of (integration point nb,
 * section point nb)-tuples. This sequence only consists of as many items
 * before the sequence repeats itself, i.e. only one period.
 * E.g. for nodal values or other non-integration-point-values this should
 * normally be only one item: (0,1); for integration point values on C3D10M
 * elements the sequence should look like this: (1,1), (2,1), (3,1), (4,1).
 * In order to get the integration point and section point number for each
 * field value you have to iterate over this integrPt/sectPt field in parallel
 * to the field-values-array. Whenever the integrPt/sectPt field is exhausted
 * you start again from its beginning.
 * The second output is an array of node or element labels.
 * The third output is an array of the field values.
 *
 * Here is the data in detail that will be written to the pipe:
 * - nbIntPtSectPt: number of integration point / section point tuples in the
 *   following array. An integer. E.g. 1 for nodal values and 4 for SDV1 on
 *   C3D10M elements.
 * - (integration point, section point) - tuples: 2*nbIntPtSectPt integers.
 * - nb of components of each field value: 1 for scalars, 3 for vectors, 6 for
 *   tensors =:nbComp
 * - nb of field values =:nbVals
 * - integer array of node or element labels (nbVals ints)
 * - float array of field values. (nbComp*nbVals floats)
 *
 * Note that the number of field values does not equal the number of e.g.
 * elements because the same element label (and corresponding value(s) might
 * occur multiple times for different integration or section points.
 */
int getFieldData() {

  /* get and check arguments */
  std::string fieldName = pipeReadString();
  int component = pipeReadInt();
  int invariant = pipeReadInt();

  debug_print("This is getFieldData, reading %s, comp %d, invar %d...\n",
              fieldName.c_str(), component, invariant);

  if (!ptOdbFrame) {
    error("getFieldData: Frame not set.\n");
  }

  /* we need the part instance for the filter sets, make sure it exists */
  if (ptInstance == NULL) {
    setPtInstance();
  }

  /* get odb_FieldOutput object from fieldName */
  debug_print("getFieldData: locating field %s in"
              " the odb ...\n", fieldName.c_str());
  odb_FieldOutput odbField;
  odb_FieldOutputRepository& fieldRepos
    = ptOdbFrame->fieldOutputs();
  if (! fieldRepos.isMember(fieldName.c_str()) ) {
    error("getFieldData: Field %s does not"
          " exist in step %d frame %d.\n",
          fieldName.c_str(), currentStep, currentFrame);
  }
  odbField = fieldRepos.get(fieldName.c_str());

  /* determine position (nodal/IP/...) */
  odb_Enum::odb_ResultPositionEnum odbFldPos
    = odbField.values()[0].position();

  /* depending on odbFldPos determine:
   * - use element or node label?
   * - filterset, possibly initialized by OdbReader_storeRemapWeights */
  debug_print("getFieldData: initializing filter set and getLabelFunc ...\n");
  int (odb_FieldValue::*getLabelFunc)() const;  // getLabelFunc: ptr to member funct
  odb_Set* filterset = NULL;

  switch (odbFldPos) {
  case odb_Enum::NODAL:
    debug_print("getFieldData: reading nodal field data from odb.\n");
    getLabelFunc = &odb_FieldValue::nodeLabel;

    if (filterNsetName && filterNsetName->length()>0) {
      filterset = &(ptInstance->nodeSets()[*(filterNsetName)]);
      debug_print("getFieldData: filter nset %s with %d nodes.\n",
                  filterset->name().cStr(), filterset->nodes().size());
    }
    break;
  case odb_Enum::INTEGRATION_POINT:
  case odb_Enum::CENTROID:
  case odb_Enum::WHOLE_ELEMENT:
    debug_print("getFieldData: reading element field data from odb.\n");
    getLabelFunc = &odb_FieldValue::elementLabel;

    if (filterElsetName && filterElsetName->length()>0) {
      filterset =&(ptInstance->elementSets()[*(filterElsetName)]);
      debug_print("getFieldData: filter elset %s with %d elements.\n",
                  filterset->name().cStr(), filterset->elements().size());
    }
    break;
  default:
    error("getFieldData: position %d for field %s not"
          " implemented. Search for odb_ResultPositionEnum in"
          " /usr/abaqus/6.13-2/code/include/odb_Enum.h",
          odbFldPos, fieldName.c_str());
    break;
  }

  /* applying the filterset if applicable */
  if (filterset) {
    debug_print("getFieldData: applying the filter set ...\n");
    odbField = odbField.getSubset(*filterset);
  }
  const odb_SequenceFieldValue& odbFieldVals = odbField.values();

  /* determine type (scalar, vector, tensor) and number of components */
  odb_Enum::odb_DataTypeEnum odbFldType
    = odbFieldVals[0].type();
  int nbComps;
  if (component>-1 || invariant>0 || odbFldType == odb_Enum::SCALAR) {
    nbComps=1;
  }
  else if (odbFldType == odb_Enum::VECTOR) {
    nbComps=3;
  }
  else if (odbFldType == odb_Enum::TENSOR_3D_FULL) {
    nbComps=6;
  }
  else {
    error("getFieldData: field type %d for field %s not"
          " implemented. Search for odb_DataTypeEnum in"
          " /usr/abaqus/6.13-2/code/include/odb_Enum.h",
          odbFldType, fieldName.c_str());
  }
  debug_print("getFieldData: determined the number of components to store"
              " for field %s: %d\n", fieldName.c_str(), nbComps);

  /* initializing the arrays to store all data */
  /* because of it's (potential) size it's on the heap not the stack */
  int nbVals = odbFieldVals.size();
  int* labels = new int[nbVals];              /* ...labels */
  float* values = new float[nbVals*nbComps];  /* ...values */
  int* intPtSectPt = new int[2*nbVals];       /* ...(integr pt, section pt) */

  /* storing all values in arrays, invariants (scalar) values */
  if (invariant>0) {

    GetInvariantFuncPtr getInvariantFunc = invariantIntToFunc[invariant];

    for (int i=0; i<nbVals; ++i) {
      const odb_FieldValue& val = odbFieldVals[i];

      labels[i] = (val.*getLabelFunc)();
      values[i] = (val.*getInvariantFunc)();
      if (odbFldPos==odb_Enum::INTEGRATION_POINT)
        intPtSectPt[i*2+0] = val.integrationPoint();
      else
        intPtSectPt[i*2+0] = 0;
      intPtSectPt[i*2+1] = val.sectionPoint().number();
    }
    debug_print("getFieldData: stored data for invariant %d for %d field"
                " values in arrays.\n", invariant, nbVals);
    //################################33 DEBUG
    for (int i=0; i<nbVals; ++i) {
      debug_print("getFieldData: label %d intpt %d, sectpt %d, value %g\n",
                  labels[i], intPtSectPt[i*2+0], intPtSectPt[i*2+1], values[i]);
    }
    //################################33 DEBUG
  }

  /* storing all values in arrays, scalar values */
  else if (nbComps==1) {

    /* just to make sure: for scalar values ignore -possibly wrong- component
     * argument */
    if (odbFldType == odb_Enum::SCALAR) {
      debug_print("getFieldData: fldType: SCALAR.\n");
      component = 0;
    }

    for (int i=0; i<nbVals; ++i) {
      const odb_FieldValue& val = odbFieldVals[i];

      labels[i] = (val.*getLabelFunc)();
      values[i] = val.data()[component];
      // //################################33 DEBUG
      // debug_print("getFieldData: label %d, value %g, %g.\n",
      //             labels[i], values[i], odbFieldVals[i].data()[0]);
      // //################################33 DEBUG
      if (odbFldPos==odb_Enum::INTEGRATION_POINT)
        intPtSectPt[i*2+0] = val.integrationPoint();
      else
        intPtSectPt[i*2+0] = 0;
      intPtSectPt[i*2+1] = val.sectionPoint().number();
    }
    debug_print("getFieldData: stored scalar data for %d field values in"
                " arrays.\n", nbVals);
    //################################33 DEBUG
    for (int i=0; i<nbVals; ++i) {
      debug_print("getFieldData: label %d intpt %d, sectpt %d, value %g\n",
                  labels[i], intPtSectPt[i*2+0], intPtSectPt[i*2+1], values[i]);
    }
    //################################33 DEBUG
  }

  /* storing all values in arrays, vector/tensor */
  else {
    for (int i=0; i<nbVals; ++i) {
      const odb_FieldValue& val = odbFieldVals[i];

      labels[i] = (val.*getLabelFunc)();
      for (int j=0; j<nbComps; ++j) {
        values[i*nbComps+j] = val.data()[j];
      }
      if (odbFldPos==odb_Enum::INTEGRATION_POINT)
        intPtSectPt[i*2+0] = val.integrationPoint();
      else
        intPtSectPt[i*2+0] = 0;
      intPtSectPt[i*2+1] = val.sectionPoint().number();
    }
    debug_print("getFieldData: stored vector data for %d field values in"
                " arrays.\n", nbVals);
    //################################33 DEBUG
    for (int i=0; i<nbVals; ++i) {
      debug_print("getFieldData: label %d intpt %d, sectpt %d, value %g\n",
                  labels[i], intPtSectPt[i*2+0], intPtSectPt[i*2+1], values[i]);
    }
    //################################33 DEBUG
  }

  /* check if intPtSectPt is periodic */
  bool isPeriodic = false;
  int nbIntPtSectPt = -1;
  for (int i=1; i<nbVals; ++i) {
    if (!isPeriodic) {
      /* periodic series not found yet: check if it starts here */
      if ((intPtSectPt[i*2+0]==intPtSectPt[0*2+0])
          and (intPtSectPt[i*2+1]==intPtSectPt[0*2+1])) {
        isPeriodic = true;
        nbIntPtSectPt = i;
      }
    }
    else {
      /* if periodic: check that it stays periodic */
      if (!((intPtSectPt[i*2+0]==intPtSectPt[(i%nbIntPtSectPt)*2+0])
            and (intPtSectPt[i*2+1]==intPtSectPt[(i%nbIntPtSectPt)*2+1]))) {
        isPeriodic = false;
        break;
      }
    }
  }
  if (!isPeriodic) {
    nbIntPtSectPt = nbVals;
  }
  debug_print("getFieldData: Checked periodicity of intPtSectPt."
              " nbIntPtSectPt=%d\n", nbIntPtSectPt);

  /* Send intPtSectPt-array through the pipe */
  debug_print("getFieldData: Will write %d intPr/SecPt tuples.\n",
              nbIntPtSectPt);
  pipeWrite(nbIntPtSectPt);
  pipeWrite(intPtSectPt, 2*nbIntPtSectPt);
  debug_print("getFieldData: Wrote intPtSectPt.\n");

  /* Send the field data through the pipe */
  debug_print("getFieldData: Will write %d labels and field values with %d"
              " components.\n", nbVals, nbComps);
  pipeWrite(nbComps);
  pipeWrite(nbVals);
  pipeWrite(labels, nbVals);
  pipeWrite(values, nbComps*nbVals);
  debug_print("getFieldData: Wrote labels and field values.\n");

  /* release data arrays */
  delete[] labels;
  delete[] values;
  delete[] intPtSectPt;

  debug_print("getFieldData: Ending.\n");
  return 0;
}


/* getMatNames
 *
 * This is to read the material assignment from the odb.
 *
 * If createFilterElset or setFilterElset has been called before only elements
 * from this set will be considered.
 *
 * input arguments: firstMatCode (int), nbSecTypeFilters (int),
 *     secTypeFilters (strings)
 *
 * firstMatCode determines the lowest (possible) material number in the result
 * corresponding to the first material name returned. It's typically =1.
 *
 * nbSecTypeFilters is the number of items in the following secTypeFilters
 * sequence of strings.
 * If you don't need this feature then pass zero for nbSecTypeFilters and skip
 * the secTypeFilters argument(s) alltogether.
 *
 * secTypeFilters are strings. Only section / material names starting with one
 * of those strings will be considered in the output. Elements whose section /
 * material name does not begin with one of these strings will be ignored.
 *
 *
 * Output:
 *
 * The first output of getMatNames is an array of element labels.
 * The second output is an array of the material numbers.
 * The third output is a sequence of material names.
 *
 * Here is the data in detail that will be written to the pipe:
 * - nb of items =:nbElems
 * - integer array of element labels (nbElems ints)
 * - integer array material numbers (nbElems ints)
 * - nb of material names =:nbMatNames
 * - a sequence of nbMatNames strings giving the material (section) names from
 *   the odb.
 */
int getMatNames() {
  // no matcodes allowed smaller or equal to this
  const int noMatCode(-9999999);

  /* get and check arguments */
  int firstMatCode = pipeReadInt();
  int nbSecTypeFilters = pipeReadInt();
  std::vector<std::string> secTypeFilters;
  for (int i=0; i<nbSecTypeFilters; ++i) {
    secTypeFilters.push_back(pipeReadString());
  }
  debug_print("This is getMatNames, firstMatCode %d, %d type filter strings.\n",
              firstMatCode, nbSecTypeFilters);
  if (firstMatCode<=noMatCode) {
    error("getMatNames: firstMatCode must be greater than %d. Got %d.\n",
          noMatCode, firstMatCode);
  }

  /* initialize/get Abaqus instance */
  if (ptInstance == NULL) {
    setPtInstance();
  }

  /* determine applicable elements / apply filterset */
  bool applyFilter = (filterElsetName && filterElsetName->length()>0);
  std::set<int> filterset;
  int nbElems;
  if (applyFilter) {
    odb_Set odbfilterset = ptInstance->elementSets()[*(filterElsetName)];
    const odb_SequenceElement& elements = odbfilterset.elements();
    for (int i=0; i<elements.size(); ++i) {
      filterset.insert(elements[i].label());
    }
    nbElems = filterset.size();
    debug_print("getMatNames: filter elset %s with %d elements.\n",
                odbfilterset.name().cStr(), nbElems);
  }
  else {
    const odb_SequenceElement& elements = ptInstance->elements();
    nbElems = elements.size();
    debug_print("getMatNames: no filterset, considering all %d elements.\n",
                nbElems);
  }

  /* examine section assignments */
  int nextNewMatCode = firstMatCode;
  int matCode;
  bool secTypeIsInFilterTuple;

  // initialize result objects: arrays to store all data
  // because of it's (potential) size it's on the heap not the stack
  int iElems = 0;
  int* labels = new int[nbElems];
  int* matCodes = new int[nbElems];
  std::vector<std::string> matNamesList;
  const char* matName = NULL;
  const char* secType = NULL;
  
  odb_SequenceSectionAssignment sectionAssignmentSeq =
    ptInstance->sectionAssignments();
  for (int i = 0; i < sectionAssignmentSeq.size(); ++i) {
    odb_SectionAssignment sa = sectionAssignmentSeq[i];
    // we get a std::bad_alloc exception if the elset is not explicitly created:
    odb_Set saElset = sa.region();
    const odb_SequenceElement& elements = saElset.elements();
    int elemsInThisSect = elements.size();
    const odb_Section& sect = sa.section();

    if (odb_isA(odb_CohesiveSection, sect)) {
      odb_CohesiveSection sect2 = odb_dynamicCast(odb_CohesiveSection, sect);
      matName = sect2.material().cStr();
      secType = "cohesive";
      debug_print("getMatNames: found cohesive section %s with material %s and"
                  " %d elements\n",
                  sa.sectionName().cStr(), matName, elemsInThisSect);
    }
    else if (odb_isA(odb_HomogeneousSolidSection, sect)) {
      odb_HomogeneousSolidSection sect2 =
        odb_dynamicCast(odb_HomogeneousSolidSection, sect);
      matName = sect2.material().cStr();
      secType = "solid";
      debug_print("getMatNames: found solid section %s with material %s and"
                  " %d elements\n",
                  sa.sectionName().cStr(), matName, elemsInThisSect);
    }
    else if (odb_isA(odb_BeamSection, sect)) {
      odb_BeamSection sect2 = odb_dynamicCast(odb_BeamSection, sect);
      matName = sect2.material().cStr();
      secType = "beam";
      debug_print("getMatNames: found beam section %s with material %s and"
                  " %d elements\n",
                  sa.sectionName().cStr(), matName, elemsInThisSect);
    }
    else if (odb_isA(odb_HomogeneousShellSection, sect)) {
      odb_HomogeneousShellSection sect2 =
        odb_dynamicCast(odb_HomogeneousShellSection, sect);
      matName = sect2.material().cStr();
      secType = "shell";
      debug_print("getMatNames: found shell section %s with material %s and"
                  " %d elements\n",
                  sa.sectionName().cStr(), matName, elemsInThisSect);
    }
    else if (odb_isA(odb_TrussSection, sect)) {
      odb_TrussSection sect2 = odb_dynamicCast(odb_TrussSection, sect);
      matName = sect2.material().cStr();
      debug_print("getMatNames: found truss section %s with material %s and"
                  " %d elements\n",
                  sa.sectionName().cStr(), matName, elemsInThisSect);
      secType = "truss";
    }
    else {
      error("getMatNames: section type not implemented for section %s with"
            " %d elements.\n", sa.sectionName().cStr(), elemsInThisSect);
    }

    /* check that secType passes secTypeFilters */

    secTypeIsInFilterTuple = false;
    for (int j=0; j<secTypeFilters.size(); ++j) {
      if (strncmp(secType, secTypeFilters[j].c_str(),
                  secTypeFilters[j].size()) == 0) {
        secTypeIsInFilterTuple = true;
        break;
      }
    }
    debug_print("getMatNames: secTypeIsInFilterTuple: %s\n",
                (secTypeIsInFilterTuple?"True":"False"));

    if (secTypeFilters.size() && !secTypeIsInFilterTuple) {
      /* if filtertuple is specified and if current section type not in the
       * filter list then ignore this section assignment */
      debug_print("getMatNames: ignoring section assignment %s.\n", secType);
      // skip storing of new section assignment
      // and element label
      continue;
    }

    /* store matName in list */
    matNamesList.push_back(matName);

    /* get new matCode */
    matCode = nextNewMatCode;
    ++nextNewMatCode;
    debug_print("getMatNames: Stored mat-name %s for mat code %d\n",
                matName, matCode);

    /* loop over elements and update  */
    int ii;
    for (ii=0; ii<elemsInThisSect; ++ii) {
      int elemNb = elements[ii].label();
      if (!applyFilter || (filterset.find(elemNb) != filterset.end())) {
        if (iElems>=nbElems) {
          error("getMatNames: Too many elements in this section %s: %d."
                "\nAlready stored: %d, total expected: %d."
		"\nThis is very likely caused by overlapping, i.e. ambiguous"
		" section assignments.",
                sa.sectionName().cStr(), elemsInThisSect, iElems, nbElems);
        }

        // store (only valid) matCode for element
        // debug_print("getMatNames: storing element %d with matCode %d.\n",
        //             elemNb, matCode);
        labels[iElems] = elemNb;
        matCodes[iElems] = matCode;
        ++iElems;
      }
    }
  }

  /* Send the result data through the pipe */
  debug_print("getMatNames: Will write %d elem numbers and matCodes.\n", iElems);
  pipeWrite(iElems);
  if (iElems) {
    pipeWrite(labels, iElems);
    pipeWrite(matCodes, iElems);
  }

  debug_print("getMatNames: Will write %d matNames.\n", matNamesList.size());
  pipeWrite((int) matNamesList.size());
  for (int i=0; i<matNamesList.size(); ++i) {
    pipeWrite(matNamesList[i]);
  }

  /* release data arrays */
  delete[] labels;
  delete[] matCodes;

  debug_print("getMatNames: Ending.\n");
  return 0;
}


/******************************************************************************
*******************************************************************************
*** map of server functions                                                 ***
*******************************************************************************
******************************************************************************/

typedef int (*FnPtr)();
std::map<std::string, FnPtr> funcMap;
void initFuncMap() {
  funcMap["test"] = test;
  funcMap["init"] = init;
  funcMap["openOdb"] = openOdb;
  funcMap["closeOdb"] = closeOdb;
  funcMap["getNodeCoordsFromOdb"] = getNodeCoordsFromOdb;
  funcMap["getElementsFromOdb"] = getElementsFromOdb;
  funcMap["getElsetsFromOdb"] = getElsetsFromOdb;
  funcMap["setFrame"] = setFrame;
  funcMap["getFieldProps"] = getFieldProps;
  funcMap["getFieldData"] = getFieldData;
  funcMap["getMatNames"] = getMatNames;
  funcMap["createFilterNset"] = createFilterNset;
  funcMap["setFilterNset"] = setFilterNset;
  funcMap["createFilterElset"] = createFilterElset;
  funcMap["setFilterElset"] = setFilterElset;
  // funcMap[""] = ;
  // funcMap[""] = ;
  // funcMap[""] = ;
  // funcMap[""] = ;
  // funcMap[""] = ;
  // funcMap[""] = ;
  // funcMap[""] = ;
  // funcMap[""] = ;

  // these are "private"...
  // funcMap["setPtInstance"] = setPtInstance;
}
