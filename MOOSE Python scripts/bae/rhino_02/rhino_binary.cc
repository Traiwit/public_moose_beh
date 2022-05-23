
/*****************************************************************************/
/*   module for OpenNurbs API / Python integration                              */
/*****************************************************************************/

/******************************************************************************
 *****************************************************************************/
static const char
module_version[] = "1.3";
/* version history:
 *
 * 1.1: first versioned
 * 1.2: unit now metre. Added Mesh.setUserText
 * 1.3: moved Mech object into separate file, added PointCloud object
 *
 ******************************************************************************
 *****************************************************************************/


#include <stdio.h>
#include <string.h>
#include <map>
#include <vector>
#include <string>

#include <stdlib.h>
// #include <assert.h>

/* Python integration */
#include <Python.h>
#include <structmember.h>   // needed for e.g. T_INT

/* OpenNurbs integration */
#include "opennurbs/opennurbs.h"
// #include "opennurbs/examples_linking_pragmas.h"
// #include "opennurbs/example_userdata/example_ud.h"


// #define DEBUG
#undef DEBUG

#ifdef DEBUG
#define debug_print(fmt, ...) fprintf(stderr, fmt, ##__VA_ARGS__)
#else
#define debug_print(fmt, ...) do {} while (0)
#endif


/******************************************************************************
*******************************************************************************
***  Base class for all Rhino Geometry objects                              ***
*******************************************************************************
******************************************************************************/

class GeomObj {
public:
  PyObject_HEAD
  ON_Geometry *m_geomObj;
  ON_3dmObjectAttributes *m_attributes;
};


/* getting and setting objectName
 */
static PyObject *
GeomObj_get_objectName(GeomObj *self, void *closure)
{

  if (self->m_attributes) {
    PyObject *objectName;
    objectName = PyUnicode_FromWideChar
      (self->m_attributes->m_name,
       (Py_ssize_t) self->m_attributes->m_name.Length());
    return objectName;
  }
  else {
    Py_RETURN_NONE;
  }
}
static int
GeomObj_set_objectName(GeomObj *self, PyObject *value, void *closure)
{
  if (!self->m_attributes) {
    PyErr_SetString
      (PyExc_TypeError,
       "The geometry object has no attribute yet. Can't assign objectName.");
    return -1;
  }
  if (value == NULL) {
    // for: del GeomObj.objectName
    self->m_attributes->m_name.Empty();
    return 0;
  }
  else if (PyString_Check(value)) {
    self->m_attributes->m_name = PyString_AsString(value);
    return 0;
  }
  else if (PyUnicode_Check(value)) {
    Py_ssize_t len = PyUnicode_GetSize(value);
    wchar_t* w = new wchar_t[len+1];
    if (!w) {
      PyErr_SetString(PyExc_TypeError, "Could not assign new objectName"
                      " attribute.");
      return -1;
    }
    int rc = PyUnicode_AsWideChar((PyUnicodeObject*) value, w, len);
    if (rc<0) {
      PyErr_SetString(PyExc_TypeError, "Could not convert new objectName"
                      " attribute to wide char.");
      delete[] w;
      return -1;
    }
    w[rc] = 0;  /* add terminating zero */
    self->m_attributes->m_name = w;
    delete[] w;
    return 0;
  }

  PyErr_SetString(PyExc_TypeError, 
                  "The objectName attribute value must be a string");
  return -1;
}

/* getting and setting objectLayerIndex
 */
static PyObject *
GeomObj_get_objectLayerIndex(GeomObj *self, void *closure)
{
  if (self->m_attributes) {
    PyObject *objectLayerIndex;
    objectLayerIndex = PyInt_FromLong
      ((long) self->m_attributes->m_layer_index);
    return objectLayerIndex;
  }
  else {
    Py_RETURN_NONE;
  }
}
static int
GeomObj_set_objectLayerIndex(GeomObj *self, PyObject *value, void *closure)
{
  if (!self->m_attributes) {
    PyErr_SetString
      (PyExc_TypeError,
       "Geometry object has no attribute yet. Can't assign objectLayerIndex.");
    return -1;
  }
  if (value == NULL) {
    // for: del GeomObj.objectLayerIndex
    self->m_attributes->m_layer_index = 0;
    return 0;
  }
  else if (PyInt_Check(value)) {
    self->m_attributes->m_layer_index = (int) PyInt_AsLong(value);
    return 0;
  }

  PyErr_SetString
    (PyExc_TypeError, 
     "GeomObj objectLayerIndex attribute value must be an integer");
  return -1;
}

/* setting UserString values
 */
static PyObject*
GeomObj_setUserString(GeomObj *self, PyObject *args, PyObject *kwds) {

  const char* key;
  const char* value;
  wchar_t* wc_key;
  wchar_t* wc_value;
  size_t len_key;
  size_t len_value;
  // const wchar_t* key;
  // const wchar_t* value;

  /* parse arguments */
  static char *kwlist[]
    = {"key", "value", NULL};
  if (! PyArg_ParseTupleAndKeywords(args, kwds, "ss:GeomObj_setUserString",
                                    kwlist, &key, &value) ) {
    return NULL;
  }
  debug_print("GeomObj_setUserString key: %s, value %s.\n", key, value);

  /* convert char to wchar */
  len_key = strlen(key);
  wc_key = new wchar_t[len_key+1];
  wc_key[len_key] = 0;
  len_value = strlen(value);
  wc_value = new wchar_t[len_value+1];
  wc_value[len_value] = 0;
  if (!wc_key || !wc_value) {
    PyErr_SetString(PyExc_TypeError, "Could not assign new UserString."
                    " (wchar_t-string init failed.)");
    return -1;
  }
  if (mbstowcs(wc_key, key, strlen(key)) != len_key){
    PyErr_SetString(PyExc_TypeError, "Could not assign new UserString."
                    " (char to wchar_t string conversion of key failed.)");
    return -1;
  }
  if (mbstowcs(wc_value, value, strlen(value)) != len_value){
    PyErr_SetString(PyExc_TypeError, "Could not assign new UserString."
                    " (char to wchar_t string conversion of value failed.)");
    return -1;
  }

  /* update the ON_3dmObjectAttributes object with provided userStrings */
  self->m_attributes->SetUserString(wc_key, wc_value);

  Py_RETURN_NONE;
}

/******************************************************************************
*******************************************************************************
*** classes for Rhino objects in separate files                             ***
*******************************************************************************
******************************************************************************/

#include "rhino_bin_mesh.cc"
#include "rhino_bin_pointcloud.cc"

/******************************************************************************
*******************************************************************************
***  RhinoModel class                                                       ***
*******************************************************************************
******************************************************************************/

/* data structure for the RhinoModel object
 */
typedef struct {
  PyObject_HEAD
  ONX_Model *m_model;

  /* list of objects in the Rhino model */
  PyObject *objects;

  /* list of layers in the Rhino model
   * each layer being a [name, colour, parent]-list, colour is an RGB tuple
   * of integers 0...255, parent is an integer index of the parent layer
   * (None for no parent) */
  PyObject *layers;
} RhinoModel;

static PyMemberDef RhinoModel_members[] = {
  {"objects", T_OBJECT_EX, offsetof(RhinoModel, objects), 0,
   "list of Rhino objects, for example of type bae.rhino_02.rhino_binary.Mesh"},
  {"layers", T_OBJECT_EX, offsetof(RhinoModel, layers), 0,
   "list of layers, each layer being a [name, colour, parent_layer_id]-list,"
   " colour is an RGB tuple, parent_layer_id an integer index of the parent"
   " layer (None for no parent)."},
  {NULL}  /* Sentinel */
};


/******************************************************************************
*******************************************************************************
***  function to write the complete Rhino file                              ***
*******************************************************************************
******************************************************************************/

/* RhinoModel_write
 *
 * write a Rhino file
 *
 */
static PyObject*
RhinoModel_write(RhinoModel* self, PyObject *args, PyObject *kwds) {

  /* 1. argument: filename */
  const char* filename = NULL;

  // The OpenNURBS toolkit will write version 2 and 3 and read
  // version 1, 2 and 3 of the 3DM file format.
  //
  // version 1 is the legacy Rhino I/O tookit format and was used by Rhino 1.x.
  // version 2 is the OpenNURBS format (released 1 July 2000) and is used by
  //           Rhino 2.x
  // version 3 is the OpenNURBS format (released 1 November 2002) and is used
  //           by Rhino 3.x
  // version 4 is the OpenNURBS format (released September 2006) and is used
  //           by Rhino 4.x
  // version 5 is the OpenNURBS format (released September 2009) and is used
  //           by Rhino 5.x

  /* second argument: version to write */
  int version = 0; // version will be ON_BinaryArchive::CurrentArchiveVersion()

  /* parse arguments */
  static const char *kwlist[] = {"filename", "version", NULL};
  debug_print("RhinoModel_write Step 1.\n");
  if (!PyArg_ParseTupleAndKeywords
      (args, kwds, "s|i:write", const_cast<char **>(kwlist),
       &filename, &version) ) {
    return NULL;
  }
  debug_print("RhinoModel_write Step 2.\n");

  ON::Begin();
  // If you want to learn to write b-rep models, first work through
  // this example paying close attention to write_trimmed_surface_example(),
  // then examime example_brep.cpp.

  // errors printed to stdout
  ON_TextLog error_log;

  // messages printed to stdout
  ON_TextLog message_log;

  // pointer to file handler
  FILE* fp;
  
  Py_BEGIN_ALLOW_THREADS
  debug_print("Opening %s for writing...\n", filename);
  fp = ON::OpenFile( filename, "wb" );

  // set revision history information
  self->m_model->m_properties.m_RevisionHistory.NewRevision();
  
  // set application information
  self->m_model->m_properties.m_Application.m_application_name
    = "rhino_binary module";
  self->m_model->m_properties.m_Application.m_application_URL
    = "http://www.opennurbs.org";
  self->m_model->m_properties.m_Application.m_application_details
    = "Rhino file writing Python API.";
  debug_print("Added some file properties.\n");

  // set unit to metres
  self->m_model->m_settings.m_ModelUnitsAndTolerances.m_unit_system
    = ON::meters;

  Py_END_ALLOW_THREADS

  // layer table
  {
    int layer_count = PyList_Size(self->layers);
    
    // Each object in the object table (written below)
    // should be on a defined layer.  There should be
    // at least one layer with layer index 0 in every file.

    // layer table indices begin at 0
    if ( layer_count <= 0 ) {
      ON_Layer default_layer;
      default_layer.SetLayerIndex(0);
      default_layer.SetLayerName("Default");
      self->m_model->m_layer_table.Append(default_layer);
      layer_count = 1;
    }
    else {

      debug_print("%d layers to be written to the Rhino file.\n", layer_count);
      ON_Layer layer;
      PyObject* layerData;
      PyObject* value;
      self->m_model->m_layer_table.Reserve(layer_count);
      debug_print("...'Reserved' %d layers.\n", layer_count);
      for (int i = 0; i < layer_count; i++) {
        layer.SetLayerIndex(i);

        layerData = PyList_GET_ITEM(self->layers, i);
        if (!PyList_Check(layerData)) {
          PyErr_Format
            (PyExc_ValueError,
             "RhinoModel_write: self.layers must only contain lists. The %d."
             " item of self.layers is not a list.", i);
          return NULL;
        }
        Py_XINCREF(layerData);

        // assign the layer name (item idx 0 in layerData)
        value = PyList_GetItem(layerData, 0);
        if (PyString_Check(value)) {
          layer.SetLayerName(PyString_AsString(value));
        }
        else if (PyUnicode_Check(value)) {
          Py_ssize_t len = PyUnicode_GetSize(value);
          wchar_t* w = new wchar_t[len+1];
          if (!w) {
            PyErr_Format
              (PyExc_TypeError, "RhinoModel_write: Could not assign layer"
               " name for layer index %d.", i);
            Py_XDECREF(layerData);
            return NULL;
          }
          int rc = PyUnicode_AsWideChar((PyUnicodeObject*) value, w, len);
          if (rc<0) {
            PyErr_Format(PyExc_TypeError, "Could not convert new layer name"
                         " to wide char for layer index %d.", i);
            delete[] w;
            Py_XDECREF(layerData);
            return NULL;
          }
          w[rc] = 0;  /* add terminating zero */
          layer.SetLayerName(w);
          delete[] w;
        }
        else {
          // if the first item is not any type of string...
          PyErr_Format
            (PyExc_TypeError, "RhinoModel_write: No layer name given for"
             " layer index %d.", i);
          Py_XDECREF(layerData);
          return NULL;
        }
        debug_print("Assigned name to layer %d.\n", i);

        // assign colour (item idx 1 in layerData)
        value = PyList_GetItem(layerData, 1);
        debug_print("Retrieved colour value for layer %d.\n", i);
        if (value == NULL || value == Py_None) {
          // if there is no second item in the list describing the particular
          // layer then assume default clour RGB=0,0,0
          debug_print("No colour value for layer %d, will assign black.\n", i);
          PyErr_Clear();
          layer.SetColor( ON_Color(0, 0, 0) );
          debug_print("Assigned default colour (0,0,0) to layer %d.\n", i);
        }
        else {
          PyObject *colVal;
          long rgb[3];
          debug_print("Interpreting colour value for layer %d...\n", i);
          Py_XINCREF(value);
          if (PyTuple_Check(value)) {
            if (PyTuple_Size(value)!=3) {
              PyErr_Format(PyExc_ValueError, "Colour value must be a RGB tuple"
                           " of 3 integers. For layer %d the colour value tuple"
                           " contains %d items.", i, PyTuple_Size(value));
              Py_XDECREF(value);
              Py_XDECREF(layerData);
              return NULL;
            }              
            for (int j=0; j<3; j++) {
              colVal = PyTuple_GetItem(value, j);
              if (colVal==NULL) {
                PyErr_Format(PyExc_ValueError, "Colour value must be a RGB"
                             " tuple of 3 integers. It's not for layer index"
			     " %d.", i);
                Py_XDECREF(value);
                Py_XDECREF(layerData);
                return NULL;
              }
              if (!PyInt_Check(colVal)) {
                PyErr_Format(PyExc_ValueError, "Colour value must be a RGB"
                             " tuple of 3 integers. It's not for layer index"
                             " %d, colour index %d.", i, j);
                Py_XDECREF(value);
                Py_XDECREF(layerData);
                return NULL;
              }
              rgb[j] = PyInt_AS_LONG(colVal);
            }  // end for (int j=0; j<3; j++) ...
          }  // end if value is tuple ...
          else if (PyList_Check(value)) {
            PyObject *colVal;
            if (PyList_Size(value)!=3) {
              PyErr_Format(PyExc_ValueError, "Colour value must be a RGB list"
                           " of 3 integers. For layer %d the colour value list"
                           " contains %d items.", i, PyList_Size(value));
              Py_XDECREF(value);
              Py_XDECREF(layerData);
              return NULL;
            }              
            for (int j=0; j<3; j++) {
              colVal = PyList_GetItem(value, j);
              if (colVal==NULL) {
                PyErr_Format(PyExc_ValueError, "Colour value must be a RGB list"
                             " of 3 integers. It's not for layer index %d.", i);
                Py_XDECREF(value);
                Py_XDECREF(layerData);
                return NULL;
              }
              if (!PyInt_Check(colVal)) {
                PyErr_Format(PyExc_ValueError, "Colour value must be a RGB list"
                             " of 3 integers. It's not for layer index %d,"
                             " colour index %d.", i, j);
                Py_XDECREF(value);
                Py_XDECREF(layerData);
                return NULL;
              }
              rgb[j] = PyInt_AS_LONG(colVal);
            }  // end for (int j=0; j<3; j++) ...
          }  // end if value is list ...
          else {
            PyErr_Format(PyExc_ValueError, "Colour value must be a RGB tuple"
                         " or list of 3 integers. It's no list or tuple for"
                         " layer index %d.", i);
            Py_XDECREF(value);
            Py_XDECREF(layerData);
            return NULL;
          }
          Py_XDECREF(value);
          debug_print("Assigning colour (%ld,%ld,%ld) to layer %d.\n",
                      rgb[0], rgb[1], rgb[2], i);
          layer.SetColor( ON_Color(rgb[0], rgb[1], rgb[2]) );
        }  // end if colour value given, or == NULL
        debug_print("Assigned RGB colour to layer %d.\n", i);


        // assign parent layer index (item idx 2 in layerData)
        value = NULL;
        if (PyList_Size(layerData)>2) {
          value = PyList_GetItem(layerData, 2);
          debug_print("Retrieved parent layer index for layer %d.\n", i);
        }
        if (value!=NULL && PyInt_Check(value)) {
          // assign parent layer index
          int parentIdx = (int) PyInt_AS_LONG(value);
          ON_Layer *parentLayer = self->m_model->m_layer_table.At(parentIdx);
          char parentLayerUuidStr[64];
          if (parentLayer==NULL) {
            PyErr_Format(PyExc_ValueError, "Parent layer index %ld does not"
                         " exist for layer index %d.", parentIdx, i);
            Py_XDECREF(layerData);
            return NULL;
          }
#ifdef DEBUG
          char uuidStr[64];
          ON_UuidToString(parentLayer->m_layer_id, uuidStr);
          debug_print
            ("Assigning parent layer %s (idx %ld, UUID %s) to layer %d.\n",
             ON_String(parentLayer->m_name).Array(),
             parentIdx, uuidStr, i);
#endif
          layer.m_parent_layer_id = parentLayer->m_layer_id;
        }
        else {
          // parent not assigned
          // ON_nil_uuid is declared in opennurbs_uuid.h
          layer.m_parent_layer_id = ON_nil_uuid;
        }

        Py_XDECREF(layerData);

        layer.SetVisible(true);
        layer.SetLocked(false);
        ON_CreateUuid(layer.m_layer_id);

        self->m_model->m_layer_table.Append(layer);
        debug_print("Stored layer %d in Rhino file.\n", i);
      }  // end for (int i = 0; i < layer_count; i++) ....
    }

    debug_print("Added %d layer.\n", layer_count);
  }

  // write objects
  for (int i=0; i<PyList_Size(self->objects); ++i) {

    Mesh* ptObj;
    ptObj = (Mesh*) PyList_GET_ITEM(self->objects, i);

    ONX_Model_Object& mo = self->m_model->m_object_table.AppendNew();
    mo.m_object = ptObj->m_geomObj;
    mo.m_bDeleteObject = false;
    mo.m_attributes = *(ptObj->m_attributes);
  }
  debug_print("Added %d mesh objects.\n", PyList_Size(self->objects));
  Py_BEGIN_ALLOW_THREADS

  // archive to write to
  ON_BinaryFile archive( ON::write3dm, fp );
  debug_print("RhinoModel_write Step 5.\n");

  // Set uuid's, indices, etc.
  self->m_model->Polish();
  debug_print("RhinoModel_write Step 6.\n");
  // writes model to archive
  bool ok = self->m_model->Write
    ( archive, version, __FILE__ " RhinoModel.write() " __DATE__, &error_log );

  ON::CloseFile( fp );
  if (ok)
    debug_print("Successfully wrote %s.\n",filename);
  else
    message_log.Print("Errors while writing %s.\n",filename);

  ON::End();

  Py_END_ALLOW_THREADS
  Py_RETURN_NONE;
}



/******************************************************************************
*******************************************************************************
***  function for initialize, open, close and such...                       ***
*******************************************************************************
******************************************************************************/

/* RhinoModel constructor: __new__()
 *
 */
static PyObject *
RhinoModel_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  RhinoModel *self;

  self = (RhinoModel *)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->m_model = NULL;
    self->objects = NULL; /* list of objects in the Rhino model */
    self->layers = NULL; /* list of layers in the Rhino model */
    debug_print("Initialized new C-RhinoModel.\n");
  }
  else {
    debug_print("RhinoModel_new failed to create new C-RhinoModel.\n");
  }
  return (PyObject *)self;
}


/* RhinoModel constructor: __init__()
 *
 */
static int
RhinoModel_init(RhinoModel *self, PyObject *args, PyObject *kwds)
{
  /* initialize ONX_Model object */
  if (self->m_model) {
    delete self->m_model;
  }
  self->m_model = new ONX_Model();
  if (self->m_model == NULL) {
    debug_print("RhinoModel_init failed to initialize new C-RhinoModel.\n");
    self->ob_type->tp_free((PyObject*)self);
    return -1;
  }

  /* initialize list of Rhino objects */
  Py_XDECREF(self->objects);
  self->objects = PyList_New(0);
  if (self->objects == NULL) {
    debug_print("Failed to initialize object list for new C-RhinoModel.\n");
    delete self->m_model;
    self->m_model = NULL;
    self->ob_type->tp_free((PyObject*)self);
    return -1;
  }

  /* initialize list of layers */
  Py_XDECREF(self->layers);
  self->layers = PyList_New(0);
  if (self->layers == NULL) {
    debug_print("Failed to initialize layer list for new C-RhinoModel.\n");
    delete self->m_model;
    self->m_model = NULL;
    self->ob_type->tp_free((PyObject*)self);
    return -1;
  }
  
  debug_print("End of RhinoModel_init.\n");
  return 0;
}


/* RhinoModel destructor
 *
 */
static void
RhinoModel_dealloc(RhinoModel* self)
{
  debug_print("In RhinoModel_dealloc (1).\n");
  if (self->m_model != NULL) {
    delete self->m_model;
    self->m_model = NULL;
  }
  Py_XDECREF(self->objects);
  Py_XDECREF(self->layers);
  self->ob_type->tp_free((PyObject*)self);
  debug_print("In RhinoModel_dealloc (2).\n");
}


/******************************************************************************
*******************************************************************************
***  Python API stuff                                                       ***
*******************************************************************************
******************************************************************************/


/* RhinoModel methods
 */
static PyMethodDef RhinoModel_methods[] = {
  {"write",
   (PyCFunction)RhinoModel_write, METH_VARARGS | METH_KEYWORDS,
   "Write a Rhino file."
   "\nArguments: filename, version (optional) of the Rhino file format."
   "\nSpecifying 0 results in the 'current version', this is the default."
   "\nVersion 1 is the legacy Rhino I/O tookit format and was used by"
   " Rhino 1.x."
   "\nVersion 2 is the OpenNURBS format (released 1 July 2000) and is used by"
   " Rhino 2.x."
   "\nVersion 3 is the OpenNURBS format (released 1 November 2002) and is used"
   " by Rhino 3.x."
   "\nVersion 4 is the OpenNURBS format (released September 2006) and is used"
   " by Rhino 4.x."
   "\nVersion 5 is the OpenNURBS format (released September 2009) and is used"
   " by Rhino 5.x."
   "\n"
  },
  {NULL}  /* Sentinel */
};


/* RhinoModel type description
 */
static PyTypeObject RhinoModelType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "rhino_binary.RhinoModel", /*tp_name*/
    sizeof(RhinoModel),         /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)RhinoModel_dealloc, /*tp_dealloc*/
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
    "RhinoModel class",        /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    RhinoModel_methods,        /* tp_methods */
    RhinoModel_members,        /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)RhinoModel_init, /* tp_init */
    0,                         /* tp_alloc */
    RhinoModel_new,            /* tp_new */
};

/*   Python module stuff
 */
static PyMethodDef ModuleMethods[] = {
  /* ... */
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initrhino_binary(void)
{
  // Initialize the module
  PyObject* m;
  m = Py_InitModule3
    ("rhino_binary", ModuleMethods,
     "Extension module for the (binary-) RhinoModel class."
     "\n"
     "\nExample"
     "\n======="
     "\n"
     "\n >>> from bae.rhino_02.rhino_binary import RhinoModel, Mesh"
     "\n >>>"
     "\n >>> model = RhinoModel()"
     "\n >>>"
     "\n >>> # Layer: name and RGB-tuple"
     "\n >>> model.layers.append( ['Default', (0,0,0)] )"
     "\n >>> model.layers.append( ['green layer', (0,1,0)] )"
     "\n >>>"
     "\n >>> # welded mesh: each vertex defined only once"
     "\n >>> vertices = [ [0,0,0], [2,0,0], [2,3,1], [0.5,3,1]]"
     "\n >>> faces = [ [0,3,1], [1, 2, 3] ]  # vertex indices"
     "\n >>> # texcoords: (element number, node number)-tuple for each vertex"
     "\n >>> texCoords = [[7,9], [7,9],[7,9],[7,9]]"
     "\n >>> objectName = 'quadrilateral made of two triangles'"
     "\n >>> layerIndex = 1  # index of 'green layer' in model.layers"
     "\n >>> obj = RhinoMesh(vertices, faces, texCoords, objectName, layerIndex)"
     "\n >>> model.objects.append(obj)"
     "\n >>>"
     "\n >>> # write result"
     "\n >>> model.write('MyRhinoFile.3dm')"
     "\n"
     );
  if (m == NULL)
    return;

  // add version string to module 
  if (PyModule_AddStringConstant(m, "__version__", module_version) != 0)
    return;

  // add Mesh class to module
  if (PyType_Ready(&MeshType) < 0)
    return;
  Py_INCREF(&MeshType);
  PyModule_AddObject(m, "Mesh", (PyObject *)&MeshType);

  // add PointCloud class to module
  if (PyType_Ready(&PointCloudType) < 0)
    return;
  Py_INCREF(&PointCloudType);
  PyModule_AddObject(m, "PointCloud", (PyObject *)&PointCloudType);

  // add RhinoModel class to module
  if (PyType_Ready(&RhinoModelType) < 0)
    return;
  Py_INCREF(&RhinoModelType);
  PyModule_AddObject(m, "RhinoModel", (PyObject *)&RhinoModelType);

}
