
// The pointcloud Python-object is of type GeomObj. Its m_geomObj attribute is
// a (ON_Geometry*)-pointer but actually points to a ON_PointCloud object.

static PyMemberDef PointCloud_members[] = {
  // {"userString",  T_OBJECT_EX, offsetof(GeomObj, m_userString), 0,
  //  "UserString-dictionary: Rhino.DocObjects.ObjectAttributes.GetUserString()"
  //  " or GetUserStrings()"},
  {NULL}  /* Sentinel */
};

/* PointCloud constructor: __new__()
 */
static PyObject *
PointCloud_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  GeomObj *self;

  self = (GeomObj *)type->tp_alloc(type, 0);
  if (self != NULL) {
    debug_print("Initialized new C-PointCloud.\n");
  }

  /* initialize mesh and attributes pointers */
  self->m_geomObj = NULL;
  self->m_attributes = NULL;

  return (PyObject *)self;
}

/* PointCloud constructor: __init__()
 */
static int
PointCloud_init(GeomObj *self, PyObject *args, PyObject *kwds)
{
  /* parse arguments */
  PyObject* points;
  PyObject* normals = NULL;
  const char *objectName = NULL;
  int layerIndex = 0;
  
  static char *kwlist[]
    = {"points", "normals", "objectName", "layerIndex", NULL};
  if (! PyArg_ParseTupleAndKeywords
      (args, kwds, "O|Osi:PointCloud_init", kwlist,
       &points, &normals, &objectName, &layerIndex)) {
    return -1;
  }

  if (!PyList_Check(points)) {
    PyErr_SetString(PyExc_ValueError,
       "PointCloud_init: points argument must be a list of float-triples.");
    return -1;
  }
  
  bool bHasNormals; // specify point normals?
  if (!normals) {
    bHasNormals = false;
  }
  else if (normals==Py_None) {
    bHasNormals = false;
  }
  else if (!PyList_Check(normals)) {
    PyErr_SetString(PyExc_ValueError,
       "PointCloud_init: normals argument must be a list of float-triples.");
    return -1;
  }
  else {
    bHasNormals = true;
  }

  Py_XINCREF(points);
  Py_XINCREF(normals);
  const int point_count = PyList_Size(points);

  /* discard old mesh object if exists */
  if (self->m_geomObj) {
    delete self->m_geomObj;
  }

  /* create a new ON_PointCloud object */
  self->m_geomObj = new ON_PointCloud(point_count);

  if (self->m_geomObj == NULL) {
    debug_print("Failed to initialize new C-PointCloud.\n");
    self->ob_type->tp_free((PyObject*)self);
    Py_XDECREF(points);
    Py_XDECREF(normals);
    return -1;
  }

  /* create new object attributes object */
  if (self->m_attributes) {
    /* discard old attribute object if exists */
    delete self->m_attributes;
  }
  self->m_attributes = new ON_3dmObjectAttributes;
  if (self->m_attributes == NULL) {
    debug_print("Failed to initialize new C-PointCloud-attributes.\n");
    delete self->m_geomObj;
    self->m_geomObj = NULL;
    self->ob_type->tp_free((PyObject*)self);
    Py_XDECREF(points);
    Py_XDECREF(normals);
    return -1;
  }

  /* update the ON_PointCloud object with provided data arguments */
  bool ok = true;

  /* create points */
  for (int i=0; ok && (i<point_count); ++i) {
    double x[3];
    PyObject* coords;
    
    coords = PyList_GET_ITEM(points, i);
    if (PyList_Check(coords)) {
      for (int j=0; j<3; ++j) {
	x[j] = PyFloat_AsDouble(PyList_GET_ITEM(coords, j));
      }
    }
    else if (PyTuple_Check(coords)) {
      for (int j=0; j<3; ++j) {
	x[j] = PyFloat_AsDouble(PyTuple_GET_ITEM(coords, j));
      }
    }
    else {
      ok = false;
      PyErr_SetString
	(PyExc_ValueError,
	 "PointCloud_init: First argument (points) must be a list of lists or"
	 " tuples.");
      break;
    }
    ((ON_PointCloud*)self->m_geomObj)->AppendPoint( ON_3dPoint(x[0], x[1], x[2]) );
  }

  /* add normals */
  if (ok && bHasNormals) {
    if (PyList_Size(normals) != point_count) {
      ok = false;
      PyErr_Format
	(PyExc_ValueError,
	 "PointCloud_init: Nb of normals (%d) not equal to nb of"
	 " points (%d).", PyList_Size(normals), point_count);
    }

    for (int i=0; ok && (i<point_count); ++i) {
      const int dim = 3;
      double x[dim];
      PyObject* coords;
    
      coords = PyList_GET_ITEM(normals, i);
      if (PyList_Check(coords)) {
	for (int j=0; j<dim; ++j) {
	  x[j] = PyFloat_AsDouble(PyList_GET_ITEM(coords, j));
	}
      }
      else if (PyTuple_Check(coords)) {
	for (int j=0; j<dim; ++j) {
	  x[j] = PyFloat_AsDouble(PyTuple_GET_ITEM(coords, j));
	}
      }
      else {
	ok = false;
	PyErr_SetString
	  (PyExc_ValueError,
	   "PointCloud_init: normals argument must be a list of lists or"
	   " tuples.");
	break;
      }
      // vert-id, 3 normal coords
      ((ON_PointCloud*)self->m_geomObj)->m_N.Append(ON_3dVector(x[0], x[1], x[2]));
    }
  }

  /* update the ON_3dmObjectAttributes object with provided data arguments */
  if (ok && objectName) {
    self->m_attributes->m_name = objectName;
  }
  if (ok) {
    self->m_attributes->m_layer_index = layerIndex;
    // assign a uuid to the attributes
    ON_CreateUuid(self->m_attributes->m_uuid);
  }
  #ifdef DEBUG
  ON_wString uuidStr;
  ON_UuidToString(self->m_attributes->m_uuid, uuidStr);
  debug_print("attrib-uuid: %s\n", (const char*) ((ON_String) uuidStr));
  #endif
  
  Py_XDECREF(points);
  Py_XDECREF(normals);
  if (!ok) {
    delete self->m_attributes;
    self->m_attributes = NULL;
    delete self->m_geomObj;
    self->m_geomObj = NULL;
    self->ob_type->tp_free((PyObject*)self);
    return -1;
  }
  return 0;
}

/* PointCloud destructor
 */
static void
PointCloud_dealloc(GeomObj* self)
{
  debug_print("In PointCloud_dealloc (1).\n");
  if (self->m_geomObj != NULL) {
    delete self->m_geomObj;
    self->m_geomObj = NULL;
  }
  debug_print("In PointCloud_dealloc (2).\n");
  if (self->m_attributes != NULL) {
    delete self->m_attributes;
    self->m_attributes = NULL;
  }
  self->ob_type->tp_free((PyObject*)self);
  debug_print("In PointCloud_dealloc (3).\n");
}


/* PointCloud methods
 */
static PyMethodDef PointCloud_methods[] = {
  {"setUserString",
   (PyCFunction)GeomObj_setUserString,  METH_VARARGS | METH_KEYWORDS,
   "Set Rhino object user string."
   "\nArguments: key, value (both strings)"
   "\n"
  },
  {NULL}  /* Sentinel */
};

/* PointCloud attribute get-setters
 */
static PyGetSetDef PointCloud_getseters[] = {
    {"objectName",
     (getter)GeomObj_get_objectName, (setter)GeomObj_set_objectName,
     "objectName attribute of the Rhino mesh object", NULL},
    {"objectLayerIndex",
     (getter)GeomObj_get_objectLayerIndex, (setter)GeomObj_set_objectLayerIndex,
     "Zero based index in the layertable of the RhinoModel object", NULL},
    {NULL}  /* Sentinel */
};

/* PointCloud type description
 */
static PyTypeObject PointCloudType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "rhino_binary.PointCloud", /*tp_name*/
    sizeof(GeomObj),           /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)PointCloud_dealloc,  /*tp_dealloc*/
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
    "PointCloud class",        /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    PointCloud_methods,        /* tp_methods */
    PointCloud_members,        /* tp_members */
    PointCloud_getseters,      /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)PointCloud_init, /* tp_init */
    0,                         /* tp_alloc */
    PointCloud_new,            /* tp_new */
};
