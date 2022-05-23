
// The mesh Python-object is of type GeomObj. Its m_geomObj attribute is
// a (ON_Geometry*)-pointer but actually points to a ON_Mesh object.

/* data structure for the Mesh object
 */
class Mesh {
public:
  PyObject_HEAD
  ON_Mesh *m_geomObj;
  ON_3dmObjectAttributes *m_attributes;
};

static PyMemberDef Mesh_members[] = {
  // {"userString",  T_OBJECT_EX, offsetof(GeomObj, m_userString), 0,
  //  "UserString-dictionary: Rhino.DocObjects.ObjectAttributes.GetUserString()"
  //  " or GetUserStrings()"},
  {NULL}  /* Sentinel */
};

/* Mesh constructor: __new__()
 */
static PyObject *
Mesh_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  GeomObj *self;

  self = (GeomObj *)type->tp_alloc(type, 0);
  if (self != NULL) {
    debug_print("Initialized new C-Mesh.\n");
  }

  /* initialize mesh and attributes pointers */
  self->m_geomObj = NULL;
  self->m_attributes = NULL;

  return (PyObject *)self;
}

/* Mesh constructor: __init__()
 */
static int
Mesh_init(GeomObj *self, PyObject *args, PyObject *kwds)
{
  /* parse arguments */
  PyObject* vertices;
  PyObject* faces;
  PyObject* texCoords = NULL;
  const char *objectName = NULL;
  int layerIndex = 0;
  
  static char *kwlist[]
    = {"vertices", "faces", "texCoords", "objectName", "layerIndex", NULL};
  if (! PyArg_ParseTupleAndKeywords
      (args, kwds, "OO|Osi:Mesh_init", kwlist,
       &vertices, &faces, &texCoords, &objectName, &layerIndex))
    return -1;

  bool bHasVertexNormals = false; // we will NOT specify vertex normals

  if (!PyList_Check(vertices)) {
    PyErr_SetString(PyExc_ValueError,
       "Mesh_init: First argument must be a list of vertex coords.");
    return -1;
  }
  if (!PyList_Check(faces)) {
    PyErr_SetString(PyExc_ValueError,
       "Mesh_init: Second argument must be a list of faces.");
    return -1;
  }

  bool bHasTexCoords; // specify texture coordinates?
  if (!texCoords) {
    bHasTexCoords = false;
  }
  else if (texCoords==Py_None) {
    bHasTexCoords = false;
  }
  else if (!PyList_Check(texCoords)) {
    PyErr_SetString(PyExc_ValueError,
       "Mesh_init: texCoords argument must be a list of 2-tuples.");
    return -1;
  }
  else {
    bHasTexCoords = true;
  }

  Py_XINCREF(vertices);
  Py_XINCREF(faces);
  Py_XINCREF(texCoords);
  const int vertex_count = PyList_Size(vertices);
  const int face_count = PyList_Size(faces);

  /* discard old mesh object if exists */
  if (self->m_geomObj) {
    delete self->m_geomObj;
  }

  /* create a new ON_Mesh object */
  self->m_geomObj = (ON_Geometry*) new ON_Mesh
    (face_count, vertex_count, bHasVertexNormals, bHasTexCoords);

  if (self->m_geomObj == NULL) {
    debug_print("Failed to initialize new C-Mesh.\n");
    self->ob_type->tp_free((PyObject*)self);
    Py_XDECREF(vertices);
    Py_XDECREF(faces);
    Py_XDECREF(texCoords);
    return -1;
  }

  /* create new object attributes object */
  if (self->m_attributes) {
    /* discard old attribute object if exists */
    delete self->m_attributes;
  }
  self->m_attributes = new ON_3dmObjectAttributes;
  if (self->m_attributes == NULL) {
    debug_print("Failed to initialize new C-mesh-attributes.\n");
    delete self->m_geomObj;
    self->m_geomObj = NULL;
    self->ob_type->tp_free((PyObject*)self);
    Py_XDECREF(vertices);
    Py_XDECREF(faces);
    Py_XDECREF(texCoords);
    return -1;
  }

  /* update the ON_Mesh object with provided data arguments */
  bool ok = true;

  /* create vertices */
  for (int i=0; ok && (i<vertex_count); ++i) {
    double x[3];
    PyObject* coords;
    
    coords = PyList_GET_ITEM(vertices, i);
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
	 "Mesh_init: First argument (vertices) must be a list of lists or"
	 " tuples.");
      break;
    }
    ok = ((ON_Mesh*)self->m_geomObj)->SetVertex( i, ON_3dPoint(x[0], x[1], x[2]) );
    if (!ok) {
      PyErr_Format
	(PyExc_ValueError,
	 "Mesh_init: Could not create vertex nb %d", i);
    }
  }

  /* create faces */
  for (int i=0; ok && (i<face_count); ++i) {
    long v[3];
    PyObject* verts;
    PyObject* oneInt;
    
    verts = PyList_GET_ITEM(faces, i);
    if (PyList_Check(verts)) {
      for (int j=0; j<3; ++j) {
	oneInt = PyList_GET_ITEM(verts, j);
	if (PyInt_Check(oneInt)) {
	  v[j] = PyInt_AS_LONG(oneInt);
	}
	else {
	  ok = false;
	  PyErr_Format
	    (PyExc_ValueError,
	     "Mesh_init: Faces should be lists or tuples of integers."
	     " Could not convert vertex id %d of face %d to int.", j, i);
	  break;
	}
      }
    }
    else if (PyTuple_Check(verts)) {
      for (int j=0; j<3; ++j) {
	oneInt = PyTuple_GET_ITEM(verts, j);
	if (PyInt_Check(oneInt)) {
	  v[j] = PyInt_AS_LONG(oneInt);
	}
	else {
	  ok = false;
	  PyErr_Format
	    (PyExc_ValueError,
	     "Mesh_init: Faces should be lists or tuples of integers."
	     " Could not convert vertex id %d of face %d to int.", j, i);
	  break;
	}
      }
    }
    else {
      ok = false;
      PyErr_SetString
	(PyExc_ValueError,
	 "Mesh_init: Second argument (faces) must be a list of lists or"
	 " tuples.");
    }
    if (ok) {
      ok = ((ON_Mesh*)self->m_geomObj)->SetTriangle( i, v[0], v[1], v[2] );
      if (!ok) {
	PyErr_Format
	  (PyExc_ValueError,
	   "Mesh_init: Could not create face nb %d", i);
      }
    }
  }

  /* add texture coordinates */
  if (ok && bHasTexCoords) {
    if (PyList_Size(texCoords) != vertex_count) {
      ok = false;
      PyErr_Format
	(PyExc_ValueError,
	 "Mesh_init: Nb of texture coordinates (%d) not equal to nb of"
	 " vertices (%d).", PyList_Size(texCoords), vertex_count);
    }

    for (int i=0; ok && (i<vertex_count); ++i) {
      double x[2];
      PyObject* coords;
    
      coords = PyList_GET_ITEM(texCoords, i);
      if (PyList_Check(coords)) {
	for (int j=0; j<2; ++j) {
	  x[j] = PyFloat_AsDouble(PyList_GET_ITEM(coords, j));
	}
      }
      else if (PyTuple_Check(coords)) {
	for (int j=0; j<2; ++j) {
	  x[j] = PyFloat_AsDouble(PyTuple_GET_ITEM(coords, j));
	}
      }
      else {
	ok = false;
	PyErr_SetString
	  (PyExc_ValueError,
	   "Mesh_init: texCoords argument must be a list of lists or"
	   " tuples.");
	break;
      }
      // vert-id, 2 texture coords
      ok = ((ON_Mesh*)self->m_geomObj)->SetTextureCoord( i, x[0], x[1] );
      if (!ok) {
	PyErr_Format
	  (PyExc_ValueError,
	   "Mesh_init: Could not create texture coords for vertex nb %d", i);
      }
    }
  }

  // Most applications expect vertex normals.
  // If they are not present, ComputeVertexNormals sets
  // them by averaging face normals.
  if ( ok && ((ON_Mesh*)self->m_geomObj)->IsValid()
       && !((ON_Mesh*)self->m_geomObj)->HasVertexNormals() ) {
    ((ON_Mesh*)self->m_geomObj)->ComputeVertexNormals();
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
  
  Py_XDECREF(vertices);
  Py_XDECREF(faces);
  Py_XDECREF(texCoords);
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

/* Mesh destructor
 */
static void
Mesh_dealloc(GeomObj* self)
{
  debug_print("In Mesh_dealloc (1).\n");
  if (self->m_geomObj != NULL) {
    delete self->m_geomObj;
    self->m_geomObj = NULL;
  }
  debug_print("In Mesh_dealloc (2).\n");
  if (self->m_attributes != NULL) {
    delete self->m_attributes;
    self->m_attributes = NULL;
  }
  self->ob_type->tp_free((PyObject*)self);
  debug_print("In Mesh_dealloc (3).\n");
}

/* Mesh methods
 */
static PyMethodDef Mesh_methods[] = {
  {"setUserString",
   (PyCFunction)GeomObj_setUserString,  METH_VARARGS | METH_KEYWORDS,
   "Set Rhino object user string."
   "\nArguments: key, value (both strings)"
   "\n"
  },
  {NULL}  /* Sentinel */
};

/* Mesh attribute get-setters
 */
static PyGetSetDef Mesh_getseters[] = {
    {"objectName",
     (getter)GeomObj_get_objectName, (setter)GeomObj_set_objectName,
     "objectName attribute of the Rhino mesh object", NULL},
    {"objectLayerIndex",
     (getter)GeomObj_get_objectLayerIndex, (setter)GeomObj_set_objectLayerIndex,
     "Zero based index in the layertable of the RhinoModel object", NULL},
    {NULL}  /* Sentinel */
};

/* Mesh type description
 */
static PyTypeObject MeshType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "rhino_binary.Mesh",       /*tp_name*/
    sizeof(GeomObj),           /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)Mesh_dealloc,  /*tp_dealloc*/
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
    "Mesh class",              /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    Mesh_methods,              /* tp_methods */
    Mesh_members,              /* tp_members */
    Mesh_getseters,            /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Mesh_init,       /* tp_init */
    0,                         /* tp_alloc */
    Mesh_new,                  /* tp_new */
};

/******************************************************************************
*******************************************************************************
***  more mesh methods to come (?).....                                     ***
*******************************************************************************
******************************************************************************/
