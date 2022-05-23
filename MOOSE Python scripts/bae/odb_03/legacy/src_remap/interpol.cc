#include "interpol.h"



// a map unsigned int -> constant length vector.
//  The vectors contain T-elements and are all of the same size N
template< class T, unsigned int N >
class IntCLVectDict {
public:

  //--- constructors and so on
  IntCLVectDict() : keys(0) {};
  virtual ~IntCLVectDict() {
    if (keys) delete keys;
  };

  //--- Initialization from two vectors one for the keys one for the values.

  // Read the key vector from inFile. nbItems ... number of items to be read.
  IntCLVectDict& loadKeysFromBinaryFile
  (std::ifstream& inFile, unsigned int nbItems) {
    if (keys) delete keys;
    keys = new std::vector<unsigned int>(nbItems);
    inFile.read(reinterpret_cast<char*>(keys->data()),
		nbItems*sizeof(unsigned int));
    return *this;
  }

  // Read the value vector from inFile. nbItems ... number of items to be read.
  IntCLVectDict& loadValuesFromBinaryFile(std::ifstream& inFile,
					  unsigned int nbItems) {
    values.clear();
    values.resize(nbItems*N);
    inFile.read(reinterpret_cast<char*>(values.data()), nbItems*N*sizeof(T));
    return *this;
  }

  // initialize the map 
  IntCLVectDict& initializeFromStore() {

    unsigned int nbItems;

    if (!keys) {
      cerr << "IntCLVectDict can't initializeFromStore: no keys." << endl;
      exit(2);
    }

    nbItems = keys->size();
    if (nbItems != values.size()/N) {
      cerr << "IntCLVectDict can't initializeFromStore:"
	   << " inconsistent size of keys (" << nbItems << ") and values ("
	   << values.size()/N << ")" << endl;
      exit(2);
    }

    for (unsigned int i=0; i<nbItems; ++i) {
      indexMap[(*keys)[i]] = N*i;
    }

    delete keys;
    keys = 0;
    
    return *this;
  }

  //--- getting, setting individual values

  // for getting a value
  T* operator[] (const unsigned int label) const {
    return values->data()+indexMap[label];
  }

  // for setting a value
  const T* operator[] (const unsigned int label) {
    
    std::map<unsigned int, unsigned int>::const_iterator it;
    unsigned int newIdx;

    it = indexMap.find(label);

    if (it == indexMap.end()) {
      // not found -> create new
      newIdx = values.size();
      indexMap[label] = newIdx;
      values.reserve(newIdx+N);
      return values.data()+newIdx;
    }
    else {
      // item found -> overwrite
      return values.data()+(it->second);
    }
  }

  unsigned int size() const {
    return indexMap.size();
  }

  // this is just for debugging ...
  unsigned int firstKey() const {
    return indexMap.begin()->first;
  }

protected:
  std::map<unsigned int, unsigned int> indexMap;
  std::vector<T> values;
  std::vector<unsigned int>* keys;
  
};



// a map unsigned int -> variable length vector.
//  The vectors contain T-elements of potentially different sizes
template< class T >
class IntVLVectDict {
public:

  //--- constructors and so on
  IntVLVectDict() : keys(0), values(0) {};
  ~IntVLVectDict() {
    if (keys) delete keys;
    if (values) delete values;
  };

  //--- Initialization from two vectors one for the keys one for the values.

  // Read the key vector from inFile. nbItems ... number of items to be read.
  IntVLVectDict& loadKeysFromBinaryFile(std::ifstream& inFile,
					unsigned int nbItems) {
    if (keys) delete keys;
    keys = new std::vector<unsigned int>(nbItems);
    inFile.read(reinterpret_cast<char*>(keys->data()),
		nbItems*sizeof(unsigned int));
    return *this;
  }

  // Read the value vector from inFile.
  // vecSize ... size of all vectors to be stored during this session 
  // nbItems ... number of items to be read.
  IntVLVectDict& loadValuesFromBinaryFile
  ( std::ifstream& inFile, unsigned int nbItems, unsigned int vecSize ) {
    if (values) delete values;
    values = new std::vector< std::vector<T> >
      (nbItems, std::vector<T>(vecSize));
    for (unsigned int i; i<nbItems && inFile.good(); ++i) {
      void* data = (*values)[i].data();
      inFile.read(reinterpret_cast<char*>(data), vecSize*sizeof(T));
    }
    return *this;
  }

  // initialize the map 
  IntVLVectDict& updateFromStore() {

    unsigned int nbItems;

    if (!keys) {
      cerr << "IntVLVectDict can't updateFromStore: no keys." << endl;
      exit(2);
    }
    if (!values) {
      cerr << "IntVLVectDict can't updateFromStore: no values." << endl;
      exit(2);
    }

    nbItems = keys->size();
    if (nbItems != values->size()) {
      cerr << "IntVLVectDict can't updateFromStore:"
	   << " inconsistent size of keys (" << nbItems << ") and values ("
	   << values->size() << ")" << endl;
      exit(2);
    }

    // might not be most efficient because vector is copied...
    // but it's possible to understand
    for (unsigned int i=0; i<nbItems; ++i) {
      map[(*keys)[i]] = (*values)[i];
    }

    delete keys;    keys = 0;
    delete values;  values = 0;
    
    return *this;
  }

  //--- getting, setting individual values

  // for getting (and setting -I think-) a value
  std::vector<T>& operator[] (const unsigned int label) {
    return map[label];
  }

  unsigned int size() const {
    return map.size();
  }

  unsigned int firstKey() const {
    return map.begin()->first;
  }

protected:
// public:
  std::map<unsigned int, std::vector<T> > map;
  std::vector< std::vector<T> > *values;
  std::vector<unsigned int> *keys;
  
};



typedef IntCLVectDict<float, 3> NodeCoordsType;
typedef IntVLVectDict<unsigned int> ElNodesType;



// ---------
// --- F --- Reading interpolation data from file
// ---------


// class MeshStructuredPoints function bodies

unsigned int MeshStructuredPoints::size() const {
  return dims[0]*dims[1]*dims[2];
}

void MeshStructuredPoints::readFromFile(std::ifstream& inFile) {

  std::string line;  // read buffer for one line
  std::string token;
  std::istringstream lineStr;
  int numGridPoints;

  // DIMENSIONS nx ny nz
  std::getline(inFile, line);
  lineStr.str(line); lineStr.clear();
  lineStr >> token >> dims[0] >> dims[1] >> dims[2];
  numGridPoints=size();
  if (! (lineStr and inFile)) {
    cerr << "Error reading DIMENSIONS line: " << line << endl;
    exit(2);
  }
  if (token != "DIMENSIONS") {
    cerr << "Expected DIMENSIONS but found " << token << endl;
    exit(2);
  }
    
  // ORIGIN x y z
  std::getline(inFile, line);
  lineStr.str(line); lineStr.clear();
  lineStr >> token >> org[0] >> org[1] >> org[2];
  if (! (lineStr and inFile)) {
    cerr << "Error reading ORIGIN line: " << line << endl;
    exit(2);
  }
  if (token != "ORIGIN") {
    cerr << "Expected ORIGIN but found " << token << endl;
    exit(2);
  }

  // SPACING sx sy sz
  std::getline(inFile, line);
  lineStr.str(line); lineStr.clear();
  lineStr >> token >> spc[0] >> spc[1] >> spc[2];
  if (! (lineStr and inFile)) {
    cerr << "Error reading SPACING line: " << line << endl;
    exit(2);
  }
  if (token != "SPACING") {
    cerr << "Expected SPACING but found " << token << endl;
    exit(2);
  }
}

class MeshUnstructuredGrid : public MeshBaseType {
public:
  NodeCoordsType nodeCoords;
  ElNodesType elNodes;

  unsigned int size() const {
    return elNodes.size();
  }
  
  void readFromFile(std::ifstream& inFile) {

    std::string line;  // read buffer for one line
    std::string token;
    std::istringstream lineStr;
    unsigned int nbNodes, nbNodes2, nbShapes, nbElems, nbElems2;

    //--- load nodeCoords

    // POINTS nbNodes
    std::getline(inFile, line);
    lineStr.str(line); lineStr.clear();
    lineStr >> token >> nbNodes;
    if (token != "POINTS") {
      cerr << "Expected POINTS but found " << token << endl;
      exit(2);
    }
    nodeCoords.loadValuesFromBinaryFile(inFile, nbNodes);

    // NODE_NUMBERS nbNodes
    std::getline(inFile, line);
    lineStr.str(line); lineStr.clear();
    lineStr >> token >> nbNodes2;
    if (token != "NODE_NUMBERS") {
      cerr << "Expected NODE_NUMBERS but found " << token << endl;
      exit(2);
    }
    if (nbNodes!=nbNodes2) {
      cerr << "NODE_NUMBERS cnt " << nbNodes2
		<< " differs from POINTS cnt " << nbNodes << endl;
      exit(2);
    }
    nodeCoords.loadKeysFromBinaryFile(inFile, nbNodes);
    nodeCoords.initializeFromStore();

    //--- load elNodes

    // UNSTRUCTURED_MESH <nbShapes> shapes
    std::getline(inFile, line);
    lineStr.str(line); lineStr.clear();
    lineStr >> token >> nbShapes;  // ignore "shapes" at the end...
    if (token != "UNSTRUCTURED_MESH") {
      cerr << "Expected UNSTRUCTURED_MESH but found " << token << endl;
      exit(2);
    }

    for (unsigned int i=0; i<nbShapes; ++i) {

      // ELSHAPE shape
      std::getline(inFile, line);
      lineStr.str(line); lineStr.clear();
      lineStr >> token;  // ignore shape at the end...
      if (token != "ELSHAPE") {
        cerr << "Expected ELSHAPE but found " << token << endl;
        exit(2);
      }

      // ELNODES <nbNodes> <nbElems>
      // nbNodes ... nodes per element
      std::getline(inFile, line);
      lineStr.str(line); lineStr.clear();
      lineStr >> token >> nbNodes >> nbElems;
      if (token != "ELNODES") {
        cerr << "Expected ELNODES but found " << token << endl;
        exit(2);
      }
      elNodes.loadValuesFromBinaryFile(inFile, nbElems, nbNodes);
      if (!inFile.good()) {
	cerr << "Error while reading ELNODES-values pos: " << inFile.tellg()
	     << endl;
	exit(2);
      }

      //ELEMENT_NUMBERS nbElems
      std::getline(inFile, line);
      lineStr.str(line); lineStr.clear();
      lineStr >> token >> nbElems2;
      if (token != "ELEMENT_NUMBERS") {
        cerr << "Expected ELEMENT_NUMBERS but found " << token << endl;
        exit(2);
      }
      elNodes.loadKeysFromBinaryFile(inFile, nbElems);
      if (!inFile.good()) {
	cerr << "Error while reading ELEMENT_NUMBERS pos: " << inFile.tellg()
	     << endl;
	exit(2);
      }
      elNodes.updateFromStore();
    }
  }

  // This function returns a vector of the values of the shape function
  // associated with the ten nodes of the quad tet at the given point.
  // The point is identified by its element coordinates (i.e. the vector
  // [1-xi-eta-zeta, xi, eta, zeta].
  // Note: the output argument must already be a vector of size 10.
  // Note: for historical reasons this formula only considers [xi,eta,zeta].
  static const float _tetQ_sf0coeff[10];
  void getElemShapeFuncs_tetQ
  ( std::vector<float>& shapeFuncValues,   // output
    float* elemCoords) {
    std::vector<float> xxx(10);

    // 1.0, xi, eta, zeta, xi**2, eta**2, zeta**2, xi*eta, xi*zeta, eta*zeta
    xxx[0] = 1.0;
    xxx[1] = elemCoords[1];
    xxx[2] = elemCoords[2];
    xxx[3] = elemCoords[3];
    xxx[4] = elemCoords[1]*elemCoords[1];
    xxx[5] = elemCoords[2]*elemCoords[2];
    xxx[6] = elemCoords[3]*elemCoords[3];
    xxx[7] = elemCoords[1]*elemCoords[2];
    xxx[8] = elemCoords[1]*elemCoords[3];
    xxx[9] = elemCoords[2]*elemCoords[3];

    shapeFuncValues[0] = 0.0;
    for (int i=0; i<10; ++i) {
      shapeFuncValues[0] += _tetQ_sf0coeff[i]*xxx[i];
    }
    shapeFuncValues[1] = 2.0*xxx[4] - xxx[1];
    shapeFuncValues[2] = 2.0*xxx[5] - xxx[2];
    shapeFuncValues[3] = 2.0*xxx[6] - xxx[3];
    shapeFuncValues[4] = 4.0*xxx[1] - 4.0*xxx[4] - 4.0*xxx[7] - 4.0*xxx[8];
    shapeFuncValues[5] = 4.0*xxx[7];
    shapeFuncValues[6] = 4.0*xxx[2] - 4.0*xxx[5] - 4.0*xxx[7] - 4.0*xxx[9];
    shapeFuncValues[7] = 4.0*xxx[3] - 4.0*xxx[6] - 4.0*xxx[8] - 4.0*xxx[9];
    shapeFuncValues[8] = 4.0*xxx[8];
    shapeFuncValues[9] = 4.0*xxx[9];
    return;
  }

  // This function returns a vector of the values of the interpolating
  // functions associated with the four integration points of the
  // quad tet element at the given point.
  // The point is identified by its element coordinatese.
  // (i.e. the vector [1-xi-eta-zeta, xi, eta, zeta])
  // Note: the output argument must already be a vector of size 4.
  static const float _gpQTetInterpolation[16];
  void getElemIPInterpolFuncs_tetQ
  ( std::vector<float>& interpolFuncValues,   // output
    float* elemCoords) {

    for (int i=0; i<4; ++i) {
      interpolFuncValues[i] = 0.0;
      for (int j=0; j<4; ++j) {
	interpolFuncValues[i] += _gpQTetInterpolation[i*4+j]*elemCoords[j];
      }
    }
    return;
  }
};
const float MeshUnstructuredGrid::_tetQ_sf0coeff[10] =
  {1., -3., -3., -3., 2., 2., 2., 4., 4., 4.};
// gpTetC = (1+3*sqrt(5))/4.
static const float _gpTetC = 1.927050983124842272306880251548457176581;
// gpTetD = (1-sqrt(5))/4.
static const float _gpTetD = -0.3090169943749474241022934171828190588603;
//
const float MeshUnstructuredGrid::_gpQTetInterpolation[16] =
  { _gpTetC,_gpTetD,_gpTetD,_gpTetD,
    _gpTetD,_gpTetC,_gpTetD,_gpTetD,
    _gpTetD,_gpTetD,_gpTetC,_gpTetD,
    _gpTetD,_gpTetD,_gpTetD,_gpTetC };


class Interpolation {
public:
  MeshStructuredPoints* grid;
  MeshUnstructuredGrid* mesh;

  unsigned int nbCoords;  // coords per element
  std::vector<unsigned int> elemNbs;
  std::vector<float> elemCoords;

  Interpolation& readFromPickleFile(const std::string &fileNameInterpol) {

    cout << "Reading interpolation data from file " << fileNameInterpol
	 << " ." << endl;

    // processing lines of file:
    std::ifstream inFile(fileNameInterpol.c_str());
    ASSERT(inFile); // make sure file exists
    std::string line;  // read buffer for one line

    // verify interpol-file format
    std::getline(inFile, line);
    if (line!="Interpolation Data v2.0: binary grid,mesh,elemCoords") {
      cerr << "Unrecognized file format: " << line << endl;
      exit(2);
    }

    // read grid data
    std::getline(inFile, line);
    if (line=="STRUCTURED_POINTS") {
      grid = new MeshStructuredPoints;
      grid -> readFromFile(inFile);
    }
    else {
      // currently we don't do unstruct points
      cerr << "Only STRUCTURED_POINTS for the grid at the moment, no "
	   << line << endl;
      exit(2);
    }
  
    // read mesh data
    mesh = new MeshUnstructuredGrid;
    mesh -> readFromFile(inFile);

    // read elemCoords
    unsigned int nbItems;
    inFile.read(reinterpret_cast<char*>(&nbItems), sizeof(unsigned int));
    inFile.read(reinterpret_cast<char*>(&nbCoords), sizeof(unsigned int));

    if (nbItems != grid->size()) {
      cerr << "ERROR: Nb of elemCoords (" << nbItems << ") inconsistent with"
	   << "nb of grid points (" << grid->size() << ")." << endl;
      exit(2);
    }

    // read element numbers
    elemNbs.clear();
    elemNbs.resize(nbItems);
    inFile.read(reinterpret_cast<char*>(elemNbs.data()),
		nbItems*sizeof(unsigned int));

    // read element coordinates
    elemCoords.clear();
    elemCoords.resize(nbItems*nbCoords);
    inFile.read(reinterpret_cast<char*>(elemCoords.data()),
		nbItems*nbCoords*sizeof(float));

    inFile.close();

    return *this;
  }


  // 1) Nodal vars, such as U or V:
  // Each gridpoint is inside exactly one element, and each
  // element has exactly ten nodes with one weight each.
  //
public: // protected:
  // vector of 10 x numGridPoints nodes:
  std::vector< std::vector<unsigned int> > node_of_pt;  
  // 10 weights, one for each of the ten nodes:
  //      first gridpoint     second gridpoint    ...   last gridpoint
  // [ [ 0.24 0.38 0.55 ...][ 0.41 0.25 0.51 ...] ... [       ...       ] ]
  std::vector< std::vector< float > > n_wght_of_pt;
public:
  // set of all nodes and elements involved in the interpolation
  std::set<unsigned int> allNodes;
  std::set<unsigned int> allElems;

public:
  void initNodalWeights() {
    // vector of 10 floats (fixed size):
    std::vector<float> n_wghts(10);
    std::vector<unsigned int> nodesEmpty(10); // will be 10 zeros

    allNodes.clear();
    for (int i=0; i<grid->size(); ++i) {      

      unsigned int elem = elemNbs[i];
      // add node list at the end, update allNodes
      if (elem>0) {
	std::vector<unsigned int>& nodes = mesh->elNodes[elem];
	node_of_pt.push_back(nodes);
	// store nodes in allNodes
	allNodes.insert(nodes.begin(), nodes.end());
      }
      else {
	node_of_pt.push_back(nodesEmpty);
      }

      mesh->getElemShapeFuncs_tetQ(n_wghts, elemCoords.data()+i*nbCoords);
      // add weight list at the end:
      n_wght_of_pt.push_back(n_wghts);
    }
  }

  // 2) Element vars, such as S or S_MIN_PRINCIPAL or S_MAX_PRINCIPAL:
  // Each gridpoint is inside exactly one element, and
  // each element has exactly 4 integration points
  // with one weight each:

public: // protected:
  // vector of numGridPoints integers:
  // [ 278, 2256, ... ]
  std::vector< int > elem_of_pt;
  // vector of numGridPoints vectors of four weight factors, one for each IP:
  // [ [0.1, 0.3, 0.2, 0.4], [0.2, 0.3, 0.2, 0.3], ... ]
  std::vector< std::vector< float > > e_wght_of_pt;

public:
  void initElemIpWeights() {

    // vector of 4 floats (fixed size):
    std::vector< float > ip_wghts(4);

    for (int i=0; i<grid->size(); ++i) {      
      unsigned int elem = elemNbs[i];
      mesh->getElemIPInterpolFuncs_tetQ(ip_wghts, elemCoords.data()+i*nbCoords);
      

      // store elem numbers and weights in elem_of_pt and e_wght_of_pt
      elem_of_pt.push_back(elem);
      e_wght_of_pt.push_back(ip_wghts);

      // store elem in allElems
      if (elem>0) {
	allElems.insert(elem);
      }
    }

  }

  // diagnostic information
public:
  
  void printDiagnosticOutput() {
    cout << "     --> Interpolation data:" << endl;
    cout << "         Grid with " << grid->size() << " points." << endl
	 << "         . origin: " << toStr(grid->org) << endl
	 << "         . resolution: " << toStr(grid->dims) << endl
	 << "         . spacing: " << toStr(grid->spc) << endl
	 << "         Touching " << allElems.size()
	 << " elements and " << allNodes.size() << " nodes." << endl;
  }

  void test() {

    cout << grid->size() << " points in grid." << endl;

    cout << mesh->nodeCoords.size() << " nodes have coordinates." << endl;
    int node = mesh->nodeCoords.firstKey();
    const float* coords = mesh->nodeCoords[node];
    cout << "E.g. node " << node << " has coords "
      	 << coords[0] << "," << coords[1] << "," << coords[2] << endl;

    cout << mesh->elNodes.size() << " elems have nodes." << endl;
    int elem = mesh->elNodes.firstKey();
    std::vector<unsigned int> nodes;

    nodes = mesh->elNodes[elem];
    cout << "E.g. elem " << elem << " has nodes ";
    for (int i=0; i<nodes.size(); ++i) {
      cout << nodes[i] << ",";
    }
    cout << endl;

  }
};

