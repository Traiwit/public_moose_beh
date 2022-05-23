/* 
doc_usage_start

start the executable either with abaqus:
$ abaqus remap ...

or set appropriate environment variable:
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/6.13-2/code/bin
$ ./remap ...


supply arguments to remap on the command line:
$ abaqus remap  test8.odb  box1.interpol S_11,U,S_22,SDV1,U_2  ascii  10,20,5,100.,0.,0.,1.,2.,4.  Step-2,99,prj_99.vtk  Step-2,100,prj_100.vtk
$ abaqus remap  test8.odb  box1.interpol SDV1,S_MIN_PRINCIPAL,U_3  ascii  10,20,5,100.,0.,0.,1.,2.,4.  Step-2,99,prj_99.vtk  Step-2,77,prj_77.vtk
$ abaqus remap  test8.odb  box1.interpol U_1,S_11,SDV12  binary  100,200,50,100.,0.,0.,1.,2.,4.  Step-2,99,prj_99.vtk  Step-2,77,prj_77.vtk


or use parameter file:
$ cat >para_myproj.txt
/mnt/boab/data3/abaqus/telfer2014/.../P201445_R10_G04_S01_Q05_M05.odb
/mnt/boab/data3/abaqus/telfer2014/.../P201445_G04_box4m5.interpol
S_MIN_PRINCIPAL,S_11,U_1,U_2,U_3,SDV4,SDV12
binary
Step-2,55,box4m5FromRemap_F055.vtk
Step-2,60,box4m5FromRemap_F060.vtk
Step-2,65,box4m5FromRemap_F065.vtk

$ abaqus remap para_myproj.txt


or use stdin (from tty):
$ abaqus remap

1: Processing input from stdin.
Path to odb: test8.odb
Path to interpolation-file: box1.interpol
Requested variables (comma separated, no spaces): U_1,U_2,U_3,SDV4
ascii/binary: binary
...


or pipe a file to stdin:
$ abaqus remap <para_myproj.txt


doc_usage_end


doc_desc_start


Program Algorithm
-----------------

 0: Declaration section.

    A: Variables for handling input data
    B: Output variables
    C: Other

 1: Argument Processing

    A: Find out how to read arguments
    B: Process/read arguments    

 2: Opening odb, Reading weights.

    A: Opening odb.
    B: Read interpolation data from file

 3: Read, interpolate, write vtk

    A: Reading from odb    
    B: Interpolating
    C: Writing vtk

doc_desc_end

doc_make_start

Compiling/Linking 
-----------------

1) Compile-command (if this file is called remap.cc): 

   working-dir$ abaqus make job=remap

bash-Environment
---------------- 

# Abaqus 6.11:
# LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/6.11-1/exec/lbr
# LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/6.11-1/External
# LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/6.11-1/External/Backbone

# Abaqus 6.12-3:
# LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/6.12-3/code/bin
# export GXX_INCLUDE=/usr/include/c++/4.4.5/

# Abaqus 6.13-2:
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/6.13-2/code/bin

doc_make_end


*/

/* TODO:

   - call from python script ...? calling procedure / installation / setup
   - csv output (unstructured points)

   - read element coordinates and calculate weights inside
   - new .interpol version completely read-/writable from C and Python
   - generalize variables (make available all odb variables)
   - parallel processing

   - C-program to generate element coordinates
 */



// doc_include_start


// header-files:

#include <stdlib.h>  // for malloc(), free()
#include <stdio.h>  // for stdin
#include <unistd.h>  // for isatty(), fileno()

#include <iostream>
#include <fstream>
#include <sstream>

#include <algorithm>   // for find()
#include <string>
#include <vector>
#include <map>
#include <set>

#include <odb_API.h>

// Some abbreviations
using std::cout;
using std::cerr;
using std::endl;

class ListOfStrings : public std::vector<std::string> {};

// for debugging....   #### DEBUG
#include <typeinfo>
//... std::cout << typeid(a).name() << '\n';


// Some helper functions:
#include "helper_functions.cc"
#include "odbField.cc"
#include "interpol.cc"

// openmp:
// #include <omp.h>

// doc_include_end

// doc_0_start

int ABQmain(int argc, char** argv) {

  // ---------
  // --- 0 --- Declaration section.
  // ---------

  // ---------------------------------------------------------
  // --- 0 A: Variables for handling input data = odb data ---
  // ---------------------------------------------------------

  // Mapping types/dictionaries:
  // 1) Nodal vars:
  typedef std::map< int, const float * > int_floatsp;
  int_floatsp node_U_ptr;

  // 2) Element vars:
  typedef std::map< int, std::vector<const float *> > int_vec_floatsp;
  int_vec_floatsp elem_S_ptr;
  
  typedef std::map< int, std::vector<float> > int_vec_floats;
  int_vec_floats elem_minP;
  int_vec_floats elem_maxP;


  // -------------------------------------------------------
  // --- 0 B: Output variables = variables on gridpoints ---
  // -------------------------------------------------------
  //
  // Scalar variables:
  //
  // 1) Nodal vars:
  // U-Components as scalar variables:
  // 3 vectors of floats combined in one array (of 3 vectors of floats):
  std::vector< float > U_cmp[3];
  // U_cmp is 3 arrays of numGridPoints floats:
  // { [ ux ux ... ][ uy uy ... ][ uz uz ... ] }
  // { [  U_cmp[0] ][  U_cmp[1] ][  U_cmp[2] ] }

  // 2) Element vars:
  // S-Components as scalar variables:
  std::vector<float> S_cmp[6];
  // S_cmp is 6 arrays of numGridPoints floats:
  // { [sxx sxx ...][syy syy ...][szz szz ...][sxy sxy ...][sxz sxz ...][syz syz ...] }
  // { [  S_cmp[0] ][  S_cmp[1] ][  S_cmp[2] ][  S_cmp[3] ][  S_cmp[4] ][  S_cmp[5] ] }
  //
  // Min./max. princ. stress as scalar variable:
  // S_minP, S_maxP are vectors of numGridPoints floats:
  std::vector<float> S_minP;
  std::vector<float> S_maxP;
  //
  //
  // Vector or tensor variables:
  //
  // 1) Nodal vars:
  // U as vector:
  std::vector<float> U_vec;
  // U_vec is one array containg 3 x numGridPoints floats:
  // { [ ux  uy  uz      ux  uy  uz      ux  uy  uz   ... }
  // { [  0   1   2       3   4   5       6   7   8   ... }

  // 2) Element vars:
  std::vector<float> S_tnsr;
  // S_tnsr is one array containg 6 x numGridPoints floats:
  // { [ xx  yy  zz  xy  xz  yz      xx  yy  zz  xy  xz  yz      ... }
  // { [  0   1   2   3   4   5       6   7   8   9  10  11      ... }


  // -------------------------------
  // --- 0 C: Other declarations ---
  // -------------------------------
  
  // For counting time:
  double wall_start, wall_end;
  
  // For reading odb-data:
  int val_size, d_size;

  // Find out, if to read U or S or both:
  bool read_U    = false;
  bool read_S    = false;
  bool needsNodesWeights = false;
  bool needsElemIpWeights = false;

  // wallclock-time:
  wall_start  = get_wall_time();

  // doc_0_end


//   // --------------
//   // --- openmp --- 
//   // --------------
//   // 

//   // n_threads:
//   omp_set_num_threads(8);

// #pragma omp parallel for
//   for (int i = 0; i < 8; ++i)
//     {
//       const int id = omp_get_thread_num();      
//       // printf("This is thread %d\n", id);
      
//       if (id == 0)
// 	/* Nur im Masterthread ausführen */
// 	printf("There are %d threads\n", omp_get_num_threads());
//       for (int j = 0; j < 3; ++ j){
// 	for ( int k = 0; k < 2000000000; ++k ) int c = 17;
//       }
//     }


//   wall_end  = get_wall_time();
//   std::cout << "All took: " << wall_end - wall_start << " s" << std::endl;
//   exit(0);




  // doc_1_start

  // ---------
  // --- 1 --- Argument processing
  // ---------
  std::cout << std::endl << std::endl << "1: Processing input ";


  // -------------------------------------------
  // --- 1 A: Find out how to read arguments ---
  // -------------------------------------------

  ListOfStrings programArgumentList;

  // Four cases:
  //
  // CASE i:
  //
  // Program is called like this: 
  // $ ./remap datei.txt
  // Note: parameter file accepts empty lines and comments
  if (argc == 2){
    std::cout << "from file " << argv[1] << " .";
    std::cout << std::endl << std::endl;
    // processing lines of file:
    std::ifstream inFile(argv[1]);
    ASSERT(inFile); // make sure file exists

    for( std::string line; getline( inFile, line ); ){
      if (line.empty() || line[0]=='#') {
	// ingore empty line and comments
	continue;
      }
      programArgumentList.push_back(line);
    }
    inFile.close();

  }
  //
  // CASE ii:
  // No command line arguments: read from stdin
  //
  else if (argc==1) {
    bool interactive = isatty(fileno(stdin));
    std::string value;

    ListOfStrings promptsList;
    if (interactive) {
      promptsList.push_back("Path to odb: ");
      promptsList.push_back("Path to interpolation-file: ");
      promptsList.push_back(
	     "Requested variables (comma separated, no spaces): ");
      promptsList.push_back("ascii/binary: ");
      promptsList.push_back("stepName,frame,vtk-filename: ");
    }

    std::cout << "from stdin." << std::endl;

    ListOfStrings::const_iterator prompt=promptsList.begin();
    // for each prompt in the list get one value
    for (; prompt!=promptsList.end(); ++prompt) {
      if (interactive) {
	std::cout << *prompt;
      }
      std::getline(std::cin, value);
      programArgumentList.push_back(value);
    }
    // now continue with the last prompt until user presses Enter twice
    prompt = promptsList.end() - 1;
    bool stop = false;
    while (!stop) {
      if (interactive) {
	std::cout << *prompt;
      }
      std::getline(std::cin, value);

      if (interactive) {
	stop = value.empty();
      }
      else {
	stop = std::cin.eof();
      }
      if (!stop) {
	programArgumentList.push_back(value);
      }
    }

    // // for DEBUGGING:
    // for (ListOfStrings::const_iterator val=programArgumentList.begin();
    // 	 val != programArgumentList.end();
    // 	 ++val) {
    //   std::cout << "*" << *val << "*" << std::endl;
    // }
    // exit(0);

  }

  //
  // CASE iii:
  // Program does *not* get enough command line arguments:
  //
  else if (argc<6){
    std::cout << "from command line arguments, but... ";
    std::cout << std::endl;
    std::cout << "   ... sorry. You must provide 5 or more arguments:";

    // print some info about the arguments,
    // which need to be provided and exit:
    print_info();
    exit(1);
  }
  //
  // CASE iv:
  // Program does get enough command line arguments:
  else{
    std::cout << "from command line args.";
    for(int i=1; i<argc; ++i){
      programArgumentList.push_back(argv[i]);
    }
  }
    

  // ------------------------------
  // --- 1 B: Process arguments ---
  // ------------------------------


  // ----------
  // --- 1. --- Path to odb 
  // ----------
  static const odb_String odb_path = programArgumentList[0].c_str();

  // ---------------
  // --- 2. --- Path to interpolation data file
  // ---------------
  const std::string fileNameInterpol = programArgumentList[1];

  // ----------
  // --- 3. --- list of variable-names
  // ----------

  // populate requestedVarNameList from second argument
  std::string vars_chars = programArgumentList[2];
  ListOfStrings requestedVarNameList;
  std::istringstream vars_buf(vars_chars);
  for(std::string token; getline(vars_buf, token, ','); ) {
    requestedVarNameList.push_back(token);
  }

  // iterate over requestedVarNameList ...
  // decide, if U and/or S have to be read:  
  // Create centering_v and vardims_v;
  std::vector<int> centering_v, vardims_v;

  // list of pointers to actual data
  std::vector<float *> interpolatedDataArrays;
  interpolatedDataArrays.resize(requestedVarNameList.size());

  ListOfStrings read_Scalars;
  for (ListOfStrings::iterator it=requestedVarNameList.begin();
       it!=requestedVarNameList.end();
       ++it) {
    std::string& varName = *it;
    if ( (varName == "U") or
	 (varName == "U_1") or
	 (varName == "U_2") or
	 (varName == "U_3") ) {
      read_U=true;
      needsNodesWeights = true;
    }
    else if ( (varName == "S")      or
	      (varName == "S_11")   or
	      (varName == "S_22")   or
	      (varName == "S_33")   or
	      (varName == "S_12")   or
	      (varName == "S_13")   or
	      (varName == "S_23")   or
	      (varName == "S_MIN_PRINCIPAL") or
	      (varName == "S_MAX_PRINCIPAL") ) {
      read_S=true;
      needsElemIpWeights = true;
    }
    else {
      // all others (i.e. SDVs) are assumed to be scalars at integration points
      read_Scalars.push_back(varName);
      needsElemIpWeights = true;
    }

    if (varName == "U") {
      centering_v.push_back(1);
      vardims_v.push_back(3);
    }
    else if (varName == "S") {
      centering_v.push_back(1);
      vardims_v.push_back(9);
    }
    else {
      centering_v.push_back(1);
      vardims_v.push_back(1);
    }
  }

  // ----------
  // --- 4. --- vtk output-format
  // ----------
  const std::string bin_asc = programArgumentList[3];
  int i_bin;
  if      (bin_asc == std::string("ascii"))  i_bin=0;
  else if (bin_asc == std::string("binary")) i_bin=1;
  else {
    cout << "Vtk output format must be ascii or binary. Instead got "
	 << bin_asc << " ." << endl;
    exit(2);
  }

  // -----------------
  // --- 5, 6, ... --- Frame-indices, vtk-filenames
  // -----------------
  // vector of step names:
  ListOfStrings stepNamesList;
  // vector of frame indices:
  std::vector<int> frame_idxs;
  // vector of vtk filenames:
  ListOfStrings vtk_fnames;

  int n_frms = programArgumentList.size()-4;
  for (int i_frm=0; i_frm<n_frms; i_frm++) {
    std::string chars = programArgumentList[4+i_frm];
    std::istringstream buf(chars);
    ListOfStrings v;
    for(std::string token; getline(buf, token, ','); ) {
      v.push_back(token);
    }
    // step name
    stepNamesList.push_back(v[0]);
    // frame index:
    int frm_idx = stringToInteger( v[1] );
    frame_idxs.push_back( frm_idx );
    // vtk filename:
    vtk_fnames.push_back( v[2] );
    v.clear();
  }
  // doc_1_end

  // doc_2_start

  // ---------
  // --- 2 --- Opening odb, Reading weights.
  // ---------

  // ---------------------
  // --- 2 A: Open odb ---
  // ---------------------

  cout << endl;
  cout << endl;
  cout << "2: Opening odb...";
  cout << endl;

  odb_initializeAPI();

  // Open odb read-only:
  static const odb_String name = "";  
  odb_Odb& myOdb = openOdb(name, odb_path, true);
  odb_StepRepository& steps = myOdb.steps();

  // instance, nsets, elsets... for all fields
  OdbFieldReaderCommonData odbCommon(myOdb);
  OdbFieldReader::m_odbCommonData = &odbCommon;

  cout << "     --> odb:  " << odb_path.CStr() << endl;

  // ------------------------------------
  // --- 2 B: Read interpolation file ---
  // ------------------------------------

  cout << endl;
  cout << endl;
  cout << "2b: Reading interpolation file...";
  cout << endl;

  Interpolation* interpol;
  interpol = new Interpolation;
  interpol -> readFromPickleFile(fileNameInterpol);
  

  // ------------------------------------
  // --- 2 C: Calculate nodes-weights ---
  // ------------------------------------

  if (needsNodesWeights){
    cout << "     --> Calculating nodes-weights." << endl;
    interpol -> initNodalWeights();
  }
  odbCommon.initNset(interpol->allNodes);

  // ------------------------------------
  // --- 2 D: Calculate elems-weights ---
  // ------------------------------------

  // elements required for the interpolation
  odb_SequenceInt e_lbls( interpol->allElems.size() );

  if (needsElemIpWeights){
    cout << "     --> Calculating element IP weights." << endl;
    interpol -> initElemIpWeights();
  }
  // store set of elems required for the interpolation
  odbCommon.initElset(interpol->allElems);

  interpol -> printDiagnosticOutput();

  // doc_2_end

  // doc_3A_start

  // ---------
  // --- 3 --- Read, interpolate, write vtk
  // ---------

  std::cout << std::endl;
  int n_frames = programArgumentList.size()-6;
  std::cout << "3: Processing " << n_frames << " frames." << std::endl ;

  // -----------------------------
  // --- 3 A: Reading from odb ---
  // -----------------------------
  
  for (int i_frm=0; i_frm<n_frames; i_frm++) { // Frame-loop-start.

    std::string stepName = stepNamesList[i_frm];
    odb_String stepNameOS = stepName.c_str();
    odb_Step& step = steps[stepNameOS];
    odb_SequenceFrame& frameSequence = step.frames();

    // frame index:
    int frm_idx = frame_idxs[i_frm];
    ASSERT(frm_idx >=0 and frm_idx < frameSequence.size());
    // vtk-filename:
    const char *vtk_fname = vtk_fnames[i_frm].c_str();

    odb_Frame& frame = frameSequence[frm_idx];

    std::cout << "\n( Step, Frame, Step-time ) = ( "\
	      << stepName << "," << frm_idx << ", "\
	      << frame.frameValue() << " )" << std::endl;

    if (read_U){
      std::cout << "  A: Reading U from odb." << std::endl;
      // Get field output U for subset:
      odb_FieldOutput&              field_op_U = frame.fieldOutputs()["U"];
      odb_FieldOutput           sub_field_op_U = field_op_U.getSubset(*(odbCommon.m_nset));
      const odb_SequenceFieldValue& seq_vals_U = sub_field_op_U.values();
      val_size = seq_vals_U.size();
      // std::cout << "   Values in sequence: " << val_size << std::endl;
      d_size = 0;
      // Iterating over all sequence values:
      for (int l=0; l<val_size; l++) {
	const odb_FieldValue val = seq_vals_U[l];
	const float* const data = val.data(d_size);
	int nl = val.nodeLabel();
	// Storing pointer to first component:
	node_U_ptr[nl]=data;
      }
    }    
    if (read_S){
      std::cout << "  A: Reading S from odb." << std::endl;
      odb_FieldOutput&              field_op_S = frame.fieldOutputs()["S"];
      odb_FieldOutput           sub_field_op_S = field_op_S.getSubset((*odbCommon.m_elset));
      const odb_SequenceFieldValue& seq_vals_S = sub_field_op_S.values();
      // const odb_SequenceInvariant invars = sub_field_op_S.validInvariants();
      val_size = seq_vals_S.size();
      d_size = 0;

      for (int l=0; l<val_size; l++) {
	const odb_FieldValue val = seq_vals_S[l];
	const float* const data = val.data(d_size);

	// val.elementLabel();
	// val.integrationPoint();
	int el = val.elementLabel();

	// 1. Pointer to components of S:
	elem_S_ptr[el].push_back(data);

	// 2. Storing minP / maxP
	float minP = val.minPrincipal();
	float maxP = val.maxPrincipal();
	elem_minP[el].push_back(minP);
	elem_maxP[el].push_back(maxP);
      }
    }

    // doc_3A_end    
    // doc_3B_start

    // --------------------------------------
    // --- 3 B: Interpolating / Remapping ---
    // --------------------------------------

    int numGridPoints = interpol->grid->size();
    if (read_U){
      cout << "  B: Interpolating U." << endl;
      // Interpolate U
      int n_components=3;

      for (int iGridPt=0; iGridPt < numGridPoints; iGridPt++) {
	// const int id = omp_get_thread_num();      
	// printf("This is thread %d\n", id);

	if (interpol->elemNbs[iGridPt] > 0) {
	  float *wght_ptr = &(interpol->n_wght_of_pt[iGridPt][0]);
	  unsigned int *node_ptr = &(interpol->node_of_pt[iGridPt][0]);
	  // U has n=3 components:
	  for (int n=0; n<n_components; n++){
	    float U_int=0.;
	    for (int i=0; i<10; i++){
	      unsigned int node = *(node_ptr+i);
	      // add displacement of i-th node times weight:
	      U_int+= *(node_U_ptr[node]+n) * *(wght_ptr+i);
	    }
	    U_vec.push_back(U_int);
	    U_cmp[n].push_back(U_int);
	  }
	}
	else {
	  for (int n=0; n<n_components; n++){
	    U_vec.push_back(0.0);
	    U_cmp[n].push_back(0.0);
	  }
	}
      }
    } // if (read_U) ends here
    
    if (read_S){
      // Interpolate S
      cout << "  B: Interpolating S." << endl;
      int n_components=6;

      for (int iGridPt=0; iGridPt < numGridPoints; iGridPt++) {
	// element of this line:
	int elem = interpol->elem_of_pt[iGridPt];

	if (elem > 0) {
	  // pointer to weights:
	  float *wght_ptr = &(interpol->e_wght_of_pt[iGridPt][0]);
	  // stress has n=6 components:
	  for (int n=0; n<n_components; n++){

	    float S_int=0.;
	    // sum over IPs:
	    for (int i=0; i<4; i++){
	      // Interpolate all components using pointers:
	      S_int+=*(elem_S_ptr[elem][i]+n) * *(wght_ptr+i);
	    }
	    // if stress tensor needs to be written:
	    S_tnsr.push_back(S_int);
	    S_cmp[n].push_back(S_int);
	  }
      
	  // min/max princ stress:
	  float S_minP_int=0.;
	  float S_maxP_int=0.;
	  for (int i=0; i<4; i++){
	    S_minP_int+=elem_minP[elem][i] * *(wght_ptr+i);
	    S_maxP_int+=elem_maxP[elem][i] * *(wght_ptr+i);
	  }
	  S_minP.push_back(S_minP_int);
	  S_maxP.push_back(S_maxP_int);
	}
	else {
	  for (int n=0; n<n_components; n++){
	    S_tnsr.push_back(0.0);
	    S_cmp[n].push_back(0.0);
	  }
	  S_minP.push_back(0.0);
	  S_maxP.push_back(0.0);
	}
      }
    }
    cout << "  B: Finished interpolating U and S." << endl;
    // doc_3B_end

    // doc_3C_start
    // --------------------------------------
    // --- 3 C: Reading and Interpolating SDVs and storing U and S
    // --------------------------------------

    // loop over all varName in requestedVarNameList
    // and iterate interpolDataIdx as index in interpolatedDataArrays
    int interpolDataIdx = 0;
    OdbFieldReader* odbFieldReader;
    for (ListOfStrings::const_iterator it=requestedVarNameList.begin();
    	 it!=requestedVarNameList.end();
    	 ++it, ++interpolDataIdx) {
      const std::string& varName = *it;

      // if current varName in read_Scalars
      if (std::find(read_Scalars.begin(), read_Scalars.end(), varName
		    ) != read_Scalars.end()) {

	// Read field data from odb
	cout << "  A: Reading " << varName << " from odb." << endl;

	odbFieldReader = new OdbFieldReader(varName);
	odbFieldReader->getFieldFromFrame(frame);

	// Interpolate field data
	std::cout << "  B: Interpolating " << varName << "." << std::endl;

	// allocate memory
	void* buffer = malloc(sizeof(float)*numGridPoints);
	if (!buffer) {
	  cerr << "Could not allocate memory for " << numGridPoints
	       << " points for " << varName << "." << endl;
	  exit(2);
	}
	interpolatedDataArrays[interpolDataIdx] = (float*) buffer;

	// calculate interpolation and store
	if (odbFieldReader->m_location == odb_Enum::INTEGRATION_POINT) {
	  OdbFieldDataMapIntegrationPt& map
	    = *((OdbFieldDataMapIntegrationPt*) odbFieldReader->m_map);
	  for (int iGridPt=0; iGridPt < numGridPoints; iGridPt++) {
	    // element of this line:
	    int elem = interpol->elem_of_pt[iGridPt];

	    if (elem>0) {
	      // pointer to weights:
	      float *weightPtr = &(interpol->e_wght_of_pt[iGridPt][0]);
	      std::vector<const float*>& valuePtrVec = map[elem];
	      float value=0.;
	      // sum over IPs:
	      for (int i=0; i<4; i++){
		// Interpolate using pointers:
		value += *(valuePtrVec[i]) * *(weightPtr+i);
	      }
	      interpolatedDataArrays[interpolDataIdx][iGridPt] = value;
	    }
	    else {
	      interpolatedDataArrays[interpolDataIdx][iGridPt] = 0.0;
	    }
	  }
	}
	delete odbFieldReader;
      }

      // else if current varName is U or S
      else {
	if      (varName ==   "S_11")
	  interpolatedDataArrays[interpolDataIdx] = &(S_cmp[0][0]);
	else if (varName ==   "S_22")
	  interpolatedDataArrays[interpolDataIdx] = &(S_cmp[1][0]);
	else if (varName ==   "S_33")
	  interpolatedDataArrays[interpolDataIdx] = &(S_cmp[2][0]);
	else if (varName ==   "S_12")
	  interpolatedDataArrays[interpolDataIdx] = &(S_cmp[3][0]);
	else if (varName ==   "S_13")
	  interpolatedDataArrays[interpolDataIdx] = &(S_cmp[4][0]);
	else if (varName ==   "S_23")
	  interpolatedDataArrays[interpolDataIdx] = &(S_cmp[5][0]);
	else if (varName == "S_MIN_PRINCIPAL")
	  interpolatedDataArrays[interpolDataIdx] = &(S_minP[0]  );
	else if (varName == "S_MAX_PRINCIPAL")
	  interpolatedDataArrays[interpolDataIdx] = &(S_maxP[0]  );
	else if (varName ==    "U_1")
	  interpolatedDataArrays[interpolDataIdx] = &(U_cmp[0][0]);
	else if (varName ==    "U_2")
	  interpolatedDataArrays[interpolDataIdx] = &(U_cmp[1][0]);
	else if (varName ==    "U_3")
	  interpolatedDataArrays[interpolDataIdx] = &(U_cmp[2][0]);
	else if (varName ==      "U")
	  interpolatedDataArrays[interpolDataIdx] = &(U_vec[0]   );
      }
    }
    // doc_3C_end


    // doc_4_start

    // ------------------------
    // --- 4: Writing vtk ---
    // ------------------------
    std::cout << "  C: " << bin_asc << " data to: " << vtk_fname << std::endl;


    // Write vtk-file:
    // parameters:
    // 1. output file
    // 2. variable descriptions
    // 3. stored variables
    // 4. grid

    // write vtk using simple writer
    write_to_vtk(vtk_fname, i_bin, *(interpol->grid),
		 requestedVarNameList.size(),
		 &(vardims_v[0]), requestedVarNameList,
		 &(interpolatedDataArrays[0]));

    // // write vtk using visit writer
    // // create varnamesArr as an (const char* const*) representation of varnames
    // int n_vars=requestedVarNameList.size();
    // const char **varnamesArr = (const char **)malloc(n_vars * sizeof(char *));
    // for(int i = 0; i < n_vars; ++i)
    //   varnamesArr[i] = requestedVarNameList[i].c_str();
  
    // write_structured_points_mesh(vtk_fname, i_bin,
    // 				 dims, n_vars,
    // 				 &(vardims_v[0]), &(centering_v[0]),
    // 				 varnamesArr, &(interpolatedDataArrays[0]),
    // 				 org, spc);
    // free(varnamesArr);
    
    // Clearing memory:
    // www.cplusplus.com/reference/vector/vector/clear/
    U_vec.clear();
    for (int i=0; i<3; i++)
      U_cmp[i].clear();
    for (int i=0; i<6; i++)
      S_cmp[i].clear();
    S_minP.clear();
    S_maxP.clear();
    // www.cplusplus.com/reference/map/map/clear/
    // node_U_ptr.clear(); not necessary, will be overwritten.
    elem_S_ptr.clear();
    elem_minP.clear();
    elem_maxP.clear();


  } // Frame-loop-end.



  // interpol -> test();
  // cout << "len(nset_U): " << nset_U.size() << endl;
  // cout << "len(elset_S): " << elset_S.size() << endl;
  // std::cout << "THE END." << std::endl;
  // exit(0);
  // /***************************************************
  // *****************************************************************/

  std::cout << std::endl;
  wall_end  = get_wall_time();
  std::cout << "All took: " << wall_end - wall_start << " s" << std::endl;

  odb_finalizeAPI();
  // doc_4_end

  return 0;
}



