/* 
doc_usage_start

abaqus remap parameter.txt

with parameter.txt looking like this:
--- snip ---
test8.odb
S_11,U,S_22,SDV1,U_2
ascii
nodes_weights_1k.bin
elems_weights_1k.bin
10,20,5,100.,0.,0.,1.,2.,4.
99,prj_99.vtk
100,prj_100.vtk
--- snap ---

or...
--- snip ---
test8.odb
SDV1,S_MIN_PRINCIPAL,U_3
ascii
nodes_weights_1k.bin
elems_weights_1k.bin
10,20,5,100.,0.,0.,1.,2.,4.
99,prj_99.vtk
77,prj_77.vtk
--- snap ---

or...
--- snip ---
test8.odb
U,S,SDV1
binary
nodes_weights_1mio.bin
elems_weights_1mio.bin
100,200,50,100.,0.,0.,1.,2.,4.
99,prj_99.vtk
77,prj_77.vtk
--- snap ---

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
    B: Read nodes-weights file
    C: Read elems-weights file

 3: Read, interpolate, write vtk

    A: Reading from odb    
    B: Interpolating
    C: Writing vtk

doc_desc_end

doc_make_start

Compiling/Linking 
-----------------

1) Compile-command (if this file is called remap.cpp): 

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

// doc_include_start

// header-files:
#include "headers.h"


// hard coded field names
// const std::string varName_U = "UEQUIB";
// const std::string varName_U_1 = "UEQUIB_1";
// const std::string varName_U_2 = "UEQUIB_2";
// const std::string varName_U_3 = "UEQUIB_3";
// const std::string varName_SDV1 = "SDV4";
// const std::string varName_SDV2 = "SDV12";
// const std::string varName_SDV3 = "SDV15";
#include "hardCodedFieldNames.cpp"

// visit_writer to write vtk, see visit_writer_mod.c for details.
#include "visit_writer.h"
#include "visit_writer_mod.c"

// Some helper functions:
#include "helper_functions.cpp"

// doc_include_end


// doc_0_start

int ABQmain(int argc, char** argv) {

  // ---------
  // --- 0 --- Declaration section.
  // ---------

  // ---------------------------------------------------------
  // --- 0 A: Variables for handling input data = odb data ---
  // ---------------------------------------------------------

  // 1) Nodal vars, such as U or V:
  // Each gridpoint is inside exactly one element, and each
  // element has exactly ten nodes with one weight each.
  //
  // vector of 10 x n_GP nodes:
  std::vector< std::vector<int> > node_of_pt;  
  // vector of 10 ints (fixed size):
  std::vector<  int  > nods(10);
  // vector of 10 floats (fixed size):
  std::vector< float > n_wghts(10);
  // 10 weights, one for each of the ten nodes:
  std::vector< std::vector< float > > n_wght_of_pt;
  //      first gridpoint     second gridpoint    ...   last gridpoint
  // [ [ 0.24 0.38 0.55 ...][ 0.41 0.25 0.51 ...] ... [       ...       ] ]

  // 2) Element vars, such as S or S_minP or S_maxP:
  // Each gridpoint is inside exactly one element, and
  // each element has exactly 4 integration points
  // with one weight each:

  // vector of n_GP integers:
  // [ 278, 2256, ... ]
  std::vector< int > elem_of_pt;
  // vector of 4 floats (fixed size):
  std::vector< float > ip_wghts(4);
  std::vector< std::vector< float > > e_wght_of_pt;

  // Mapping types/dictionaries:
  // 1) Nodal vars:
  typedef std::map< int, const float * > int_floatsp;
  int_floatsp node_U_ptr;

  // 2) Element vars:
  typedef std::map< int, std::vector<const float *> > int_vec_floatsp;
  int_vec_floatsp elem_S_ptr;
  int_vec_floatsp elem_SDV1_ptr;
  int_vec_floatsp elem_SDV2_ptr;
  int_vec_floatsp elem_SDV3_ptr;
  
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
  // U_cmp is 3 arrays of n_GP floats:
  // { [ ux ux ... ][ uy uy ... ][ uz uz ... ] }
  // { [  U_cmp[0] ][  U_cmp[1] ][  U_cmp[2] ] }

  // 2) Element vars:
  // S-Components as scalar variables:
  std::vector<float> S_cmp[6];
  // S_cmp is 6 arrays of n_GP floats:
  // { [sxx sxx ...][syy syy ...][szz szz ...][sxy sxy ...][sxz sxz ...][syz syz ...] }
  // { [  S_cmp[0] ][  S_cmp[1] ][  S_cmp[2] ][  S_cmp[3] ][  S_cmp[4] ][  S_cmp[5] ] }
  //
  // Min./max. princ. stress as scalar variable:
  // S_minP, S_maxP are vectors of n_GP floats:
  std::vector<float> S_minP;
  std::vector<float> S_maxP;
  std::vector<float> SDV1;
  std::vector<float> SDV2;
  std::vector<float> SDV3;
  //
  //
  // Vector or tensor variables:
  //
  // 1) Nodal vars:
  // U as vector:
  std::vector<float> U_vec;
  // U_vec is one array containg 3 x n_GP floats:
  // { [ ux  uy  uz      ux  uy  uz      ux  uy  uz   ... }
  // { [  0   1   2       3   4   5       6   7   8   ... }

  // 2) Element vars:
  std::vector<float> S_tnsr;
  // S_tnsr is one array containg 6 x n_GP floats:
  // { [ xx  yy  zz  xy  xz  yz      xx  yy  zz  xy  xz  yz      ... }
  // { [  0   1   2   3   4   5       6   7   8   9  10  11      ... }


  // -------------------------------
  // --- 0 C: Other declarations ---
  // -------------------------------
  bool check;
  check = false;
  check = true;
  
  // For counting time:
  double wall_start, wall_end;
  
  // For reading BINARY files:
  FILE *f_bin;
  int f_bin_size;
  int node, elem;
  float wght;

  // For reading odb-data:
  int val_size, d_size;

  // Find out, if to read U or S or both:
  bool read_U    = false;
  bool read_S    = false;
  bool read_SDV1 = false;
  bool read_SDV2 = false;
  bool read_SDV3 = false;

  // wallclock-time:
  wall_start  = get_wall_time();

  // doc_0_end


  // doc_1_start

  // ---------
  // --- 1 --- Argument processing
  // ---------
  std::cout << std::endl << std::endl << "1: Processing input ";
  if (argc != 2){
    std::cout << std::endl;
    std::cout << "   ... sorry. You must provide parameter file name "
	      << "as single argument." << std::endl;

    // print some info about the arguments,
    // which need to be provided and exit:
    print_info();
    exit(0);
  }


  std::vector<std::string> parameterlist;

  {
    std::cout << "from file " << argv[1] << " .";
    std::cout << std::endl << std::endl;
    // processing lines of file:
    std::ifstream inFile(argv[1]);
    ASSERT(inFile); // make sure file exists

    for( std::string line; getline( inFile, line ); ){
      ASSERT(!line.empty()); // make sure line is not empty
      parameterlist.push_back(line);
    }
    inFile.close();

  }
    

  // ------------------------------
  // --- 1 B: Process arguments ---
  // ------------------------------


  // ----------
  // --- 1. --- Path to odb 
  // ----------
  const char *odb_path = parameterlist[0].c_str();

  // ----------
  // --- 2. --- list of variable-names
  // ----------
  // push variable names into vector of strings:
  std::string vars_chars = parameterlist[1];
  std::vector<std::string> vars_v;
  std::istringstream vars_buf(vars_chars);
  for(std::string token; getline(vars_buf, token, ','); )
    vars_v.push_back(token);
  int n_vars = vars_v.size();

  // decide, if U and/or S have to be read:  
  if (
      (std::find(vars_v.begin(), vars_v.end(), varName_U)   != vars_v.end()) or
      (std::find(vars_v.begin(), vars_v.end(), varName_U_1) != vars_v.end()) or
      (std::find(vars_v.begin(), vars_v.end(), varName_U_2) != vars_v.end()) or
      (std::find(vars_v.begin(), vars_v.end(), varName_U_3) != vars_v.end())
      )
    read_U=true;
  if (
      (std::find(vars_v.begin(), vars_v.end(), "S")      != vars_v.end()) or
      (std::find(vars_v.begin(), vars_v.end(), "S_11")   != vars_v.end()) or
      (std::find(vars_v.begin(), vars_v.end(), "S_22")   != vars_v.end()) or
      (std::find(vars_v.begin(), vars_v.end(), "S_33")   != vars_v.end()) or
      (std::find(vars_v.begin(), vars_v.end(), "S_12")   != vars_v.end()) or
      (std::find(vars_v.begin(), vars_v.end(), "S_13")   != vars_v.end()) or
      (std::find(vars_v.begin(), vars_v.end(), "S_23")   != vars_v.end()) or
      (std::find(vars_v.begin(), vars_v.end(), "S_MIN_PRINCIPAL") != vars_v.end()) or
      (std::find(vars_v.begin(), vars_v.end(), "S_MAX_PRINCIPAL") != vars_v.end())
      )
    read_S=true;
  if (
      (std::find(vars_v.begin(), vars_v.end(), varName_SDV1)      != vars_v.end())
      )
    read_SDV1=true;
  if (
      (std::find(vars_v.begin(), vars_v.end(), varName_SDV2)      != vars_v.end())
      )
    read_SDV2=true;
  if (
      (std::find(vars_v.begin(), vars_v.end(), varName_SDV3)      != vars_v.end())
      )
    read_SDV3=true;

  /* Create list of var names from the vector of strings, 
     add centering_v and vardims_v;
     so that we have: n_vars, centering_v, vardims_v, varnames */

  std::vector<int> centering_v, vardims_v;
  const char **varnames = NULL;
  varnames = (const char **)malloc(n_vars * sizeof(char *));
  for(int i = 0; i < n_vars; ++i){
      char tmp[100];
      sprintf(tmp, vars_v[i].c_str());
      varnames[i] = strdup(tmp);
      if (vars_v[i]==varName_U){
	centering_v.push_back(1);
	vardims_v.push_back(3);
      }
      else if (vars_v[i]=="S"){
	centering_v.push_back(1);
	vardims_v.push_back(9);
      }
      else{
	centering_v.push_back(1);
	vardims_v.push_back(1);
      }
  }


  // ----------
  // --- 3. --- vtk output-format
  // ----------
  const std::string bin_asc = parameterlist[2];
  int i_bin;
  ASSERT(  bin_asc == std::string("ascii") or \
	   bin_asc == std::string("binary") );
  if      (bin_asc == std::string("ascii"))  i_bin=0;
  else if (bin_asc == std::string("binary")) i_bin=1;

  // ---------------
  // --- 4. / 5. --- Path to nodes- / elems-weights-file
  // ---------------
  const std::string f_nds_wghts = parameterlist[3];
  const std::string f_els_wghts = parameterlist[4];


  // ----------
  // --- 6. --- Grid-data
  // ----------
  std::string grid_chars = parameterlist[5];
  std::vector<std::string> grid_v;
  std::istringstream grid_buf(grid_chars);
  for(std::string token; getline(grid_buf, token, ','); )
    grid_v.push_back(token);

  // resolution:
  const int   nx    = atoi(grid_v[0].c_str());
  const int   ny    = atoi(grid_v[1].c_str());
  const int   nz    = atoi(grid_v[2].c_str());
  // compute number of gridpoints:
  const int n_GP=nx*ny*nz;

  // origin:
  const float org_x = atof(grid_v[3].c_str());
  const float org_y = atof(grid_v[4].c_str());
  const float org_z = atof(grid_v[5].c_str());
  float org[]={org_x,org_y,org_z};

  // spacing:
  const float spc_x = atof(grid_v[6].c_str());
  const float spc_y = atof(grid_v[7].c_str());
  const float spc_z = atof(grid_v[8].c_str());
  float spc[]={spc_x,spc_y,spc_z};

  // -----------------
  // --- 7, 8, ... --- Frame-indices, vtk-filenames
  // -----------------
  // vector of frame indices:
  std::vector<int> frame_idxs;
  // vector of vtk filenames:
  std::vector<std::string> vtk_fnames;

  int n_frms = parameterlist.size()-6;
  for (int i_frm=0; i_frm<n_frms; i_frm++) {
    std::string chars = parameterlist[6+i_frm];
    std::vector<std::string> v;
    std::istringstream buf(chars);
    for(std::string token; getline(buf, token, ','); )
      v.push_back(token);
    // frame index:
    int frm_idx=atoi(v[0].c_str());
    frame_idxs.push_back( frm_idx );
    // vtk filename:
    vtk_fnames.push_back( v[1].c_str() );
    v.clear();
  }

  // some diagnostic output
  std::cout << std::endl;
  std::cout << "Nodes weights file: " << f_nds_wghts << std::endl;
  std::cout << "Elem weights file: " << f_els_wghts << std::endl;
  std::cout << "Nb of grid points: " << n_GP << std::endl;
  std::cout << std::endl;

  // doc_1_end

  // doc_2_start

  // ---------
  // --- 2 --- Opening odb, Reading weights.
  // ---------
  
  // ---------------------
  // --- 2 A: Open odb ---
  // ---------------------

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "2: Opening odb, Reading binary weights-file(s)...";
  std::cout << std::endl;

  odb_initializeAPI();
  // nodes & elements as appearing in weights file:
  odb_SequenceInt n_lbls, e_lbls;

  // Specify odb, step, part, instance, frame:
  static const odb_String name         = "";  
  static const odb_String path         = odb_path;
  static const odb_String stepName     = "Step-2";
  static const odb_String instanceName = "PART-1-1";
  // Open odb read-only:
  std::cout << "     ... opening odb " << odb_path << " ..."; std::cout.flush();
  odb_Odb& myOdb = openOdb(name, path, true);
  std::cout << " ok." << std::endl;
  odb_Assembly rootAssy = myOdb.rootAssembly();
  odb_Instance instance = rootAssy.instances()[instanceName];
  odb_StepRepository& steps = myOdb.steps();
  odb_Step& step = steps[stepName];
  odb_SequenceFrame& frameSequence = step.frames();
  std::cout << "     ... step <" << stepName.cStr() << "> contains "
	    << step.frames().size() << " frames." << std::endl;

  // ------------------------------------
  // --- 2 B: Read nodes-weights file ---
  // ------------------------------------

  if (read_U){
    // Open nodes-weights-file:
    std::cout << "     --> Reading nodes/weights-file " << f_nds_wghts
	      << " for " << n_GP << " grid points." << std::endl;
    f_bin=fopen(f_nds_wghts.c_str(),"rb");
    f_bin_size = filesize(f_nds_wghts.c_str());
    std::cout << "     ... filesize " << f_bin_size
	      << ", expected " << n_GP*80 << std::endl;

    // check filesize:
    // int = 4 bytes und float = 4 bytes. 
    // Each gridpoint: 10 x (4+4) bytes = 80 bytes
    ASSERT(f_bin_size == n_GP * 80);

    // read nodes and weights from binary file:
    for (int GP=0; GP<n_GP; GP++){      
      // read nodes:
      for ( int cnt=0; cnt<10; cnt++){
    	fread(&node,sizeof(int),1,f_bin);
    	nods[cnt]=node;
    	n_lbls.append(node);
    	fread(&wght,sizeof(float),1,f_bin);
    	n_wghts[cnt]=wght;
      }
      // add node list at the end:
      node_of_pt.push_back(nods);
      // add weight list at the end:
      n_wght_of_pt.push_back(n_wghts);
    }
    // close
    fclose(f_bin);
    // 
  }
  // store set:
  static const odb_String nodes_U  = "nodes_U";
  instance.NodeSet(nodes_U, n_lbls);
  odb_Set& nset_U = instance.nodeSets()[nodes_U];

  // ------------------------------------
  // --- 2 C: Read elems-weights file ---
  // ------------------------------------

  if (read_S or read_SDV1 or read_SDV2 or read_SDV3){

    std::cout << "     --> Reading elems/weights-file " << f_els_wghts
	      << " for " << n_GP << " grid points." << std::endl;
    f_bin=fopen(f_els_wghts.c_str(),"rb");
    f_bin_size = filesize(f_els_wghts.c_str());
    std::cout << "     ... filesize " << f_bin_size
	      << ", expected " << n_GP*20 << std::endl;

    // Each gridpoint: 1x4 bytes + 4x4 bytes = 20 bytes.
    ASSERT(f_bin_size == n_GP * 20);

    // read elems and weights from binary file:
    for (int GP=0; GP<n_GP; GP++){      
      // read elems:
      fread(&elem,sizeof(int),1,f_bin);
      e_lbls.append(elem);
      elem_of_pt.push_back(elem);
      // // read weights one by one:
      // for ( int cnt=0; cnt<4; cnt++){
      // 	fread(&wght,sizeof(float),1,f_bin);
      // 	ip_wghts[cnt]=wght;
      // }
      // read 4 weights all at once:
      fread(&ip_wghts[0],sizeof(float),4,f_bin);
      // add weight list at the end:
      e_wght_of_pt.push_back(ip_wghts);
    }
    // close
    fclose(f_bin);
  }
  // store set:
  static const odb_String elems_S  = "elems_S";
  instance.ElementSet(elems_S, e_lbls);
  odb_Set& elset_S = instance.elementSets()[elems_S];

  // doc_2_end

  // doc_3A_start

  // ---------
  // --- 3 --- Read, interpolate, write vtk
  // ---------

  std::cout << std::endl;
  int n_frames = parameterlist.size()-6;
  std::cout << "3: Processing " << n_frames << " frames." << std::endl ;

  // -----------------------------
  // --- 3 A: Reading from odb ---
  // -----------------------------
  
  for (int i_frm=0; i_frm<n_frames; i_frm++) { // Frame-loop-start.

    int frm_idx = frame_idxs[i_frm];
    ASSERT(frm_idx >=0 and frm_idx < frameSequence.size());
    // frame index:
    if (frm_idx==0){
      std::cout << "*** Warning... atoi returns zero." << std::endl;
      std::cout << "    This can mean, that no integer was provided" << std::endl;
      std::cout << "    as the frame index. However: If you intended" << std::endl;
      std::cout << "    to write results for frame 0, this is ok." << std::endl;
    }
    // vtk-filename:
    const char *vtk_fname = vtk_fnames[i_frm].c_str();

    odb_Frame& frame = frameSequence[frm_idx];

    std::cout << "\n( Frame-index, Step-time ) = ( " << frm_idx;
    std::cout << ", " << frame.frameValue() << " )" << std::endl;

    odb_FieldOutputRepository& fieldRepos = frame.fieldOutputs();
    
    std::cout << "  A: Reading U and/or S and/or SDV from odb." << std::endl;
    // ********************* some lines below not shown ******************
    // doc_3A_hide_block_1_start
    if (read_U){
      // Get field output U for subset:
      // odb_FieldOutput&              field_op_U = fieldRepos[varName_U.c_str()];
      odb_FieldOutput& field_op_U = checkGetField(fieldRepos, varName_U.c_str(), frm_idx);
      odb_FieldOutput sub_field_op_U = field_op_U.getSubset(nset_U);
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
      // odb_FieldOutput&              field_op_S = fieldRepos["S"];
      odb_FieldOutput& field_op_S = checkGetField(fieldRepos, "S", frm_idx);
      odb_FieldOutput           sub_field_op_S = field_op_S.getSubset(elset_S);
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

    if (read_SDV1){
      // odb_FieldOutput& field_op_SDV1 = fieldRepos[varName_SDV1.c_str()];
      odb_FieldOutput& field_op_SDV1 = checkGetField(fieldRepos, varName_SDV1.c_str(), frm_idx);
      odb_FieldOutput sub_field_op_SDV1 = field_op_SDV1.getSubset(elset_S);
      const odb_SequenceFieldValue& seq_vals_SDV1 = sub_field_op_SDV1.values();

      val_size = seq_vals_SDV1.size();
      d_size = 0;

      std::cout << "  ... reading " << varName_SDV1 << ". There are "
		<< val_size << " values available." << std::endl;
      for (int l=0; l<val_size; l++) {
	const odb_FieldValue val = seq_vals_SDV1[l];
	const float* const data = val.data(d_size);

	// val.elementLabel();
	// val.integrationPoint();
	int el = val.elementLabel();

	// Pointer to values of SDV1:
	elem_SDV1_ptr[el].push_back(data);
      }
      std::cout << "  ... finished reading " << varName_SDV1
        << ", got values for " << elem_SDV1_ptr.size() << " elements."
        << std::endl;
      if (elem_SDV1_ptr.size() == 0) {
        std::cout << "ERROR: No values found. Stopping now." << std::endl;
        exit(2);
      }
    }

    if (read_SDV2){
      // odb_FieldOutput& field_op_SDV2 = fieldRepos[varName_SDV2.c_str()];
      odb_FieldOutput& field_op_SDV2 = checkGetField(fieldRepos, varName_SDV2.c_str(), frm_idx);
      odb_FieldOutput sub_field_op_SDV2 = field_op_SDV2.getSubset(elset_S);
      const odb_SequenceFieldValue& seq_vals_SDV2 = sub_field_op_SDV2.values();

      val_size = seq_vals_SDV2.size();
      d_size = 0;

      std::cout << "  ... reading " << varName_SDV2 << ". There are "
		<< val_size << " values available." << std::endl;
      for (int l=0; l<val_size; l++) {
	const odb_FieldValue val = seq_vals_SDV2[l];
	const float* const data = val.data(d_size);

	// val.elementLabel();
	// val.integrationPoint();
	int el = val.elementLabel();

	// Pointer to values of SDV2:
	elem_SDV2_ptr[el].push_back(data);
      }
      std::cout << "  ... finished reading " << varName_SDV2
        << ", got values for " << elem_SDV2_ptr.size() << " elements."
        << std::endl;
      if (elem_SDV2_ptr.size() == 0) {
        std::cout << "ERROR: No values found. Stopping now." << std::endl;
        exit(2);
      }
    }

    if (read_SDV3){
      // odb_FieldOutput& field_op_SDV3 = fieldRepos[varName_SDV3.c_str()];
      odb_FieldOutput& field_op_SDV3 = checkGetField(fieldRepos, varName_SDV3.c_str(), frm_idx);
      odb_FieldOutput sub_field_op_SDV3 = field_op_SDV3.getSubset(elset_S);
      const odb_SequenceFieldValue& seq_vals_SDV3 = sub_field_op_SDV3.values();

      val_size = seq_vals_SDV3.size();
      d_size = 0;

      std::cout << "  ... reading " << varName_SDV3 << ". There are "
		<< val_size << " values available." << std::endl;
      for (int l=0; l<val_size; l++) {
	const odb_FieldValue val = seq_vals_SDV3[l];
	const float* const data = val.data(d_size);

	// val.elementLabel();
	// val.integrationPoint();
	int el = val.elementLabel();

	// Pointer to values of SDV3:
	elem_SDV3_ptr[el].push_back(data);
      }
      std::cout << "  ... finished reading " << varName_SDV3
        << ", got values for " << elem_SDV3_ptr.size() << " elements."
        << std::endl;
      if (elem_SDV3_ptr.size() == 0) {
        std::cout << "ERROR: No values found. Stopping now." << std::endl;
        exit(2);
      }
    }

    // doc_3A_hide_block_1_end
    // ********************* some lines above not shown ******************
    // doc_3A_end    
    // doc_3B_start

    // --------------------------------------
    // --- 3 B: Interpolating / Remapping ---
    // --------------------------------------

    std::cout << "  B: Interpolating." << std::endl;

    if (read_U){
      // Interpolate U
      std::cout << "  ... interpolating U." << std::endl;
      int n_components=3;

      for (int ln=0; ln < n_GP; ln++) {
	// const int id = omp_get_thread_num();      
	// printf("This is thread %d\n", id);

	float *wght_ptr = &(n_wght_of_pt[ln][0]);
	int   *node_ptr = &(node_of_pt[ln][0]);
	// U has n=3 components:
	for (int n=0; n<n_components; n++){
	  float U_int=0.;
	  for (int i=0; i<10; i++){
	    int node = *(node_ptr+i);
	    // add displacement of i-th node times weight:
	    U_int+= *(node_U_ptr[node]+n) * *(wght_ptr+i);
	  }
	  U_vec.push_back(U_int);
	  U_cmp[n].push_back(U_int);
	}
      }
      std::cout << "  ... finished interpolating U." << std::endl;
    } // if (read_U) ends here
    
    if (read_S){
      // Interpolate S
      std::cout << "  ... interpolating S." << std::endl;
      int n_components=6;

      for (int ln=0; ln < n_GP; ln++) {
	// element of this line:
	int elem = elem_of_pt[ln];
	// pointer to weights:
	float *wght_ptr = &(e_wght_of_pt[ln][0]);
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
      std::cout << "  ... finished interpolating S." << std::endl;
    }
    if (read_SDV1){
      // Interpolate SDV1
      std::cout << "  ... interpolating " << varName_SDV1 << "." << std::endl;
      for (int ln=0; ln < n_GP; ln++) {
	// element of this line:
	int elem = elem_of_pt[ln];
	// pointer to weights:
	float *wght_ptr = &(e_wght_of_pt[ln][0]);
	float SDV1_int=0.;
	// sum over IPs:
	for (int i=0; i<4; i++){
	  // Interpolate using pointers:
	  SDV1_int+=*(elem_SDV1_ptr[elem][i]) * *(wght_ptr+i);
	}
	SDV1.push_back(SDV1_int);
      }
      std::cout << "  ... finished interpolating " << varName_SDV1 << "." << std::endl;
    }
    if (read_SDV2){
      // Interpolate SDV2
      std::cout << "  ... interpolating " << varName_SDV2 << "." << std::endl;
      for (int ln=0; ln < n_GP; ln++) {
	// element of this line:
	int elem = elem_of_pt[ln];
	// pointer to weights:
	float *wght_ptr = &(e_wght_of_pt[ln][0]);
	float SDV2_int=0.;
	// sum over IPs:
	for (int i=0; i<4; i++){
	  // Interpolate using pointers:
	  SDV2_int+=*(elem_SDV2_ptr[elem][i]) * *(wght_ptr+i);
	}
	SDV2.push_back(SDV2_int);
      }
      std::cout << "  ... finished interpolating " << varName_SDV2 << "." << std::endl;
    }
    if (read_SDV3){
      // Interpolate SDV3
      std::cout << "  ... interpolating " << varName_SDV3 << "." << std::endl;
      for (int ln=0; ln < n_GP; ln++) {
	// element of this line:
	int elem = elem_of_pt[ln];
	// pointer to weights:
	float *wght_ptr = &(e_wght_of_pt[ln][0]);
	float SDV3_int=0.;
	// sum over IPs:
	for (int i=0; i<4; i++){
	  // Interpolate using pointers:
	  SDV3_int+=*(elem_SDV3_ptr[elem][i]) * *(wght_ptr+i);
	}
	SDV3.push_back(SDV3_int);
      }
      std::cout << "  ... finished interpolating " << varName_SDV3 << "." << std::endl;
    }
    // doc_3B_end
    // doc_3C_start

    // ------------------------
    // --- 3 C: Writing vtk ---
    // ------------------------
    std::cout << "  C: " << bin_asc << " data to: " << vtk_fname << std::endl;


    // Write vtk-file:
    // parameters:
    // 1. output file
    // 2. variable descriptions
    // 3. stored variables
    // 4. grid
    int dims[] = {nx, ny, nz};
    write_to_vtk( vtk_fname, i_bin,
		  dims, n_vars, &(vardims_v[0]), &(centering_v[0]), varnames,
		  S_cmp, S_minP, S_maxP, SDV1, SDV2, SDV3, U_cmp, U_vec,
		  org, spc );
    
    // Clearing memory:
    // www.cplusplus.com/reference/vector/vector/clear/
    U_vec.clear();
    for (int i=0; i<3; i++)
      U_cmp[i].clear();
    for (int i=0; i<6; i++)
      S_cmp[i].clear();
    S_minP.clear();
    S_maxP.clear();
    SDV1.clear();
    SDV2.clear();
    SDV3.clear();
    // www.cplusplus.com/reference/map/map/clear/
    // node_U_ptr.clear(); not necessary, will be overwritten.
    elem_S_ptr.clear();
    elem_minP.clear();
    elem_maxP.clear();
    elem_SDV1_ptr.clear();
    elem_SDV2_ptr.clear();
    elem_SDV3_ptr.clear();


  } // Frame-loop-end.



  std::cout << std::endl;
  wall_end  = get_wall_time();
  std::cout << "All took: " << wall_end - wall_start << " s" << std::endl;


  odb_finalizeAPI();
  // doc_3C_end

  return 0;
}



