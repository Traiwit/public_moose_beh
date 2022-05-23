// contains functions for:

// Z: Assert
// A: Float Swap
// B: Write to vtk using visit_writer
// C: Measuring time
// D: Checking filesize
// E: Printing some info to screen


// ---------
// --- Z --- Assert
// ---------

#define ASSERT(condition) {						\
    if(!(condition)){							\
      std::cerr << std::endl << std::endl				\
		<< "ASSERT FAILED: "					\
		<< #condition << " @ "					\
		<< __FILE__ << " ("					\
		<< __LINE__ << ")"					\
		<< std::endl << std::endl;				\
      exit(0);								\
    }									\
  }


// check that there is a field "fieldName" in the repository and return
// a reference to it.
odb_FieldOutput& checkGetField(
			       odb_FieldOutputRepository& fieldRepos,
			       const char* fieldName,
			       int frameIdx )
{
  if (! fieldRepos.isMember(fieldName) ) {
	std::cerr << std::endl << std::endl
		  << "ERROR: Field " << fieldName << " does not exist in frame "
		  << frameIdx
		  << std::endl << std::endl;
	exit(11);
  }
  return fieldRepos[fieldName];
}

// ---------
// --- A --- Float Swap
// ---------

// endianness
float FloatSwap( float f )
{    
   union
   {
      float f;
      unsigned char b[4];
   } dat1, dat2;

   dat1.f = f;
   dat2.b[0] = dat1.b[3];
   dat2.b[1] = dat1.b[2];
   dat2.b[2] = dat1.b[1];
   dat2.b[3] = dat1.b[0];
   return dat2.f;
}

// ---------
// --- B --- Write to vtk using visit_writer
// ---------

int write_to_vtk(
		 const char * vtk_fname, int i_bin, 
		 int *dims, 
		 int n_vars, int *vardim, int *centering,
		 const char * const *varnames, 
		 std::vector<float> S_cmp[6],
		 std::vector<float> S_minP,
		 std::vector<float> S_maxP,
		 std::vector<float> SDV1,
		 std::vector<float> SDV2,
		 std::vector<float> SDV3,
		 std::vector<float> U_cmp[3],
		 std::vector<float> U_vec,
		 // float **vars,
		 float *org, float *spc
		 ){

  // Write stress tensor to vtk?
  bool write_S_manually=false;

  // print info about output file:
  // std::cout << "     --> Output file: \"" << vtk_fname << "\"." << std::endl;


  // Picking variables according to variable names:
  std::vector<float *> vars;
  for (int i=0; i<n_vars; i++){

    if      (strcmp(varnames[i],   "S_11") == 0) vars.push_back(&(S_cmp[0][0]));
    else if (strcmp(varnames[i],   "S_22") == 0) vars.push_back(&(S_cmp[1][0]));
    else if (strcmp(varnames[i],   "S_33") == 0) vars.push_back(&(S_cmp[2][0]));
    else if (strcmp(varnames[i],   "S_12") == 0) vars.push_back(&(S_cmp[3][0]));
    else if (strcmp(varnames[i],   "S_13") == 0) vars.push_back(&(S_cmp[4][0]));
    else if (strcmp(varnames[i],   "S_23") == 0) vars.push_back(&(S_cmp[5][0]));
    else if (strcmp(varnames[i], "S_MIN_PRINCIPAL") == 0) vars.push_back(&(S_minP[0]  ));
    else if (strcmp(varnames[i], "S_MAX_PRINCIPAL") == 0) vars.push_back(&(S_maxP[0]  ));
    else if (strcmp(varnames[i], varName_SDV1.c_str()) == 0) vars.push_back(&(SDV1[0]  ));
    else if (strcmp(varnames[i], varName_SDV2.c_str()) == 0) vars.push_back(&(SDV2[0]  ));
    else if (strcmp(varnames[i], varName_SDV3.c_str()) == 0) vars.push_back(&(SDV3[0]  ));
    else if (strcmp(varnames[i], varName_U_1.c_str()) == 0) vars.push_back(&(U_cmp[0][0]));
    else if (strcmp(varnames[i], varName_U_2.c_str()) == 0) vars.push_back(&(U_cmp[1][0]));
    else if (strcmp(varnames[i], varName_U_3.c_str()) == 0) vars.push_back(&(U_cmp[2][0]));
    else if (strcmp(varnames[i], varName_U.c_str()) == 0) vars.push_back(&(U_vec[0]   ));
    else if (strcmp(varnames[i],      "S") == 0){
      vars.push_back(NULL);
      write_S_manually=true;
    }
    else{
      std::cout << std::endl;
      std::cout << "Unknown variable: " << varnames[i];
      std::cout << std::endl;
      std::cout << "Exiting.";
      std::cout << std::endl;
      exit(0);
    }
  }

  
  write_structured_points_mesh(vtk_fname, i_bin, 
			       dims, n_vars, &(vardim[0]), centering, varnames, 
			       &(vars[0]), org, spc);    

  if (write_S_manually){
    std::cout << "     --> Writing S manually instead..." << std::endl;
    if (i_bin==0){
      std::ofstream vtk;
      vtk.open(vtk_fname, std::ios_base::app); 
      vtk << "TENSORS S float" << "\n";
      for (int ln=0; ln < dims[0]*dims[1]*dims[2]; ln++) {
	vtk << S_cmp[0][ln] << " " << S_cmp[3][ln] << " " << S_cmp[4][ln] << " " << "\n";
	vtk << S_cmp[3][ln] << " " << S_cmp[1][ln] << " " << S_cmp[5][ln] << " " << "\n";
	vtk << S_cmp[4][ln] << " " << S_cmp[5][ln] << " " << S_cmp[3][ln] << " " << "\n";
      }
      vtk.close();
    }
    else{
      float sxx,syy,szz,sxy,sxz,syz;
      char s[256]; 
      FILE *fp;
      fp = fopen(vtk_fname,"a+w"); 
      sprintf(s, "TENSORS S float\n");              
      fwrite(s, sizeof(char), strlen(s), fp); 
      for (int ln=0; ln < dims[0]*dims[1]*dims[2]; ln++) {
	sxx = FloatSwap(S_cmp[0][ln]);
	syy = FloatSwap(S_cmp[1][ln]);
	szz = FloatSwap(S_cmp[2][ln]);
	sxy = FloatSwap(S_cmp[3][ln]);
	sxz = FloatSwap(S_cmp[4][ln]);
	syz = FloatSwap(S_cmp[5][ln]);
	fwrite((void *)&sxx, sizeof(float), 1, fp);
	fwrite((void *)&sxy, sizeof(float), 1, fp);
	fwrite((void *)&sxz, sizeof(float), 1, fp);
	fwrite((void *)&sxy, sizeof(float), 1, fp);
	fwrite((void *)&syy, sizeof(float), 1, fp);
	fwrite((void *)&syz, sizeof(float), 1, fp);
	fwrite((void *)&sxz, sizeof(float), 1, fp);
	fwrite((void *)&syz, sizeof(float), 1, fp);
	fwrite((void *)&szz, sizeof(float), 1, fp);
      }
      
      fclose(fp);
      std::cout << "     --> Done.\n";
    }
  }

  // Clear memory:
  vars.clear();
  return 0;
}

// ----------------------------------------------
// --- Writing vtk (ascii or binary) manually ---
// ----------------------------------------------
//
// if (read_U){
//   const char *debug_vtk_fname = "debug_U.vtk";
//   std::cout << "  D: debug data to: " << debug_vtk_fname << std::endl;
//   if (i_bin){ //binary
// 	float val;
// 	char s[256]; 
// 	FILE *fp;
// 	fp = fopen(debug_vtk_fname,"w");
// 	// BINARY:
// 	sprintf(s,"# vtk DataFile Version 2.0\n");
// 	fwrite(s,sizeof(char),strlen(s), fp);
// 	sprintf(s,"written manually with byteswap\n");
// 	fwrite(s,sizeof(char),strlen(s), fp);
// 	sprintf(s,"BINARY\n");
// 	fwrite(s,sizeof(char),strlen(s), fp);
// 	sprintf(s,"DATASET STRUCTURED_POINTS\n");
// 	fwrite(s,sizeof(char),strlen(s), fp); 
// 	sprintf(s,"DIMENSIONS %d %d %d\n",nx,ny,nz);     
// 	fwrite(s,sizeof(char),strlen(s), fp); 
// 	sprintf(s,"ORIGIN %f %f %f\n",org_x,org_y,org_z);
// 	fwrite(s,sizeof(char),strlen(s), fp);
// 	sprintf(s,"SPACING %f %f %f\n",spc_x,spc_y,spc_z);
// 	fwrite(s,sizeof(char),strlen(s), fp);
// 	sprintf(s,"CELL_DATA %d\n",(nx-1)*(ny-1)*(nz-1));
// 	fwrite(s,sizeof(char),strlen(s), fp); 
// 	sprintf(s,"POINT_DATA %d\n",nx*ny*nz);
// 	fwrite(s,sizeof(char),strlen(s), fp);
// 	//
// 	sprintf(s,"SCALARS U_X float\n");
// 	fwrite(s,sizeof(char),strlen(s), fp); 
// 	sprintf(s,"LOOKUP_TABLE default\n");
// 	fwrite(s,sizeof(char),strlen(s), fp); 
// 	for (int ln=0; ln < n_GP; ln++) {
// 	  val = FloatSwap(U_cmp[0][ln]);
// 	  fwrite((void *)&val, sizeof(float),1,fp);
// 	}
// 	sprintf(s,"SCALARS U_Y float\n");
// 	fwrite(s,sizeof(char),strlen(s), fp); 
// 	sprintf(s,"LOOKUP_TABLE default\n");
// 	fwrite(s,sizeof(char),strlen(s), fp); 
// 	for (int ln=0; ln < n_GP; ln++) {
// 	  val = FloatSwap(U_cmp[1][ln]);
// 	  fwrite((void *)&val, sizeof(float),1,fp);
// 	}
// 	sprintf(s,"SCALARS U_Z float\n");
// 	fwrite(s,sizeof(char),strlen(s), fp);
// 	sprintf(s,"LOOKUP_TABLE default\n");
// 	fwrite(s,sizeof(char),strlen(s), fp); 
// 	for (int ln=0; ln < n_GP; ln++) {
// 	  val = FloatSwap(U_cmp[2][ln]);
// 	  fwrite((void *)&val, sizeof(float),1,fp);
// 	}
// 	fclose(fp);
//   }
//   else{
// 	std::ofstream vtkstream(debug_vtk_fname);
// 	vtkstream << "# vtk DataFile Version 2.0" << "\n";
// 	vtkstream << "written manually with ofstream" << "\n";
// 	vtkstream << "ASCII" << "\n";
// 	vtkstream << "DATASET STRUCTURED_POINTS" << "\n";
// 	vtkstream << "DIMENSIONS "  << nx << " " <<    ny << " " <<    nz << "\n";
// 	vtkstream << "ORIGIN "   << org_x << " " << org_y << " " << org_z << "\n";
// 	vtkstream << "SPACING "  << spc_x << " " << spc_y << " " << spc_z << "\n";
// 	vtkstream << "CELL_DATA " << (nx-1)*(ny-1)*(nz-1) << "\n";
// 	vtkstream << "POINT_DATA "<<       nx*ny*nz       << "\n";
// 	// write U_X:
// 	vtkstream << "SCALARS U_X float"    << "\n";
// 	vtkstream << "LOOKUP_TABLE default" << "\n";
// 	for (int ln=0; ln < n_GP; ln++) {
// 	  vtkstream << U_cmp[0][ln] << " \n";
// 	}
// 	// write U_Y:
// 	vtkstream << "SCALARS U_Y float"    << "\n";
// 	vtkstream << "LOOKUP_TABLE default" << "\n";
// 	for (int ln=0; ln < n_GP; ln++) {
// 	  vtkstream << U_cmp[1][ln] << " \n";
// 	}
// 	// write U_Z:
// 	vtkstream << "SCALARS U_Z float"    << "\n";
// 	vtkstream << "LOOKUP_TABLE default" << "\n";
// 	for (int ln=0; ln < n_GP; ln++) {
// 	  vtkstream << U_cmp[2][ln] << " \n";
// 	}
// 	vtkstream.close();
//   }
// }






// ---------
// --- C --- Measuring time
// ---------


//  Posix/Linux
#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
  return (double)clock() / CLOCKS_PER_SEC;
}

// ---------
// --- D --- Checking filesize
// ---------

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}


// ---------
// --- E --- Printing some info to screen
// ---------

void print_info(){
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "  Contents of the parameter file:";
  std::cout << std::endl;
  std::cout << "  1:         Path to odb";
  std::cout << std::endl;
  std::cout << "  2:         Comma-separated list, e.g. U,U_Z,S_XX,S_YZ";
  std::cout << std::endl;
  std::cout << "  3:         Word binary or ascii - specifying vtk-file-format";
  std::cout << std::endl;
  std::cout << "  4:         Path to binary file containing nodes-weights";
  std::cout << std::endl;
  std::cout << "  5:         Path to binary file containing elems-weights";
  std::cout << std::endl;
  std::cout << "  6:         Comma-separated list  nx,ny,nz,x,y,z,sx,sy,sz";
  std::cout << "  - specifying grid dimensions, origin and spacing";
  std::cout << std::endl;
  std::cout << "  7:         Comma-separated tupel, e.g. 99,prj_xxx_99.vtk";
  std::cout << "  - specifying frame-index and vtk-filename";
  std::cout << std::endl;
  std::cout << "  8 etc.:    similar to 7, next tuple";
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "Examples:    http://be.fiziko.de/read_data/docs.html#benutzung ";
  std::cout << std::endl;
  std::cout << std::endl;

}
