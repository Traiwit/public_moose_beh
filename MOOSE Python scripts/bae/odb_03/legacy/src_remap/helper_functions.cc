// contains functions for:

// Z: Assert
// A: converting strings to values
// B: Write to vtk
// C: Measuring time
// D: Checking filesize
// E: Printing some info to screen

#include "helper_functions.h"
#include <vector>
#include <algorithm>
#include <iterator>

// ---------
// --- Z --- Assert
// ---------

#define ASSERT(condition) {						\
    if(!(condition)){							\
      cerr << endl << endl						\
	   << "ASSERT FAILED: "						\
	   << #condition << " @ "					\
	   << __FILE__ << " ("						\
	   << __LINE__ << ")"						\
	   << endl << endl;						\
      exit(100);							\
    }									\
  }

// ---------
// --- A --- converting strings to values
// ---------

int stringToInteger(std::string& s) {

  int number;
  std::stringstream str(s);
  std::string overflow;
  str >> number;
  if (str && !(str >> overflow)) {
    return number;
  }
  else {
    cerr << "Error: can't convert " << s
	 << " to an integer." << endl;
    exit(10);
  }
}

// ---------
// --- A1 --- pretty print vector
// ---------

// pretty print std::vector objects
// code from the following source, slightly modified
// http://stackoverflow.com/questions/15435313/pretty-print-a-stdvector-in-c
template<typename T>
std::ostream & operator<<(std::ostream & os, std::vector<T> vec) {
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, ","));
  return os;
}

template<class T, size_t N>
std::string toStr(T (& a)[N]) {
  std::stringstream res;
  for (int i=0; i<N; ++i) {
    if (i>0) res << ",";
    res << a[i];
  }
  return res.str();
}

// ---------
// --- B --- Write to vtk
// ---------

// ----------------------------------------------
// --- Writing vtk (ascii or binary) manually ---
// ----------------------------------------------
//

/* ****************************************************************************
 *
 *  Purpose:
 *      Determines if the machine is little-endian.  If so, then,
 *      convert 32bit float and integer vectors to be big-endian.
 *
 *  Note:   This assumes a vector of 4 bytes words.
 *
 *  derived from force_big_endian in visit_writer.c
 *  from Hank Childs, September 3, 2004
 * 
 * ************************************************************************* */

static void convert_to_big_endian(void* data, unsigned int N)
{
    static int doneTest = 0;
    static int shouldSwap = 0;
    if (!doneTest)
    {
        int tmp1 = 1;
        unsigned char *tmp2 = (unsigned char *) &tmp1;
        if (*tmp2 != 0)
            shouldSwap = 1;
        doneTest = 1;
    }

    if (shouldSwap) {
      unsigned char* bytes = (unsigned char*) data;
      for (int i=0; i<N; ++i)
	{
	  unsigned char tmp = bytes[0];
	  bytes[0] = bytes[3];
	  bytes[3] = tmp;
	  tmp = bytes[1];
	  bytes[1] = bytes[2];
	  bytes[2] = tmp;

	  // next word
	  bytes += 4;
	}
    }
}

#include <stdio.h>
int write_to_vtk(
    const char * vtk_fname, // file name
    int binary_flag,        // binary flag
    const MeshStructuredPoints& grid,
    int n_vars,             // # of requested vars
    int *vardim,            // per var: 1=scalar, 3=vec, 9=tensor
    ListOfStrings& varNames, // ...
    float** dataArrays
    ){

  int numGridPoints = grid.size();
  
  if (binary_flag) {
    // convert to big endian
    for (int i=0; i<n_vars; ++i)
      convert_to_big_endian(dataArrays[i], numGridPoints);
  }

  FILE *fp;
  fp = fopen(vtk_fname, "w");

  // header
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"some comment here\n");
  if (binary_flag)
    fprintf(fp,"BINARY\n");
  else
    fprintf(fp,"ASCII\n");
  fprintf(fp,"DATASET STRUCTURED_POINTS\n");
  fprintf(fp,"DIMENSIONS %d %d %d\n",
	  grid.dims[0],grid.dims[1],grid.dims[2]);     
  fprintf(fp,"ORIGIN %f %f %f\n", grid.org[0],grid.org[1],grid.org[2]);
  fprintf(fp,"SPACING %f %f %f\n", grid.spc[0],grid.spc[1],grid.spc[2]);
  fprintf(fp,"CELL_DATA %d\n",
	  (grid.dims[0]-1)*(grid.dims[1]-1)*(grid.dims[2]-1));
  fprintf(fp,"POINT_DATA %d\n", numGridPoints);
  // write all data fields
  for (int i=0; i<n_vars; ++i) {
    fprintf(fp,"SCALARS %s float\n", varNames[i].c_str());
    fprintf(fp,"LOOKUP_TABLE default\n");
    if (binary_flag) {
      fwrite((void*)dataArrays[i], sizeof(float), numGridPoints,fp);
    }
    else {
      for (int j=0; j<numGridPoints; ++j) {
	fprintf(fp, "%g\n", dataArrays[i][j]);
      }
    }
  }
  fclose(fp);
}


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
  cout << endl << endl
       << "  1:         Path to odb" << endl
       << "  2:         Path to interpolation file" << endl
       << "  3:         Comma-separated list, e.g. U,U_3,S_11,S_23" << endl
       << "  4:         Word binary or ascii - specifying vtk-file-format"<<endl
       << "  5:         Comma-separated tupel, e.g. 99,prj_xxx_99.vtk"
       << "  - specifying frame-index and vtk-filename" << endl
       << "  6 etc.:    similar to 5, next tuple"
       << endl << endl
       << "Examples:    http://be.fiziko.de/read_data/docs.html#benutzung "
       << endl << endl;
}
