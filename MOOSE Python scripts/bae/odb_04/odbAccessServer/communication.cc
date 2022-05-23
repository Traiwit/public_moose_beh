/*****************************************************************************
 * communication through pipes: read and write                               *
 *****************************************************************************/

// for basic I/O using integer file handles
// needed to read from and write to pipes
#include <unistd.h>
#include <string>

#include <stdio.h>
#include <stdarg.h>  // for pipeException - varargs

class ExceptionPipeIO {
 public:
  int label;
  ExceptionPipeIO(int label_arg) : label(label_arg) {}
};

/**********************
 * read from the pipe *
 *********************/
// Note: A templated function might be more elegant but is less convenient
// to type.
int pipeReadInt() {
  int v;
  int rc = read(0, &v, sizeof(int));
  if (rc!=sizeof(int)) throw ExceptionPipeIO(1);
  return v;
}

float pipeReadFloat() {
  float v;
  int rc = read(0, &v, sizeof(float));
  if (rc!=sizeof(float)) throw ExceptionPipeIO(1);
  return v;
}

std::string pipeReadString() {

  int length = pipeReadInt();

  char str[length+1];
  
  int rc = read(0, str, length);
  if (rc!=length) throw ExceptionPipeIO(1);
  str[length] = 0;  // append end-marker
  return str;
}

/********************
* write to the pipe *
********************/

/*--- write a string */
void pipeWrite (const std::string& val) {

  int length = val.size();
  int rc = write(1, &length, sizeof(int));
  if (rc!=sizeof(int)) throw ExceptionPipeIO(2);
  // debug_print("pipeWrite string wrote size: %d\n", length);
  
  rc = write(1, val.c_str(), length);
  if (rc!=length) throw ExceptionPipeIO(2);
  // debug_print("pipeWrite string wrote string: %s\n", val.c_str());
}

// template <typename T>
// void pipeWrite(const T val) {
//   int rc = write(1, &val, sizeof(T));
//   if (rc!=sizeof(T)) throw ExceptionPipeIO(2);
//   // debug_print("pipeWrite <T> wrote something...\n");
// }

/*--- write an integer */
void pipeWrite(const int val) {
  int rc = write(1, &val, sizeof(int));
  if (rc!=sizeof(int)) throw ExceptionPipeIO(2);
  // debug_print("pipeWrite <int> wrote something...\n");
}

/*--- write a float */
void pipeWrite(const float val) {
  int rc = write(1, &val, sizeof(float));
  if (rc!=sizeof(float)) throw ExceptionPipeIO(2);
  // debug_print("pipeWrite <float> wrote something...\n");
}

/*--- write an array of floats */
// note this does not write the size but only the actual array data
// size must be transferred separately
void pipeWrite(const float* const val, const int size) {
  int rc = write(1, val, sizeof(float)*size);
  if (rc!=sizeof(float)*size) throw ExceptionPipeIO(2);
  // debug_print("pipeWrite <float*> wrote something...\n");
}

/*--- write an array of integers */
// note this does not write the size but only the actual array data
// size must be transferred separately
void pipeWrite(const int* const val, const int size) {
  int rc = write(1, val, sizeof(int)*size);
  if (rc!=sizeof(int)*size) throw ExceptionPipeIO(2);
  // debug_print("pipeWrite <int*> wrote something...\n");
}

/*--- write an array of bytes/chars */
// note this does not write the size but only the actual array data
// size must be transferred separately
void pipeWrite(const char* const val, const int size) {
  int rc = write(1, val, size);
  if (rc!=size) throw ExceptionPipeIO(2);
  // debug_print("pipeWrite <char*> wrote something...\n");
}

/*--- send an exception to the pipe *
 *
 * An Exception will be sent through the pipe as three strings:
 * 1. "#EXCEPTION"
 * 2. the Python type name of the exception, e.g. "ValueError"
 * 3. the error message to be passed to the constructor of the exception.
 */
void pipeException(const char* const excString, const char* const fmt, ...) {
  char buffer[2048];
  va_list args;
  va_start (args, fmt);
  vsnprintf (buffer, 2047, fmt, args);
  buffer[2047] = 0;
  va_end (args);

  /* send #ERROR to the pipe */
  pipeWrite("#EXCEPTION");
  pipeWrite(excString);
  pipeWrite(buffer);
}

