/*****************************************************************************/
/* odbAccessServer serves requests to access the Abaqus odb. Communication   */
/* is through stdin, stdout, stderr and separate data pipes.                 */
/* The rationale for separate data pipes is that the data processing on the  */
/* side of the requester my be in a separate parallel process.               */
/*****************************************************************************/

/******************************************************************************
 *****************************************************************************/
static const char
version[] = "1.01";
/* version history:
 *
 * 1.01: first version
 *
 ******************************************************************************
 *****************************************************************************/


#include <stdio.h>

#include <map>
#include <vector>
#include <string>

/* Abaqus integration */
#include <odb_API.h>

/******************************************************************************
*******************************************************************************
*** Debugging                                                               ***
*******************************************************************************
******************************************************************************/

/* Debugging */
// #define DEBUG
#undef DEBUG

// #define DEBUG_MEMORYLEAK
#undef DEBUG_MEMORYLEAK

#ifdef DEBUG
#define debug_print(fmt, ...) fprintf(stderr, fmt, ##__VA_ARGS__)
#else
#define debug_print(fmt, ...) do {} while (0)
#endif

#ifdef DEBUG_MEMORYLEAK
/* prints /proc/self/statm
 * from man 5 proc:
              Provides information about memory usage, measured in pages.  The
              columns are:

                  size       (1) total program size
                             (same as VmSize in /proc/[pid]/status)
                  resident   (2) resident set size
                             (same as VmRSS in /proc/[pid]/status)
                  share      (3) shared pages (i.e., backed by a file)
                  text       (4) text (code)
                  lib        (5) library (unused in Linux 2.6)
                  data       (6) data + stack
                  dt         (7) dirty pages (unused in Linux 2.6)
 */
void debug_print_mem(char* where) {
  FILE* f;
  char status[1000];
  int i;

  f=fopen("/proc/self/statm", "r");
  fgets(status, 1000, f);
  fclose(f);

  for (i=0; i<1000; ++i) {
    if (status[i] == '\n') {
      status[i] = 0;
    }
  }
  fprintf(stderr, "%s: current /proc/self/statm: <%s>\n", where, status);
}
#else
#define debug_print_mem(where) do {} while (0)
#endif


/******************************************************************************
*******************************************************************************
*** various service routines                                                ***
*******************************************************************************
******************************************************************************/

/* pipe communication service routines */
#include "communication.cc"

/* error handling for fatal errors */
void error(const char * fmt, ...) {
  char buffer[2048];
  va_list args;
  va_start (args, fmt);
  vsnprintf (buffer, 2047, fmt, args);
  buffer[2047] = 0;
  va_end (args);

  /* send #ERROR to the pipe */
  pipeWrite("#ERROR");
  pipeWrite(buffer);

  /* throw an exception to stop the odbAccessServer */
  throw odb_Exception(odb_TEXT_MESSAGE, atr_raw(buffer));
}


/* server functions to access the odb */
#include "odbFunctions.cc"


/******************************************************************************
*******************************************************************************
*** main server loop                                                        ***
*******************************************************************************
******************************************************************************/

int ABQmain(int argc, char** argv) {

  std::string commandStr;
  FnPtr func;
  int rc;

  initFuncMap();
  odb_initializeAPI();

  try {
    do {

      // get command from the ctrl pipe
      debug_print("odbAccessServer waiting for command\n");
      commandStr = pipeReadString();
      debug_print("odbAccessServer received command %s\n", commandStr.c_str());

      if (commandStr=="#STOP") {
        debug_print("odbAccessServer: stopping\n");
        break;
      }

      // ... std::map: if commandStr not found: insert and return NULL-ptr
      func = funcMap[commandStr];
      if (func) {
        debug_print("odbAccessServer running command %s\n", commandStr.c_str());
        rc = func();
        if (rc==0) {
          pipeWrite("#RETURN");
          debug_print("odbAccessServer finished command %s\n", commandStr.c_str());
        }
        else if (rc==1) {
          // in this case pipeException has already been used to send the
          // #EXCEPTION phrase through the pipe.
          debug_print("odbAccessServer finished command %s with an"
            " exception.\n", commandStr.c_str());
        }
        else {
          // rc==2 indicates an unrecoverable error
          error("ERROR: odbAccessServer command %s caused an unrecoverable"
                " error.\n", commandStr.c_str());
          debug_print("odbAccessServer command %s caused an unrecoverable"
                      " error.\n", commandStr.c_str());
        }          
      }
      else {
        error("odbAccessServer unknown command %s\n", commandStr.c_str());
        break;
      }
      {
        char msg[200];
        sprintf(msg, "after %s", commandStr.c_str());
        debug_print_mem(msg);
      }

    } while(1);
  }
  catch(ExceptionPipeIO& exc) {
    debug_print("Exception from %d\n", exc.label);
  }
  
  odb_finalizeAPI();
  debug_print("odbAccessServer: ending\n");
  return 0;
}
