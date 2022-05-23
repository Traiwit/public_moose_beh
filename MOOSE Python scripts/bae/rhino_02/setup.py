"""Description for the rhino_02/rhino_binary extension module
to be imported by the global setup.py script


TODO:

Debian: apt-get install uuid-dev
SuSE: yast ... install  libuuid-devel
CentOS: [if not allready done: yum install epel-release] yum install libuuid-devel

implement check in setup.py: uuid installed? if not issue error message with
tip how to install uuid-dev...

// The uuid_generate function creates a new universally unique identifier (UUID). The uuid will be generated based on high-quality randomness from /dev/urandom, if available. If it is not available, then uuid_generate will use an alternative algorithm which uses the current time, the local ethernet MAC address (if available), and random data generated using a pseudo-random generator. 
// The newly created UUID is returned in the memory location pointed to by out.

#include <uuid/uuid.h>
void uuid_generate(uuid_t out);


"""

import os.path
from glob import iglob
from distutils.core import Extension

#-------------------------------------------
#--- create a list of extension modules

opennurbsBaseDir = "bae/rhino_02/opennurbs"
opennurbsSources = ['bae/rhino_02/rhino_binary.cc', ]
opennurbsSources.extend(
    x for x in iglob(os.path.join(opennurbsBaseDir, "*.cpp"))
    if os.path.basename(x) not in "opennurbs/opennurbs_gl.cpp")
opennurbsSources.extend(
    x for x in iglob(os.path.join(opennurbsBaseDir, "*.c")))
opennurbsSources.extend(
    x for x in iglob(os.path.join(opennurbsBaseDir, "zlib", "*.c")))

# # for the example
# opennurbsSources.append(os.path.join(
#     opennurbsBaseDir, "example_userdata", "example_ud.cpp"))

ext_modules = [
    Extension(
        name="bae.rhino_02.rhino_binary",
        sources=opennurbsSources,
        include_dirs=[opennurbsBaseDir, os.path.join(opennurbsBaseDir, 'zlib')],
        extra_compile_args=["-fpermissive", "-w",],
        libraries=["uuid",],
        )
    ]
