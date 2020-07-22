# This file contains plateform-specifics settings.

#------------------------------ general options ---------------------
CONFIG_MK_VERSION = 0.2.0
MODE = release
OS_NAME = Linux
NPROC = 1
WITH_EXCEPTION = no
LOGGED = no
PRINT_MK_DEBUG = no

#------------------------------ compiler, linker (etc), flags------------------------------
#compiler
CXX =  /your_path/mpic++_compiler

#linker
LDD = ldd
#LDD = otool -L

#compiler flags
# -m64 -arch x86_64  -stdlib=libc++
CXXFLAGS_ARCH = -fPIC -fno-common
CXXFLAGS_DEBUG = -O0 -g -Wall -Wextra
CXXFLAGS_RELEASE = -O3

#preprocessor flags
CPPFLAGS_ARCH = -DCARDIOXCOMP_PETSC_NEXT # uncomment for petsc 3.4
# CPPFLAGS_ARCH = # uncomment for petsc 3.3
CPPFLAGS_DEBUG = 
CPPFLAGS_RELEASE = -DNDEBUG

#linker flags
#Note that LIBRARY_DIRS and LIBRARIES are defined in rules.mk

# redfine LIBRARIES, and invert petsc and slepc order.
LIBRARIES_12_04 = \
	$(CARDIOXCOMP_LIB) \
	$(if $(findstring $(CARDIOXCOMP_WITH_SUNDIALS),yes),$(SUNDIALS_LIB),) \
	$(PETSC_LIB) \
	$(if $(findstring $(CARDIOXCOMP_WITH_SLEPC),yes),$(SLEPC_LIB),) \
	$(SCALAPACK_LIB) \
	$(BLACS_LIB) \
	$(BLAS_LIB) \
	$(LAPACK_LIB) \
	$(SUPERLU_LIB) \
	$(MUMPS_LIB) \
	$(PARMETIS_LIB) \
	$(METIS_LIB) \
	$(LM5_LIB) \
	$(MPI_LIB) \
	$(X11_LIB) \
	$(BOOST_LIB) \
	$(LUA_LIB) \
	$(ADDITIONAL_LIBS) \

LDFLAGS = \
  $(foreach dir,$(LIBRARY_DIRS),-L$(dir)) \
  $(foreach dir,$(LIBRARY_DIRS),-Wl,-R$(dir)) \
  $(foreach lib,$(LIBRARIES_12_04),-l$(lib)) \

#command to create a library
#This is a function, that take as arguments:
#  - $1, the library name, for example 'libfoo.so'
#  - $2, files objects, for example: 'foo.o bar.o baz.o'

CREATELIB = \
  $(CXX) \
  -shared \
  -Wl,-undefined,dynamic_lookup,-soname,$(abspath $(1)) \
  -o $(1) $(2)

#CREATELIB = \
#  $(CXX) \
#  -dynamiclib \
#  -Wl,-undefined,dynamic_lookup,-install_name,$(abspath $(1)) \
#  -o $(1) $(2)

#header dependencies flags
DEPENDFLAGS = -MM
CXXDEPEND = mpicxx
#CXXDEPEND = g++


#include flags
#Note that INCLUDE_DIRS is defined in rules.mk
INCLUDES = $(foreach dir,$(INCLUDE_DIRS),-I$(dir))

#library extension
LIBEXT = so
#LIBEXT = dylib
#LIBEXT = o

#---------------------- directory of local libraries  ------------------------------
LOCAL_LIB_DIR = /your_lib_path
PETSC_DIR = $(LOCAL_LIB_DIR)/ubuntu-16/petsc/petsc-3.4.5
PETSC_LIB_DIR = $(PETSC_DIR)/$(PETSC_ARCH)/lib

#PETSC_DIR = $(LOCAL_LIB_DIR)/petsc-ubuntu16/3.4.5
#PETSC_ARCH = ubuntu-16.04-release
PETSC_LIB_DIR = $(PETSC_DIR)/ubuntu-16.04-release/lib

#------------------------------ framework ---------------------
FRAMEWORK = 

#------------------------------ blacs ---------------------
BLACS_LIB_DIR = 
BLACS_LIB = 

#------------------------------ blas ---------------------
BLAS_LIB_DIR = 
#BLAS_LIB_DIR = /usr/lib
BLAS_LIB = blas



#------------------------------ CardioXcomp ---------------------
CARDIOXCOMP_LIB_DIR = ./build_Linux_release/LibCardioXcomp
CARDIOXCOMP_LIB = cardioxcomp


#------------------------------ lua ---------------------
LUA_INC_DIR = /path_lua_include
LUA_LIB_DIR = /path_lua_lib
LUA_LIB = lua53


#------------------------------ boost ------------------------------
BOOST_INC_DIR = /usr/include
BOOST_LIB_DIR = 
BOOST_LIB = boost_regex

#------------------------------ getpot ------------------------------
GETPOT_INC_DIR = $(LOCAL_LIB_DIR)/getpot-c++

#------------------------------ lapack ---------------------
LAPACK_LIB_DIR = 
#LAPACK_LIB_DIR = /usr/lib
LAPACK_LIB = lapack

#------------------------------ lm5 ------------------------------
LM5_INC_DIR = $(LOCAL_LIB_DIR)/lm5-10.04
LM5_LIB_DIR = $(LOCAL_LIB_DIR)/lm5-10.04
LM5_LIB = mesh5

#------------------------------ metis ---------------------
#METIS_LIB_DIR = $(PETSC_LIB_DIR)
METIS_LIB_DIR = $(PETSC_LIB_DIR)/externalpackages/metis-5.0.2-p3/ubuntu-16.04-release/libmetis
METIS_INC_DIR = $(PETSC_LIB_DIR)/externalpackages/metis-5.0.2-p3/ubuntu-16.04-release/include
METIS_LIB = metis

#------------------------------ mpi ---------------------
MPI_INC_DIR = $(LOCAL_LIB_DIR)/ubuntu-16/mpich/3.0.4/include
MPI_LIB_DIR = 
MPI_LIB = 

MPI_EXEC = $(LOCAL_LIB_DIR)/ubuntu-16/mpich/3.0.4/bin/mpiexec

#------------------------------ mumps ---------------------
MUMPS_INC_DIR = $(PETSC_DIR)/externalpackages/MUMPS/include
MUMPS_LIB_DIR = $(PETSC_LIB_DIR)
MUMPS_LIB = cmumps dmumps smumps zmumps mumps_common

#------------------------------ parmetis ---------------------
PARMETIS_LIB_DIR = $(PETSC_LIB_DIR)/externalpackages/parmetis-4.0.2-p5/ubuntu-16.04-release/libparmetis
PARMETIS_INC_DIR = $(PETSC_LIB_DIR)/externalpackages/parmetis-4.0.2-p5/ubuntu-16.04-release/include


PARMETIS_LIB = parmetis

#------------------------------ petsc ------------------------------
PETSC_INC_DIR = $(PETSC_DIR)/include $(PETSC_DIR)/ubuntu-16.04-release/include
PETSC_LIB = petsc

#------------------------------ slepc ------------------------------
CARDIOXCOMP_WITH_SLEPC = no
SLEPC_DIR = $(LOCAL_LIB_DIR)/slepc/3.4.3
SLEPC_ARCH = $(PETSC_ARCH)
SLEPC_INC_DIR = $(SLEPC_DIR)/include $(SLEPC_DIR)/$(SLEPC_ARCH)/include
SLEPC_LIB_DIR = $(SLEPC_DIR)/$(SLEPC_ARCH)/lib
SLEPC_LIB = slepc

#------------------------------ sundials ------------------------------
CARDIOXCOMP_WITH_SUNDIALS = yes
SUNDIALS_DIR = $(LOCAL_LIB_DIR)/sundials-2.5.0
SUNDIALS_INC_DIR = $(SUNDIALS_DIR)/include
SUNDIALS_LIB_DIR = $(SUNDIALS_DIR)/lib
SUNDIALS_LIB = sundials_cvode sundials_nvecserial


#------------------------------ python wrapping ---------------------
WITH_PYTHON_USE_WRAPPING=yes
WITH_PYTHON_GENERATE_WRAPPING=yes
PYTHON = /usr/bin/python
SWIG = /usr/bin/swig
SWIG_OPT = -python -c++ -outcurrentdir -I$(CARDIOXCOMP_ROOT)/Source -I$(CARDIOXCOMP_ROOT)/Python/swig_lib
CXXFLAGS_PYWRAP = -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -Wno-write-strings
LDFLAGS_PYWRAP = -g -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions
PYTHON_INC_DIR = /usr/include/python2.7
INCLUDES_PYWRAP = $(foreach dir,$(INCLUDE_PYWRAP_DIRS),-I$(dir))

#------------------------------ scalapack ---------------------
SCALAPACK_LIB_DIR = $(PETSC_LIB_DIR)
SCALAPACK_LIB = scalapack

#------------------------------ superlu ---------------------
SUPERLU_LIB_DIR = $(PETSC_LIB_DIR)
SUPERLU_LIB = superlu_dist_3.3

#------------------------------ X11 ---------------------
X11_LIB_DIR = /usr/lib/x86_64-linux-gnu/
X11_LIB = X11

#--------------------------- additional_libs ---------------------
ADDITIONAL_LIBS = \
    m \
    gfortran \
    pord

ADDITIONAL_LIBS_DIR = ./build_Linux_release/LibCardioXcomp/

#--------------------------- others stuffs ---------------------
# Absolute path of the directory containing this file (config.mk).
# Works even if the file is included by another/path/Makefile.
# (in which case, $(CURDIR) is incorrect because it returns another/path)
CARDIOXCOMP_ROOT := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
BUILDDIR_LOCATION = here

#for people using git-svn
#CARDIOXCOMP_ARCH_HOOK = _$(shell git branch --no-color 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/\1/')

