include ../../config.mk

BUILDDIR_LOCATION = src_builddir # or: here

SRC =  \
 stdafx.cpp \
 luafile.cpp \
 mesh.cpp \
 dof.cpp \
 csrmatrixpattern.cpp \
 boundaryCondition.cpp \
 problem.cpp \
 matElem.cpp \
 vecElem.cpp \
 ionicModel.cpp \
 bdf.cpp \
 stimulation.cpp \
 output.cpp \
 heterogeneity.cpp \
 drugCompound.cpp \

LIB = libcardioxcomp.$(LIBEXT)

default: create_lib
clean: clean_libfiles


MODELFLAGS = bidomain

ionicmodel: ionicModel.o
	g++ -o ionicmodel.so ionicmodel.o

ionicmodel.o: ionicmodel.cpp
	g++ -I/path_sundials/include -o ionicModel.o ionicModel.cpp -L/path_sundials/lib -lsundials_cvode -lsundials_nvecserial



include ../../rules.mk
