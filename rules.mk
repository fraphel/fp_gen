# This file contains all rules to be used in Makefiles.
#
#     some rules
#     include ../../rules.mk
#
#   It would be a bad idea to change to:
#     some variables
#     some rules
#     include ../../rules.mk  #rules.mk include config.mk
# 
# Indeed, it's common and convenient that targets and prerequisites of 'some rules'
# are constructed using variables defined in config.mk, and 
# contrary to variables, targets and prerequisites are evaluated immediately,
# i.e before config.mk is included.
SHELL:=/bin/bash
#----------------------------------- string utility functions -------------------------------------
# test if 2 stripped string are equal
# Example:
# RES = $(call eqstring,foo,foo)  # RES is yes
# RES = $(call eqstring,foo,bar)  # RES is empty string
# RES = $(call eqstring,foo,foo ) # RES is yes
eqstring = $(if $(and $(findstring $(strip $1),$(strip $2)),$(findstring $(strip $2),$(strip $1))),yes,)

# test if a list contains a word
# Example:
# LIST = one two three four
# RES = $(call list_contains_word,$(LIST),one)  # RES is yes
# RES = $(call list_contains_word,$(LIST),two)  # RES is yes
# RES = $(call list_contains_word,$(LIST),five) # RES is empty string
# RES = $(call list_contains_word,$(LIST),hree) # RES is empty string
# LIST = one two three four one one onenotmatch
# RES = $(call list_contains_word,$(LIST),hree) # RES is yes yes yes
list_contains_word = $(strip $(foreach w,$1,$(call eqstring,$(w),$2)))

# test if a list of goals contains $(MAKEGOAL)
# Example:
# GOALOPT_FOO_YES_FOR_GOALS = one two three four
# GOALOPT_FOO = $(call yes_if_makegoal_in_list,$(GOALOPT_FOO_YES_FOR_GOALS))#no trailing space
# if make is invoked with 'make one', GOALOPT_FOO is 'yes'.
# if make is invoked with 'make five', GOALOPT_FOO is 'no'.
yes_if_makegoal_in_list = $(if $(call list_contains_word,$1,$(MAKEGOAL)),yes,no)

# test if a list of goals does not contain $(MAKEGOAL)
# Example:
# GOALOPT_FOO_NO_FOR_GOALS = one two three four
# GOALOPT_FOO = $(call yes_if_makegoal_in_list,$(GOALOPT_FOO_NO_FOR_GOALS))#no trailing space
# if make is invoked with 'make one', GOALOPT_FOO is 'no'.
# if make is invoked with 'make five', GOALOPT_FOO is 'yes'.
no_if_makegoal_in_list = $(if $(call list_contains_word,$1,$(MAKEGOAL)),no,yes)

# test if a variable is defined (ie, does not contain empty string)
# Example 1:
#    FOO=foo
#    $(call assert-variable-defined,FOO) # no error
# Example 2:
#    $(call assert-variable-defined,FOO) # an error is raised
define assert-variable-defined
$(if $(strip $($1)),,$(error Variable '$1' must be defined))
endef

# Raise an error is a variable is not an integer
# Example:
#     FOO=123
#     $(call assert-is-integer,FOO)
# do nothing, but:
#     FOO=a
#     $(call assert-is-integer,FOO)
# stop the program and print an error message
define assert-is-integer
$(if $(findstring True,$(shell python -c "import re; print re.compile('^\d+$$').match('$($(1))')!=None")),,$(error Variable $1 is '$($1)', but expected an integer. ))
endef

# Raise an error is config.mk version is not >= to give value
# Example:
#     $(call assert-config-mk-version-ge,3,5,2,2,8,4,$(ERROR_MSG))
# do nothing, because 3.5.2 >= 2.8.4, but:
#     $(call assert-config-mk-version-ge,3,5,2,4,8,4,$(ERROR_MSG))
# stop the program and print $(ERROR_MSG), because 3.5.2 >= 4.8.4 is false.
define assert-config-mk-version-ge
$(if $(findstring True,$(shell python -c "print ($1,$2,$3) >= ($4,$5,$6)")),,$(error $7))
endef

debug_flow = $(if $(call eqstring,yes,$(PRINT_MK_DEBUG)), $(info MAKE:DEBUG: -> $1),)
debug_var = $(if $(call eqstring,yes,$(PRINT_MK_DEBUG)), $(info MAKE:DEBUG: @@ $1 is '$($1)'),)

info_msg = $(info MAKE:INFO: $1)
info_var = $(info MAKE:INFO: $1 is '$($1)')

error_var = $(error $1 has invalid value '$($1)')

$(call debug_flow,--------------------------------- Entering 'rules.mk' file. ---------------------------------)
$(call debug_var,MAKECMDGOALS)
$(call debug_var,CURDIR)
$(call debug_var,MODE)
$(call debug_var,CARDIOXCOMP_ROOT)

# limit number of goals to 0 or 1
ifneq "$(words $(MAKECMDGOALS))" "0"
ifneq "$(words $(MAKECMDGOALS))" "1"
$(error Zero or one goal supported, got $(words $(MAKECMDGOALS)))
endif
endif
MAKEGOAL = $(MAKECMDGOALS)

################################ check config.mk version ###############################
REQUIRED_CONFIG_MK_VERSION_MAJOR = 0
REQUIRED_CONFIG_MK_VERSION_MINOR = 2
REQUIRED_CONFIG_MK_VERSION_PATCH = 0

CONFIG_MK_VERSION_MAJOR = $(word 1,$(subst ., ,$(CONFIG_MK_VERSION)))
CONFIG_MK_VERSION_MINOR = $(word 2,$(subst ., ,$(CONFIG_MK_VERSION)))
CONFIG_MK_VERSION_PATCH = $(word 3,$(subst ., ,$(CONFIG_MK_VERSION)))

$(call assert-is-integer,CONFIG_MK_VERSION_MAJOR)
$(call assert-is-integer,CONFIG_MK_VERSION_MINOR)
$(call assert-is-integer,CONFIG_MK_VERSION_PATCH)

define MSG_OLD_VERSION_ERROR 
Sorry, your config.mk version is too old, and you will have to update it.

0. archive your old config.mk: 
       cp config.mk config_old.mk

1. choose an example of config.mk among $(CARDIOXCOMP_ROOT)/config/config-XXX.mk

2. create a new config.mk:
       cp $(CARDIOXCOMP_ROOT)/config/config-XXX.mk config.mk

3. edit the new config.mk and set values to fit your configuration (you
   will find the good values in the config_old.mk)

   some indications:
       EXT_PACK is now LOCAL_LIB_DIR
       XXX_INC_DIR is the directory containing headers xxx.h ...
       XXX_LIB_DIR is the directory containing libraries libxxx.a, libxxx.so ...
       XXX_LIB is the part xxx of the library libxxx.a or libxxx.so or libxxx.dylib
			 for example: LAPACK_LIB = lapack

4. check the config.mk with the command:
       make checkconfig

5. fix possible errors

6. If you have trouble, please ask help at felisce-user@lists.gforge.inria.fr

Thanks!

David Froger
endef

define CONFIG_MK_VERSION_UPDATE_MSG_ERROR
config.mk need update.

==============================================================================
Your config.mk file needs to be updated, please follow instructions in file:
$(CARDIOXCOMP_ROOT)/config/CHANGELOG

Your config.mk version    : $(CONFIG_MK_VERSION_MAJOR).$(CONFIG_MK_VERSION_MINOR).$(CONFIG_MK_VERSION_PATCH)
Current config.mk version : $(REQUIRED_CONFIG_MK_VERSION_MAJOR).$(REQUIRED_CONFIG_MK_VERSION_MINOR).$(REQUIRED_CONFIG_MK_VERSION_PATCH)

If you have trouble, please ask help at felisce-user@lists.gforge.inria.fr
Thanks.
==============================================================================


endef

ifeq '$(CONFIG_MK_VERSION_MAJOR)' ''
$(info $(MSG_OLD_VERSION_ERROR))
$(error config.mk needs update)
endif

$(call assert-config-mk-version-ge,$(CONFIG_MK_VERSION_MAJOR),$(CONFIG_MK_VERSION_MINOR),$(CONFIG_MK_VERSION_PATCH),$(REQUIRED_CONFIG_MK_VERSION_MAJOR),$(REQUIRED_CONFIG_MK_VERSION_MINOR),$(REQUIRED_CONFIG_MK_VERSION_PATCH),$(CONFIG_MK_VERSION_UPDATE_MSG_ERROR))

$(call debug_flow,Passed config.mk version checks)

################################ release/debug ###############################
ifeq "$(MODE)" "debug"
PETSC_ARCH = $(PETSC_ARCH_DEBUG)
CXXFLAGS_ALL =  $(CXXFLAGS_ARCH) $(CXXFLAGS_DEBUG) $(CXXFLAGS)
CPPFLAGS_ALL =  $(CPPFLAGS_ARCH) $(CPPFLAGS_DEBUG) $(CPPFLAGS)

else ifeq "$(MODE)" "release"
PETSC_ARCH = $(PETSC_ARCH_RELEASE)
CXXFLAGS_ALL =  $(CXXFLAGS_ARCH) $(CXXFLAGS_RELEASE) $(CXXFLAGS)
CPPFLAGS_ALL =  $(CPPFLAGS_ARCH) $(CPPFLAGS_RELEASE) $(CPPFLAGS)

else
$(call error_var,MODE)

endif

$(call debug_var,PETSC_ARCH)
$(call debug_var,CXXFLAGS_ALL)
$(call debug_var,CPPFLAGS_ALL)
$(call debug_var,BUILDDIR_LOCATION)

$(call debug_flow,Set variables according to release/debug)

################################ optional libraries ###############################

ifeq "$(CARDIOXCOMP_WITH_SLEPC)" "yes"
CPPFLAGS_ALL +=  -DCARDIOXCOMP_WITH_SLEPC
endif

ifeq "$(CARDIOXCOMP_WITH_SUNDIALS)" "yes" 	 
CPPFLAGS_ALL +=  -DCARDIOXCOMP_WITH_SUNDIALS 	 
endif

ifeq "$(CARDIOXCOMP_WITH_LIBXFM)" "yes" 	 
CPPFLAGS_ALL +=  -DCARDIOXCOMP_WITH_LIBXFM 	 
endif

############################## build base directory #############################
ifeq "$(WITH_EXCEPTION)" "yes"

CARDIOXCOMP_ARCH = $(OS_NAME)_$(MODE)$(CARDIOXCOMP_ARCH_HOOK)_with_exception

else
CARDIOXCOMP_ARCH = $(OS_NAME)_$(MODE)$(CARDIOXCOMP_ARCH_HOOK)

endif

print_felisce_arch:
	@echo $(CARDIOXCOMP_ARCH)

$(call debug_var,CARDIOXCOMP_ARCH)

$(call debug_flow,Evaluated CARDIOXCOMP_ARCH)

################################## goalopt #################################
# Goal options are activated depending on which goal is invoked on command
# line. (we guarantee that make is invoked with only 0 or 1 goal)
# For example, if goal is 'clean', that is if make is invoked like this:
# make clean
# we don't want the c++ header dependencies to be generated.
# So we will set GOALOPT_HEADERDEPEND_MK to 'no'.
# Goal option variable names starts with 'GOALOPT_', have values 'yes' or
# 'no', are all set here, to be used later.
# For convenience, a list of goals is first defined, for which the goal option
# must be set to 'yes' or must be set to 'no'. It's easy and not error prone to add
# new goal to this list or remove exstining ones. Then, the function yes_if_makegoal_in_list
# or no_if_makegoal_in_list is used to set the goal option according to this
# list and to which goal make is invoked with.

# Whether we will run felisce application
GOALOPT_RUN_YES_FOR_GOALS = run
GOALOPT_RUN = $(call yes_if_makegoal_in_list,$(GOALOPT_RUN_YES_FOR_GOALS))#no trailing space

# Whether C++ header dependencies need to be generated
GOALOPT_HEADERDEPEND_MK_NO_FOR_GOALS = checkconfig clean clean_libfiles clean_exefiles testdiff_clean testdiff_pyclean testdiff_run testdiff_check pytestdiff_run pytestdiff_check lsbuild pywrap pywrap_clean
GOALOPT_HEADERDEPEND_MK = $(call no_if_makegoal_in_list,$(GOALOPT_HEADERDEPEND_MK_NO_FOR_GOALS))#no trailing space

# Whether Swig dependencies need to be generated
GOALOPT_PY_HEADERDEPEND_MK_YES_FOR_GOALS = pywrap
GOALOPT_PY_HEADERDEPEND_MK = $(call yes_if_makegoal_in_list,$(GOALOPT_PY_HEADERDEPEND_MK_YES_FOR_GOALS))#no trailing space

$(call debug_var,GOALOPT_RUN)
$(call debug_var,GOALOPT_HEADERDEPEND_MK)
$(call debug_var,GOALOPT_PY_HEADERDEPEND_MK)
$(call debug_flow,set goal options)

############################## build directory and executable path #############################
# WARNING:
#
# Makefile is first invoked, for example, from:
# - trunk/Source/Model/Makefile
# - this Makefile 'include ../../../config.mk'
# - then we move in directory trunk/build_dir/Model and reinvoke Makefile from here.
# - 'include ../../../conifg.mk' must work, relatively to 'trunk/build_dir/Model'
# - it means that directories 'trunk/Source/Model' and 'trunk/build_dir/Model'
#   must have the same depth relatively to trunk/ (depth 2).
#
# By constrast, trunk/Tests/Examples/XXXX have depth 3.

EXEDIR = .

ifeq ($(strip $(BUILDDIR_LOCATION)),here)
BUILDDIR = $(CURDIR)

else ifeq ($(strip $(BUILDDIR_LOCATION)),src_builddir)
BUILDDIR = $(CARDIOXCOMP_ROOT)/build_$(CARDIOXCOMP_ARCH)/$(notdir $(CURDIR))

else ifeq ($(strip $(BUILDDIR_LOCATION)),pywrap_builddir)
BUILDDIR = $(CARDIOXCOMP_ROOT)/build_$(CARDIOXCOMP_ARCH)/$(notdir $(CURDIR)_pywrap)

else ifeq ($(strip $(BUILDDIR_LOCATION)),example_builddir)

#we change to builddir if building lib or main.
#if executing main, we stay in place. this require $(EXEPATH) to exists:
#we can not build it since we do not move in the build directory.
ifeq ($(GOALOPT_RUN),yes)
EXEDIR = $(CARDIOXCOMP_ROOT)/build_$(CARDIOXCOMP_ARCH)/Examples/$(notdir $(CURDIR))
BUILDDIR = $(CURDIR)
else ifeq ($(GOALOPT_RUN),no)
BUILDDIR = $(CARDIOXCOMP_ROOT)/build_$(CARDIOXCOMP_ARCH)/Examples/$(notdir $(CURDIR))
else
$(call error_var,GOALOPT_RUN)
endif

else
$(call error_var,BUILDDIR_LOCATION)

endif

EXEPATH = $(EXEDIR)/$(EXE)

$(call debug_var,BUILDDIR)
$(call debug_var,EXEDIR)
$(call debug_var,EXEPATH)

lsbuild:
	@echo $(BUILDDIR):
	@ls $(BUILDDIR) || true

$(call debug_flow,Evaluated BUILDDIR)

############################## moving to build directory #############################
ifeq (,$(filter %$(BUILDDIR),$(CURDIR)))
#this is to cd to the build directory and re run make from there
#see: http://mad-scientist.net/make/multi-arch.html

$(call debug_flow,Moving to BUILDDIR and reinvoking make from there)

.SUFFIXES:

MAKETARGET = $(MAKE) --no-print-directory -C $@ -f $(CURDIR)/Makefile SRCDIR=$(CURDIR) BUILDDIR=$(BUILDDIR) $(MAKECMDGOALS) 

.PHONY: $(BUILDDIR)
$(BUILDDIR):
	+@[ -d $@ ] || mkdir -p $@ 
	+@$(MAKETARGET) 

Makefile: ;
%.mk :: ; 

% :: $(BUILDDIR) ;

#we are in builddir
else

$(call debug_flow,We are in the build directory)

ifeq ($(GOALOPT_RUN),yes)

vpath main $(SRCDIR)

else

vpath %.c   $(SRCDIR)
vpath %.h   $(SRCDIR)
vpath %.hpp $(SRCDIR)
vpath %.cpp $(SRCDIR)
vpath %.hxx $(SRCDIR)
vpath %.cxx $(SRCDIR)
vpath %.i   $(SRCDIR)

endif

################################## libraires #################################
LIBRARIES = \
	cardioxcomp \
	$(if $(findstring $(CARDIOXCOMP_WITH_SLEPC),yes),$(SLEPC_LIB),) \
	$(if $(findstring $(CARDIOXCOMP_WITH_SUNDIALS),yes),$(SUNDIALS_LIB),) \
  $(if $(findstring $(CARDIOXCOMP_WITH_LIBXFM),yes),$(LIBXFM_LIB),) \
	$(PETSC_LIB) \
	$(SCALAPACK_LIB) \
	$(BLACS_LIB) \
	$(BLAS_LIB) \
	$(LAPACK_LIB) \
	$(SUPERLU_LIB) \
	$(MUMPS_LIB) \
	$(PARAMETIS_LIB) \
	$(METIS_LIB) \
	$(LM5_LIB) \
	$(LUA_LIB) \
	$(MPI_LIB) \
	$(X11_LIB) \
	$(BOOST_LIB) \
	$(GSL_LIB) \
	$(SOBOL_LIB) \
	$(ADDITIONAL_LIBS)

############################# library directories  ###########################
LIBRARY_DIRS = \
	$(CARDIOXCOMP_ROOT)/build_$(CARDIOXCOMP_ARCH)/LibCardioXcomp \
	$(BLACS_LIB_DIR) \
	$(BLAS_LIB_DIR) \
	$(LAPACK_LIB_DIR) \
	$(LM5_LIB_DIR) \
	$(LUA_LIB_DIR) \
	$(MPI_LIB_DIR) \
	$(MUMPS_LIB_DIR) \
	$(METIS_LIB_DIR) \
	$(PARMETIS_LIB_DIR) \
	$(PETSC_LIB_DIR) \
	$(if $(findstring $(CARDIOXCOMP_WITH_SLEPC),yes),$(SLEPC_LIB_DIR),) \
	$(if $(findstring $(CARDIOXCOMP_WITH_SUNDIALS),yes),$(SUNDIALS_LIB_DIR),) \
	$(if $(findstring $(CARDIOXCOMP_WITH_LIBXFM),yes),$(LIBXFM_LIB_DIR),) \
	$(SCALAPACK_LIB_DIR) \
	$(SUPERLU_LIB_DIR) \
	$(X11_LIB_DIR) \
	$(BOOST_LIB_DIR) \
	$(GSL_LIB_DIR) \
	$(SOBOL_LIB_DIR) \
	$(ADDITIONAL_LIBS_DIR)

# Lists dynamic librairies that the executable is linked to.
CARDIOXCOMP_DYNLIB_FILENAMES = \

CARDIOXCOMP_PY_DYNLIB_FILENAMES = \

print_dynlib:
	@echo $(CARDIOXCOMP_DYNLIB_FILENAMES)

print_py_dynlib:
	@echo $(CARDIOXCOMP_PY_DYNLIB_FILENAMES)

############################# includes directories ###########################
INCLUDE_DIRS = \
	. \
	$(SRCDIR) \
	$(BOOST_INC_DIR) \
	$(LM5_INC_DIR) \
	$(LUA_INC_DIR) \
	$(MPI_INC_DIR) \
	$(METIS_INC_DIR) \
	$(PARMETIS_INC_DIR) \
	$(PETSC_INC_DIR) \
	$(if $(findstring $(CARDIOXCOMP_WITH_SLEPC),yes),$(SLEPC_INC_DIR),) \
	$(if $(findstring $(CARDIOXCOMP_WITH_SUNDIALS),yes),$(SUNDIALS_INC_DIR),) \
	$(if $(findstring $(CARDIOXCOMP_WITH_LIBXFM),yes),$(LIBXFM_INC_DIR),) \
	$(GSL_INC_DIR) \
	$(SOBOL_INC_DIR) \
	$(ADDITIONAL_INC_DIR)

INCLUDE_PYWRAP_DIRS = \
  $(PYTHON_INC_DIR) \
	$(CARDIOXCOMP_ROOT)/Python/Addition

$(call debug_flow,Evaluated compiler flags)
############################## header dependencies #############################

# When compiling an example, headerdepend.mk depends not only on $(SRC), the
# example source code, but on all felisce source code (so that for example,
# if main.cpp include felisce.hpp and felisce.hpp changes, main.cpp is
# rebuild). It would be hard to tell make to rebuild headerdepend.mk if
# one of the felisce source file changes, it simple to rebuild it every
# time.
.PHONY: headerdepend.mk

HEADERDEPEND_MK = headerdepend.mk

ifeq ($(GOALOPT_HEADERDEPEND_MK),yes)

$(HEADERDEPEND_MK): $(SRC)
	@echo "# DO NOT DELETE -- make depend depends on it." > $@
	@echo >> $@
	$(CXXDEPEND) $(DEPENDFLAGS) $(INCLUDES) $^ >> $@

-include $(HEADERDEPEND_MK)

$(call debug_flow,Included headerdepend.mk)
endif

############################### python wrapping ##############################
PYWRAP_I   =  $(PYWRAP_MODULE).i
PYWRAP_CXX =  $(PYWRAP_MODULE)_wrap.cxx
PYWRAP_H   =  $(PYWRAP_MODULE)_wrap.h
PYWRAP_PY  =  $(PYWRAP_MODULE).py
PYWRAP_SO  = _$(PYWRAP_MODULE).so
PYWRAP_O   =  $(PYWRAP_MODULE)_wrap.o

LDFLAGS_WITH_EXCEPTION=$(subst $(BUILDBASE)/Core,$(BUILDBASE)_with_exception/Core,$(LDFLAGS))

.PRECIOUS: $(PYWRAP_O) $(PYWRAP_CXX) $(PYWRAP_H) $(PYWRAP_I)

############################## swig header dependencies #############################

.PHONY: py_headerdepend_mk

PY_HEADERDEPEND_MK = py_headerdepend.mk

ifeq ($(GOALOPT_PY_HEADERDEPEND_MK),yes)

$(PY_HEADERDEPEND_MK): $(PYWRAP_I)
	@echo "# DO NOT DELETE -- make depend depends on it." > $@
	@echo >> $@
	$(SWIG) $(SWIG_OPT) -MM $< > $@

ifeq ($(MAKECMDGOALS),pywrap)
-include $(PY_HEADERDEPEND_MK)
endif

endif

################################## some check ################################
# function signature: $(call assert-no-whitespace,DIRECTORY)
define assert-no-whitespace
  $(if $(word 2,[$($1)]),$(error There is a whitespace inside variable '$(1)', please correct it),)
endef

# If for example $(CARDIOXCOMP_ROOT) contains trailing whitepace(s),
# $(CARDIOXCOMP_ROOT)/Source would be two words, 
# resulting is two wrong directories.
CHECK_FOR_WHITESPACE = \
    CARDIOXCOMP_ROOT \
    PETSC_DIR \
    PETSC_ARCH \
    EXT_PACK \
    X11_DIR
$(foreach dir,$(CHECK_FOR_WHITESPACE),$(call assert-no-whitespace,$(dir)))


$(call debug_flow,Performed some check)
############################### exceptions ##############################
ifeq "$(WITH_EXCEPTION)" "yes"
CPPFLAGS_ALL += -DCARDIOXCOMP_WITH_EXCEPTION
endif


pyfelisce_addition: $(PYCARDIOXCOMP_ADDITION) install_py_felisce_addition

clean_pyfelisce_addition:
	rm -f felisce_additionmodule.o  felisce_addition.so  headerdepend.mk

install_py_felisce_addition: $(PYCARDIOXCOMP_ADDITION)
	+@[ -d $(CARDIOXCOMP_ROOT)/Python/lib_$(MODE)/wrap/ ] || mkdir -p $(CARDIOXCOMP_ROOT)/Python/lib_$(MODE)/wrap/
	cp $(PYCARDIOXCOMP_ADDITION) $(CARDIOXCOMP_ROOT)/Python/lib_$(MODE)/wrap/
$(call debug_flow,Set exception flags)
ifneq "$(LOGGED)" "yes"

############################### usual build rules ##############################

OBJ = $(SRC:.cpp=.o)

#----------------------------------- compile -----------------------------------
%.o: %.cpp
	$(CXX) $(CXXFLAGS_ALL) $(CPPFLAGS_ALL) $(INCLUDES) -c $<

#------------------------------- create library --------------------------------
create_lib: $(LIB)

$(LIB): $(OBJ)
	$(call CREATELIB,$@,$^)

#------------------------------ create executable ------------------------------
create_exe: $(EXEPATH)

ifeq ($(GOALOPT_RUN),yes)
$(EXEPATH):
	$(error '$(EXEPATH)' does not exsits or is not up to date. Please build/rebuild it)

else ifeq ($(GOALOPT_RUN),no)
$(EXEPATH): $(OBJ) $(CARDIOXCOMP_DYNLIB_FILENAMES)
	$(CXX) $(OBJ) $(LDFLAGS) -o $@

else
$(call error_var,GOALOPT_RUN)
endif

#------------------------------ python wrapping ------------------------------
pywrap: $(PY_HEADERDEPEND_MK) $(PYWRAP_SO) $(PYWRAP_PY) pywrap_install

felisce_%_wrap.cxx felisce_%_wrap.h felisce_%.py: felisce_%.i
	$(SWIG) $(SWIG_OPT) $<

felisce_%_wrap.o: felisce_%_wrap.cxx felisce_%_wrap.h
	$(CXX) $(CXXFLAGS_ALL) $(CPPFLAGS_ALL) -DCARDIOXCOMP_WITH_EXCEPTION $(CXXFLAGS_PYWRAP) $(INCLUDES) $(INCLUDES_PYWRAP) -c $< -o $@

_felisce_%.so: felisce_%_wrap.o
	$(CXX) $(LDFLAGS_PYWRAP) $< $(LDFLAGS_WITH_EXCEPTION) -o $@

pywrap_clean:
	@rm -f *.o $(PY_HEADERDEPEND_MK) $(PYWRAP_CXX) $(PYWRAP_H) $(PYWRAP_PY) *.so *.pyc headerdepend.mk *.log *.failure_log

pywrap_install: $(PYWRAP_PY) $(PYWRAP_SO)
	+@[ -d $(CARDIOXCOMP_ROOT)/Python/lib_$(MODE)/wrap/ ] || mkdir -p $(CARDIOXCOMP_ROOT)/Python/lib_$(MODE)/wrap/
	cp $(PYWRAP_PY) $(PYWRAP_SO) $(CARDIOXCOMP_ROOT)/Python/lib_$(MODE)/wrap/

else # ifneq "$(LOGGED)" "yes"

######################### build rules for logged mode ##########################

# When LOGGED is true, if an object is up to date, print log of last compilation
# instead of silently pass.  If a  object is not up to date,  compilation output
#  is outputed  as usually.  Usefull  for  keeping  traces  of  warnings without
# recompiling.

#----------------------------------  output  -----------------------------------
# If  the file target.failure_log exits,  it  means that the  compiling has just
# been performed,  so  the  output  has  already  been  sent to standard output.
#  (target  can   be   file.o   or   file.so   or   file.dylib  or  main...)  If
# target.failure_log  does not exists,  it means  that target  is already  up to
#  date,  so  we  sent  to  standard  output  the  last  output  compiling  log.
#
# Note: %.output depend on %, so for example main.output forces main to be
# created.
#                                                                               
# Warning:  an %.output target is never  up-to-date (there is no %.output file),
# so be carefull of targets that depend on %.output:  the rule they depend on is
#  always executed.  That's  why %.output  target  should  appear  only  on main
#  targets (e.g.  create_lib,  create_exe,  testdiff)  and  not  on intermediate
# targets.
# Is there a solution to execute a %.output without depend on it? It should be
# worth to search in GNU/Make manual.
#
.PRECIOUS: %.o %.(LIBEXT) %.log %.failure_log $(EXEPATH)

%.output: %
	@if [ -f $*.failure_log ]; \
  then \
    mv $*.failure_log $*.log; \
  elif [ -f $*.log ]; \
  then	\
		sed 's/^/logged: /' $*.log; \
	else \
	  echo -n "make:WARNING: '$*' is up-to-date, but build log can't be found. "; \
	  echo "You may want to force regeneration of build log by cleaning."; \
  fi

#--------------------------------- loggedcmd -----------------------------------
define loggedcmd
@echo "$(1)" 2>&1 | tee $@.failure_log
@set -o pipefail ; $(1) 2>&1 | tee -a $@.failure_log
endef

#----------------------------------- compile -----------------------------------
OBJ = $(SRC:.cpp=.o)
OBJ_OUTPUT = $(foreach obj,$(OBJ),$(obj).output)

%.o: %.cpp
	$(call loggedcmd,$(CXX) $(CXXFLAGS_ALL) $(CPPFLAGS_ALL) $(INCLUDES) -c $<)

#------------------------------- create library --------------------------------
LIB_OUTPUT = $(LIB).output

create_lib: $(OBJ_OUTPUT) $(LIB_OUTPUT)

$(LIB): $(OBJ)
	$(call loggedcmd,$(call CREATELIB,$@,$^))

#------------------------------ create executable ------------------------------
EXEPATH_OUTPUT = $(EXEPATH).output

create_exe: $(OBJ_OUTPUT) $(EXEPATH_OUTPUT)

ifeq ($(GOALOPT_RUN),yes)
$(EXEPATH):
	$(error '$(EXEPATH)' does not exsits or is not up to date. Please build/rebuild it)

else ifeq ($(GOALOPT_RUN),no)
$(EXEPATH): $(OBJ) $(CARDIOXCOMP_DYNLIB_FILENAMES)
	$(call loggedcmd,$(CXX) $(OBJ) $(LDFLAGS) -o $@)

else
$(call error_var,GOALOPT_RUN)
endif

$(call debug_flow,Defined rules for logged MODE)

#------------------------------ python wrapping ------------------------------
PYWRAP_CXX_OUTPUT = $(PYWRAP_CXX).output
PYWRAP_O_OUTPUT = $(PYWRAP_O).output
PYWRAP_SO_OUTPUT = $(PYWRAP_SO).output

pywrap: $(PY_HEADERDEPEND_MK) $(PYWRAP_CXX_OUTPUT) $(PYWRAP_O_OUTPUT) $(PYWRAP_SO_OUTPUT) pywrap_install

felisce_%_wrap.cxx: felisce_%.i
	$(call loggedcmd,$(SWIG) $(SWIG_OPT) $<)

felisce_%_wrap.o: felisce_%_wrap.cxx 
	$(call loggedcmd,$(CXX) $(CXXFLAGS_ALL) $(CPPFLAGS_ALL) -DCARDIOXCOMP_WITH_EXCEPTION $(CXXFLAGS_PYWRAP) $(INCLUDES) $(INCLUDES_PYWRAP) -c $< -o $@,)

_felisce_%.so: felisce_%_wrap.o
	$(call loggedcmd,$(CXX) $(LDFLAGS_PYWRAP) $< $(LDFLAGS_WITH_EXCEPTION) -o $@,)

pywrap_clean:
	@rm -f *.o $(PY_HEADERDEPEND_MK) $(PYWRAP_CXX) $(PYWRAP_H) $(PYWRAP_PY) *.so *.pyc headerdepend.mk *.log *.failure_log

pywrap_install: $(PYWRAP_SO)
	+@[ -d $(CARDIOXCOMP_ROOT)/Python/lib_$(MODE)/wrap/ ] || mkdir -p $(CARDIOXCOMP_ROOT)/Python/lib_$(MODE)/wrap/
	cp $(PYWRAP_PY) $(PYWRAP_SO) $(CARDIOXCOMP_ROOT)/Python/lib_$(MODE)/wrap/

endif # ifneq "$(LOGGED)" "yes"

##################################### clean ####################################
.PHONY: clean cleanfiles

LIBFILES = *.o *.$(LIBEXT) *.failure_log *.log $(HEADERDEPEND_MK)

clean_libfiles:
	@find $(BUILDDIR) \( $(foreach pat,$(LIBFILES),-name "$(pat)" -o ) -name "*.o" \) -exec rm {} \;

clean_exefiles:
	@rm -f $(HEADERDEPEND_MK) $(OBJ) $(EXEPATH) *.failure_log *.log

$(call debug_flow,Defined clean rules)

#################################### diff ######################################

# Simply diff of files, to be always run on the same machine.
#
# Usage:
# ------
#
# In Makefile, declare  FILES_TO_DIFF, for example:
#
#     FILES_TO_DIFF = SolutionDirectory/file1.ext SolutionDirectory/file2.ext
#
#	After that (only after), include rules.mk:
#    include ../../../rules.mk
#
# Use the included 'testdiff' target:
#
#     SHELL$ make testdiff
#
# WARNING: Do not delete Solution files without  cleanin executable (see
# note-A).
#
# You may need/want to override default variables:
#
#     SHELL$ make testdiff CARDIOXCOMP_RESULTDIR_REF=/path/to/solution/ref/
#     SHELL$ make testdiff DIFF_VERBOSITY=200
#
# CARDIOXCOMP_RESULTDIR_REF is a directory where the reference solution files reside.
# It must be must organised like this:
#     CARDIOXCOMP_RESULTDIR_REF/
#     |-- TestCase0
#     |   `-- SolutionDirectory
#     |       |-- file1
#     |       `-- file2
#     `-- TestCase0
#         `-- SolutionDirectory
#             |-- file1
#             `-- file2
#

DIFF_VERBOSITY = 100 # Number of diff lines to output

# Notes about implementation:
# ---------------------------
#
# note-A: FILES_TO_DIFF and run_testdiff dependence.
#
# The FILES_TO_DIFF are created by running $(EXE).  However, It is not guaranted
# a  priori that  running $(EXE)  will really  create FILES_TO_DIFF.  Therefore,
# it's a  bad idea to make FILES_TO_DIFF  depends on run_testdiff.  For example,
# make will re-run the  test case for each missing file.  Instead,  run_testdiff
# is executed  one  time  at  the  first  beginning,  and then FILES_TO_DIFF are
# assumed to  exists (a proper error  is raised if not).  The  risk is that,  if
#  FILES_TO_DIFF are  inadvertensly  deleted,  and  dependences  of run_testdiff
# didn't changed,  the test will failed while it should succeed.  However, there
# is no robust solution to detect that run_testdiff must be reexecuted.
#
# note-B: Use of '%.output: %' rule
#
#  This rule  is (by  definition) always  NOT up-to-date  (there is  no %.output
# file).  It results that  a rule depending on a '%.output:  %'  rule is as well
# always not  up-to-date,  so the recipe is always  executed.  That's the reason
# why,  for example,  create_exe  is  a  dependence  of  testdiff  rather than a
# dependence of run_testdiff.

# note-C: Default target VS testdiff/pytestdiff
#
# We  must force user  to explicity call  'make testdiff' or  'make pytestdiff',
# because  some variable are evaluated only  if command line
#  goal is  testdiff or  pytestdiff.  A rule  like this  for example:  'default:
# pytestdiff' will miss to evaluated some variables.

#------------------------------------- run -------------------------------------
define clean_result_dir
    @if [ -d $1 ]; \
		then \
				if [[ "$1" =~ ^"felisceresult" ]]; \
				then \
						rm -rf $1; \
						echo "MAKE:INFO: clean_result_dir: deleted directory $1"; \
				else \
				    echo "MAKE:ERROR: clean_result_dir: For safety, not delete directory '$1' because is does not starts with 'felisceresult'"; \
						exit 1; \
				fi; \
    fi
endef

# * run  is special:  is it not a  file to create but a  command to execute.  To
# avoid re-run the  command each time,  on need to remember date  of last run to
# compare it  with  dependencies.  We  can't  know  surely  a  file that will be
# created, so let's create of file called run_testdiff.

.PHONY: run pyrun

WITH_RUN_LOCK = no

$(call debug_var,WITH_RUN_LOCK)

#=====================================
ifneq "$(WITH_RUN_LOCK)" "yes"

# note: when with_run_lock is not yes, we reexecute the run everytime.
# so logged mode as no sense.

run: $(EXEPATH)
	$(MPI_EXEC) -n $(NPROC) $(EXEPATH) -f $(DATAFILE) $(CARDIOXCOMP_OPT)

pyrun: $(PYMAIN)
	$(MPI_EXEC) -n $(NPROC) $(PYTHON) $(PYMAIN) -f $(DATAFILE) $(CARDIOXCOMP_OPT)

#=====================================
else #with_run_lock

# note: when with_run_lock is yes, we reexecute the run only is a dependence of 
# the results of the run has changed.

# note: when we reexecute the run, we delete the solution directory, so it
# requires to the variable RESULT_DIR to be defined, and it will override the
# value defined in the date file.

#.....................................
ifneq "$(LOGGED)" "yes"

run: lock_$(RESULT_DIR)

pyrun: pylock_$(RESULT_DIR)

lock_$(RESULT_DIR): Makefile $(DATAFILE) $(EXEPATH) $(CARDIOXCOMP_DYNLIB_FILENAMES) $(ADDITIONNAL_RUN_DEPENDENCIES)
	$(call clean_result_dir,$(RESULT_DIR))
	$(MPI_EXEC) -n $(NPROC) $(EXEPATH) -f $(DATAFILE) $(CARDIOXCOMP_OPT) --felisce-mesh-resultDir=$(RESULT_DIR)/
	@echo "Used by make to store date of last run" > $@

pylock_$(RESULT_DIR): Makefile $(DATAFILE) $(PYSRC) $(CARDIOXCOMP_PY_DYNLIB_FILENAMES) $(ADDITIONNAL_RUN_DEPENDENCIES)
	@echo "Used by make to store date of last run" > $@
	$(MPI_EXEC) -n $(NPROC) $(PYTHON) $(PYMAIN) -f $(DATAFILE) $(CARDIOXCOMP_OPT) --felisce-mesh-resultDir=$(RESULT_DIR)/

#.....................................
else #logged

run: lock_$(RESULT_DIR).output

pyrun: pylock_$(RESULT_DIR).output

lock_$(RESULT_DIR): Makefile $(DATAFILE) $(EXEPATH) $(CARDIOXCOMP_DYNLIB_FILENAMES) $(ADDITIONNAL_RUN_DEPENDENCIES)
	$(call clean_result_dir,$(RESULT_DIR))
	$(call loggedcmd,$(MPI_EXEC) -n $(NPROC) $(EXEPATH) -f $(DATAFILE) $(CARDIOXCOMP_OPT) --felisce-mesh-resultDir=$(RESULT_DIR)/)
	@echo "Used by make to store date of last run" > $@

pylock_$(RESULT_DIR): Makefile $(DATAFILE) $(PYSRC) $(CARDIOXCOMP_PY_DYNLIB_FILENAMES) $(ADDITIONNAL_RUN_DEPENDENCIES)
	@echo "Used by make to store date of last run" > $@
	$(call loggedcmd,$(MPI_EXEC) -n $(NPROC) $(PYTHON) $(PYMAIN) -f $(DATAFILE) $(CARDIOXCOMP_OPT) --felisce-mesh-resultDir=$(RESULT_DIR)/)

#.....................................
endif # logged

#=====================================
endif # with_run_lock

$(call debug_flow,Defined run rules)

#------------------------------------ diff -------------------------------------
DIFF_FILES = $(foreach file,$(FILES_TO_DIFF),$(STORE_DIFF_DIR)/$(file).diff)
DIFF_CHECK = $(foreach file,$(FILES_TO_DIFF),$(STORE_DIFF_DIR)/$(file).check)

#check variables required for testdiff are defined
check_diff_files_required_variables:
	  $(call assert-variable-defined,FILES_TO_DIFF)
	  $(call assert-variable-defined,DIFF_DIR_1)
	  $(call assert-variable-defined,DIFF_DIR_2)
	  $(call assert-variable-defined,STORE_DIFF_DIR)
	  $(call assert-variable-defined,DIFF_VERBOSITY)

diff_files: check_diff_files_required_variables $(DIFF_CHECK)

# * Make diff file beetween current results and reference results.
# * This is a static rule,  so if a file is missing the rule is applied however,
#    so the missing file is correctly reported.
$(DIFF_FILES): $(STORE_DIFF_DIR)/%.diff: $(DIFF_DIR_1)/% $(DIFF_DIR_2)/%
	@diff $^ > $@ || :

$(call debug_flow,Defined diff files rules)

# * Check that number of lines in diff file is 0.
%.check: %.diff
	@ echo "$< ($(DIFF_VERBOSITY) first lines)"; \
   head -n $(DIFF_VERBOSITY) $< ; \
   test $$(wc -l $< | awk '{print $$1}') -eq 0

#----------------------------------- help -------------------------------------
testdiff_pyhelp:
	@echo "usage:"
	@echo "python main.py      # run test case"
	@echo "make pytestdiff_run     # run and diff results"
	@echo "make clean          # clean everything"

$(call debug_flow,Defined testdiff help rules)

#----------------------------------- checkconfig -------------------------------------
# It will check that directory BLACS_LIB_DIR contains libraries listed in
# BLACS_LIB, and so on.
PROJECT_WITH_LIB_TO_CHECK = \
    BLACS \
    BLAS \
    BOOST \
    LAPACK \
    LM5 \
    LUA \
    METIS \
    MUMPS \
    MPI \
    PARMETIS \
    PETSC \
    SCALAPACK \
    SUPERLU \
    X11

# It will check that petsc.h and petscconf.h are each one contained by
# one of the directory listed by PETSC_DIR, and so on.
PROJECT_WITH_HEADER_TO_CHECK = \
    BOOST:boost/numeric/ublas/vector.hpp \
    GETPOT:GetPot \
    LM5:libmesh5.h \
    LUA:lua.h \
    MPI:mpi.h \
    MUMPS:cmumps_c.h \
    PARMETIS:parmetis.h \
    PETSC:petsc.h:petscconf.h

ifeq "$(CARDIOXCOMP_WITH_SLEPC)" "yes" 
PROJECT_WITH_LIB_TO_CHECK += SLEPC
PROJECT_WITH_HEADER_TO_CHECK += SLEPC:slepc.h:slepcconf.h
endif

ifeq "$(CARDIOXCOMP_WITH_SUNDIALS)" "yes" 	 
PROJECT_WITH_LIB_TO_CHECK += SUNDIALS 	 
#PROJECT_WITH_HEADER_TO_CHECK += SUNDIALS: 	 
endif

ifeq "$(CARDIOXCOMP_WITH_LIBXFM)" "yes" 	 
PROJECT_WITH_LIB_TO_CHECK += LIBXFM
PROJECT_WITH_HEADER_TO_CHECK += LIBXFM:xfm_library.h
endif

ifeq '$(WITH_PYTHON_USE_WRAPPING)' 'yes'

PROJECT_WITH_HEADER_TO_CHECK += \
  PYTHON:Python.h

EXECUTABLE_PATH_TO_CHECK += \
    $(PYTHON) \
    $(CXX) \
    $(CXXDEPEND) \
    $(MPI_EXEC)

PYTHON_SYS_PATH = $(shell python -c "import sys; print ' '.join(sys.path)")

CHECK_PYTHON_SYS_PATH_CARDIOXCOMP = felisce.py
CHECK_PYTHON_SYS_PATH_CARDIOXCOMP_HELP = "in felisce SVN repository in trunk/Python/lib"

CHECK_PYTHON_SYS_PATH_PETSC4PY = petsc4py/PETSc.py
CHECK_PYTHON_SYS_PATH_PETSC4PY_HELP = "at http://code.google.com/p/petsc4py/"

PYTHON_FILES_IN_SYS_PATH_TO_CHECK = \
    CARDIOXCOMP \
    PETSC4PY

#Test existance of a Python file in sys.path
# $1 is for example: felisce.py
#                or: petsc4py/PETSc.py
# $2 is a help message on how to install the missing file.
define python-sys-path-contains-file
    $(shell \
		for dir in $(PYTHON_SYS_PATH); \
		do \
		    [ -f $$dir/$1 ] && exit 0 || true ; \
	    done; \
		echo "WARNING: Python sys.path does not contained file '$1' (can be found $2).";\
	)
endef

endif

ifeq '$(WITH_PYTHON_GENERATE_WRAPPING)' 'yes'

EXECUTABLE_PATH_TO_CHECK += \
  $(SWIG)

endif

CHECKCONFIG_COUNTER = 0

#Test existance of a library
# $1 is for example: MPI_LIB_DIR
# $2 is for example: MPI_LIB
# $3 if for example: mpi_f77
define dir-contains-lib
    $(shell \
				for dir in $($1); \
				do \
				    [ -f $$dir/lib$3.a -o -f $$dir/lib$3.so -o -f $$dir/lib$3.dylib ] && exit 0 || true ; \
				done;
				echo "WARNING: Directory(ies) '$($1)' ('$1') does not contain library '$3' (listed in '$2')."; \
		 )
endef

#Test existance of a header file
# $1 is for example: PETSC_INC_DIR
# $2 is for example: petsc.h
define dir-contains-header
    $(shell \
				[ -z "$($1)" ] && exit 0 || true; \
				for dir in $($1); \
				do \
						[ -f $$dir/$2 ] && exit 0 || true; \
					done; \
				echo "WARNING: Directory(ies) '$($1)' ('$1') does not contained required header '$2'.";\
	)
endef

#Test existance of a executable in the PATH
define executable-path-exists
    $(shell \
		[ -z "$1" ] && exit 0 || true; \
		which $1 &> /dev/null && exit 0 || true; \
		echo "WARNING: Executable path '$1' does not exists";\
	)
endef

#Print Warning if parameter is a non-empty string
#Usage:
# $(call warns-if-not-empty-string $foo)
define warns-if-not-empty-string
    CHECKCONFIG_COUNTER=$(shell echo $$(($(CHECKCONFIG_COUNTER)+1))$$ )
    $(if $(strip $1),$(warning $1),)
endef

check_project_with_lib:
	$(foreach PROJECT_WITH_LIB, $(PROJECT_WITH_LIB_TO_CHECK), \
	    $(foreach LIBRARY, $($(PROJECT_WITH_LIB)_LIB), \
		    $(eval \
	            $(call \
                     warns-if-not-empty-string, \
                     $(call dir-contains-lib,$(PROJECT_WITH_LIB)_LIB_DIR,$(PROJECT_WITH_LIB)_LIB,$(LIBRARY)) \
                 ) \
		    ) \
      ) \
    )

check_project_with_header:
	$(foreach PROJECT_WITH_HEADER, $(PROJECT_WITH_HEADER_TO_CHECK), \
		$(foreach HEADER, $(wordlist 2, $(words $(subst :, ,$(PROJECT_WITH_HEADER))), $(subst :, ,$(PROJECT_WITH_HEADER))), \
		$(eval \
			$(call \
				warns-if-not-empty-string, \
				$(call dir-contains-header,$(firstword $(subst :, ,$(PROJECT_WITH_HEADER)))_INC_DIR,$(HEADER)) \
			) \
			) \
		) \
	)

check_python_files_in_sys_path:
	$(foreach PYTHON_FILE_IN_SYS_PATH, $(PYTHON_FILES_IN_SYS_PATH_TO_CHECK), \
			$(eval \
			  $(call \
				warns-if-not-empty-string, \
				$(call python-sys-path-contains-file,$(CHECK_PYTHON_SYS_PATH_$(PYTHON_FILE_IN_SYS_PATH)),$(CHECK_PYTHON_SYS_PATH_$(PYTHON_FILE_IN_SYS_PATH)_HELP)) \
			) \
	) \
	)

check_executables:
	$(foreach EXECUTALBE_PATH, $(EXECUTABLE_PATH_TO_CHECK), \
	  $(eval \
		  $(call \
	  warns-if-not-empty-string, \
				$(call executable-path-exists,$(EXECUTALBE_PATH)) \
		  ) \
		) \
  )

checkconfig: check_project_with_lib check_project_with_header check_executables check_python_files_in_sys_path
	@echo "checkconfig completed ($(CHECKCONFIG_COUNTER) tests executed)"

$(call debug_flow,Defined checkconfig rules)

#we are in builddir
endif

#----------------------------------- safeshell function -------------------------------------

# The safeshell function does the same thing as $(shell ...),  but cause make to
# stop if an error has occured and print informative report of the error,
# what $(shell ...) does not.
# 
# example:
# --------
#     CMD = cat
#     ARGS = file1.txt file2.txt file3.txt
#     CONTENT = $(call safeshell,$(CMD),$(ARGS))
# or in one line:
#     CONTENT = $(call safeshell,cat,file1.txt file2.txt file3.txt)
# Be carefull to the two commas separating the make function name,
# the shell command, and the blank separeted shell arguments list.
#
# errors management:
# ------------------
# An error is raised '$(CMD)' is not found (i.e. 'which $(CMD)' return nothing).
# An error is raised if '$(CMD) $(ARGS)' fails (i.e. $? is not zero).
#
# when to use it, when to not use it:
# -----------------------------------
#
# $(CMD) $(ARGS) is called multiple times, so:
#    - safeshell is suitable for commands that run extremelly fast,
#    - $(CMD) $(ARGS) returned value must be the same for the 1st call, the 2nd
#    call, the 3rd call etc. (typically, the result must not depend on a
#    incremented variable...).
#    - it is usefull for example to run a command that extract a value
#    from a file.
#
# Note about implementation:
# --------------------------
# The line is pretty unreadable,  but blanks (whitespaces,  newlines, tabs) does
# matter in Makefile syntax.
# In peusdo-code, it reads like this ($1 is the command, $2 are the arguments):
#     cmdpath = shell( which CMD )
#     if cmdpath == '':
#         raise the error: 'No such shell command: CMD'
#     else:
#         errcode = shell( CMD ARGS > /dev/null 2>&1; echo $? )
#         if errcode == 0:
#             result = shell( CMD ARGS )
#         else:
#             raise the error: 'Shell command failed: CMD ARSG'
#     return result
safeshell = $(if $(strip $(shell which $1)),$(if $(filter 0,$(shell $1 $2 > /dev/null 2>&1; echo $$?)),$(shell $1 $2),$(shell $1 $2 1>&2) $(error Shell command failed: '$1 $2')),$(error No such shell command: '$1'))
