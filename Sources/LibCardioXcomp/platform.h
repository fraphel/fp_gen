#pragma once
#include <stdint.h>
#include <petscsys.h>

namespace cardioxcomp{

#ifdef _WIN32
#ifdef _WIN64
	typedef int64_t platformInt;
#else
	typedef int platformInt;
#endif
#else //ubuntu, mac
	typedef PetscInt platformInt;
#endif

}
