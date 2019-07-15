
/* Hack: including a C file to wrap in the appropriate namespace */

#include "CodeGenWrapper.h"
#include <cmath>

namespace PLMD {
namespace curvature_codegen {

#include "sympy_codegen/curvature_codegen.c"

}
}
