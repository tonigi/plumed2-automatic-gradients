
/* Hack: including a C file to wrap in the appropriate namespace */

#include "Wrapper.h"
#include <cmath>

namespace PLMD {
namespace curvature {

#include "sympy_codegen/curvature_codegen.c"

}
}
