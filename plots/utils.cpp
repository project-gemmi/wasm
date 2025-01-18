//#include <gemmi/version.hpp>   // for GEMMI_VERSION
#include "describe.h"
#include <emscripten/emscripten.h>

extern "C" {

const char* EMSCRIPTEN_KEEPALIVE get_version() {
  return GEMMI_DESCRIBE " (compiled on " __DATE__ ")";
}

} // extern "C"
