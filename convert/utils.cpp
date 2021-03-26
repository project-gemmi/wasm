#include <string>
#include <gemmi/version.hpp>   // for GEMMI_VERSION
#include <emscripten/emscripten.h>

std::string global_str;
std::string global_str2;

extern "C" {

const char* EMSCRIPTEN_KEEPALIVE get_version() {
  return GEMMI_VERSION;
}

const char* EMSCRIPTEN_KEEPALIVE get_str2() {
  return global_str2.c_str();
}

void EMSCRIPTEN_KEEPALIVE clear_string() {
  global_str.clear();
  global_str2.clear();
}

} // extern "C"
