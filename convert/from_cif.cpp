
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/to_pdb.hpp>
#include <gemmi/cif2mtz.hpp>
#include <emscripten/emscripten.h>

extern std::string global_str;
extern std::string global_str2;

extern "C" {

const char* EMSCRIPTEN_KEEPALIVE cif2pdb(char* data, size_t size) {
  try {
    gemmi::cif::Document doc = gemmi::cif::read_memory(data, size, "input");
    std::free(data);
    gemmi::Structure st = gemmi::make_structure(doc);
    if (!st.models.empty() && !st.models[0].chains.empty()) {
      std::ostringstream os;
      gemmi::write_pdb(st, os);
      global_str = os.str();
    } else {
      global_str = "ERROR: probably not an mmCIF file.";
    }
  } catch (std::exception& e) {
    global_str = "ERROR: ";
    global_str += e.what();
  }
  return global_str.c_str();
}

const char* EMSCRIPTEN_KEEPALIVE cif2mtz(char* data, size_t size) {
  global_str.clear();
  global_str2.clear();
  try {
    auto rblocks = gemmi::as_refln_blocks(
        gemmi::cif::read_memory(data, size, "input").blocks);
    std::free(data);
    if (!rblocks.empty()) {
      std::ostringstream out;
      gemmi::CifToMtz c2m;
      c2m.verbose = true;
      gemmi::Mtz mtz = c2m.convert_block_to_mtz(rblocks[0], out);
      global_str2 = out.str();
      mtz.write_to_string(global_str);
    } else {
      global_str = "ERROR: probably not a SF-mmCIF file.";
    }
  } catch (std::exception& e) {
    global_str = "ERROR: ";
    global_str += e.what();
  }
  return global_str.c_str();
}

} // extern "C"
