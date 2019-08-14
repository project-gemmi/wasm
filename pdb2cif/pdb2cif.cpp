// Convert pdb file from BUSTER to mmCIF file for deposition.

#include <cstdlib>
#include <sstream>
#define GEMMI_WRITE_IMPLEMENTATION
#include <gemmi/pdb.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/to_mmcif.hpp>  // for update_cif_block
#include <gemmi/remarks.hpp>   // for read_metadata_from_remarks
#include <gemmi/version.hpp>   // for GEMMI_VERSION

std::string global_str;

extern "C" {

const char* get_version() {
  return GEMMI_VERSION;
}

const char* pdb2cif(char* data, size_t size) {
  try {
    gemmi::Structure st = gemmi::read_pdb_from_memory(data, size, "input.pdb");
    std::free(data);
    if (!st.models.empty() && !st.models[0].chains.empty()) {
      gemmi::read_metadata_from_remarks(st);
      setup_entities(st);
      gemmi::cif::Document doc;
      doc.blocks.resize(1);
      update_cif_block(st, doc.blocks[0]);
      std::ostringstream os;
      write_cif_to_stream(os, doc, gemmi::cif::Style::PreferPairs);
      global_str = os.str();
    } else {
      global_str = "ERROR: probably it is not a PDB file.";
    }
  } catch (std::runtime_error& e) {
    global_str = "ERROR: ";
    global_str += e.what();
  }
  return global_str.c_str();
}

void clear_string() {
  global_str.clear();
}

} // extern "C"
