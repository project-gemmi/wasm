// Convert either pdb file (in particular BUSTER output)
// or MTZ file (merged or unmerged) to mmCIF.

#include <cstdlib>
#include <sstream>
#define GEMMI_WRITE_IMPLEMENTATION
#include <gemmi/pdb.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/to_mmcif.hpp>  // for update_cif_block
#include <gemmi/remarks.hpp>   // for read_metadata_from_remarks
#include <gemmi/version.hpp>   // for GEMMI_VERSION
#include <gemmi/mtz.hpp>       // for Mtz
#include <gemmi/mtz2cif.hpp>   // for MtzToCif
#include <gemmi/align.hpp>     // for assign_label_seq_id

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
      gemmi::setup_entities(st);
      gemmi::assign_label_seq_id(st, false);
      std::ostringstream os;
      gemmi::cif::write_cif_to_stream(os, gemmi::make_mmcif_document(st),
                                      gemmi::cif::Style::PreferPairs);
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

const char* mtz2cif(char* data, size_t size) {
  try {
    gemmi::Mtz mtz;
    mtz.read_stream(gemmi::MemoryStream(data, data + size), true);
    std::free(data);
    mtz.switch_to_original_hkl();
    std::ostringstream os;
    gemmi::MtzToCif mtz_to_cif;
    mtz_to_cif.write_cif(mtz, nullptr, os);
    global_str = os.str();
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
