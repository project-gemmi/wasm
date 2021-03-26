// Convert either pdb file (in particular BUSTER output)
// or MTZ file (merged or unmerged) to mmCIF.

#include <cstdlib>
#include <sstream>
#include <gemmi/pdb.hpp>
#include <gemmi/to_cif.hpp>    // for write_cif_to_stream
#include <gemmi/to_mmcif.hpp>  // for update_cif_block
#include <gemmi/remarks.hpp>   // for read_metadata_from_remarks
#include <gemmi/mtz.hpp>       // for Mtz
#include <gemmi/mtz2cif.hpp>   // for MtzToCif
#include <gemmi/align.hpp>     // for assign_label_seq_id
#include <emscripten/emscripten.h>

extern std::string global_str;
extern std::string global_str2;

extern "C" {

const char* EMSCRIPTEN_KEEPALIVE pdb2cif(char* data, size_t size) {
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

const char* EMSCRIPTEN_KEEPALIVE mtz2cif(char* data, size_t size) {
  try {
    gemmi::Mtz mtz;
    mtz.read_stream(gemmi::MemoryStream(data, size), true);
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

const char* EMSCRIPTEN_KEEPALIVE mtzpair2cif(char* data1, size_t size1,
                                             char* data2, size_t size2) {
  try {
    gemmi::Mtz mtz1;
    mtz1.read_stream(gemmi::MemoryStream(data1, size1), true);
    std::free(data1);

    gemmi::Mtz mtz2;
    mtz2.read_stream(gemmi::MemoryStream(data2, size2), true);
    std::free(data2);

    if (mtz1.is_merged() && mtz2.is_merged())
      gemmi::fail("two merged MTZ files given");
    if (!mtz1.is_merged() && !mtz2.is_merged())
      gemmi::fail("two unmerged MTZ files given");
    gemmi::Mtz& merged = mtz1.is_merged() ? mtz1 : mtz2;
    gemmi::Mtz& unmerged = mtz1.is_merged() ? mtz2 : mtz1;

    std::ostringstream validate_out;
    gemmi::Intensities mi = gemmi::read_mean_intensities_from_mtz(merged);
    validate_merged_intensities(unmerged, mi, validate_out);
    global_str2 = validate_out.str();

    unmerged.switch_to_original_hkl();

    std::ostringstream os;
    gemmi::MtzToCif mtz_to_cif;
    mtz_to_cif.write_cif(merged, &unmerged, os);
    global_str = os.str();
  } catch (std::runtime_error& e) {
    global_str = "ERROR: ";
    global_str += e.what();
  }
  return global_str.c_str();
}

} // extern "C"
