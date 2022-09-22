// Convert either pdb file (in particular BUSTER output)
// or MTZ file (merged or unmerged) to mmCIF.

#include <cstdlib>
#include <sstream>
#include <gemmi/pdb.hpp>
#include <gemmi/cif.hpp>       // for cif::read_memory
#include <gemmi/to_cif.hpp>    // for write_cif_to_stream
#include <gemmi/to_mmcif.hpp>  // for update_cif_block
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
      gemmi::setup_entities(st);
      gemmi::assign_label_seq_id(st, false);
      std::ostringstream os;
      gemmi::cif::write_cif_to_stream(os, gemmi::make_mmcif_document(st),
                                      gemmi::cif::Style::PreferPairs);
      global_str = os.str();
    } else {
      global_str = "ERROR: probably it is not a PDB file.";
    }
  } catch (std::exception& e) {
    global_str = "ERROR: ";
    global_str += e.what();
  }
  return global_str.c_str();
}

const char* EMSCRIPTEN_KEEPALIVE mtz2cif(char* data, size_t size) {
  try {
    std::ostringstream os;
    gemmi::MtzToCif mtz_to_cif;
    if (size > 17 && strncmp(data, "!FORMAT=XDS_ASCII", 17) == 0) {
      gemmi::XdsAscii xds_ascii;
      xds_ascii.read_stream(gemmi::MemoryStream(data, size), "<input>");
      std::free(data);
      mtz_to_cif.write_cif_from_xds(xds_ascii, os);
    } else {
      gemmi::Mtz mtz;
      mtz.read_stream(gemmi::MemoryStream(data, size), true);
      std::free(data);
      mtz.switch_to_original_hkl();
      gemmi::SMat33<double> b;
      mtz_to_cif.staraniso_version = gemmi::read_staraniso_b_from_mtz(mtz, b);
      gemmi::SMat33<double>* staraniso_b = b.all_zero() ? nullptr : &b;
      mtz_to_cif.write_cif(mtz, nullptr, staraniso_b, os);
    }
    global_str = os.str();
  } catch (std::exception& e) {
    global_str = "ERROR: ";
    global_str += e.what();
  }
  return global_str.c_str();
}

const char* EMSCRIPTEN_KEEPALIVE mxdepo(char* data1, size_t size1,
                                        char* data2, size_t size2) {
  try {
    std::unique_ptr<gemmi::Mtz> mtz1, mtz2;
    std::unique_ptr<gemmi::XdsAscii> xds_ascii;

    if (data1[0] == 'M') {
      mtz1.reset(new gemmi::Mtz);
      mtz1->read_stream(gemmi::MemoryStream(data1, size1), true);
      std::free(data1);
    }

    if (data2[0] == 'M') {
      mtz2.reset(new gemmi::Mtz);
      mtz2->read_stream(gemmi::MemoryStream(data2, size2), true);
      if (mtz1 && !mtz1->is_merged() && mtz2->is_merged())
        mtz1.swap(mtz2);
    } else if (data2[0] == '!') {
      xds_ascii.reset(new gemmi::XdsAscii);
      xds_ascii->read_stream(gemmi::MemoryStream(data2, size2), "<input>");
      xds_ascii->gather_iset_statistics();
    } else {
      gemmi::fail("the second file is neither MTZ nor XDS_ASCII");
    }
    std::free(data2);

    if (mtz1 && !mtz1->is_merged())
      gemmi::fail("the first file is not merged");
    if (mtz2 && mtz2->is_merged())
      gemmi::fail("the second file should not be merged");

    bool ok = true;
    global_str2.clear();
    if (mtz1) {
      std::ostringstream out;
      gemmi::remove_appendix_from_column_names(*mtz1, out);
      ok = gemmi::validate_merged_mtz_deposition_columns(*mtz1, out);
      global_str2 = out.str();
    }
    gemmi::MtzToCif mtz_to_cif;
    gemmi::Intensities mi;
    try {
      std::ostringstream validate_out;
      if (mtz1) {
        mtz_to_cif.staraniso_version = mi.take_staraniso_b_from_mtz(*mtz1);
        if (!mtz_to_cif.staraniso_version.empty()) {
          validate_out << "Merged MTZ went through STARANISO "
                       << mtz_to_cif.staraniso_version;
          if (mi.staraniso_b.ok())
            validate_out << ". Taking into account anisotropic scaling.\n";
          else
            validate_out << ". B tensor is unknown. Intensities won't be checked.\n";
        }
        mi.read_merged_intensities_from_mtz(*mtz1);
      } else {
        gemmi::ReflnBlock rblock = gemmi::get_refln_block(
            gemmi::cif::read_memory(data1, size1, "<input>").blocks, {});
        if (!rblock.entry_id.empty() && rblock.entry_id != mtz_to_cif.entry_id) {
          global_str2 += "Note: using _entry.id from mmCIF: ";
          global_str2 += rblock.entry_id;
          global_str2 += ".\n";
          mtz_to_cif.entry_id = rblock.entry_id;
        }
        mi.take_staraniso_b_from_mmcif(rblock.block);
        mi.read_merged_intensities_from_mmcif(rblock);
      }
      gemmi::Intensities ui;
      if (mtz2)
        ui.read_unmerged_intensities_from_mtz(*mtz2);
      else if (xds_ascii)
        ui.read_unmerged_intensities_from_xds(*xds_ascii);

       // If an old StarAniso version was used that doesn't store B tensor,
       // allow intensities to differ.
       bool relaxed_check = !mtz_to_cif.staraniso_version.empty() && !mi.staraniso_b.ok();

       if (!gemmi::validate_merged_intensities(mi, ui, relaxed_check, validate_out))
        ok = false;
      global_str2 += validate_out.str();
    } catch (std::runtime_error& e) {
      global_str2 += "Intensity merging not validated: ";
      global_str2 += e.what();
      ok = false;
    }
    if (!ok)
      return "Error. If you think the files are correct, contact us.";

    if (mtz2)
      mtz2->switch_to_original_hkl();

    std::ostringstream os;
    mtz_to_cif.write_special_marker_for_pdb = true;
    mtz_to_cif.with_history = false;
    mtz_to_cif.less_anomalous = 1;
    mtz_to_cif.gemmi_run_from = "mxdepo.html";
    gemmi::SMat33<double>* staraniso_b = mi.staraniso_b.ok() ? &mi.staraniso_b.b : nullptr;
    if (mtz1)
      mtz_to_cif.write_cif(*mtz1, nullptr, staraniso_b, os);
    else
      os.write(data1, size1);
    os << "\n\n";
    mtz_to_cif.block_name = "unmerged";
    if (mtz2)
      mtz_to_cif.write_cif(*mtz2, nullptr, nullptr, os);
    else if (xds_ascii)
      mtz_to_cif.write_cif_from_xds(*xds_ascii, os);
    global_str = os.str();
  } catch (std::exception& e) {
    global_str = "ERROR: ";
    global_str += e.what();
  }
  return global_str.c_str();
}

} // extern "C"
