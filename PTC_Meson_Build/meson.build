# ext/PTC/meson.build

# Include the source directory relative to the current directory
inc_dirs = include_directories('source')

# Define the Fortran source files with relative paths
ptc_sources = files(
  'source/m_real_polymorph.f90',
  'source/n_complex_polymorph.f90',
  'source/o_tree_element.f90',
  'source/c_dabnew.f90',
  'source/d_lielib.f90',
  'source/h_definition.f90',
  'source/i_tpsa.f90',
  'source/j_tpsalie.f90',
  'source/k_tpsalie_analysis.f90',
  'source/l_complex_taylor.f90',
  'source/a_scratch_size.f90',
  'source/b_da_arrays_all.f90',
  'source/Sq_orbit_ptc.f90',
  'source/Sra_fitting.f90',
  'source/Sr_spin.f90',
  'source/Ss_fake_mad.f90',
  'source/St_pointers.f90',
  'source/Sk_link_list.f90',
  'source/Sl_family.f90',
  'source/Sma0_beam_beam_ptc.f90',
  'source/Sma_multiparticle.f90',
  'source/Sm_tracking.f90',
  'source/Sn_mad_like.f90',
  'source/So_fitting.f90',
  'source/Spb_fake_gino_sub.f90',
  'source/Sp_keywords.f90',
  'source/Se_status.f90',
  'source/Sf_def_all_kinds.f90',
  'source/Sg_sagan_wiggler.f90',
  'source/Sh_def_kind.f90',
  'source/Si_def_element.f90',
  'source/Sa_extend_poly.f90',
  'source/Sb_sagan_pol_arbitrary.f90',
  'source/Sc_euclidean.f90',
  'source/Sd_frame.f90'
)

# Build the PTC static library
ptc_lib = static_library(
  'ptc',
  ptc_sources,
  include_directories: inc_dirs,
  install: false
)
