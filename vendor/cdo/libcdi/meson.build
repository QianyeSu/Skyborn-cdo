# libcdi - Climate Data Interface library
# Simplified Meson build

# This is a simplified version focusing on core functionality
# Full CDI library has many more features

cdi_sources = files(
  # Core CDI sources - listing key files
  # Note: Full list would be extracted from libcdi/src/Makefile.am
  'src/basetime.c',
  'src/cdi.c',
  'src/cdi_att.c',
  'src/cdi_cksum.c',
  'src/cdi_error.c',
  'src/cdi_int.c',
  'src/cdi_key.c',
  'src/cdi_limits.c',
  'src/cdi_util.c',
  'src/cdi_uuid.c',
  'src/dmemory.c',
  'src/file.c',
  'src/resource_handle.c',
  'src/resource_unpack.c',
  'src/serialize.c',
  'src/stream_cdf.c',
  'src/stream_cgribex.c',
  'src/stream_echam.c',
  'src/stream_ext.c',
  'src/stream_fcommon.c',
  'src/stream_ieg.c',
  'src/stream_record.c',
  'src/stream_scan.c',
  'src/stream_srv.c',
  'src/stream_var.c',
  'src/stream_grb.c',
  'src/stream_gribapi.c',
  'src/stream_history.c',
  'src/subtype.c',
  'src/vlist.c',
  'src/vlist_var.c',
  'src/zaxis.c',
  'src/table.c',
  'src/taxis.c',
  'src/timebase.c',
  'src/varscan.c',
  # Add more source files as needed
)

# CDI include directory
cdi_inc = include_directories('src', 'interfaces')

# Build libcdi as static library
libcdi = static_library('cdi',
  cdi_sources,
  include_directories: cdi_inc,
  dependencies: [netcdf_dep, thread_dep],
  install: false,
)

# CDI dependency for use in main CDO
libcdi_dep = declare_dependency(
  link_with: libcdi,
  include_directories: cdi_inc,
)
