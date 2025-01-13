#project('PyORBIT', ['cpp', 'fortran'], version: run_command('python', '-m',  'setuptools_scm').stdout().strip())
project('PyORBIT', ['cpp', 'fortran'], 
        version: run_command('python', '-m', 'setuptools_scm', check: true).stdout().strip(),
        meson_version: '>=1.4.1',
        default_options: ['buildtype=release'])
# Add Fortran compiler
fc = meson.get_compiler('fortran')

subdir('src')
subdir('py/orbit')

# Add PTC directory
subdir('ext/PTC')

#~ pyorbit_sources = files(
#~     *glob(pyorbit_src_dir, '*.cc', '*.cpp')  # Adjust extensions as needed
#~ )

#~ pyorbit_lib = shared_library(
#~     'pyorbit',
#~     pyorbit_sources,
#~     link_with: [ptc_lib],  # Link with the built PTC library
#~     include_directories: include_directories('ext/PTC/source'),
#~     install: true  # Install pyorbit library via pip
#~ )

