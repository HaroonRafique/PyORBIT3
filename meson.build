#project('PyORBIT', ['cpp', 'fortran'], version: run_command('python', '-m',  'setuptools_scm').stdout().strip())
project('PyORBIT', ['cpp', 'fortran'], 
        version: run_command('python', '-m', 'setuptools_scm', check: true).stdout().strip(),
        meson_version: '>=1.4.1',
        default_options: ['buildtype=release'])
# Add Fortran compiler
#fc = meson.get_compiler('fortran')

# Include subdirectories
subdir('ext/PTC')  # Add the PTC library subdirectory
subdir('src')
subdir('py/orbit')
