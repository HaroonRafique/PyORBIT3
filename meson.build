project('PyORBIT', ['cpp', 'fortran'], 
        version: run_command('python', '-m', 'setuptools_scm', check: true).stdout().strip(),
        meson_version: '>=1.4.1',
        default_options: ['buildtype=release'])


# Add libptc as a subproject
libptc_proj = subproject('libptc')

libptc_dep = dependency('libptc', fallback:['libptc', 'libptc_dep'])


libptc_dep = libptc_proj.get_variable('libptc_dep')

# Include subdirectories
subdir('src')
subdir('py/orbit')
