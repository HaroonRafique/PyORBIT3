project('PyORBIT', ['cpp'], version: run_command('python', '-m',  'setuptools_scm').stdout().strip())

subdir('src')
subdir('py/orbit')

if get_option('enable_ptc')
    subdir('py/ext')
endif
