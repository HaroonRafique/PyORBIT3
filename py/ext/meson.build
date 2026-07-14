python = import('python').find_installation('python3', pure: false)

subdir('ptc_orbit')

python.install_sources(
    files('__init__.py'),
    subdir : 'ext',
)
