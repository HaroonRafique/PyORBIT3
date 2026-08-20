# PyORBIT3 PTC branch installation

This README documents the `PTC` branch of Haroon Rafique's PyORBIT3 fork. It
covers the integrated PyORBIT3 + PTC build path and the PyORBIT3-only fallback.

The PTC branch adds:

- optional Meson PTC integration through `-Denable_ptc=true`;
- the compiled `pylibptc_orbit` extension;
- the Python wrapper package `ext.ptc_orbit`;
- selected CERN-era compatibility updates for 1D longitudinal space charge,
  tune analysis, analytical Gaussian space charge, and foil hit/loss counters.

The PTC source is not vendored in this repository. It is consumed as a Meson
subproject from a sibling checkout of the PTC repository.

## Repositories

Clone both repositories into the same workspace:

```bash
git clone --branch PTC git@github.com:HaroonRafique/PyORBIT3.git
git clone git@github.com:HaroonRafique/PTC.git
```

Expected layout:

```text
workspace/
  PyORBIT3/
  PTC/
```

Link the PTC checkout into PyORBIT3 as the Meson subproject:

```bash
cd workspace/PyORBIT3
mkdir -p subprojects
ln -s ../../PTC subprojects/libptc_orbit
```

If `subprojects/libptc_orbit` already exists and is not the intended PTC
checkout, fix that before building. The PTC-enabled build expects that path to
provide the `libptc_orbit` Meson subproject variables used by
`src/meson.build`.

## System requirements

Install the normal PyORBIT3 build tools plus the PTC requirements:

- `python3`
- `git`
- `pkg-config`
- `gcc`
- `g++`
- `gfortran`
- FFTW development files
- MPICH or OpenMPI, preferably with `mpicc` and `mpirun`

Ubuntu example:

```bash
sudo apt-get update
sudo apt-get install -y build-essential gfortran python3 python3-venv \
    libpython3-dev libfftw3-dev pkg-config git
```

Red Hat/Fedora example:

```bash
sudo dnf group install -y "Development Tools"
sudo dnf install -y gcc-c++ gcc-gfortran python3-devel fftw-devel pkgconf-pkg-config git
```

macOS example:

```bash
brew install pkg-config fftw gcc
```

MPI is part of the validated PTC-PyORBIT3 path. The build still supports a
PyORBIT3-only/no-MPI fallback, but use MPI for PTC validation unless you are
debugging a specific build problem.

## Python environment

Create and activate a local virtual environment:

```bash
cd workspace/PyORBIT3
python3 -m venv ../venv
source ../venv/bin/activate
python -m pip install -U pip
python -m pip install -r requirements.txt
python -m pip install -U setuptools
python -m pip install numpy scipy matplotlib pytest
```

For notebook workflows, also install:

```bash
python -m pip install notebook ipykernel nbconvert nbformat
python -m ipykernel install --sys-prefix \
    --name pyorbit3-ptc \
    --display-name "PyORBIT3 + PTC (venv)"
```

The `ptc_pyorbit/build_instructions_PTC_PyORBIT3` script uses the same policy:
all Python packages are installed into a local virtual environment, not into the
system Python.

## PTC-enabled build

Build the branch in editable mode with PTC enabled:

```bash
python -m pip install --no-build-isolation --editable . \
    --config-settings=setup-args="-DUSE_MPI=mpich" \
    --config-settings=setup-args="-Denable_ptc=true"
```

Use `-DUSE_MPI=ompi` for OpenMPI. Use `-DUSE_MPI=auto` to let Meson choose the
first detected MPI implementation. The supported values are:

| MPI choice | Meaning |
| --- | --- |
| `mpich` | Require MPICH through `pkg-config` package `mpich`. |
| `ompi` | Require OpenMPI through `pkg-config` package `ompi`. |
| `auto` | Try MPICH, then OpenMPI, then build without MPI if neither is found. |
| `none` | Build without MPI support. |

Meson discovers MPI through `pkg-config`. If your MPI installation provides
`mpicc` but not a `.pc` file, generate a local pkg-config file and prepend its
directory to `PKG_CONFIG_PATH`. The integrated handoff script writes this under
`integrated_build_artifacts/pkgconfig/<mpi>.pc` from `mpicc -show`.

## PyORBIT3-only fallback

To build without PTC, disable the PTC option:

```bash
python -m pip install --no-build-isolation --editable . \
    --config-settings=setup-args="-DUSE_MPI=auto" \
    --config-settings=setup-args="-Denable_ptc=false"
```

With `-Denable_ptc=false`, `pylibptc_orbit` and `ext.ptc_orbit` are not
installed. Importing `ext.ptc_orbit` will fail by design.

## Reproducing the handoff build script

The PTC-PyORBIT3 handoff script is in the sibling handoff directory:

```text
ptc_pyorbit/build_instructions_PTC_PyORBIT3
```

It is intended for a fresh workspace. To use that script directly for this fork
and branch, run it with PyORBIT3 overridden to Haroon Rafique's fork:

```bash
PYORBIT_REPO_URL=git@github.com:HaroonRafique/PyORBIT3.git \
PYORBIT_REF=PTC \
PTC_REPO_URL=git@github.com:HaroonRafique/PTC.git \
./build_instructions_PTC_PyORBIT3
```

The important defaults are:

| Variable | Default | Purpose |
| --- | --- | --- |
| `PYORBIT_REPO_URL` | `git@github.com:PyORBIT-Collaboration/PyORBIT3.git` | PyORBIT3 repo used by the original handoff script. For this branch, use `git@github.com:HaroonRafique/PyORBIT3.git` and `PYORBIT_REF=PTC` or clone this branch directly. |
| `PTC_REPO_URL` | `git@github.com:HaroonRafique/PTC.git` | PTC source repository. |
| `PYORBIT_DIR` | `./PyORBIT3` | PyORBIT3 checkout path. |
| `PTC_DIR` | `./PTC` | PTC checkout path. |
| `PYORBIT_REF` | empty | Optional PyORBIT3 branch, tag, or commit to check out after clone. Use `PTC` for this branch. |
| `PTC_REF` | empty | Optional PTC branch, tag, or commit to check out after clone. |
| `CLONE_PYORBIT` | `1` | Clone PyORBIT3 if `PYORBIT_DIR` is empty. Set to `0` to require an existing checkout. |
| `CLONE_PTC` | `1` | Clone PTC if `PTC_DIR` is empty. Set to `0` to require an existing checkout. |
| `VENV_DIR` | `./venv` | Local virtual environment. |
| `ARTIFACT_DIR` | `./integrated_build_artifacts` | Generated artifacts, including MPI pkg-config fallback files. |
| `LOG_DIR` | `./integrated_build_outputs` | Build logs. |
| `PYTHON_BOOTSTRAP` | `python3` | Python used to create the virtual environment. |
| `MPI_CHOICE` | `auto` | MPI choice passed to `-DUSE_MPI`. |
| `ENABLE_PTC` | `1` | Build PTC by default. |
| `RUN_PTC_LATTICE_TEST` | `1` | Read the PTC flat file after building. |
| `RUN_FULL_TESTS` | `0` | Skip full pytest unless explicitly enabled. |
| `FLAT_FILE` | `./PTC-PyORBIT_flat_file.madx.flt` | First PTC smoke-test lattice. |
| `VERIFY_REQUIREMENTS` | `numpy scipy matplotlib pytest` | Extra verification/runtime packages installed into the venv. |
| `JUPYTER_REQUIREMENTS` | `notebook ipykernel nbconvert nbformat` | Notebook packages installed into the venv. |
| `JUPYTER_KERNEL_NAME` | `pyorbit3-ptc` | Venv-backed Jupyter kernel name. |
| `JUPYTER_KERNEL_DISPLAY_NAME` | `PyORBIT3 + PTC (venv)` | Jupyter display name for the venv kernel. |

The script also installs these verification/notebook packages into the venv:

```text
numpy scipy matplotlib pytest notebook ipykernel nbconvert nbformat
```

For this already-integrated `PTC` branch, the script's patching steps are
already present in the source tree: `enable_ptc`, `py/ext/ptc_orbit`, and
`pylibptc_orbit` wiring in `src/meson.build`.

## PTC smoke test

Use the MAD-X-generated flat file from the handoff bundle for the first PTC
runtime check:

```text
PTC-PyORBIT_flat_file.madx.flt
```

After building, verify the integrated extension and wrapper:

```bash
python -c "import orbit; print('orbit import ok', orbit.__file__)"
python -c "from orbit.core.bunch import Bunch; b = Bunch(); print('Bunch ok', b.getSize())"
python -c "import pylibptc_orbit; print('pylibptc_orbit import ok', pylibptc_orbit.__file__)"
python -c "from ext.ptc_orbit import PTC_Lattice; print('PTC_Lattice import ok')"
```

Then read the PTC flat file:

```bash
python - /path/to/PTC-PyORBIT_flat_file.madx.flt <<'PY'
from pathlib import Path
from ext.ptc_orbit import PTC_Lattice

flat_file = Path(__import__("sys").argv[1]).resolve()
lattice = PTC_Lattice("PTC_SMOKE")
lattice.readPTC(str(flat_file))

print("n_nodes", len(lattice.getNodes()))
print("lattice_length", lattice.getLength())
print("n_harm", lattice.nHarm)
print("gamma_t", lattice.gammaT)

if len(lattice.getNodes()) <= 0:
    raise SystemExit("PTC lattice read produced no nodes")
if lattice.getLength() <= 0.0:
    raise SystemExit("PTC lattice length is not positive")
PY
```

The recorded baseline for the accepted MAD-X flat file is:

```text
n_nodes 404
lattice_length 163.36282
n_harm 2
gamma_t 5.069309599302147
```

Run the MPI smoke test after the basic imports:

```bash
mpirun -n 2 python examples/MPI_Tests/mpi_initialization_test.py
```

If `mpirun` is unavailable, run the script directly as a limited local check:

```bash
python examples/MPI_Tests/mpi_initialization_test.py
```

## Verification ladder

Use this order when validating a fresh PTC-PyORBIT3 build:

1. Build editable PyORBIT3 with PTC enabled.
2. `import orbit`.
3. Import the compiled `orbit.core.*` modules.
4. Instantiate `orbit.core.bunch.Bunch`.
5. Run the minimal MPI-aware PyORBIT3 example.
6. `import pylibptc_orbit`.
7. `from ext.ptc_orbit import PTC_Lattice`.
8. Read `PTC-PyORBIT_flat_file.madx.flt` with `PTC_Lattice.readPTC(...)`.
9. Confirm positive node count and lattice length.
10. Run full PyORBIT3 tests only after the import/lattice-read smoke test
    passes:

```bash
python -m pytest tests
```

The first PTC smoke test proves build/import/flat-file-read readiness. It does
not replace the standalone `ptc_pyorbit3_examples` acceptance suite or later
production physics validation.

## Run a standard PyORBIT3 example

After building, the normal PyORBIT3 examples remain available. For example:

```bash
cd examples/SNS_Linac/pyorbit3_linac_model/
python pyorbit3_sns_linac_mebt_hebt2.py
```
