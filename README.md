# PyORBIT3 PTC branch installation

This README documents the `PTC` branch of Haroon Rafique's PyORBIT3 fork. It
covers the integrated PyORBIT3 + PTC build path and the PyORBIT3-only fallback.

The PTC branch adds:

- optional Meson PTC integration through `-Denable_ptc=true`;
- the compiled `pylibptc_orbit` extension;
- the Python wrapper package `ext.ptc_orbit`;
- selected CERN-era compatibility updates for 1D longitudinal space charge,
  tune analysis, analytical Gaussian space charge, and foil hit/loss counters.

The PTC source is not vendored in this repository. The PTC checkout can live
beside PyORBIT3, but Meson must be able to find it through the subproject path
`subprojects/libptc_orbit`.

## Repositories

Clone PyORBIT3 and PTC side by side:

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
cd PyORBIT3
mkdir -p subprojects
ln -s ../../PTC subprojects/libptc_orbit
```

This side-by-side layout is the recommended developer setup when editing both
repositories. It keeps each repository's Git state separate while still giving
Meson the `subprojects/libptc_orbit` path required by
`subproject('libptc_orbit')`.

For a one-off build where you do not plan to edit PTC, you can instead clone
PTC directly into `subprojects/libptc_orbit`:

```bash
cd PyORBIT3
mkdir -p subprojects
git clone git@github.com:HaroonRafique/PTC.git subprojects/libptc_orbit
```

That avoids the symlink, but it creates a nested Git checkout under PyORBIT3,
which can make branch/status handling less clear when developing both projects.

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
cd PyORBIT3
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

Install Python packages into a virtual environment, not into the system Python.

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
directory to `PKG_CONFIG_PATH`.

## PyORBIT3-only fallback

To build without PTC, disable the PTC option:

```bash
python -m pip install --no-build-isolation --editable . \
    --config-settings=setup-args="-DUSE_MPI=auto" \
    --config-settings=setup-args="-Denable_ptc=false"
```

With `-Denable_ptc=false`, `pylibptc_orbit` and `ext.ptc_orbit` are not
installed. Importing `ext.ptc_orbit` will fail by design.

## PTC smoke test

The repository does not currently include a validated PTC flat file. To test
`PTC_Lattice.readPTC(...)`, provide a MAD-X-generated PTC flat file from your
own validation inputs.

After building, verify the integrated extension and wrapper:

```bash
python -c "import orbit; print('orbit import ok', orbit.__file__)"
python -c "from orbit.core.bunch import Bunch; b = Bunch(); print('Bunch ok', b.getSize())"
python -c "import pylibptc_orbit; print('pylibptc_orbit import ok', pylibptc_orbit.__file__)"
python -c "from ext.ptc_orbit import PTC_Lattice; print('PTC_Lattice import ok')"
```

Then read the PTC flat file:

```bash
PTC_FLAT_FILE=valid_ptc_flat_file.flt
python - "$PTC_FLAT_FILE" <<'PY'
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
8. Read a validated PTC flat file with `PTC_Lattice.readPTC(...)`.
9. Confirm positive node count and lattice length.
10. Run full PyORBIT3 tests only after the import/lattice-read smoke test
    passes:

```bash
python -m pytest tests
```

The first PTC smoke test proves build/import/flat-file-read readiness. It does
not replace later production physics validation.

## Run a standard PyORBIT3 example

After building, the normal PyORBIT3 examples remain available. For example:

```bash
cd examples/SNS_Linac/pyorbit3_linac_model/
python pyorbit3_sns_linac_mebt_hebt2.py
```
