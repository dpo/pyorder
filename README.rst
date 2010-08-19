PyOrder --- A Library of Sparse Matrix Reordering Methods in Python

Requirements
------------

1. NumPy
2. The Harwell Subroutine Library subroutines MC21 and MC60. More interfaces
   may be added as needed. Head to
   http://hsl.rl.ac.uk/hsl2007/hsl20074researchers.html for license
   terms and for obtaining the software. You need to download the double
   precision real version. The source files must be (re)named mc21ad.f and
   mc60ad.f. Edit site.cfg to specify the location of the source files.

Installation
------------

For now, just::

    python setup.py build
    python setup.py install

To select another Fortran compiler::

    python setup.py config_fc --fcompiler=<Fortran compiler name> build

For instance::

    python setup.py config_fc --fcompiler=gfortran build

assuming `gfortran` is in your `PATH`.

To see a list of available Fortran compilers and their names::

    python setup.py config_fc --help-fcompiler
