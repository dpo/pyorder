Installation
------------

For now, just::

    python setup.py build
    python setup.py install

To select another Fortran compiler::

    python setup.py config_fc --fcompiler=<Fortran compiler name> build

To see a list of available Fortran compilers and their names::

    python setup.py config_fc --help-fcompiler
