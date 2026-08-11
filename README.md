# nbody-iv-conditioner

A N-body initial values conditioner. Given state vectors and desired first integrals, this lib conditionates these state vectors to obtain these first integrals (as long as it's physically possible).

The theory behind all this thing is in the [documentation](doc/documentation.pdf).

## Fortran interface

To compile the Fortran interface, just go the directory `fortran` and run

```bash
make lib
```

It's possible to change the directory where the library goes. The default is `/lib` (at the root of this directory). To change, just use the parameter `INSTALL_DIR`.

To link the compiled library to your script if you installed it on the directory "/abobrinha", just add the flags to your `gfortran` command:

```bash
-I/abobrinha -L/abobrinha -lconditioners_mod -lutils
```

Example:
```bash
gfortran -o my_program my_script.f90 -I/abobrinha -L/abobrinha -lconditioners_mod -lutils
```

To use the subroutines, just import the modules:

```Fortran
use conditioners_mod
use utils_mod
```

### Tests/examples

There are some tests to guarantee that all the things works. To run it:

```bash
make test
./test
```