# nbody-iv-conditioner

A N-body initial values conditioner. Given state vectors and desired first integrals, this lib conditionates these state vectors to obtain these first integrals (as long as it's physically possible).

The theory behind all this thing is in the [documentation](doc/documentation.pdf).

## Python interface

The python interface is available on a [PyPI repository](https://pypi.org/project/nbody-iv-conditioner/) and to install it just run

```bash
pip install nbody-iv-conditioner
```

The interface is almost the same as the Fortran, and it have two big classes: `Utils` and `Conditioners`. Although the modules are available too, the parameters needs to be in Fortran format in many functions (because in Fortran they may be subroutines with `intent(inout)`), and using the classes `Utils` and `Conditioners` all the functions have the automatic conversor.

```python
import nbody_iv_conditioner as nivc
import numpy as np

N = 10
ms = np.ones(N)/N
qs = 2.0 * np.random.random((3,N)) - 1.0

# move the com to origin
ms, qs, ps = nivc.Conditioners.comori(ms, qs)

com = nivc.Utils.center_of_mass(ms, qs)
```

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