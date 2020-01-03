SyNumSeS
========

SyNumSeS is a `Python <https://www.python.org>` package for simulating various semiconductor devices
in one dimension.
It solves the *Van Roosbroeck equations*.

It can simulate:
- pn-diodes
- LEDs
- MOS-diodes
- bipolar transistors
- solar cells
- hetero junctions


About SyNumSeS
--------------

SyNumSeS was developed to get a better
understanding on how semiconductor devices work and
visualize their parameters like
current density, charge concentration, quasi-Fermi-levels ...

**SyNumSeS** stands for
**Sy**mbolically, **Num**erically
**Se**miconductor **S**imulation.
It uses `SymPy <https://www.sympy.org/en/index.html>` to derive
Jacobian matrix from the the discretized *Van Roosbroeck equations*
to solve for the vector containing the
potential and the quasi-Fermi-level by the Newton method.
The linear solver uses
`SciPy <https://www.scipy.org>` further
`NumPy <https://numpy.org/>` is used for fast operations.


Documentation
-------------

The best way to get started is to have a look at the
`jupyter lab examples <https://github.com/pabele/synumses/tree/master/examples/jupyter>`.

Or for documentation of the modules used in this package
you can clone this repository
```sh
git clone https://github.com/pabele/synumses.git
```
change to the **doc** directory
```sh
cd synumses/doc
```
and build the documentation
```sh
make html
```
Open **synumses/docs/build/html/index.html** with a browser.
Installation
------------
Installation with python 3.7:
```sh
python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps synumses-pkg-pabele
```
Or clone this repository.
```sh
git clone https://github.com/pabele/synumses.git
```
