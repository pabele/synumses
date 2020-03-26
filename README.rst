SyNumSeS
========

**SyNumSeS** is a Python https://www.python.org package for simulating various semiconductor devices
in one dimension.
It solves the *Van Roosbroeck equations* using the *Scharfetter Gummel algorithm* [CIT01]_.



It can simulate:

- pn-diodes
- LEDs
- MOS-diodes
- bipolar transistors
- solar cells
- hetero junctions


About SyNumSeS
--------------

**SyNumSeS** was developed to get a better
understanding on how semiconductor devices work and
visualize their parameters like:

- current density,
- charge concentration,
- quasi-Fermi-levels ...

**SyNumSeS** stands for
Symbolicall, Numerically
Semiconductor Simulation.

It uses **SymPy** https://www.sympy.org/en/index.html to derive
the Jacobian matrix from the the discretized *Van Roosbroeck equations*
using the *Scharfetter Gummel algorithm*
to solve for the vector containing the
potential and the quasi-Fermi-levels by the Newton method.
This way the *Van Roosbroeck equations* can be simply extended by,
for example, recombination terms. 
The linear solver uses
**SciPy** https://www.scipy.org further
**NumPy** https://numpy.org/ is used for fast operations.


Documentation
-------------

For documentation of the modules used in this package
you can clone this repository::

  git clone https://github.com/pabele/synumses.git

change to the *doc* directory::

  cd synumses/doc

and build the documentation (not needed)::

  make html

Open the file *synumses/docs/build/html/index.html* in a browser.



Installation
------------

It is recommended to use **Anaconda  Python 3.7 version** https://www.anaconda.com/distribution,
use this link for download,
and install **JupyterLab** https://github.com/jupyterlab/jupyterlab with::
  
  conda install -c conda-forge jupyterlab


Installation of **SyNumSeS** with python 3.7::

  python3 -m pip install synumses-pkg-pabele


Or clone this repository if you want to do changes::

  git clone https://github.com/pabele/synumses.git

and::

  pip install -e synumses


Examples
--------

The best way to get started is to have a look at the
**JupyterLab** examples https://github.com/pabele/synumses/tree/master/examples/jupyter.
or in the folder *synumses/examples/jupyter/*.

Some Python scripts can be found in *synumses/examples/*

Known Issues
------------

- Quasi-Fermi-levels are not constant in thermodynamic equilibrium.
- The simulation of wide band gap material shows convergence problems.

.. [CIT01] D. Scharfetter, H. Gummel, Large-signal analyses of a
	   silicon read diode oscillator,
	   IEEE Trans. Electron Devices 16(1) (1969) 64-77
