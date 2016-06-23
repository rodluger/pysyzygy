.. image:: https://coveralls.io/repos/github/rodluger/pysyzygy/badge.svg?branch=master :target: https://coveralls.io/github/rodluger/pysyzygy?branch=master
.. image:: https://travis-ci.org/rodluger/pysyzygy.svg?branch=master :target: https://travis-ci.org/rodluger/pysyzygy
    
py·sy·zy·gy
-----------
/ˈpīsizijē/ |speaker|

.. |speaker| image:: img/speaker.png?raw=True
             :target: http://www.astro.washington.edu/users/rodluger/pysyzygy.mp3

*noun*

**1.** A fast and general planet transit (`syzygy <http://en.wikipedia.org/wiki/Syzygy_%28astronomy%29>`_) code written in C and in Python.

**2.** ``pysyzygy`` computes **fast** lightcurves for the most general case of a *massive*, *eccentric* planet orbiting a limb-darkened star. Here's a sample output image of an assymetric transit:

.. image:: img/transit.png?raw=True
   :alt: pysyzygy
   :align: center

Installation
============
Clone the repository and run

>>> python setup.py install


Calling pysyzygy...
===================

... is super easy.

.. code-block:: python
  
  import numpy as np
  
  # Orbital elements
  kwargs = {'rhos': 1.0,          # Stellar density in g/cm^3
            'MpMs': 0.001,        # Planet-star mass ratio
            'esw': 0.1,           # Eccentricity vectors
            'ecw': 0.1, 
            'per': 1.0,           # Period in days
            'RpRs': 0.1,          # Planet-star radius ratio
            't0': 0.0,            # Time of first transit in days
            'q1': 1.0,            # Kipping (2013) quadratic limb darkening coefficients
            'q2': 0.0,
            'b': 0.5}             # Circular impact parameter
  
  # Instantiate a transit object
  trn = ps.Transit(**kwargs) 
  
  # Now evaluate the lightcurve on a grid of observation times
  t = np.arange(0., 10., ps.transit.KEPLONGCAD)
  flux = trn(t)
        

Notes
=====

Feel free to change, adapt, or incorporate this code into your project, but please make sure to cite this repository, as well as `Mandel and Agol (2002) <http://adsabs.harvard.edu/abs/2002ApJ...580L.171M>`_, the transit model on which ``pysyzygy`` is based.
