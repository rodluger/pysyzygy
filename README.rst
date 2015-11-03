py·sy·zy·gy
-----------
/ˈpīsizijē/ |speaker|

.. |speaker| image:: /../img/speaker.png?raw=True
             :target: http://www.astro.washington.edu/users/rodluger/pysyzygy.mp3

*noun*

**1.** A fast and general planet transit (`syzygy <http://en.wikipedia.org/wiki/Syzygy_%28astronomy%29>`_) code written in C and in Python.

**2.** ``pysyzygy`` computes **fast** lightcurves for the most general case of a *massive*, *eccentric* planet orbiting a limb-darkened star. Here's a sample output image of an assymetric transit:

.. image:: /../img/transit.png?raw=True
   :alt: pysyzygy
   :scale: 50 %
   :align: center

**3.** ``pysyzygy`` can also generate sweet animations, such as this one (*WARNING: not to scale!*):

.. image:: /../img/transit.gif?raw=True
   :alt: pysyzygy
   :align: center

Installation
============
Clone the repository and run

>>> make

to build the transit module. Add the top-level folder to your ``$PATH``, then call or 
import ``pysyzygy`` to begin plotting transits!

Calling pysyzygy...
===================

... is super easy. If you want to plot stuff, try the following examples:

.. code-block:: python
  
  import pysyzygy as ps
  ps.PlotTransit(per = 2.0, RpRs = 0.1, ecc = 0, rhos = 1.0, 
                 b = 0.75, q1 = 1., q2 = 0., w = 0.)
 
.. code-block:: python  
  
  import matplotlib.pyplot as pl
  fig, ax = ps.PlotImage(M = 0., per = 1.0, RpRs = 0.5, ecc = 0, w = 0, rhos = 1.0,
                         b = 1.5, q1 = 1., q2 = 0., bkgimage = 'stars')
  pl.show() 
  
.. code-block:: python 
 
  ps.AnimateImage(per = 0.5, RpRs = 0.45, ecc = 0, rhos = 1.0, w = 0,
                  b = 0.3, u1 = 1., u2 = 0., delay = 0,
                  bkgimage = 'stars', nsteps = 100,
                  image_map = 'earth', size_inches = (12, 9),
                  lightcurve = True)

Or, if you're interested in the lightcurve model itself, the ``Transit`` class is
what you want:

.. code-block:: python
  
  import numpy as np
  
  # Orbital elements
  kwargs = {'rhos': 1.0, 'MpMs': 0.001, 'esw': 0.1, 
            'ecw': 0.1, 'per': 1.0, 'RpRs': 0.1, 
            't0': 0.0, 'q1': 1.0, 'q2': 0.0,
            'bcirc': 0.5}
  
  # Instantiate a transit object
  trn = ps.Transit(**kwargs) 
  
  # Compute the orbital solution and the lightcurve
  trn.Compute()
  
  # Bin the lightcurve to simulate an observation with finite exposure time
  trn.Bin()
  
  # Now interpolate to get the lightcurve on a grid of observation times
  t = np.arange(0., 10., ps.transit.KEPLONGCAD)
  flux = trn.Interpolate(t, 'binned')
        
Stay tuned; detailed documentation is coming soon!

Notes
=====

Feel free to change, adapt, or incorporate this code into your project, but please make sure to cite this repository, as well as `Mandel and Agol (2002) <http://adsabs.harvard.edu/abs/2002ApJ...580L.171M>`_, the transit model on which ``pysyzygy`` is based.
