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
   :scale: 50 %
   :align: center

**3.** ``pysyzygy`` can also generate sweet animations, such as this one (*WARNING: not to scale!*):

.. image:: img/transit.gif?raw=True
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

... is super easy.

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
  flux = trn.Interpolate(t)
        
Stay tuned; detailed documentation is coming soon!

Notes
=====

Feel free to change, adapt, or incorporate this code into your project, but please make sure to cite this repository, as well as `Mandel and Agol (2002) <http://adsabs.harvard.edu/abs/2002ApJ...580L.171M>`_, the transit model on which ``pysyzygy`` is based.
