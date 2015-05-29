py·sy·zy·gy
-----------
/ˈpīsizijē/ `|speaker|</../img/pysyzygy.mp3?raw=true>`_

.. |speaker| image:: /../img/speaker.png?raw=true
   [height=14 width=14]

*noun*

**1.** Simple planet transit visualizations, coded in Python.

**2.** The ``pysyzygy`` output looks like this:

.. image:: https://raw.githubusercontent.com/rodluger/pysyzygy/master/transit.png
   :alt: pysyzygy
   :scale: 50 %
   :align: center

Installation
============
Clone the repository and run

>>> f2py -c transit.f -m transit

to build the transit module. Then call or import ``pysyzygy`` to begin plotting transits!

Calling pysyzygy...
===================

... is super easy.

.. code-block:: python
  
  import pysyzygy
  pysyzygy.plot(**kwargs)

where ``kwargs`` can be any of:

* **e** (*float*) -
  The eccentricity of the orbit, in the range ``(0., 1.]``. Default ``0.``

* **i** (*float*) -
  The inclination of the orbit in degrees, in the range ``(0., 90.)``. Default ``90.``

* **MpMs** (*float*) -
  The ratio of the planet mass to the stellar mass. Default ``0.001``

* **per** (*float*) -
  The period of the planet in days. Default ``2.``

* **rhos** (*float*) -
  The density of the star in g/cm^3. Default ``1.4``
  
* **RpRs** (*float*) -
  The ratio of the planet radius to the stellar radius, between ``[0., 0.5]``. Default ``0.1``
  
* **u1** (*float*) -
  The linear stellar limb darkening coefficient, in the range ``(0., 1.)``. Default ``0.8``

* **u2** (*float*) -
  The quadratic LD coefficient. ``u1 + u2`` must be in the range ``(0., 1.)``. Default ``-0.4``

* **w** (*float*) -
  The argument of periapsis in degrees. Default ``270.``


* **exptime** (*float*) -
  The exposure time in days. Default ``0.020434`` (Kepler long cadence)

* **exp_pts** (*int*) -
  The number of calls to the transit module in the exposure window. Default ``10``

* **plot_name** (*str*) -
  The name of the file to save the plot to. Default ``transit.png``

* **plot_title** (*str*) -
  The plot title. Default ""
