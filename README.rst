py·sy·zy·gy
-----------
/ˈpīsizijē/ |speaker|

.. |speaker| image:: /../img/speaker.png?raw=True
             :target: http://www.astro.washington.edu/users/rodluger/pysyzygy.mp3

*noun*

**1.** Simple planet transit (`syzygy <http://en.wikipedia.org/wiki/Syzygy_%28astronomy%29>`_) visualizations, coded in Python.

**2.** The ``pysyzygy`` output looks like this:

.. image:: /../master/transit.png?raw=True
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
    The inclination of the orbit in degrees, in the range ``(0., 90.)``. 
    Default ``90.``

  * **MpMs** (*float*) -
    The ratio of the planet mass to the stellar mass. Default ``0.001``

  * **per** (*float*) -
    The period of the planet in days. Default ``2.``

  * **rhos** (*float*) -
    The density of the star in g/cm^3. Default ``1.4``
  
  * **RpRs** (*float*) -
    The ratio of the planet radius to the stellar radius, in the range ``[0., 0.5]``. 
    Default ``0.1``
  
  * **u1** (*float*) -
    The linear stellar limb darkening coefficient, in the range ``(0., 1.)``. 
    Default ``0.8``

  * **u2** (*float*) -
    The quadratic LD coefficient. Note that the sum ``u1 + u2`` must be in the 
    range ``(0., 1.)``. Default ``-0.4``

  * **w** (*float*) -
    The argument of periapsis in degrees. Default ``270.``

  * **exptime** (*float*) -
    The exposure time in days. Default ``0.020434`` (Kepler long cadence)

  * **exp_pts** (*int*) -
    The number of calls to the transit module in the exposure window. Default ``10``

  * **lightcurve** (*str*) -
    Options are ``"ideal"`` (plots only the ideal lightcurve), ``"observed"`` (plots
    only the observed lightcurve), or ``"both"``. Default ``"both"``

  * **ldplot** (*bool*) -
    Plot the limb darkening profile inset? Default ``True``

  * **plot_name** (*str*) -
    The name of the file to save the plot to. Default ``transit.png``

  * **plot_title** (*str*) -
    The plot title. Default ``""``

  * **show_params** (*bool*) -
    Show the planet parameters on the top plot? Default ``True``
  
  * **xypts** (*int*) -
    The number of points to use when plotting the orbit. Default ``1000``
    
Notes
=====

By default, ``pysyzygy`` outputs a PNG image file with two subplots. The top plot is the stellar lightcurve, centered at the planet transit. The light blue line is the **ideal** lightcurve, which is what a telescope would observe with an infinitely small exposure time. The dark blue line is the **observed** lightcurve, which is what you actually measure with a telescope like Kepler, whose default exposure time is about 30 minutes. The observed lightcurve is smoother and wider because of a blurring effect; it is what you get when you average the ideal lightcurve with a window of size equal to the exposure time.

The bottom plot is the orbit of the planet, as seen by us. Everything is to scale, normalized to the stellar radius. The star is colored according to the assigned limb darkening profile, which is plotted at the top left. The inset on the bottom left is a top view of the orbit, aligned so that the observer is at the bottom.

Feel free to change, adapt, or incorporate this code into your project, but please make sure to cite this repository, as well as `Mandel and Agol (2002) <http://adsabs.harvard.edu/abs/2002ApJ...580L.171M>`_, the transit model on which ``pysyzygy`` is based.
