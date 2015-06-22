#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
main.py
-------
A really simple implementation of the Mandel & Agol (2002) transit model
for planet orbits of arbitrary eccentricity. Calculates both the transit
lightcurve and the orbital path as seen by the observer.

Copyright Rodrigo Luger, May 2015.

"""

import matplotlib.pyplot as pl
import numpy as np
from scipy.optimize import newton
import os, sys, subprocess
G = 6.672e-8
DAYSEC = 86400

try:
  import transit
except ImportError:
  raise ImportError('Failed to import the transit module. Did you compile the '
                    'transit.f file? To do so, cd to pysyzygy/pysyzygy and '
                    'run:\n\n    >>> f2py -c transit.f -m transit\n\n')

__all__ = ['lightcurve', 'xy', 'plot', 'I', 'animate']

def transit_times(tstart, tstop, t0, per, tdur):
  '''
  Calculates all transit times between ``tstart`` and ``tstop``.

  '''
  n, r = divmod(tstart - t0, per)
  if r < tdur:
    t0 = t0 + n*per
  else:
    t0 = t0 + (n + 1)*per
  
  return np.arange(t0, tstop, per)

def diff(E, e, M):
  '''
  Kepler's equation relating the eccentric anomaly ``E``, the eccentricity
  ``e``, and the mean anomaly ``M``.
  
  '''
  return E - e*np.sin(E) - M

def der1(E, e, M):
  '''
  First derivative of Kepler's equation.
  
  '''
  return 1. - e*np.cos(E)

def der2(E, e, M):
  '''
  Second derivative of Kepler's equation.
  
  '''
  return e*np.sin(E)

def I(r, u1, u2):
  '''
  The standard quadratic limb darkening law.
  
  '''
  return (1-u1*(1-np.sqrt(1-r**2))-u2*(1-np.sqrt(1-r**2))**2)/(1-u1/3-u2/6)/np.pi

def lightcurve(t, **kwargs):
  '''
  Returns an array of flux values for each point in the array ``t`` corresponding
  to the transit lightcurve of the planet.
  
  :Keyword Arguments:
  
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
  
  * **t0** (*float*) -
    The time of transit center for the first transit. Default ``0.0``
  
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
    
  .. note:: To get the ideal lightcurve (unblurred by the exposure time), \
  simply set ``exp_pts`` to 1.
  
  '''

  per = kwargs.get('per', 2.)
  assert (0 < per), 'Invalid value for the period.'
  i = kwargs.get('i', 90.)*np.pi/180.
  assert (0 <= i < 2*np.pi), 'Invalid value for the inclination.'
  rhos = kwargs.get('rhos', 1.4)
  assert (0 < rhos), 'Invalid value for the stellar density.'
  MpMs = kwargs.get('MpMs', 0.001)
  assert (0 <= MpMs), 'Invalid value for the planet mass.'
  e = kwargs.get('e', 0.0)
  assert (0 <= e < 1.0), 'Invalid value for the eccentricity.'
  w = kwargs.get('w', 270.)*np.pi/180.
  assert (0 <= w < 2*np.pi), 'Invalid value for omega.'
  u1 = kwargs.get('u1', 0.8)
  u2 = kwargs.get('u2', -0.4)
  assert (0. <= u1 <= 1.) and (0. <= u1 + u2 <= 1.), \
         'Invalid values for the limb darkening coefficients.'
  RpRs = kwargs.get('RpRs', 0.1)
  assert (0 < RpRs < 0.5), 'Invalid value for the planet radius.'
  exptime = kwargs.get('exptime', 1765.5/86400)
  assert (0 < exptime), 'Invalid value for the exposure time.'
  exp_pts = kwargs.get('exp_pts', 10)
  assert (1 <= exp_pts), 'Invalid value for the number of exposure points.'
  t0 = kwargs.get('t0', 0.)
         
  # Derived stuff
  flux = np.ones_like(t, dtype=float)
  aRs = ((G*rhos*(1. + MpMs)*(per*DAYSEC)**2)/(3*np.pi))**(1./3)
  if aRs*(1-e) <= 1.:
    raise Exception('Error: star-crossing orbit!')
  bcirc = aRs*np.cos(i)
  esw = e*np.sin(w)
  ecw = e*np.cos(w)
  becc = bcirc*(1.-e**2)/(1.+e*np.sin(w-np.pi))
  with np.errstate(invalid='ignore'):
    tdur = per/(2*np.pi)*np.arcsin(((1+RpRs)**2-becc**2)**0.5/(np.sin(i)*aRs))        # (half) transit duration
    tdur *= np.sqrt(1.-e**2.)/(1.+e*np.sin(w-np.pi))                                  # Correct for eccentricity
  tdur *= 1.5                                                                         # The correction may be off for high eccentricity, so let's do this for safety
  if np.isnan(tdur):
    return flux                                                                       # No transits!
  
  tN = transit_times(t[0], t[-1], t0, per, tdur)
  ntrans = len(tN)
  try:
    err = transit.transit(t,flux,bcirc,rhos,MpMs,esw,ecw,per,u1,u2,RpRs,
                          exptime,tN,exp_pts,ntrans,len(t))
  except:
    raise Exception('Error: something went wrong while computing the transit.')
  
  return flux

def xyz(t, t0, rhos, MpMs, per, bcirc, esw, ecw, mask_star = True):
  '''
  Compute the sky-projected coordinates (x,y) of the orbit given an array of times 
  ``t``, the transit center time ``t0``, the stellar density ``rhos``, the mass ratio
  ``MpMs``, the period ``per``, the circular impact parameter ``bcirc``, 
  e*sin(omega) (``esw``) and e*cos(omega) (``ecw``).
  
  :returns: A tuple ``(x, y)`` of the sky-projected coordinates over the array of \
  times ``t``. If ``mask_star`` is ``True``, the part of the orbit that is obscured \
  by the star is masked with ``np.nan``
  
  '''
  aRs = ((G*rhos*(1. + MpMs)*(per*DAYSEC)**2)/(3*np.pi))**(1./3)                      # semi-major axis / stellar radius
  sini = np.sqrt(1. - (bcirc/aRs)**2 )
  e = np.sqrt( esw**2 + ecw**2 )
  w = np.arctan2(esw, ecw)
    
  fi = 3.*np.pi/2 - w                                                                 # True anomaly at transit center (approximate)
  t_peri = t0 + per*np.sqrt(1 - e**2)/(2*np.pi)*(e*np.sin(fi)/(1+e*np.cos(fi)) - 
    2./np.sqrt(1-e**2)*np.arctan2(np.sqrt(1-e**2)*np.tan(fi/2.), 1+e))                # Time of pericenter

  x = np.zeros_like(t)                                                                # Sky-projected coordinates
  y = np.zeros_like(t)
  z = np.zeros_like(t)                                                                # perp. to sky plane
  for j, ti in enumerate(t):
    M = 2*np.pi/per*(ti - t_peri)                                                     # Mean anomaly
    E = newton(diff, M, args = (e, M), fprime = der1, fprime2 = der2)                 # Eccentric anomaly
    f = 2*np.arctan(((1. + e)/(1. - e))**(1./2)*np.tan(E/2.))                         # True anomaly
    rRs = aRs*(1. - e**2)/(1. + e*np.cos(f))                                          # r/Rs
    b = rRs*np.sqrt(1. - (np.sin(w + f)*sini)**2)                                     # Impact parameter
    
    x[j] = rRs*np.cos(w + f)
    z[j] = rRs*np.sin(w + f)
    if (b**2 - x[j]**2) < 1e-10:                                                      # Prevent numerical error
      y[j] = 0.
    else:
      if (0 < (f + w) % (2*np.pi) < np.pi):
        y[j] = np.sqrt(b**2 - x[j]**2)
      else:
        y[j] = -np.sqrt(b**2 - x[j]**2)
    if (mask_star and (((x[j]**2 + y[j]**2) < 1.) and (y[j] > 0))):
      x[j] = np.nan                                                                   # These get masked when plotting
      y[j] = np.nan
      z[j] = np.nan
  return x, y, z

def plot(**kwargs):
  '''
  The main plotting routine. Plots the lightcurve and the sky-projected orbit, and
  saves the image to disk.
  
  :Keyword Arguments:
  
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

  * **lc** (*str*) -
    Options are ``"ideal"`` (plots only the ideal lightcurve), ``"observed"`` (plots
    only the observed lightcurve), or ``"both"``. Default ``"both"``

  * **ldplot** (*bool*) -
    Plot the limb darkening profile inset? Default ``True``

  * **plot_name** (*str*) -
    The name of the file to save the plot to. Default ``transit.png``

  * **plot_title** (*str*) -
    The plot title. Default ``""``

  * **compact** (*bool*) -
    Compact plotting? Default ``False``
  
  * **xypts** (*int*) -
    The number of points to use when plotting the orbit. Default ``1000``
    
  '''
  
  # Keywords
  per = kwargs.get('per', 2.)
  assert (0 < per), 'Invalid value for the period.'
  i = kwargs.get('i', 90.)*np.pi/180.
  assert (0 <= i < 2*np.pi), 'Invalid value for the inclination.'
  rhos = kwargs.get('rhos', 1.4)
  assert (0 < rhos), 'Invalid value for the stellar density.'
  MpMs = kwargs.get('MpMs', 0.001)
  assert (0 <= MpMs), 'Invalid value for the planet mass.'
  e = kwargs.get('e', 0.0)
  assert (0 <= e < 1.0), 'Invalid value for the eccentricity.'
  w = kwargs.get('w', 270.)*np.pi/180.
  assert (0 <= w < 2*np.pi), 'Invalid value for omega.'
  u1 = kwargs.get('u1', 0.8)
  u2 = kwargs.get('u2', -0.4)
  assert (0. <= u1 <= 1.) and (0. <= u1 + u2 <= 1.), \
         'Invalid values for the limb darkening coefficients.'
  RpRs = kwargs.get('RpRs', 0.1)
  assert (0 < RpRs < 0.5), 'Invalid value for the planet radius.'
  exptime = kwargs.get('exptime', 1765.5/86400)
  assert (0 < exptime), 'Invalid value for the exposure time.'
  exp_pts = kwargs.get('exp_pts', 10)
  assert (1 <= exp_pts), 'Invalid value for the number of exposure points.'
  plot_name = kwargs.get('plot_name', 'transit')
  plot_title = kwargs.get('plot_title', '')
  compact = kwargs.get('compact', False)
  xypts = kwargs.get('xypts', 1000)
  assert (0 < xypts), 'Invalid value for xypts.'
  lc = kwargs.get('lc', 'both')
  assert (lc == 'both') or (lc == 'ideal') or \
         (lc == 'observed'), 'Invalid option for lc.'
  ldplot = kwargs.get('ldplot', True)
  WFAC = 3.0
  
  # Derived stuff
  aRs = ((G*rhos*(1. + MpMs)*(per*DAYSEC)**2)/(3*np.pi))**(1./3)
  if aRs*(1-e) <= 1.:
    raise Exception('Error: star-crossing orbit!')
  bcirc = aRs*np.cos(i)
  esw = e*np.sin(w)
  ecw = e*np.cos(w)
  becc = bcirc*(1.-e**2)/(1.+e*np.sin(w-np.pi))
  with np.errstate(invalid='ignore'):
    window = per/(2*np.pi)*np.arcsin(((1+RpRs)**2-becc**2)**0.5/(np.sin(i)*aRs))
    window *= np.sqrt(1.-e**2.)/(1.+e*np.sin(w-np.pi))
  t = np.arange(-WFAC*window,WFAC*window,1.e-4)
  tN = np.array([0.])                                                                 # Array of transit times
  ntrans = len(tN)                                                                    # Number of transits

  # Plotting
  fig = pl.figure()
  fig.set_size_inches(12,8)
  fig.subplots_adjust(hspace=0.3)
  ax1, ax2 = pl.subplot(211), pl.subplot(212)

  if not np.isnan(window):
    try:
      # The ideal lc
      if (lc == 'both') or (lc == 'ideal'):
        flux = np.ones_like(t, dtype=float)
        err = np.array([0], dtype=int)                                                # TODO: Add error handling
        transit.transit(t,flux,bcirc,rhos,MpMs,esw,ecw,per,u1,u2,RpRs,
                        exptime,tN,1,err,ntrans,len(t))
        if lc == 'both':
          ax1.plot(t,flux, '-', color='DarkBlue', alpha = 0.25, label='Ideal')
        elif lc == 'ideal':
          ax1.plot(t,flux, '-', color='DarkBlue', alpha = 1.0, label='Ideal')
  
      # The blurred lc due to the finite exposure time
      if (lc == 'both') or (lc == 'observed'):
        flux = np.ones_like(t, dtype=float)
        err = np.array([0], dtype=int)
        transit.transit(t,flux,bcirc,rhos,MpMs,esw,ecw,per,u1,u2,RpRs,
                        exptime,tN,exp_pts,err,ntrans,len(t))
        ax1.plot(t,flux, '-', color='DarkBlue', label='Observed')
    
      rng = np.max(flux) - np.min(flux)
      ax1.set_ylim(np.min(flux) - 0.1*rng, np.max(flux) + 0.1*rng)
      ax1.set_xlim(-WFAC*window,WFAC*window)
      ax1.set_xlabel('Time (Days)', fontweight='bold')
      ax1.set_ylabel('Normalized Flux', fontweight='bold')
      ax1.legend(loc='lower left', fontsize=10)
      
    except:
      ax1.text(0.5,0.5,'No transit',ha='center',va='center')
      pass
  else:
    ax1.text(0.5,0.5,'No transit',ha='center',va='center')

  # Sky-projected motion
  t = np.linspace(-per/2,per/2,xypts)
  x, y, _ = xyz(t, tN[0], rhos, MpMs, per, bcirc, esw, ecw)
  
  # The star
  r = np.linspace(0,1,100)
  Ir = I(r,u1,u2)/I(0,u1,u2)
  
  for ri,Iri in zip(r[::-1],Ir[::-1]):
    star = pl.Circle((0, 0), ri, color=str(0.95*Iri), alpha=1.)
    ax2.add_artist(star)

  # Inset: Limb darkening
  if ldplot:
    if compact:
      inset1 = pl.axes([0.145, 0.32, .09, .1])
    else:
      inset1 = fig.add_axes([0.925,0.3,0.2,0.15])
    inset1.plot(r,Ir,'k-')
    pl.setp(inset1, xlim=(-0.1,1.1), ylim=(-0.1,1.1), xticks=[0,1], yticks=[0,1])
    for tick in inset1.xaxis.get_major_ticks() + inset1.yaxis.get_major_ticks():
      tick.label.set_fontsize(8)
    inset1.set_ylabel(r'I/I$_0$', fontsize=8, labelpad=-8)
    inset1.set_xlabel(r'r/R$_\star$', fontsize=8, labelpad=-8)
    inset1.set_title('Limb Darkening', fontweight='bold', fontsize=9)
    
  # Inset: Top view of orbit
  if compact:
    inset2 = pl.axes([0.135, 0.115, .1, .1])
  else:
    inset2 = fig.add_axes([0.925,0.1,0.2,0.15])
  pl.setp(inset2, xticks=[], yticks=[])
  t = np.linspace(-per/2,per/2,1000)
  xp, yp, _ = xyz(t, tN[0], rhos, MpMs, per, aRs, esw, ecw)
  inset2.plot(xp, yp, '-', color='DarkBlue', alpha=0.5)
  # Draw some invisible dots at the corners to set the window size
  xmin, xmax, ymin, ymax = np.nanmin(xp), np.nanmax(xp), np.nanmin(yp), np.nanmax(yp)
  xrng = xmax - xmin
  yrng = ymax - ymin
  xmin -= 0.1*xrng; xmax += 0.1*xrng;
  ymin -= 0.1*yrng; ymax += 0.1*yrng;
  inset2.scatter([xmin,xmin,xmax,xmax], [ymin,ymax,ymin,ymax], alpha = 0.)
  # Plot the star
  for ri,Iri in zip(r[::-10],Ir[::-10]):
    star = pl.Circle((0, 0), ri, color=str(0.95*Iri), alpha=1.)
    inset2.add_artist(star)
  # Plot the planet
  ycenter = yp[np.where(np.abs(xp) == np.nanmin(np.abs(xp)))][0]
  while ycenter > 0:
    xp[np.where(np.abs(xp) == np.nanmin(np.abs(xp)))] = np.nan
    ycenter = yp[np.where(np.abs(xp) == np.nanmin(np.abs(xp)))][0]
  planet = pl.Circle((0, ycenter), RpRs, color='DarkBlue', alpha=1.)
  inset2.add_artist(planet)
  inset2.set_title('Top View', fontweight='bold', fontsize=9)
  inset2.set_aspect('equal','datalim')
  
  # The orbit itself
  with np.errstate(invalid='ignore'):
    ax2.plot(x, y, '-', color='DarkBlue')
  
  # The planet
  ycenter = y[np.where(np.abs(x) == np.nanmin(np.abs(x)))][0]
  while ycenter > 0:
    x[np.where(np.abs(x) == np.nanmin(np.abs(x)))] = np.nan
    ycenter = y[np.where(np.abs(x) == np.nanmin(np.abs(x)))][0]
  planet = pl.Circle((0, ycenter), RpRs, color='DarkBlue', alpha=1.)
  ax2.add_artist(planet)
  
  # Force aspect
  ax2.set_ylim(-2.5,2.5)
  ax2.set_xlim(-8,8)
  
  ax2.set_xlabel(r'X (R$_\star$)', fontweight='bold')
  ax2.set_ylabel(r'Y (R$_\star$)', fontweight='bold')
  ax1.set_title(plot_title, fontsize=12)
  
  if not compact:
    rect = 0.925,0.55,0.2,0.35
    ax3 = fig.add_axes(rect)
    ax3.xaxis.set_visible(False)
    ax3.yaxis.set_visible(False)

    # Table of parameters
    ltable = [ r'$P:$',
               r'$e:$',
               r'$i:$',
               r'$\omega:$',
               r'$\rho_\star:$',
               r'$M_p:$',
               r'$R_p:$',
               r'$u_1:$',
               r'$u_2:$']
    rtable = [ r'$%.4f\ \mathrm{days}$' % per,
               r'$%.5f$' % e,
               r'$%.4f^\circ$' % (i*180./np.pi),
               r'$%.3f^\circ$' % (w*180./np.pi),
               r'$%.5f\ \mathrm{g/cm^3}$' % rhos,
               r'$%.5f\ M_\star$' % MpMs,
               r'$%.5f\ R_\star$' % RpRs,
               r'$%.5f$' % u1,
               r'$%.5f$' % u2]
    yt = 0.875
    for l,r in zip(ltable, rtable):
      ax3.annotate(l, xy=(0.25, yt), xycoords="axes fraction", ha='right', fontsize=16)
      ax3.annotate(r, xy=(0.35, yt), xycoords="axes fraction", fontsize=16)
      yt -= 0.1

  fig.savefig(plot_name, bbox_inches='tight')
  pl.close()

def animate(name, vals, **kwargs):
  """
  Generates an animated gif by looping over a given parameter.
  
  :param str name: The name of the parameter to loop over. Must be one of the \
  keyword arguments to ``plot()``
  
  :param array_like vals: The array of values of the parameter
  
  :Keyword Arguments:
  
  * **delay** (*int*) -
  The delay between frames in ms. Default ``20``
  
  * **loop** (*bool*) -
  Loop the gif? Default ``True``
  
  * **plot_name** (*str*) -
  The name of the file to save the animation to. Default ``anim.gif``
  
  """
  
  plot_name = kwargs.get('plot_name', 'anim')
  if not plot_name.endswith('.gif'): plot_name += '.gif'
  delay = kwargs.get('delay', 20)
  delay = str(delay)
  loop = kwargs.get('loop', True)
  if loop: loop = '-1'
  else: loop = '0'
  
  # Let's keep things realistic
  assert (len(vals) < 1000), 'The size of the array must be less than 1000.'
  
  # Create a temporary directory
  if not os.path.exists('tmp'):
    os.makedirs('tmp')
  
  for i,v in enumerate(vals):
    kwargs.update({'plot_name': 'tmp/%03d.png' % i})
    kwargs.update({name: v})
    sys.stdout.write("\rPlotting frame %03d/%03d..." % (i + 1, len(vals)))
    sys.stdout.flush()
    plot(**kwargs)
  
  # Now clear the line
  sys.stdout.write("\r "*30 + "\r")
  
  # Make gif
  subprocess.call(['convert', '-delay', delay, '-loop', loop, 'tmp/*.png', plot_name])
  subprocess.call(['rm', '-r', 'tmp'])

if __name__ == '__main__':
  # Produce a sample plot
  plot(e = 0.7, per = 1., w = 350., xypts=10000, 
       rhos=1., u1 = 0.75, u2 = -0.1, lc = 'ideal', i = 65)
  
  # Produce a sample animation
  #w = np.linspace(170., 375., 100) % 360.
  #w = np.concatenate((w, w[::-1]))
  #animate('w', w, e=0.7, per=1., rhos=1., u1=0.5, u2=0., lc = 'ideal', i = 65)