#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
plot.py
-------

'''

import numpy as np
import matplotlib.pyplot as pl
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap, colorConverter
import subprocess
from PIL import Image
import planet
from transit import Transit, QUADRATIC, KIPPING, NONLINEAR

__all__ = ['PlotTransit', 'PlotImage', 'AnimateImage']

def I(r, limbdark):
  '''
  The standard quadratic limb darkening law.
  
  '''
  
  if limbdark.ldmodel == QUADRATIC:
    u1 = limbdark.u1
    u2 = limbdark.u2
    return (1-u1*(1-np.sqrt(1-r**2))-u2*(1-np.sqrt(1-r**2))**2)/(1-u1/3-u2/6)/np.pi
  elif limbdark.ldmodel == KIPPING:
    a = np.sqrt(limbdark.q1)
    b = 2*limbdark.q2
    u1 = a*b
    u2 = a*(1 - b)
    return (1-u1*(1-np.sqrt(1-r**2))-u2*(1-np.sqrt(1-r**2))**2)/(1-u1/3-u2/6)/np.pi
  elif limbdark.ldmodel == NONLINEAR:
    raise Exception('Nonlinear model not yet implemented!')                           # TODO!
  else:
    raise Exception('Invalid limb darkening model.')
  
def Star(ax, u1=1, u2=0, x=0, y=0, r=1, n=100, color=(1.0, 0.85, 0.1), zorder=-1):
  '''
  
  '''
    
  # Ensure RGB
  color = colorConverter.to_rgb(color)
  
  # Create a simple gradient colormap (0 = center, bright; 1 = limb, dark)
  dictylbk = {'red':  ((0.0, color[0], color[0]),
                       (1.0, 0.0, 0.0)),

             'green': ((0.0, color[1], color[1]),
                       (1.0, 0.0, 0.0)),

             'blue':  ((0.0, color[2], color[2]),
                       (1.0, 0.0, 0.0))
             }
  cmap = LinearSegmentedColormap('YlBk', dictylbk)
  
  # Limb darkening profile
  rad = np.linspace(0,1,n)[::-1]
  Ir = I(rad,u1,u2)/I(0,u1,u2)
  for ri,Iri in zip(rad,Ir):
    lightness = 0.95*Iri  
    color = cmap(1 - lightness)
    star = pl.Circle((x, y), ri*r, color=color, alpha=1., zorder = zorder)
    ax.add_artist(star)

def Planet(ax, x = 0, y = 0, z = 0, r = 0.25, long0 = 0., image_map='maps/earth.jpg'):
  '''
  
  '''
  
  p = planet.Planet(RpRs = r, image_map=image_map)
  p.gen_image(xyz = (x,y,z))
  p.plot_image(ax=ax, extent=(x-r,x+r,y-r,y+r), long0 = long0)

def Trail(ax, x, y, z, f, r = 0.05, color = "#4682b4", ndots = None, RpRs = 0.25):
  '''
  
  '''
  
  if ndots is None: 
    ndots = int(len(x)/2.)

  # x,y must span exactly one orbit
  for i in range(f - ndots, f):
    # Don't draw trails inside the planet
    if np.sqrt((x[i] - x[f])**2 + (y[i] - y[f])**2 + (z[i] - z[f])**2) <= 0.95*RpRs: 
      continue

    alpha = 0.2*((i - (f - ndots))/(1.*ndots))**2
    
    if z[i] < z[f]:
      if z[i] < 0:
        # In front of planet and star
        zorder = 3
      else:
        # In front of planet, behind star
        zorder = 1
    else:
      if z[i] < 0:
        # Behind planet, in front of star
        zorder = -1
      else:
        # Behind planet and star
        zorder = -3
    trail = pl.Circle((x[i], y[i]), r, color=color, zorder = zorder, alpha = alpha)
    ax.add_artist(trail)

def PlotTransit(**kwargs):
  '''
    
  '''
  
  ldplot = kwargs.pop('ldplot', True)
  compact = kwargs.pop('compact', True)
  plottitle = kwargs.pop('plottitle', "")
  plotname = kwargs.pop('plotname', "transit")
  xlim = kwargs.pop('xlim', None)
  t0 = kwargs.pop('t0', 0.)
  
  # Plotting
  fig = pl.figure()
  fig.set_size_inches(12,8)
  fig.subplots_adjust(hspace=0.3)
  ax1, ax2 = pl.subplot(211), pl.subplot(212)

  kwargs.update({'fullorbit': True})
  trn = Transit(**kwargs)
  trn.Compute()

  time = trn.arrays.time + t0
  flux = trn.arrays.flux

  ax1.plot(time, flux, '-', color='DarkBlue')
  rng = np.max(flux) - np.min(flux)
  if rng > 0:
    ax1.set_ylim(np.min(flux) - 0.1*rng, np.max(flux) + 0.1*rng)
    left = np.argmax(flux < 1.)
    right = np.argmax(flux[left:] == 1.) + left
    rng = time[right] - time[left]
    ax1.set_xlim(time[left] - rng, time[right] + rng)

  ax1.set_xlabel('Time (Days)', fontweight='bold')
  ax1.set_ylabel('Normalized Flux', fontweight='bold')

  # Sky-projected motion
  x = trn.arrays.x
  y = trn.arrays.y
  z = trn.arrays.z
  
  # Mask the star
  for j in range(trn.arrays.npts):
    if (x[j]**2 + y[j]**2) < 1. and (z[j] > 0):
      x[j] = np.nan
      y[j] = np.nan
  
  # The star
  r = np.linspace(0,1,100)
  Ir = I(r,trn.limbdark)/I(0,trn.limbdark)
  
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
  trn.transit.bcirc = trn.transit.aRs                                                 # This ensures we are face-on
  trn.Compute()
  xp = trn.arrays.x
  yp = trn.arrays.y
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
  planet = pl.Circle((0, ycenter), trn.transit.RpRs, color='DarkBlue', alpha=1.)
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
  planet = pl.Circle((0, ycenter), trn.transit.RpRs, color='DarkBlue', alpha=1.)
  ax2.add_artist(planet)
  
  # Force aspect
  if xlim is None:
    xlim = 1.1 * max(np.nanmax(x), np.nanmax(-x))
  ax2.set_ylim(-xlim/3.2,xlim/3.2)
  ax2.set_xlim(-xlim,xlim)
  
  ax2.set_xlabel(r'X (R$_\star$)', fontweight='bold')
  ax2.set_ylabel(r'Y (R$_\star$)', fontweight='bold')
  ax1.set_title(plottitle, fontsize=12)
  
  if not compact:
    rect = 0.925,0.55,0.2,0.35
    ax3 = fig.add_axes(rect)
    ax3.xaxis.set_visible(False)
    ax3.yaxis.set_visible(False)

    # TODO!!!!!!!

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

  fig.savefig(plotname, bbox_inches='tight')
  pl.close()

def PlotImage(t=0., t0=0., rhos=0.5, RpRs=0.25, MpMs=0.01, per=1., bcirc=0.5, 
         esw=0., ecw=0., u1=1, u2=0, rot=0., bkgcolor = 'white', bkgimage = None,
         long0 = 0., image_map = 'maps/earth.jpg'):
  '''
  
  '''

  x, y, z = xyz([t], t0, rhos, MpMs, per, bcirc, esw, ecw, mask_star = False)

  fig = pl.figure()
  ax = pl.subplot(111)
  ax.set_xlim(-3.5,3.5)
  ax.set_ylim(-3,3)
  ax.xaxis.set_visible(False)
  ax.yaxis.set_visible(False)
  ax.set_aspect('equal')
  ax.patch.set_facecolor(bkgcolor)
  if z[0] < 0: 
    zorder = -2
  else: 
    zorder = 2
  Star(ax, u1 = u1, u2 = u2, zorder = zorder)
  Planet(ax, x = x[0], y = y[0], z = z[0], r = RpRs, long0 = long0, image_map=image_map)
  
  return fig, ax

def AnimateImage(t0=0., rhos=0.5, RpRs=0.25, MpMs=0.01, per=1., bcirc=0.5, 
         esw=0., ecw=0., u1=1, u2=0, rot=0., bkgcolor = 'white', bkgimage = None,
         nsteps = 1000, dpy = 4, trail = True):
  '''
  Note that dpy (= days_per_year) can be set negative for retrograde rotation
  
  '''
  subprocess.call(['mkdir', '-p', 'tmp'])
  frames = range(nsteps)
  rotation = (np.linspace(0, -dpy, nsteps) % 1)*360 - 180
  
  x, y, z = xyz(np.linspace(-per/2.,per/2.,nsteps), t0, rhos, MpMs, per, bcirc, esw, ecw, mask_star = False)
  
  for f, long0 in zip(frames, rotation):
    fig = pl.figure()
    ax = pl.subplot(111)
    ax.set_xlim(-3.5,3.5)
    ax.set_ylim(-3,3)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.set_aspect('equal')
    if z[f] < 0: 
      zorder = -2
    else: 
      zorder = 2
    Star(ax, u1 = u1, u2 = u2, zorder = zorder)
    Planet(ax, x = x[f], y = y[f], z = z[f], r = RpRs, long0 = long0)
    if trail: Trail(ax, x, y, z, f, RpRs = RpRs)

    if bkgimage is not None:
      im = Image.open(bkgimage)
      ax.imshow(im, extent=(-3.5,3.5,-3,3), zorder=-99)
    else:
      ax.patch.set_facecolor(bkgcolor)

    fig.savefig('tmp/%03d.png' % f, bbox_inches = 'tight')
    pl.close()
  
  # Make gif
  subprocess.call(['convert', '-delay', '5', '-loop', '-1', 'tmp/*.png', 'transit.gif'])

if __name__ == '__main__':
  PlotTransit(per = 1., RpRs = 0.1, ecw = 0.3, esw = 0.3, rhos = 1.0,
              bcirc = 0.3, t0 = 1., u1 = 1., u2 = 0.)