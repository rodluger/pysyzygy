import numpy as np
import matplotlib.pyplot as pl
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap, colorConverter
import subprocess
import planet
from main import xyz
from PIL import Image

def I(r, u1, u2):
  '''
  The standard quadratic limb darkening law.
  
  '''
  
  return (1-u1*(1-np.sqrt(1-r**2))-u2*(1-np.sqrt(1-r**2))**2)/(1-u1/3-u2/6)/np.pi
  
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

def Planet(ax, x = 0, y = 0, z = 0, r = 0.25, rot = 0.5):
  '''
  
  '''

  # TODO: Add rotation

  p = planet.Planet(RpRs = r)
  p.gen_image(xyz = (x,y,z))
  p.plot_image(ax=ax, image='maps/earth.jpg', extent=(x-r,x+r,y-r,y+r))

def Trail(ax, x, y, z, f, r = 0.05, ndots = 200, color = "#4682b4"):
  '''
  
  '''
  
  if ndots >= len(x): 
    ndots = len(x)
  
  # x,y must span exactly one orbit
  for i in range(f - ndots, f):
    alpha = min(1.0, 2.0/(f - i))
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
    trail = pl.Circle((x[i], y[i]), r, color=color, zorder = -2, alpha = alpha)
    ax.add_artist(trail)

def Plot(t=0., t0=0., rhos=0.5, RpRs=0.25, MpMs=0.01, per=1., bcirc=0.5, 
         esw=0., ecw=0., u1=1, u2=0, rot=0., bkgcolor = 'white', bkgimage = None):
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
  Planet(ax, x = x[0], y = y[0], z = z[0], r = RpRs, rot = rot)
  
  return fig, ax

def Animate(t0=0., rhos=0.5, RpRs=0.25, MpMs=0.01, per=1., bcirc=0.5, 
         esw=0., ecw=0., u1=1, u2=0, rot=0., bkgcolor = 'white', bkgimage = None,
         nsteps = 1000, dpy = 4, trail = True):
  '''
  Note that dpy (= days_per_year) can be set negative for retrograde rotation
  
  '''
  subprocess.call(['mkdir', '-p', 'tmp'])
  frames = range(nsteps)
  rotation = np.linspace(0, -dpy, nsteps) % 1
  x, y, z = xyz(np.linspace(-per/2.,per/2.,nsteps), t0, rhos, MpMs, per, bcirc, esw, ecw, mask_star = False)
  
  for f, rot in zip(frames, rotation):
    
    fig = pl.figure()
    ax = pl.subplot(111)
    ax.set_xlim(-3.5,3.5)
    ax.set_ylim(-3,3)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.set_aspect('equal')
    if z[f] < 0: 
      zorder = -1
    else: 
      zorder = 1
    Star(ax, u1 = u1, u2 = u2, zorder = zorder)
    Planet(ax, x = x[f], y = y[f], z = z[f], r = RpRs, rot = rot)
    if trail: Trail(ax, x, y, z, f)

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
  Animate(bkgcolor='k', nsteps=1000, bkgimage='maps/stars.jpg')
  subprocess.call(['rm', '-r', 'tmp'])