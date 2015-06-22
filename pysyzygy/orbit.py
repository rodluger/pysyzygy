import numpy as np
import matplotlib.pyplot as pl
from PIL import Image, ImageOps, ImageDraw
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
import pysyzygy as pz
import subprocess

def I(r, u1, u2):
  '''
  The standard quadratic limb darkening law.
  
  '''
  
  return (1-u1*(1-np.sqrt(1-r**2))-u2*(1-np.sqrt(1-r**2))**2)/(1-u1/3-u2/6)/np.pi

def L(x, y, alpha):
  '''
  The intensity of a Lambertian sphere given its sky-projected x and y coordinates
  
  '''
  

  # TODO!!!
  #theta = np.arctan2(x, y)
  #r = np.sqrt(x**2 + y**2)
  #phi = 2*np.arcsin(r/2.)
  #return np.sin(theta)*np.cos(phi)
  
  # DEBUG
  return 1 - (0.5)*(x**2 + y**2)
  
def Star(ax, u1 = 1, u2 = 0, x = 0, y = 0, r = 1, n = 100, cmap=None, behind=True):
  '''
  
  '''
  
  if behind:
    zorder = 0
  else:
    zorder = 99
  
  # Limb darkening profile
  rad = np.linspace(0,1,n)[::-1]
  Ir = I(rad,u1,u2)/I(0,u1,u2)
  for ri,Iri in zip(rad,Ir):
    lightness = 0.95*Iri
    
    # Create a simple yellow-black colormap (0 = center, bright; 1 = limb, dark)
    if cmap is None:
      dictylbk = {'red':  ((0.0, 1.0, 1.0),
                           (1.0, 0.0, 0.0)),

                 'green': ((0.0, 0.85, 0.85),
                           (1.0, 0.0, 0.0)),

                 'blue':  ((0.0, 0.1, 0.1),
                           (1.0, 0.0, 0.0))
                 }
      cmap = LinearSegmentedColormap('YlBk', dictylbk)
    
    # Draw the star
    color = cmap(1 - lightness)
    star = pl.Circle((x, y), ri*r, color=color, alpha=1., zorder = zorder)
    ax.add_artist(star)

def Planet(ax, x = 0, y = 0, r = 0.25, sz = 128, file = 'jupiter3.jpg', 
           alpha = 0., rot = 0.5):
  '''
  
  '''

  # A black background with a mask
  bkg = Image.new('RGBA', (sz, sz), 'white')
  bmk = Image.new('L', (sz, sz), 0)
  
  # A Lambertian mask
  mask = Image.new('L', size, 0)
  for xi in np.arange(0.,sz):
    ym = np.sqrt((sz/2.)**2 - (xi -  sz/2.)**2)
    for yi in np.arange( sz/2. - ym, sz/2. + ym):
      mask.putpixel((int(xi),int(yi)), 255*L(2.*xi/sz - 1, 2.*yi/sz - 1, alpha))
      bkg.putpixel((int(xi),int(yi)), 0)
      bmk.putpixel((int(xi),int(yi)), 255)
  
  # Plot the background
  bkg.putalpha(bmk)
  ax.imshow(bkg, extent=(x-r,x+r,y-r,y+r))
  
  # Overplot the alpha-corrected image
  im = Image.open(file)
  im = ImageOps.fit(im, mask.size, centering=(rot, 0))
  im.putalpha(mask)
  ax.imshow(im, extent=(x-r,x+r,y-r,y+r))

def Trail(ax, x, y, f, r = 0.05):
  '''
  
  '''
  
  # x,y must span exactly one orbit
  for i in range(f - 200, f):
    alpha = min(1.0, 2.0/(f - i))
    trail = pl.Circle((x[i], y[i]), r, color="#4682b4", zorder = 0, alpha = alpha)
    ax.add_artist(trail)

def Animate(nsteps = 1000, dpy = 4):
  '''
  Note that dpy (= days_per_year) can be set negative for retrograde rotation
  
  '''

  frames = range(nsteps)
  rotation = np.linspace(0, -dpy, nsteps) % 1
  x, y = pz.xy(np.linspace(-0.5,0.5,nsteps), 0, 0.5, 0.01, 1., 0.5, 0., 0., 
               mask_star = False)

  for f, rot in zip(frames, rotation):
    fig = pl.figure()
    ax = pl.subplot(111)
    ax.set_xlim(-3.5,3.5)
    ax.set_ylim(-3,3)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.set_aspect('equal')
    if (nsteps/4. < f < 3*nsteps/4.):
      Star(ax, u1 = 0.5, u2 = 0.25, behind = False)
    else:
      Star(ax, u1 = 0.5, u2 = 0.25, behind = True)
    
    Trail(ax, x, y, f)
    Planet(ax, x = x[f], y = y[f], rot = rot)

    fig.savefig('%03d.png' % f, bbox_inches = 'tight')
    pl.close()
  
  # Make gif
  subprocess.call(['convert', '-delay', '5', '-loop', '-1', '*.png', 'transit.gif'])