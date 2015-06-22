# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 12:05:15 2014 by Brett Morris (bmorris3)

To run a quick example, try:
    from pysyzygy.spheres import Sphere
    s = Sphere()
    for phase in range(0,6):
        s.gen_image(phase)
        s.plot_image(show=True)

or

    from pysyzygy.spheres import Sphere
    s = Sphere()
    s.grid()

"""
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from PIL import Image, ImageOps, ImageDraw

class Sphere(object):
    def __init__(self, albedo=1, resolution=1000, radius=450, orbitplane='xz'):
        '''
        Parameters
        ----------

        albedo : float
            The fraction of the flux that is reflected off of the surface
            of the sphere

        Returns
        -------

        An instance of the sphere class
        '''
        self.r_planet = radius
        self.dims = (resolution, resolution)
        self.plane = orbitplane
        self.albedo = albedo
        self.image_generated = False

    def vector_length(self, vector):
        '''
        the square root of the sum of the squares
        '''
        return np.sqrt(np.sum(np.array(vector) ** 2))

    def normalize(self, vector):
        '''
        Normalize a vector by the vector length
        '''
        return np.array(vector) / self.vector_length(vector)

    def pr(self, a, b):
        '''
        The projection of vector `b` onto `a` is the vector 
    
        .. math::    
        
            \mathrm{pr}_a (b) = (a \cdot b / ||a||^2) a
    
        '''
        return a * np.dot(a, b) / self.normalize(a) ** 2

    def gen_image(self, phaseangle=0):
        '''
        Generate the image of the sphere in the given phase angle. This method
        does most of the work.

        Parameters
        ----------
        phaseangle : float
            The phase angle between the star and the observer. Choose zero for
            a "full" phase sphere, pi for a "new" phase sphere.
        '''
        self.phase_angle = phaseangle + np.pi / 2
        self.image = np.zeros(self.dims)

        a = 1e6 * self.r_planet  # semi-major axis

        if self.plane == 'xy':
            s_centroid = [np.cos(self.phase_angle) * a,
                          np.sin(self.phase_angle) * a, 0]
        elif self.plane == 'yz':
            s_centroid = [0, np.cos(self.phase_angle) * a,
                          np.sin(self.phase_angle) * a]
        elif self.plane == 'xz':
            s_centroid = [np.cos(self.phase_angle) * a, 0,
                          np.sin(self.phase_angle) * a]

        self.observer = [0.5 * self.dims[0], 0.5 * self.dims[1],
                         a]  # Observer must be in z direction

        x, y = np.meshgrid(np.arange(self.dims[0]), np.arange(self.dims[1]))

        p_centroid = [0.5 * self.dims[0], 0.5 * self.dims[1], 0.0]
        with np.errstate(invalid='ignore'):
            z = np.sqrt(self.r_planet ** 2 - (x - p_centroid[0]) ** 2
                        - (y - p_centroid[1]) ** 2)

        z[np.isnan(z)] = 0  # Convert NANs to zeros

        planet_disk = ((x - p_centroid[0]) ** 2 + (
        y - p_centroid[1]) ** 2 <= self.r_planet ** 2)

        d_star = np.zeros(self.dims)  # Distance to the star
        d_star[planet_disk] += np.sqrt((x[planet_disk] - s_centroid[0]) ** 2 +
                                       (y[planet_disk] - s_centroid[1]) ** 2 +
                                       (z[planet_disk] - s_centroid[2]) ** 2)
        d_star_centroid = np.sqrt((s_centroid[0] - p_centroid[0]) ** 2 +
                                  (s_centroid[1] - p_centroid[1]) ** 2 +
                                  (s_centroid[2] - p_centroid[2]) ** 2)

        distance_observer = np.sqrt((x[planet_disk] - self.observer[0]) ** 2 +
                                    (y[planet_disk] - self.observer[1]) ** 2 +
                                    (z[planet_disk] - self.observer[2]) ** 2)
        dir_observer_x = (self.observer[0] - x[planet_disk]) / distance_observer
        dir_observer_y = (self.observer[1] - y[planet_disk]) / distance_observer
        dir_observer_z = (self.observer[2] - z[planet_disk]) / distance_observer

        self.image[(d_star <= d_star_centroid) * (d_star != 0)] += self.albedo

        dir_star_x = (s_centroid[0] - x[planet_disk]) / d_star[planet_disk]
        dir_star_y = (s_centroid[1] - y[planet_disk]) / d_star[planet_disk]
        dir_star_z = (s_centroid[2] - z[planet_disk]) / d_star[planet_disk]

        # Compute surface normals to the sphere in each dimension
        surface_norms_x = (x[planet_disk] - p_centroid[0]) / self.r_planet
        surface_norms_y = (y[planet_disk] - p_centroid[1]) / self.r_planet
        surface_norms_z = (z[planet_disk] - p_centroid[2]) / self.r_planet

        # Compute Lambertian reflections
        LdotN = (dir_star_x * surface_norms_x + dir_star_y * surface_norms_y +
                dir_star_z * surface_norms_z)
        self.image[planet_disk] *= LdotN

        self.planet_disk = planet_disk
        LdotNdotObs = (dir_star_x * dir_observer_x * surface_norms_x +
                      dir_star_y * dir_observer_y * surface_norms_y +
                      dir_star_z * dir_observer_z * surface_norms_z)
        LdotNdotObs[LdotNdotObs < 0] = 0

        self.image[~planet_disk] = -1

        self.image_generated = True
        return self.image

    def plot_image(self, fig=None, ax=None, show=False, cmap='Greys_r', image=None):
        '''
        Plot the image of the sphere. Return the figure and axis.
        '''
        if not self.image_generated:
            self.gen_image()

        colormap = cm.get_cmap(cmap)
        # Set off-planet background to transparent
        colormap.set_under(alpha=0)

        if (fig is None) or (ax is None):
            fig, ax = plt.subplots(1, 1, frameon=False)
            ax.patch.set_alpha(0)
            ax.set_aspect('equal')
            ax.set(xticks=[], yticks=[])
            ax.axis('off')
        
        if image is None:
          ax.imshow(self.image, origin='lower', cmap=colormap, vmin=0, vmax=1)
        else:
        
            # A black background with a mask
            sz = self.dims[0]
            dr = sz/2. - self.r_planet
            bkg = Image.new('RGBA', (sz, sz), 'white')
            bmk = Image.new('L', (sz, sz), 0)
  
            # The Lambertian mask
            self.image[self.image < 0] = 0
            mask = Image.fromarray((self.image*255).astype('uint8'))
            for xi in np.arange(dr, sz-dr):
              ym = np.sqrt((sz/2. - dr)**2 - (xi -  sz/2.)**2)
              for yi in np.arange( sz/2. - ym, sz/2. + ym):
                bkg.putpixel((int(xi),int(yi)), 0)
                bmk.putpixel((int(xi),int(yi)), 255)
  
            # Plot the background
            bkg.putalpha(bmk)
            ax.imshow(bkg, origin='lower')
  
            # Overplot the alpha-corrected image
            im = Image.open(image)
            im = ImageOps.fit(im, mask.size)
            im.putalpha(mask)
            ax.imshow(im)
        
        
        if show:
            def format_coord(x, y):
                '''
                Function to give data value on mouse over with imshow.
                '''
                try:
                    return ('x={0:.1f}, y={1:.1f}, '
                            'Flux={2:.2f}').format(x, y, self.image[y, x])
                except:
                    return 'x={0:d}, y={1:d}'.format(x, y)

            ax.format_coord = format_coord
            plt.show()

        return fig, ax
    
    def grid(self, n = 3, cmap='Greys_r', image=None):
        '''
        Plot a grid of all the phases of size n x n.
        '''
        fig, axes = plt.subplots(n, n, frameon=False)
        try:
          axes = axes.flatten()
        except AttributeError:
          axes = [axes]
        phases = np.linspace(0,2*np.pi,n**2)
        for ax, phase in zip(axes, phases):
          ax.patch.set_alpha(0)
          ax.set_aspect('equal')
          ax.set(xticks=[], yticks=[])
          ax.axis('off')
          self.gen_image(phase)
          self.plot_image(fig=fig, ax=ax, show=False, cmap=cmap, image=image)
        plt.show()
        
        return fig, axes

if __name__ == '__main__':
  s = Sphere()
  s.grid(n = 3) #, image='/Users/rod/Desktop/transit/earth.jpg')