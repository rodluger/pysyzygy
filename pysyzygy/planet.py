# -*- coding: utf-8 -*-
"""
Adapted from ``spheres.py`` by Brett Morris (bmorris3)

"""

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from PIL import Image, ImageOps
from . import proj
from . import PSZGPATH

class Planet(object):
    def __init__(self, albedo=1, resolution=1000, RpRs = 0.1, image_map='earth', night_alpha = 0.35):
        '''
        Parameters
        ----------

        albedo : float
            The fraction of the flux that is reflected off of the surface
            of the sphere

        Returns
        -------

        An instance of the Planet class
        '''
        self.r_planet = int(0.45*resolution)
        self.RpRs = RpRs
        self.dims = (resolution, resolution)
        self.albedo = albedo
        self.image_generated = False
        self.image_map = PSZGPATH + '/pysyzygy/maps/' + image_map + '.jpg'
        self.night_alpha = night_alpha

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

    def gen_image(self, xyz):
        '''
        Generate the image of the sphere at the given position. This method
        does most of the work.

        Parameters
        ----------
        xyz : tuple
            The position of the planet (x, y, z) relative to the star, in units
            of the stellar radius. TODO: Resolve some sign errors here.
            
        '''
        self.image = np.zeros(self.dims)
        
        # Position of the star. BUG: Why does a have to be so big?
        a = 1e6*self.RpRs*self.r_planet
        s_centroid = [-xyz[0]*a, xyz[1]*a, xyz[2]*a]
        
        # Observer must be in z direction, at infinity ( = one million! )
        self.observer = [0.5 * self.dims[0], 0.5 * self.dims[1], 1e6*a]  

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

    def plot_image(self, fig=None, ax=None, show=False, cmap='Greys_r', extent=None, long0=0.):
        '''
        Plot the image of the sphere. Return the figure and axis.
        '''
        if not self.image_generated:
            self.gen_image()

        colormap = cm.get_cmap(cmap)
        # Set off-planet background to transparent
        colormap.set_under(alpha=0)

        if ax is None:
            fig, ax = plt.subplots(1, 1, frameon=False)
            ax.patch.set_alpha(0)
            ax.set_aspect('equal')
            ax.set(xticks=[], yticks=[])
            ax.axis('off')
        
        if self.image_map is None:
            ax.imshow(self.image, origin='lower', cmap=colormap, vmin=0, vmax=1, extent=extent)
        else:
            # A black background
            bkg = -np.ones(self.dims)
            bkg[self.planet_disk] = 0
            ax.imshow(bkg, origin='lower', cmap=colormap, vmin=0, vmax=1, extent=extent)
              
            self.image[self.image < 0] = np.nan
            self.image = (1. - self.night_alpha)*self.image + self.night_alpha
            self.image[np.isnan(self.image)] = 0
            
            mask = Image.fromarray((self.image*255).astype('uint8'))
            im = Image.open(self.image_map)
            
            # Project onto a sphere
            sp = proj.SphericalProjection(im)
            im = sp.render(long0=long0)
            im = ImageOps.fit(im, mask.size)
            
            # Add the alpha mask
            im.putalpha(mask)
            ax.imshow(im, extent=extent)
        
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

        if fig is None:
            return ax
        else:
            return fig, ax