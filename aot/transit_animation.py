import pysyzygy as ps

# Main animation
ps.AnimateImage(per = 0.5, RpRs = 0.5, ecc = 0, rhos = 1.0,
               b = 0.8, u1 = 1., u2 = 0., delay = 20,
               bkgimage = 'stars', nsteps = 10,
               image_map = 'earth', size_inches = (12,9),
               lightcurve = True)