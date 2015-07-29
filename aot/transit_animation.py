import pysyzygy as ps

# Main animation
ps.AnimateImage(per = 0.5, RpRs = 0.45, ecc = 0, rhos = 1.0, w = 0,
               b = 0.3, u1 = 1., u2 = 0., delay = 0,
               bkgimage = 'stars', nsteps = 250,
               image_map = 'earth', size_inches = (12,9),
               lightcurve = True, delete = False)