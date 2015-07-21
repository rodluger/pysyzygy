import kplr
import numpy as np
import matplotlib.pyplot as pl
from pysyzygy import Transit

class Chunk(object):
  '''
  A dummy class.
  
  '''
  def __init__(self):
    self.time = []
    self.fsum = []

def GetKoi(koi_number):
    '''
    A wrapper around kplr.API().koi(), with the additional command
    ``select=*`` to query **all** columns in the Exoplanet database.
  
    :param float koi_number: The full KOI number of the target (``XX.XX``)
  
    :returns: A ``kplr`` ``koi`` object
  
    '''
    client = kplr.API()
    kois = client.kois(where="kepoi_name+like+'K{0:08.2f}'"
                     .format(float(koi_number)), select="*")
    if not len(kois):
      raise ValueError("No KOI found with the number: '{0}'".format(koi_number))
    return kois[0]

planet = GetKoi(17.01)
per = planet.koi_period
t0 = planet.koi_time0bk
tdur = planet.koi_duration/24.
tpf = planet.get_target_pixel_files(short_cadence = False)
data = [Chunk() for i in range(18)]

for fnum in range(len(tpf)):
  with tpf[fnum].open(clobber = False) as f:  
    quarter = f[0].header['QUARTER']
    aperture = f[2].data
    idx = np.where(aperture & 2)
    qdata = f[1].data
    crowding = f[1].header['CROWDSAP']

  data[quarter].time = np.array(qdata.field('TIME'), dtype='float64')
  fpix = np.array([f[idx] for f in qdata.field('FLUX')], dtype='float64')
  data[quarter].fsum = np.sum(fpix, axis = 1)

# Expected?
trn = Transit(per = 5., t0 = 0., q1 = 0.5, 
              q2 = 0.5, MpMs = 0.01, RpRs = 0.2, 
              bcirc = 0.2, ecw = 0., 
              esw = 0., rhos = 1.4)
time = np.linspace(0, 50, 1000)
tmod = trn(time, 'binned')
fig = pl.figure()
fig.set_size_inches(12,4)
pl.plot(time, tmod, 'b.')
pl.ylabel('Flux (counts)', fontsize = 24)
pl.xlabel('Time (days)', fontsize = 24)
pl.title('The ideal lightcurve', fontsize = 28)
fig.savefig('output/ideal_lc.png', bbox_inches = 'tight')

# Actual
fig = pl.figure()
fig.set_size_inches(12,4)
for chunk in data:
  pl.plot(chunk.time, chunk.fsum, 'b.')
pl.ylim(64000,72000)
pl.ylabel('Flux (counts)', fontsize = 24)
pl.xlabel('Time (days)', fontsize = 24)
pl.title('What it actually looks like', fontsize = 28)
fig.savefig('output/real_lc_17.01.png', bbox_inches = 'tight')