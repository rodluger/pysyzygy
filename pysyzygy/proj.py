from PIL import Image
import numpy as np
import matplotlib.pyplot as pl

class SphericalProjection(object):
    """
    BASED ON http://svn.navi.cx/misc/trunk/pycaptcha/Captcha/Visual/Distortions.py
    
    """
    

    def __init__(self, image, filtering = Image.BICUBIC, res = 10):
        self.filtering = filtering
        self.image = image
        self.res = res

    def c2s(self, x, z, sz):
        """
        Converts cartesian (x,y) to spherical (theta, phi) coordinates.
        Conventions as in the image at http://keisan.casio.com/exec/system/1359533867
  
        """
  
        # Normalize and shift origin
        xn = 0.5 - x/(1.*sz[0])
        zn = 0.5 - z/(1.*sz[1])
        rn = 0.5
        with np.errstate(invalid='ignore'):
          yn = np.sqrt(rn**2 - xn**2 - zn**2)
  
        # Convert to spherical
        phi = np.arctan2(np.sqrt(xn**2 + yn**2), zn)
        theta = np.arctan2(yn,xn)
  
        # Change range to [0, np.pi]
        try:
            phi[phi < 0] += 2*np.pi
        except IndexError:
            if phi < 0: phi += 2*np.pi
        try:
            theta[theta < 0] += 2*np.pi
        except IndexError:
            if theta < 0: theta += 2*np.pi
  
        # Normalize to source image
        phi *= (sz[1]/np.pi)
        theta *= (sz[0]/np.pi)
  
        return theta, phi

    def spherical(self, image):
        return lambda x, z: self.c2s(x, z, image.size)

    def render(self, long0 = 0.):
        
        phase = (long0 + 180.)/360.
        dx = int(self.image.size[0]*(phase - 0.25))
        
        szx = self.image.size[0]
        szy = self.image.size[1]
        
        # Tile and crop
        tmp = Image.new("RGBA", (3*szx, szy))
        tmp.paste(self.image, (0,0)) 
        tmp.paste(self.image, (szx,0))  
        tmp.paste(self.image, (2*szx,0))
        image = tmp.crop((szx + dx, 0, szx + szx/2 + dx, szy))
    
        r = self.res
        xPoints = image.size[0] / r + 2
        yPoints = image.size[1] / r + 2
        f = self.spherical(image)

        # Create a list of arrays with transformed points
        xRows = []
        yRows = []
        for j in xrange(yPoints):
            xRow = []
            yRow = []
            for i in xrange(xPoints):
                x, y = f(i*r, j*r)
                
                # Clamp the edges so we don't get black undefined areas
                x = max(0, min(image.size[0]-1, x))
                y = max(0, min(image.size[1]-1, y))
                
                xRow.append(x)
                yRow.append(y)
            xRows.append(xRow)
            yRows.append(yRow)

        # Create the mesh list, with a transformation for
        # each square between points on the grid
        mesh = []
        for j in xrange(yPoints-1):
            for i in xrange(xPoints-1):
                mesh.append((
                    # Destination rectangle
                    (i*r, j*r,
                     (i+1)*r, (j+1)*r),
                    # Source quadrilateral
                    (xRows[j  ][i  ], yRows[j  ][i  ],
                     xRows[j+1][i  ], yRows[j+1][i  ],
                     xRows[j+1][i+1], yRows[j+1][i+1],
                     xRows[j  ][i+1], yRows[j  ][i+1]),
                    ))

        return image.transform((image.size[0], image.size[1]), Image.MESH, mesh, self.filtering)
        
if __name__ == '__main__': 
    # An example            
    im = Image.open('maps/earth.jpg')
    sp = SphericalProjection(im)
    im = sp.render(long0=0)
    ax = pl.subplot(111)
    ax.imshow(im)
    pl.show()