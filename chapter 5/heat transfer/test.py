import numpy as np
from scipy import fftpack

poly = np.array([[1, 3],
                 [3, 4]])
fpoly = fftpack.fft2(poly) 
ifpoly = fftpack.ifft2(fpoly)
for i, j in zip(ifpoly, fpoly):
  print(i, j)