#!/usr/bin/env python

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# Normalization factor of the probability density
F = integrate.quad( lambda x: x**2, 0, 1)

# Cosine distribution, cos^2(theta)
def x2(x):
    return x**2/F[0]

# CDF of the x^2 distribution
def cdf(costh):
    return (1.-costh**3)/(3*F[0])

# Inverse CDF of x^2
def inv_cdf(c):
    return np.cbrt(1-3*c*F[0])

def findLineIntersectionPointWithPlane(p0, v, p1, n):
    t = np.dot((p1-p0), n)/np.dot(v, n)
    return p0 + v*t

def findTraversing(x, face):
    
    if x[0] >= face[0][0] and x[0] <= face[0][1] and \
       x[1] >= face[1][0] and x[1] <= face[1][1] and \
       x[2] >= face[2][0] and x[2] <= face[2][1]:
        return True
    
    return False

def genSolidAngle(nSamples):
    uniformCosth = np.random.default_rng().uniform(0, 1, nSamples)
    SolidAngle = np.array([
        [costh, phi]
        for costh, phi in zip(inv_cdf(uniformCosth), np.random.default_rng().uniform(0, 2*np.pi, nSamples))])
    return SolidAngle

def getSinFromCos(costh):
    return np.sqrt(1 - costh**2)

def genXYPosition(xyBound, zPosition, nSteps):
    iStep = np.linspace(0, xyBound, nSteps+1)
    XPosYMin = np.array([[x, -xyBound, zPosition] for x in iStep ])
    XNegYMin = np.array([[x, -xyBound, zPosition] for x in -iStep])
    XPosYMax = np.array([[x, xyBound, zPosition]  for x in iStep ])
    XNegYMax = np.array([[x, xyBound, zPosition]  for x in -iStep])
    XMinYPos = np.array([[-xyBound, y, zPosition] for y in iStep ])
    XMinYNeg = np.array([[-xyBound, y, zPosition] for y in -iStep])
    XMaxYPos = np.array([[xyBound, y, zPosition]  for y in iStep ])
    XMaxYNeg = np.array([[xyBound, y, zPosition]  for y in -iStep])
    XYPosition = np.concatenate(( XPosYMin, XNegYMin, XPosYMax, XNegYMin, XMinYPos, XMinYNeg, XMaxYPos, XMaxYNeg))
    unique_check = np.unique(XYPosition, axis = 1)
    if len(XYPosition) != len(unique_check):
        print(f'{len(unique_check)} unique positions out of {len(XYPosition)} positions!')
    # print(f'{nXYPositions} points of cosmic rays being checked')
    
    return XYPosition

if __name__ == "__main__":

    # Unit: cm
    xHalf = 25.
    yHalf = 30.
    zFull = 60.

    # Define the planes corresponding to the 6 faces by a point and a normal vector
    p1s = np.array([[0., 0., 0.], [0., -yHalf,0.], [0., yHalf, 0.], [ -xHalf, 0., 0.], [xHalf, 0., 0.], [0., 0., zFull]])
    norms = np.array([[0., 0., -1,], [0., -1., 0.], [0., 1., 0.], [-1., 0., 0.], [1., 0., 0.], [0., 0., 1.]])
    # Define the borders of the 6 faces
    eps = 1e-6
    faces = np.array([[[-xHalf, xHalf], [-yHalf, yHalf], [-eps, eps]], 
                      [[-xHalf, xHalf], [-yHalf-eps, -yHalf+eps], [0., zFull]],
                      [[-xHalf, xHalf], [yHalf-eps, yHalf+eps], [0., zFull]],
                      [[-xHalf-eps, -xHalf+eps], [-yHalf, yHalf], [0., zFull]],
                      [[xHalf-eps, xHalf+eps], [-yHalf, yHalf], [0., zFull]],
                      [[-xHalf, xHalf], [-yHalf, yHalf], [zFull-eps, zFull+eps]]])

    # Unit: cm
    # Define where to sample the cosmic muon flux
    zPosition = -800
    # zPosition = -10
    nSolidAngles = 10000

    outfile = 'nCosmicsInTPC3.npy'
    # xyBounds = np.array([10., 20., 30., 40., 50., 100., 500., 1000., 5000. ])
    xyBounds = np.array([7500, 8000., 9000., 10000.])
    nBounds = len(xyBounds)
    nTraverse = np.zeros(nBounds)
    nCounter = np.zeros(nBounds)

    for ixyBound in range(nBounds):
    
        xyBound = xyBounds[ixyBound]
        XYPosition = genXYPosition(xyBound, zPosition, round(xyBound))

        for p0 in XYPosition:
        
            SolidAngles = genSolidAngle(nSolidAngles)
    
            for costh, phi in SolidAngles:
    
                # print(f'Cosmic ray {nCounter}')
    
                sinth = getSinFromCos(costh)
                v = np.array([sinth*np.cos(phi), sinth*np.sin(phi), costh])
                # print(f'v = {v}')
    
                for p1, norm, face in zip(p1s, norms, faces):
                    x = findLineIntersectionPointWithPlane(p0, v, p1, norm)
                    # print(f'x = {x}')
                    doTraverse = findTraversing(x, face)
                    if (doTraverse):
                        # print(f'It traverses! v = ({v})')
                        nTraverse[ixyBound] += 1
                        break

                nCounter[ixyBound] += 1

    result = np.array([ [xyBound, nT, nC, float(nT)/float(nC)] for xyBound, nT, nC in zip(xyBounds, nTraverse, nCounter)])

    with open(outfile, 'wb') as f:
        np.save(f, result)
