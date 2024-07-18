#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 14:36:11 2022

@author: sangeeth
"""

from pythonImports import * 
from problemConstants import *

'''Function to rotate a vector'''
'''Ideally, it can be used to return a vector that is parallel & opposite 
to the thrust direction = acceleration direction of the thrust''' 
def rotate(x, y, angle):
    x_new = math.cos(angle)*x - math.sin(angle)*y 
    y_new = math.sin(angle)*x + math.cos(angle)*y;
    return(x_new, y_new)

'''Function to compute angle between 2 vectors'''
def angle2D(ux, uy, vx, vy):
    ptot2 = (ux*ux + uy*uy)*(vx*vx + vy*vy)
    if(ptot2 <= 0):
            return 0.0
    a = math.acos(min(1.0, max(-1.0, (ux*vx + uy*vy)/math.sqrt(ptot2))))
    s = ux*vy - uy*vx
    if (s < 0):
            a = -a
    return a
    
'''Function to plot the location of the ship'''
def plotLocationOfShip(x,y,u,v):
    fig, ax = plt.subplots()
    ax.axis("equal")
    circle = plt.Circle((0,0),radiusOfMoon,color='black', fill=False)
    circle2 = plt.Circle((0,0), radiusOfMoon+ hmin,color='red', fill=False)
    ax.set(xlim=(-2.1*radiusOfMoon+hmin,2.1*radiusOfMoon+hmin), ylim=(-2.1*radiusOfMoon+hmin,2.1*radiusOfMoon+hmin))
    ax.add_artist(circle)
    ax.add_artist(circle2)
    plt.plot(x,y, marker="x",markersize=20)
    '''To plot the elipse, we assume a radius vector that is oriented +y and 
    rotate it by the initial ship angle'''
    radiusVec_x, radiusVec_y = rotate(0.,semiMajorAxisOfTransfer-radiusOfMoon,initialShipAngle)
    ells = patches.Ellipse(xy=(radiusVec_x,radiusVec_y), width=2*semiMajorAxisOfTransfer, height=2*semiMinorAxisOfTransfer, fill=False, color='purple')
    ax.add_patch(ells)
    '''To plot the velocity vectors'''
    ax.quiver(x,y,u,v)


    
def plotVectors(u1,u2,v1,v2,w1,w2):
    uresultant = math.sqrt(u1**2+u2**2)
    vresultant = math.sqrt(v1**2+v2**2)
    wresultant = math.sqrt(w1**2+w2**2)

    alist = []
    alist.append([u1/uresultant,u2/uresultant])
    alist.append([v1/vresultant,v2/vresultant])
    alist.append([w1/wresultant,w2/wresultant])
    V = np.array(alist)
    origin = np.array([[0,0,0],[0,0,0] ]  ) # origin point
    plt.quiver(*origin, V[:,0], V[:,1], color=['r','b', 'c'], scale=10)
    plt.show()
