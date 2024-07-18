#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 17:52:41 2022

@author: sangeeth
"""

from pythonImports import *


'''Define constants'''
radiusOfMoon = 1737.4e3
massOfMoon = 0.07346e24
gravitationalConstant = 6.6743e-11
gravitationalParameter = gravitationalConstant*massOfMoon
hmin = 200e3
radiusOfInitialOrbit = radiusOfMoon+hmin
landingVel = 2.
maxThrust = 5e3
minThrust = 1e3
Isp = 450
gravityOfEarth = 9.81
massPropellant = 1000.0
massDry = 500. 
heightOfLowerHohmannOrbit = 40000
radiusAtPoweredDescentInitiation = radiusOfMoon + heightOfLowerHohmannOrbit
radiansToDegreeFactor = 180/math.pi

'''Define user inputs'''
initialShipAngle = math.pi
thrustAngle = 0.
x0 = 0.
y0 = radiusOfMoon + hmin
U_comp_vel = math.sqrt(gravitationalConstant*massOfMoon/(radiusOfInitialOrbit)) 
V_comp_vel = 0.
a_gravity  = gravitationalConstant*massOfMoon/(radiusOfInitialOrbit*radiusOfInitialOrbit)


'''Deorbit parameters'''
thrustForceForDeorbit = 1e3
thrustAngleForDeorbit = math.pi/2.
dTForDeorbit = 14 #CONTROL THIS

'''Hohmann Orbit parameters'''
nIterForHohmann = 360000
dt_elipse = 0.01
r1 = radiusOfInitialOrbit
r2 = radiusAtPoweredDescentInitiation
semiMajorAxisOfTransfer = 0.5*(r1+r2)
semiMinorAxisOfTransfer = math.sqrt(r1*r2)
eccentricityOfTransfer = math.sqrt(1. - (semiMinorAxisOfTransfer**2/semiMajorAxisOfTransfer**2))
semiLactusRectum = semiMinorAxisOfTransfer**2/semiMajorAxisOfTransfer
transferTime = math.pi*math.sqrt(pow(semiMajorAxisOfTransfer,3)/gravitationalParameter)


'''PDI phase 1 parameters'''
nIterForPDIPhase1 = 440
fullThrottleForPDIPhase1 = 1
thrustAngleForPDIPhase1 = math.pi/2.
dTForPDIPhase1 = 1
deltaVForPDIPhase1 = 3.3   #CONTROL THIS


'''PDI phase 2 parameters'''
nIterForPDIPhase2 = 38
dTForPDIPhase2 = 1

#If full thrust
fullThrottleForPDIPhase2 = 0
thrustAngleForPDIPhase2 = math.pi/2.

#If constant velocity change
constantdeltaVForPDIPhase2 = 0
deltaVForPDIPhase2 = 4.2   #CONTROL THIS

#If constant altitude change
constantdeltaAltForPDIPhase2 = 1
deltaAltForPDIPhase2 = 46
descentFactor = 10.

'''TouchDown phase parameters'''
nIterForTouchDown = 10
dTForTouchDown = 1

thrustAngleForTouchDown = math.pi

#If constant altitude change
constantdeltaAltForTouchDown = 1
deltaAltForTouchDown = 10
descentFactorForTouchDown = 5

