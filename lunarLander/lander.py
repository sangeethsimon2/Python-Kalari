#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 14:35:05 2022

@author: sangeeth
"""

from util import * 

'''Driver code'''
'''Initiate total time'''
totalTime = 0.

'''Initialize the ship to a position in the orbit''' 
shipAngle = initialShipAngle

'''Compute the position vector components'''
x1, y1 = rotate(x0, y0, shipAngle)
magOfPositionVector = math.sqrt(x1**2 + y1**2)

print("The ship is now at location %f, %f " %(x1,y1)) 

'''Compute the current altitude of the ship'''
currentAltitude = math.sqrt(x1**2+y1**2) - radiusOfMoon

print("The current altitude of the ship is %f" %currentAltitude)


'''Compute the velocity components of the ship'''
U_comp_vel, V_comp_vel = rotate(U_comp_vel, V_comp_vel, shipAngle)
resultantVelocity = math.sqrt(U_comp_vel*U_comp_vel+V_comp_vel*V_comp_vel)

print("The ship has a velocity of %f m/s" %resultantVelocity)

'''PLot the location and velocity'''
plotLocationOfShip(x1,y1,U_comp_vel, V_comp_vel)

'''Resolve the gravitational acceleration into its x,y components'''
ax_gravity = -a_gravity*x1/magOfPositionVector;
ay_gravity = -a_gravity*y1/magOfPositionVector;

#######################################################################
'''Decelerate phase: Enter Hohmann transfer orbit'''
#########################################################################

print("Phase 1: Decelerate the ship and Enter Hohmann transfer orbit \n")

'''Compute the del v required to exit orbit r1 and enter orbit r2'''
''' r1 = radiusOfMoon+hmin ; r2 = radiusOfMoon'''
print(" The properties of the elipse are, semimajoraxis=%f \n, semiminoraxis=%f \n, eccent=%f \n, semilactusrectum=%f \n" %(semiMajorAxisOfTransfer,semiMinorAxisOfTransfer, eccentricityOfTransfer,semiLactusRectum))
print("The transfer time required in Hohmann orbit is %f" %transferTime)
del_v = math.sqrt(gravitationalParameter/(r1))*( 1. - math.sqrt(r2/semiMajorAxisOfTransfer)) 
print("The velocity change required to deOrbit is %f" %del_v)

'''Compute the final mass of the ship after the delta V maneuvre'''
massAfterDeltaV = (massDry+massPropellant)/math.exp(del_v/(Isp*gravityOfEarth))

'''Compute the mass of the fuel required to achieve this delta v'''
massOfFuelForDeltaV = (massDry+massPropellant) - massAfterDeltaV
assert(massOfFuelForDeltaV>0)

'''Compute mass flow from massOfFuelForDeltaV'''
massFlowOfPropellant = abs(massOfFuelForDeltaV)/dTForDeorbit
assert(massFlowOfPropellant>0.)


'''Propellant mass remaining after Deorbit thrust'''
massPropellant-=massOfFuelForDeltaV
print("The propellant remaining after Deorbit thrust is %f" %massPropellant)


'''Set thrust angle to be parallel to surface and opposite to vel dir'''
thrustAngle = thrustAngleForDeorbit
'''find thrust from the massflowofpropellant'''
thrustForce = Isp*gravityOfEarth*massFlowOfPropellant

assert(thrustForce<maxThrust)
assert(thrustForce>minThrust)

print(" Thrust force required to effect the deltV change is %f" %thrustForce)
print(" The engine is switched on for %f s (USER INPUT) to achieve this thrust " %dTForDeorbit)
'''Determine the direction of the velocity vector'''
'''-ThrustAngle since this method rotates anti-clockwise which is negative by definition'''
thx, thy = rotate(x1/magOfPositionVector, y1/magOfPositionVector, -thrustAngle);


'''New resultant velocity is'''
resultantVelocity -=del_v
'''Compute velocity components after the deorbit manuever'''
U_comp_vel_afterDeorbit = thx*resultantVelocity
V_comp_vel_afterDeorbit = thy*resultantVelocity


print("The new velocity comp acquired by the ship are %f %f "%(U_comp_vel_afterDeorbit, V_comp_vel_afterDeorbit))
     
'''Compute the resultant acceleration acquired by the ship due to the thrust'''
a_thrust = thrustForce/(massDry+massPropellant)


'''Determine the direction of the acceleration vector'''
thx, thy = rotate(x1/magOfPositionVector, y1/magOfPositionVector, thrustAngle);
 
'''Resolve the resultant acceleration due to thrust into its components'''
ax_thrust = thx*a_thrust
ay_thrust = thy*a_thrust

'''Compute net x-dir and y-dir acceleration'''
ax_net = ax_thrust + ax_gravity
ay_net = ay_thrust + ay_gravity

'''Compute the instantaneous change in position experienced by the ship due 
to the acceleration'''
x2 = x1 + U_comp_vel*dTForDeorbit + 0.5*ax_net*dTForDeorbit*dTForDeorbit
y2 = y1 + V_comp_vel*dTForDeorbit + 0.5*ay_net*dTForDeorbit*dTForDeorbit


'''Compute the magnitude of the position vector'''
magOfPositionVector = math.sqrt(x2*x2+ y2*y2)

'''Compute the current altitude of the ship'''
currentAltitude = math.sqrt(x2*x2+y2*y2) - radiusOfMoon

print("The current altitude of the ship is %f" %currentAltitude)


'''Compute the acceleration due to gravity at the new position'''
a_gravity  = gravitationalParameter/(magOfPositionVector**2)

'''Compute the components of this acceleration due to gravity'''
ax_gravity = -a_gravity*x2/magOfPositionVector;
ay_gravity = -a_gravity*y2/magOfPositionVector;

'''Compute net acceleration'''
ax_net = ax_gravity
ay_net = ay_gravity

# plotLocationOfShip(x2, y2, ax_net, ay_net)

'''Store the initial parameters before iteration'''
xInitial = x2
yInitial = y2
uInitial = U_comp_vel_afterDeorbit
vInitial = V_comp_vel_afterDeorbit


'''Set clock to 0'''
time = 0.

print("Initiating Hohmann transfer\n")
print("\n ")
    
for i in range(0,nIterForHohmann):
       
       time+=dt_elipse
       '''Compute change in velocity in dt_elipse time due to this acceleration'''
       U_comp_vel_i = uInitial + ax_net*dt_elipse        
       V_comp_vel_i = vInitial + ay_net*dt_elipse
       
       '''Compute the change in position due to this acceleration'''
       x_i = xInitial + U_comp_vel_i*dt_elipse
       y_i = yInitial + V_comp_vel_i*dt_elipse 

       '''Compute the magnitude of the position vector'''
       magOfPositionVector = math.sqrt(x_i*x_i+ y_i*y_i)

       '''Compute the current altitude of the ship'''
       currentAltitude = math.sqrt(x_i*x_i+y_i*y_i) - radiusOfMoon       
       assert(currentAltitude>0.)
       
       '''Compute the acceleration due to gravity at the new position'''
       a_gravity  = gravitationalParameter/(magOfPositionVector**2)
       '''Compute the components of this acceleration due to gravity'''
       ax_gravity = -a_gravity*x_i/magOfPositionVector;
       ay_gravity = -a_gravity*y_i/magOfPositionVector;

       ax_net = ax_gravity
       ay_net = ay_gravity
       
       xInitial = x_i
       yInitial = y_i
       uInitial = U_comp_vel_i
       vInitial = V_comp_vel_i
       
print("The Hohmann transfer is complete\n")

'''Store time taken'''
totalTime+= time

'''Store the final positions and velocities when reaching the altitude at which we 
want to initiate powered descent'''
       
x3 = x_i
y3 = y_i

'''Compute the magnitude of the position vector'''
magOfPositionVector = math.sqrt(x3*x3+ y3*y3)

'''Store current altitude'''
currentAltitudeAtPDI = currentAltitude

'''Store current velocity vectors'''
U_comp_vel_atPDI = U_comp_vel_i
V_comp_vel_atPDI = V_comp_vel_i
resultantVelocity = math.sqrt(U_comp_vel_atPDI**2+V_comp_vel_atPDI**2)


'''Plot ship location''' 
plotLocationOfShip(x3,y3,U_comp_vel_atPDI, V_comp_vel_atPDI)

print("The current altitude of the ship is %f" %currentAltitudeAtPDI)

print(" The position vector is %f %f" %(x3, y3))            

print(" The vel comp are %f %f and the resultant vel is %f " %(U_comp_vel_atPDI, V_comp_vel_atPDI, resultantVelocity))
       
print("The time taken to reach this altitude is %f s" %totalTime)       

print("\n \n")
       
print("Initiating powered descent Phase 1\n")


################################################################
'''Powered Descent Phase 1 - Breaking phase'''
#################################################################


'''Initialize time'''
time = 0
'''Set initial position as the last known position of the ship after 
the end of Hohmann transfer'''
xInitial = x3
yInitial = y3
'''Compute the magnitude of the position vector'''
magOfPositionVector = math.sqrt(xInitial**2+ yInitial**2)
'''Compute the current altitude of the ship'''
currentAltitude = math.sqrt(xInitial**2+yInitial**2) - radiusOfMoon

'''Set initial velocity as the last known velocity of the ship after 
the end of Hohmann transfer'''
uInitial = U_comp_vel_atPDI
vInitial = V_comp_vel_atPDI
resultantVelocity = math.sqrt(uInitial**2+vInitial**2)

'''Compute the acceleration due to gravity at the new position'''
a_gravity  = gravitationalParameter/(magOfPositionVector**2)
'''Compute the components of this acceleration due to gravity'''
ax_gravity = -a_gravity*x3/magOfPositionVector;
ay_gravity = -a_gravity*y3/magOfPositionVector;


if(fullThrottleForPDIPhase1):
    print("PDIPhase1 in fullThrottle retrograde\n")
    for i in range(0,nIterForPDIPhase1):
        '''Increment time'''
        time+=dTForPDIPhase1
        '''Compute the angle between the current vel vect and the pos vec'''
        Angle = angle2D(uInitial,vInitial,xInitial/magOfPositionVector,yInitial/magOfPositionVector)        
        
        '''Set thrust angle to be in the direction of the velocity'''
        thrustAngle = math.pi - Angle
        
        '''Set thrust force to be the max that the engine can provide'''
        thrustForce = maxThrust
        
        '''Compute the initial total mass'''
        m0 = massDry + massPropellant  
         
        '''Find mass flow of propellant needed to sustain the thrust'''
        massFlowOfPropellantForPDIPhase1 = thrustForce/(Isp*gravityOfEarth)
        
        '''Compute mass flow from massOfFuelForDeltaVForPDIPhase1'''
        massOfFuelForDeltaVForPDIPhase1 = massFlowOfPropellantForPDIPhase1*dTForPDIPhase1
        assert(massOfFuelForDeltaVForPDIPhase1>0)
        '''Propellant mass remaining after dV change in velocity is'''
        massPropellant-=massOfFuelForDeltaVForPDIPhase1
        
        '''COmpute the final mass after reduction in fuel'''
        mf = massDry + massPropellant

        '''Compute the deltaV achieved by this propellant flow'''
        deltaV = Isp*gravityOfEarth*math.log(m0/mf)
    
        '''Rotate the position vector to the direction of the velocity'''        
        Thx, Thy = rotate(xInitial/magOfPositionVector, yInitial/magOfPositionVector, -Angle)
          
        '''New resultant velocity is'''
        resultantVelocity -=deltaV
        
        '''Compute velocity components'''
        u_ii = Thx*resultantVelocity
        v_ii = Thy*resultantVelocity

        
        '''Compute the resultant acceleration acquired by the ship due to the thrust'''
        a_thrust = thrustForce/(massDry+massPropellant)
        
        '''Determine the direction of the acceleration vector'''
        thx, thy = rotate(xInitial/magOfPositionVector, yInitial/magOfPositionVector, thrustAngle);
        
        '''Plot the acceleration, vel and position vec '''
        # plotVectors(thx, thy, uInitial, vInitial, xInitial, yInitial)
        
        '''Resolve the resultant acceleration due to thrust into its components'''
        ax_thrust = thx*a_thrust
        ay_thrust = thy*a_thrust
        
        '''Compute the gravitational acceleration at this point'''
        a_gravity  = gravitationalParameter/(magOfPositionVector**2)
    
        '''Compute the components of this acceleration due to gravity'''
        ax_gravity = -a_gravity*xInitial/magOfPositionVector;
        ay_gravity = -a_gravity*yInitial/magOfPositionVector;
    
        
        '''Compute net x-dir and y-dir acceleration'''
        ax_net = ax_thrust + ax_gravity
        ay_net = ay_thrust + ay_gravity
        
        
        '''Compute change in velocity in dt_elipse time due to this acceleration'''
        u_i = uInitial + ax_net*dTForPDIPhase1       
        v_i = vInitial + ay_net*dTForPDIPhase1
        
        '''Compute the change in position due to this acceleration'''
        x_i = xInitial + u_i*dTForPDIPhase1
        y_i = yInitial + v_i*dTForPDIPhase1 
        
        '''Compute the current altitude of the ship'''
        currentAltitude = math.sqrt(x_i*x_i+y_i*y_i) - radiusOfMoon       
        assert(currentAltitude>0.)
        
        
        '''Set'''
        xInitial = x_i
        yInitial = y_i
        '''Compute the magnitude of the position vector'''
        magOfPositionVector = math.sqrt(xInitial**2+ yInitial**2)
    
        uInitial = u_i
        vInitial = v_i
        resultantVelocity = math.sqrt(uInitial**2+vInitial**2)
        
        # if((i%1000)==0):
             # plotLocationOfShip(xInitial,yInitial,ax_net, ay_net)
             # plotVelOfShip(xInitial, yInitial, uInitial, vInitial)  
        
else:
    for i in range(0,nIterForPDIPhase1):
        '''Increment time'''
        time+=dTForPDIPhase1
        
        '''Decrease the resultant velocity by deltaV value which is predefined'''
        '''Compute the final mass of the ship after the delta V maneuvre'''
        massAfterDeltaVForPDIPhase1 = (massDry+massPropellant)/math.exp(deltaVForPDIPhase1/(Isp*gravityOfEarth))
    
        '''Compute the mass of the fuel required to achieve this delta v'''
        massOfFuelForDeltaVForPDIPhase1 = (massDry+massPropellant)-massAfterDeltaVForPDIPhase1
        assert(massOfFuelForDeltaVForPDIPhase1>0)
    
        '''Compute mass flow from massOfFuelForDeltaVForPDIPhase1'''
        massFlowOfPropellantForPDIPhase1 = abs(massOfFuelForDeltaVForPDIPhase1)/dTForPDIPhase1
        assert(massFlowOfPropellantForPDIPhase1>0.)
    
    
        '''Propellant mass remaining after dV change in velocity is'''
        massPropellant-=massOfFuelForDeltaVForPDIPhase1    
    
        '''Set thrust angle to be in the direction of the velocity'''
        '''Compute the angle between the current vel vect and the pos vec'''       
        Angle = angle2D(uInitial,vInitial,xInitial/magOfPositionVector,yInitial/magOfPositionVector)        
    
        '''Rotate the position vector to the direction of the velocity'''        
        Thx, Thy = rotate(xInitial/magOfPositionVector, yInitial/magOfPositionVector, -Angle)

        '''New resultant velocity is'''
        resultantVelocity -=deltaVForPDIPhase1
        
        '''Compute velocity components'''
        u_ii = Thx*resultantVelocity
        v_ii = Thy*resultantVelocity
        
        
        '''Set thrust angle to be in the direction of the velocity'''
        thrustAngle = math.pi - Angle
 
        '''find thrust from the massflowofpropellant'''
        thrustForce = Isp*gravityOfEarth*massFlowOfPropellantForPDIPhase1
        
        assert(thrustForce<maxThrust)
        assert(thrustForce>minThrust)
    
        '''Compute the resultant acceleration acquired by the ship due to the thrust'''
        a_thrust = thrustForce/(massDry+massPropellant)
        # print(a_thrust)
        
        '''Determine the direction of the acceleration vector'''
        thx, thy = rotate(xInitial/magOfPositionVector, yInitial/magOfPositionVector, thrustAngle);
        
        # print(thx, thy)
        
        '''Resolve the resultant acceleration due to thrust into its components'''
        ax_thrust = thx*a_thrust
        ay_thrust = thy*a_thrust
        
        '''Compute the gravitational acceleration at this point'''
        a_gravity  = gravitationalParameter/(magOfPositionVector**2)
    
        '''Compute the components of this acceleration due to gravity'''
        ax_gravity = -a_gravity*xInitial/magOfPositionVector;
        ay_gravity = -a_gravity*yInitial/magOfPositionVector;
    
        
        '''Compute net x-dir and y-dir acceleration'''
        ax_net = ax_thrust + ax_gravity
        ay_net = ay_thrust + ay_gravity
    
        '''Compute change in velocity in dt_elipse time due to this acceleration'''
        u_i = uInitial + ax_net*dTForPDIPhase1       
        v_i = vInitial + ay_net*dTForPDIPhase1
        
        
        '''Compute the change in position due to this acceleration'''    
        x_i = xInitial + u_i*dTForPDIPhase1
        y_i = yInitial + v_i*dTForPDIPhase1 
        
        
        '''Compute the current altitude of the ship'''
        currentAltitude = math.sqrt(x_i*x_i+y_i*y_i) - radiusOfMoon       
        
        assert(currentAltitude>0.)
    
        
        '''Set'''
        xInitial = x_i
        yInitial = y_i
        '''Compute the magnitude of the position vector'''
        magOfPositionVector = math.sqrt(xInitial**2+ yInitial**2)
    
        uInitial = u_i
        vInitial = v_i
        resultantVelocity = math.sqrt(uInitial**2+vInitial**2)
        
        # if((i%1000)==0):
             # plotLocationOfShip(xInitial,yInitial,ax_net, ay_net)
             # plotVelOfShip(xInitial, yInitial, uInitial, vInitial)  
        


print("Powered Descent phase 1 complete\n")

'''Store the position, velocity and time taken'''
totalTime+=time
x4 = x_i
y4 = y_i
currentAltitudeAtEndOfPDI = currentAltitude
U_comp_vel_atEndOfPDIPhase1 = uInitial
V_comp_vel_atEndOfPDIPhase1 = vInitial
resultantVelocity = math.sqrt(U_comp_vel_atEndOfPDIPhase1**2+V_comp_vel_atEndOfPDIPhase1**2)

plotLocationOfShip(x4,y4,U_comp_vel_atEndOfPDIPhase1,V_comp_vel_atEndOfPDIPhase1)

print("The current altitude of the ship is %f" %currentAltitudeAtEndOfPDI)

print(" The position vector is %f %f" %(x4, y4))            

print(" The vel comp are %f %f and resultant vel is %f" %(U_comp_vel_atEndOfPDIPhase1, V_comp_vel_atEndOfPDIPhase1, resultantVelocity))
    
print("The propellant remaining after PDI phase 1 is %f" %massPropellant)
   
print("The time taken to reach this altitude is %f s" %totalTime)       


print("\n")
       
print("Initiating powered descent Phase 2\n")


################################################################
'''Powered Descent Phase 2 '''
#################################################################


'''Initialize time'''
time = 0
'''Set initial position as the last known position of the ship after 
the end of Hohmann transfer'''
xInitial = x4
yInitial = y4
'''Compute the magnitude of the position vector'''
magOfPositionVector = math.sqrt(xInitial**2+ yInitial**2)

'''Compute current altitude'''
currentAltitude = currentAltitudeAtEndOfPDI

'''Set initial velocity as the last known velocity of the ship after 
the end of Hohmann transfer'''
uInitial = U_comp_vel_atEndOfPDIPhase1
vInitial = V_comp_vel_atEndOfPDIPhase1
resultantVelocity = math.sqrt(uInitial**2+vInitial**2)


'''If we choose to full throttle in this phase'''
if(fullThrottleForPDIPhase2):
    print("PDIPhase2 in fullThrottle\n")
    for i in range(0,nIterForPDIPhase2):
        time+=dTForPDIPhase2
        '''Compute the angle between the current vel vect and the pos vec'''
        Angle = angle2D(uInitial,vInitial,xInitial/magOfPositionVector,yInitial/magOfPositionVector)
       
        '''Set thrust angle to be in the direction of the velocity'''
        thrustAngle = math.pi - Angle
        
        '''Set thrust force to be the max that the engine can provide'''
        thrustForce = maxThrust
        
        '''Compute the initial total mass'''
        m0 = massDry + massPropellant  
        
        '''Find mass flow of propellant needed to sustain the thrust'''
        massFlowOfPropellantForPDIPhase2 = thrustForce/(Isp*gravityOfEarth)
        
        '''Compute mass flow from massOfFuelForDeltaVForPDIPhase1'''
        massOfFuelForDeltaVForPDIPhase2 = massFlowOfPropellantForPDIPhase2*dTForPDIPhase2
        assert(massOfFuelForDeltaVForPDIPhase2>0)
        '''Propellant mass remaining after dV change in velocity is'''
        massPropellant-=massOfFuelForDeltaVForPDIPhase2
                
        '''COmpute the final mass after reduction in fuel'''
        mf = massDry + massPropellant

        '''Compute the deltaV achieved by this propellant flow'''
        deltaV = Isp*gravityOfEarth*math.log(m0/mf)
       
        '''Rotate the position vector to the direction of the velocity'''        
        Thx, Thy = rotate(xInitial/magOfPositionVector, yInitial/magOfPositionVector, -Angle)

        '''New resultant velocity is'''
        resultantVelocity -=deltaV
        
        '''Compute velocity components'''
        u_ii = Thx*resultantVelocity
        v_ii = Thy*resultantVelocity

        
        '''Compute the resultant acceleration acquired by the ship due to the thrust'''
        a_thrust = thrustForce/(massDry+massPropellant)
        
        '''Determine the direction of the acceleration vector'''
        thx, thy = rotate(xInitial/magOfPositionVector, yInitial/magOfPositionVector, thrustAngle);
        
        '''Resolve the resultant acceleration due to thrust into its components'''
        ax_thrust = thx*a_thrust
        ay_thrust = thy*a_thrust
        
        '''Compute the gravitational acceleration at this point'''
        a_gravity  = gravitationalParameter/(magOfPositionVector**2)
    
        '''Compute the components of this acceleration due to gravity'''
        ax_gravity = -a_gravity*xInitial/magOfPositionVector;
        ay_gravity = -a_gravity*yInitial/magOfPositionVector;
    
        
        '''Compute net x-dir and y-dir acceleration'''
        ax_net = ax_thrust + ax_gravity
        ay_net = ay_thrust + ay_gravity
        
        
        '''Compute change in velocity in time due to this acceleration'''
        u_i = uInitial + ax_net*dTForPDIPhase2       
        v_i = vInitial + ay_net*dTForPDIPhase2
        
      
        '''Compute the change in position due to this acceleration'''
        x_i = xInitial + u_i*dTForPDIPhase2
        y_i = yInitial + v_i*dTForPDIPhase2 
        
        '''Compute the current altitude of the ship'''
        currentAltitude = math.sqrt(x_i*x_i+y_i*y_i) - radiusOfMoon       
        
        assert(currentAltitude>0.)
    
        
        '''Set'''
        xInitial = x_i
        yInitial = y_i
        '''Compute the magnitude of the position vector'''
        magOfPositionVector = math.sqrt(xInitial**2+ yInitial**2)
    
        uInitial = u_i
        vInitial = v_i
        resultantVelocity = math.sqrt(uInitial**2+vInitial**2)
        
        # if((i%1000)==0):
             # plotLocationOfShip(xInitial,yInitial,ax_net, ay_net)
             # plotVelOfShip(xInitial, yInitial, uInitial, vInitial)  
        
elif(constantdeltaVForPDIPhase2==1):
    print("Initiating constant deltaV descent\n")
    for i in range(0,nIterForPDIPhase2):
        time+=dTForPDIPhase2
        
        '''Decrease the resultant velocity by deltaV value which is predefined'''
        '''Compute the final mass of the ship after the delta V maneuvre'''
        massAfterDeltaVForPDIPhase2 = (massDry+massPropellant)/math.exp(deltaVForPDIPhase2/(Isp*gravityOfEarth))
    
        '''Compute the mass of the fuel required to achieve this delta v'''
        massOfFuelForDeltaVForPDIPhase2 = (massDry+massPropellant)- massAfterDeltaVForPDIPhase2
        assert(massOfFuelForDeltaVForPDIPhase2>0)
    
        '''Compute mass flow from massOfFuelForDeltaVForPDIPhase2'''
        massFlowOfPropellantForPDIPhase2 = abs(massOfFuelForDeltaVForPDIPhase2)/dTForPDIPhase2
        assert(massFlowOfPropellantForPDIPhase2>0.)
    
    
        '''Propellant mass remaining after dV change in velocity is'''
        massPropellant-=massOfFuelForDeltaVForPDIPhase2    
    
        '''Set thrust angle to be in the direction of the velocity'''
        '''Compute the angle between the current vel vect and the pos vec'''
        Angle = angle2D(uInitial,vInitial,xInitial/magOfPositionVector,yInitial/magOfPositionVector)
        
        '''Rotate the position vector to the direction of the velocity'''        
        Thx, Thy = rotate(xInitial/magOfPositionVector, yInitial/magOfPositionVector, -Angle)

        '''New resultant velocity is'''
        resultantVelocity -=deltaVForPDIPhase2
        
        '''Compute velocity components'''
        u_ii = Thx*resultantVelocity
        v_ii = Thy*resultantVelocity
        
        '''Set thrust angle to be in the direction of the velocity'''
        thrustAngle = math.pi - Angle
        
        '''Or find thrust from the massflowofpropellant'''
        thrustForce = Isp*gravityOfEarth*massFlowOfPropellantForPDIPhase2
    
        assert(thrustForce<maxThrust)
        assert(thrustForce>minThrust)
    
        '''Compute the resultant acceleration acquired by the ship due to the thrust'''
        a_thrust = thrustForce/(massDry+massPropellant)
        
        '''Determine the direction of the acceleration vector'''
        thx, thy = rotate(xInitial/magOfPositionVector, yInitial/magOfPositionVector, thrustAngle);
        
        # plotVectors(thx, thy, uInitial, vInitial, xInitial, yInitial) 
        
        '''Resolve the resultant acceleration due to thrust into its components'''
        ax_thrust = thx*a_thrust
        ay_thrust = thy*a_thrust
        
        '''Compute the gravitational acceleration at this point'''
        a_gravity  = gravitationalParameter/(magOfPositionVector**2)
    
        '''Compute the components of this acceleration due to gravity'''
        ax_gravity = -a_gravity*xInitial/magOfPositionVector;
        ay_gravity = -a_gravity*yInitial/magOfPositionVector;
    
        
        '''Compute net x-dir and y-dir acceleration'''
        ax_net = ax_thrust + ax_gravity
        ay_net = ay_thrust + ay_gravity
            
        '''Compute change in velocity in dt_elipse time due to this acceleration'''
        u_i = uInitial + ax_net*dTForPDIPhase2       
        v_i = vInitial + ay_net*dTForPDIPhase2

        
        '''Compute the change in position due to this acceleration'''        
        x_i = xInitial + u_i*dTForPDIPhase2
        y_i = yInitial + v_i*dTForPDIPhase2 
        
        
        '''Compute the current altitude of the ship'''
        currentAltitude = math.sqrt(x_i*x_i+y_i*y_i) - radiusOfMoon       
        
        assert(currentAltitude>0.)
    
        
        '''Set'''
        xInitial = x_i
        yInitial = y_i
        '''Compute the magnitude of the position vector'''
        magOfPositionVector = math.sqrt(xInitial**2+ yInitial**2)
    
        uInitial = u_i
        vInitial = v_i
        resultantVelocity = math.sqrt(uInitial**2+vInitial**2)
        
        # if((i%1000)==0):
             # plotLocationOfShip(xInitial,yInitial,ax_net, ay_net)
             # plotVelOfShip(xInitial, yInitial, uInitial, vInitial)  
elif(constantdeltaAltForPDIPhase2==1):
    print("Initiating constant delta A descent\n")
    for i in range(0,nIterForPDIPhase2):
        time+=dTForPDIPhase2
        '''Check if the velocity is some dropfactor*Altitude'''
        if( abs(currentAltitude/resultantVelocity)>=descentFactor):
            thrustForce = 0
            deltaVForPDIPhase2 = 0
        else:
            '''Set deltaV as a function of deltaA'''
            deltaVForPDIPhase2 = (1./descentFactor)*deltaAltForPDIPhase2
        
            '''Decrease the resultant velocity by deltaV'''
            '''Compute the final mass of the ship after the delta V maneuvre'''
            massAfterDeltaVForPDIPhase2 = (massDry+massPropellant)/math.exp(deltaVForPDIPhase2/(Isp*gravityOfEarth))
    
            '''Compute the mass of the fuel required to achieve this delta v'''
            massOfFuelForDeltaVForPDIPhase2 = (massDry+massPropellant)- massAfterDeltaVForPDIPhase2
            assert(massOfFuelForDeltaVForPDIPhase2>0)
    
            '''Compute mass flow from massOfFuelForDeltaVForPDIPhase2'''
            massFlowOfPropellantForPDIPhase2 = abs(massOfFuelForDeltaVForPDIPhase2)/dTForPDIPhase2
            assert(massFlowOfPropellantForPDIPhase2>0.)
    
            '''Propellant mass remaining after dV change in velocity is'''
            massPropellant-=massOfFuelForDeltaVForPDIPhase2    
             
            '''find thrust from the massflowofpropellant'''
            thrustForce = Isp*gravityOfEarth*massFlowOfPropellantForPDIPhase2
            
            if(thrustForce>maxThrust):
                thrustForce = maxThrust
            elif(thrustForce<minThrust):
                thrustForce = minThrust
    
        '''Set thrust angle to be in the direction of the velocity'''
        '''Compute the angle between the current vel vect and the pos vec'''
        Angle = angle2D(uInitial,vInitial,xInitial/magOfPositionVector,yInitial/magOfPositionVector)
        
        '''Rotate the position vector to the direction of the velocity'''        
        Thx, Thy = rotate(xInitial/magOfPositionVector, yInitial/magOfPositionVector, -Angle)

        '''New resultant velocity is'''
        resultantVelocity -=deltaVForPDIPhase2
        
        '''Compute velocity components'''
        u_ii = Thx*resultantVelocity
        v_ii = Thy*resultantVelocity
        
        '''Set thrust angle to be in the direction of the velocity'''
        thrustAngle = math.pi - Angle
        
        '''Compute the resultant acceleration acquired by the ship due to the thrust'''
        a_thrust = thrustForce/(massDry+massPropellant)
        
        '''Determine the direction of the acceleration vector'''
        thx, thy = rotate(xInitial/magOfPositionVector, yInitial/magOfPositionVector, thrustAngle);
        
        # plotVectors(thx, thy, uInitial, vInitial, xInitial, yInitial) 
        
        '''Resolve the resultant acceleration due to thrust into its components'''
        ax_thrust = thx*a_thrust
        ay_thrust = thy*a_thrust
        
        '''Compute the gravitational acceleration at this point'''
        a_gravity  = gravitationalParameter/(magOfPositionVector**2)
    
        '''Compute the components of this acceleration due to gravity'''
        ax_gravity = -a_gravity*xInitial/magOfPositionVector;
        ay_gravity = -a_gravity*yInitial/magOfPositionVector;
    
        
        '''Compute net x-dir and y-dir acceleration'''
        ax_net = ax_thrust + ax_gravity
        ay_net = ay_thrust + ay_gravity
        
        '''Compute change in velocity in dt_elipse time due to this acceleration'''
        u_i = uInitial + ax_net*dTForPDIPhase2       
        v_i = vInitial + ay_net*dTForPDIPhase2

        
        '''Compute the change in position due to this acceleration'''        
        x_i = xInitial + u_i*dTForPDIPhase2
        y_i = yInitial + v_i*dTForPDIPhase2 
        
        
        '''Compute the current altitude of the ship'''
        previousAltitude = currentAltitude
        currentAltitude = math.sqrt(x_i*x_i+y_i*y_i) - radiusOfMoon       
        assert(currentAltitude>0.)
        
    
        
        '''Set'''
        xInitial = x_i
        yInitial = y_i
        '''Compute the magnitude of the position vector'''
        magOfPositionVector = math.sqrt(xInitial**2+ yInitial**2)
    
        uInitial = u_i
        vInitial = v_i
        resultantVelocity = math.sqrt(uInitial**2+vInitial**2)
        
        # if((i%1000)==0):
             # plotLocationOfShip(xInitial,yInitial,ax_net, ay_net)
             # plotVelOfShip(xInitial, yInitial, uInitial, vInitial)              
            

print("Powered Descent phase 2 complete\n")

'''Store position velocity and time'''
totalTime+=time
x5 = x_i
y5 = y_i
'''Compute the magnitude of the position vector'''
magOfPositionVector = math.sqrt(x5**2+ y5**2)


currentAltitudeAtEndOfPDIPhase2 = currentAltitude
U_comp_vel_atEndOfPDIPhase2 = uInitial
V_comp_vel_atEndOfPDIPhase2 = vInitial
resultantVelocity = math.sqrt(U_comp_vel_atEndOfPDIPhase2**2+V_comp_vel_atEndOfPDIPhase2**2)

plotLocationOfShip(x5,y5, U_comp_vel_atEndOfPDIPhase2, V_comp_vel_atEndOfPDIPhase2)

print("The current altitude of the ship is %f" %currentAltitudeAtEndOfPDIPhase2)

print(" The position vector is %f %f" %(x4, y4))            

print(" The vel comp are %f %f and resultant vel is %f" %(U_comp_vel_atEndOfPDIPhase2, V_comp_vel_atEndOfPDIPhase2, resultantVelocity))
    
print("The propellant remaining after PDI phase 1 is %f" %massPropellant)
   
print("The time taken to reach this altitude is %f s" %totalTime)       


print("\n\n")


################################################################
'''TouchDown Phase '''
#################################################################


print(" Initiating touchdown maneuver... \n")



'''Initialize time'''
time = 0
'''Set initial position as the last known position of the ship after 
the end of Hohmann transfer'''
xInitial = x5
yInitial = y5
'''Compute the magnitude of the position vector'''
magOfPositionVector = math.sqrt(xInitial**2+ yInitial**2)

'''Compute current altitude'''
currentAltitude = currentAltitudeAtEndOfPDIPhase2

'''Set initial velocity as the last known velocity of the ship after 
the end of Hohmann transfer'''
uInitial = U_comp_vel_atEndOfPDIPhase2
vInitial = V_comp_vel_atEndOfPDIPhase2
resultantVelocity = math.sqrt(uInitial**2+vInitial**2)

print("Initiating constant delta A descent\n")
for i in range(0,nIterForTouchDown):
    time+=dTForTouchDown
    '''Check if the velocity is some dropfactor*Altitude'''
    if( abs(currentAltitude/resultantVelocity)>=descentFactorForTouchDown):
        thrustForce = 0
        deltaVForTouchDown = 0
    else:
        '''Set deltaV as a function of deltaA'''
        deltaVForTouchDown = (1./descentFactorForTouchDown)*deltaAltForTouchDown
    
        '''Decrease the resultant velocity by deltaV'''
        '''Compute the final mass of the ship after the delta V maneuvre'''
        massAfterDeltaVForTouchDown = (massDry+massPropellant)/math.exp(deltaVForTouchDown/(Isp*gravityOfEarth))

        '''Compute the mass of the fuel required to achieve this delta v'''
        massOfFuelForDeltaVForTouchDown = (massDry+massPropellant)- massAfterDeltaVForTouchDown
        assert(massOfFuelForDeltaVForTouchDown>0)

        '''Compute mass flow from massOfFuelForDeltaVForTouchDown'''
        massFlowOfPropellantForTouchDown = abs(massOfFuelForDeltaVForTouchDown)/dTForTouchDown
        assert(massFlowOfPropellantForTouchDown>0.)

        '''Propellant mass remaining after dV change in velocity is'''
        massPropellant-=massOfFuelForDeltaVForTouchDown    
        
        '''find thrust from the massflowofpropellant'''
        thrustForce = Isp*gravityOfEarth*massFlowOfPropellantForTouchDown
        
        if(thrustForce>maxThrust):
            thrustForce = maxThrust
        elif(thrustForce<minThrust):
            thrustForce = minThrust

    '''Set thrust angle to be in the direction of the velocity'''
    '''Compute the angle between the current vel vect and the pos vec'''
    Angle = angle2D(uInitial,vInitial,xInitial/magOfPositionVector,yInitial/magOfPositionVector)
    
    '''Rotate the position vector to the direction of the velocity'''        
    Thx, Thy = rotate(xInitial/magOfPositionVector, yInitial/magOfPositionVector, -Angle)

    '''New resultant velocity is'''
    resultantVelocity -=deltaVForTouchDown
    
    '''Compute velocity components'''
    u_ii = Thx*resultantVelocity
    v_ii = Thy*resultantVelocity
    
    '''Set thrust angle to be in the direction of the velocity'''
    thrustAngle = math.pi - Angle
    
    '''Compute the resultant acceleration acquired by the ship due to the thrust'''
    a_thrust = thrustForce/(massDry+massPropellant)
    
    '''Determine the direction of the acceleration vector'''
    thx, thy = rotate(xInitial/magOfPositionVector, yInitial/magOfPositionVector, thrustAngle);
    
    # plotVectors(thx, thy, uInitial, vInitial, xInitial, yInitial) 
    
    '''Resolve the resultant acceleration due to thrust into its components'''
    ax_thrust = thx*a_thrust
    ay_thrust = thy*a_thrust
    
    '''Compute the gravitational acceleration at this point'''
    a_gravity  = gravitationalParameter/(magOfPositionVector**2)

    '''Compute the components of this acceleration due to gravity'''
    ax_gravity = -a_gravity*xInitial/magOfPositionVector;
    ay_gravity = -a_gravity*yInitial/magOfPositionVector;

    
    '''Compute net x-dir and y-dir acceleration'''
    ax_net = ax_thrust + ax_gravity
    ay_net = ay_thrust + ay_gravity
    
    '''Compute change in velocity in dt_elipse time due to this acceleration'''
    u_i = uInitial + ax_net*dTForTouchDown       
    v_i = vInitial + ay_net*dTForTouchDown

    '''Compute the change in position due to this acceleration'''        
    x_i = xInitial + u_i*dTForTouchDown
    y_i = yInitial + v_i*dTForTouchDown 
    
    
    '''Compute the current altitude of the ship'''
    previousAltitude = currentAltitude
    currentAltitude = math.sqrt(x_i*x_i+y_i*y_i) - radiusOfMoon       
    
    
    '''Set'''
    xInitial = x_i
    yInitial = y_i
    '''Compute the magnitude of the position vector'''
    magOfPositionVector = math.sqrt(xInitial**2+ yInitial**2)

    uInitial = u_i
    vInitial = v_i
    resultantVelocity = math.sqrt(uInitial**2+vInitial**2)
    
    # if((i%1000)==0):
         # plotLocationOfShip(xInitial,yInitial,ax_net, ay_net)
         # plotVelOfShip(xInitial, yInitial, uInitial, vInitial)              

print("Landing Maneuver is complete\n")

'''Store position, velocity and time'''
totalTime+=time
x6 = x_i
y6 = y_i
'''Compute the magnitude of the position vector'''
magOfPositionVector = math.sqrt(x6**2+ y6**2)


currentAltitudeAtTouchDown = currentAltitude
U_comp_vel_atTouchDown = uInitial
V_comp_vel_atTouchDown = vInitial
resultantVelocity = math.sqrt(U_comp_vel_atTouchDown**2+V_comp_vel_atTouchDown**2)

plotLocationOfShip(x6,y6, U_comp_vel_atEndOfPDIPhase2, V_comp_vel_atEndOfPDIPhase2)

print("The current altitude of the ship is %f" %currentAltitudeAtTouchDown)

print(" The position vector is %f %f" %(x6, y6))            

print(" The vel comp are %f %f and resultant vel is %f" %(U_comp_vel_atTouchDown, V_comp_vel_atTouchDown, resultantVelocity))
    
print("The propellant remaining after PDI phase 1 is %f" %massPropellant)
   
print("The time taken to reach this altitude is %f s" %totalTime)       

if(currentAltitude<=0):    
    if((resultantVelocity-2.)<1e-3):       
        print("THE EAGLE HAS LANDED\n")
    else:
        print("The lander has Crashed!")

