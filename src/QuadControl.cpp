#include "Common.h"
#include "QuadControl.h"
#include <iostream>
#include <string.h>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <cmath>

#include "Utility/SimpleConfig.h"
#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

void QuadControl::Init()
{
  BaseController::Init();

  // variables needed for integral control
  integratedAltitudeError = 0;
    
#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();
   
  // Load parameters (default to 0)
  kpPosXY = config->Get(_config+".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  KiPosZ = config->Get(_config + ".KiPosZ", 0);
     
  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);

  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);

  kpPQR = config->Get(_config + ".kpPQR", V3F());

  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);
#else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
  // Convert a desired 3-axis moment and collective thrust command to 
  //   individual motor thrust commands
  // INPUTS: 
  //   desCollectiveThrust: desired collective thrust [N]
  //   desMoment: desired rotation moment about each axis [N m]
  // OUTPUT:
  //   set class member variable cmd (class variable for graphing) where
  //   cmd.desiredThrustsN[0..3]: motor commands, in [N]

  // HINTS: 
  // - you can access parts of desMoment via e.g. desMoment.x
  // You'll need the arm length parameter L, and the drag/thrust ratio kappa

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

    
  cmd.desiredThrustsN[0] = mass * 9.81f / 4.f; // front left
  cmd.desiredThrustsN[1] = mass * 9.81f / 4.f; // front right
  cmd.desiredThrustsN[2] = mass * 9.81f / 4.f; // rear left
  cmd.desiredThrustsN[3] = mass * 9.81f / 4.f; // rear right
    
    
    float Mx = momentCmd.x;
    float My = momentCmd.y;
    float Mz = -momentCmd.z;
    float Ftotal = collThrustCmd ;
    
    cout << " Mx: "<< Mx <<"\n";
    cout << " My:"<< My <<"\n";
    cout << " Mz:"<< Mz <<"\n";
    cout << "f_total:"<< Ftotal <<"\n";
    
    float l = L/sqrt(2);
    
    
    float f1 = (Mx/l  + My/l  - Mz/kappa + Ftotal) / 4.f;
    float f2 = (-Mx/l  + My/l  + Mz/kappa + Ftotal) / 4.f;
    float f3 = (Mx/l  - My/l  + Mz/kappa + Ftotal) / 4.f;
    float f4 = (-Mx/l  - My/l  - Mz/kappa + Ftotal) / 4.f;
     
    
    
    cmd.desiredThrustsN[0] = f1; // front left
    cmd.desiredThrustsN[1] = f2; // front right
    cmd.desiredThrustsN[2] = f3; // rear left
    cmd.desiredThrustsN[3] = f4; // rear right
    
    cout<<"cmd.desiredThrustsN[0]:"<< cmd.desiredThrustsN[0]<<"\n";
    cout<<"cmd.desiredThrustsN[1]:"<< cmd.desiredThrustsN[1]<<"\n";
    cout<<"cmd.desiredThrustsN[2]:"<< cmd.desiredThrustsN[2]<<"\n";
    cout<<"cmd.desiredThrustsN[3]:"<< cmd.desiredThrustsN[3]<<"\n";
    
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return cmd;
}


V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  // Calculate a desired 3-axis moment given a desired and current body rate
  // INPUTS: 
  //   pqrCmd: desired body rates [rad/s]
  //   pqr: current or estimated body rates [rad/s]
  // OUTPUT:
  //   return a V3F containing the desired moments for each of the 3 axes

  // HINTS:
  //  - you can use V3Fs just like scalars: V3F a(1,1,1), b(2,3,4), c; c=a-b;
  //  - you'll need parameters for moments of inertia Ixx, Iyy, Izz
  //  - you'll also need the gain parameter kpPQR (it's a V3F)

  V3F momentCmd;

  ////////////////////////////// BEGIN STUDENT CODE //////////////////////////
    
    
    V3F delta = pqrCmd - pqr;
    V3F u_pqr = kpPQR*delta;
    V3F I_xyz = V3F(Ixx, Iyy, -Izz);
    momentCmd = I_xyz * u_pqr;
     
    return  momentCmd; //u_pqr;
}

// returns a desired roll and pitch rate 
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  // Calculate a desired pitch and roll angle rates based on a desired global
  //   lateral acceleration, the current attitude of the quad, and desired
  //   collective thrust command
  // INPUTS: 
  //   accelCmd: desired acceleration in global XY coordinates [m/s2]
  //   attitude: current or estimated attitude of the vehicle
  //   collThrustCmd: desired collective thrust of the quad [N]
  // OUTPUT:
  //   return a V3F containing the desired pitch and roll rates. The Z
  //     element of the V3F should be left at its default value (0)

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the roll/pitch gain kpBank
  //  - collThrustCmd is a force in Newtons! You'll likely want to convert it to acceleration first

  V3F pqrCmd;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

    
 
    V3F bc = accelCmd/collThrustCmd;
    
    float acceleration = collThrustCmd/mass;
    float target02 = CONSTRAIN(accelCmd.x/acceleration, -maxTiltAngle, maxTiltAngle);
    float target12 = CONSTRAIN(accelCmd.y/acceleration, -maxTiltAngle, maxTiltAngle);
    bc = V3F(target02, target12, bc[2]);
    
    V3F ba = V3F(R(0,2) , R(1,2), 0);
    V3F b_dot_c = kpBank * (bc - ba);
    V3F c0 = V3F(R(1,0), -R(0,0), 0);
    V3F c1 = V3F(R(1,1), -R(0,1), 0);
    float pitch = c1.dot(b_dot_c)/R(2,2);
    
    pqrCmd = V3F(c0.dot(b_dot_c)/R(2,2), pitch, 0);
    cout<<"pqrCmd pictch:"<< pitch <<", "<<c1.dot(b_dot_c)/R(2,2)<< "\n";
 
  /////////////////////////////// END STUDENT CODE ////////////////////////////
 
  return pqrCmd;
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude,
                                   float accelZCmd, float dt)
{
  // Calculate desired quad thrust based on altitude setpoint, actual altitude,
  //   vertical velocity setpoint, actual vertical velocity, and a vertical 
  //   acceleration feed-forward command
  // INPUTS: 
  //   posZCmd, velZCmd: desired vertical position and velocity in NED [m]
  //   posZ, velZ: current vertical position and velocity in NED [m]
  //   accelZCmd: feed-forward vertical acceleration in NED [m/s2]
  //   dt: the time step of the measurements [seconds]
  // OUTPUT:
  //   return a collective thrust command in [N]

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the gain parameters kpPosZ and kpVelZ
  //  - maxAscentRate and maxDescentRate are maximum vertical speeds. Note they're both >=0!
  //  - make sure to return a force, not an acceleration
  //  - remember that for an upright quad in NED, thrust should be HIGHER if the desired Z acceleration is LOWER

  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float thrust = 0;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
 
    
    cout<<">>>posZCmd:"<<posZCmd<<" posZ:"<<posZ<<"\n";
    cout<<">>>velZCmd:"<<velZCmd<<" velZ:"<<velZ<<"\n";
    
    if (velZ < 0){
        velZ = - fmodf(-velZ, maxAscentRate);
    }else{
        velZ = - fmodf(-velZ, maxDescentRate);
    }
    
     
   cout<<"kpPosZ:"<<kpPosZ<<"\n";
   cout<<"kpVelZ:"<<kpVelZ<<"\n";
    V3F zk = V3F(kpPosZ, kpVelZ, 1);
    V3F delta = V3F(posZCmd-posZ, (velZCmd-velZ), accelZCmd);
    
   cout<<"posZCmd:"<<posZCmd<<" posZ:"<<posZ<<"\n";
   cout<<"velZCmd:"<<velZCmd<<" velZ:"<<velZ<<"\n";
    float u_bar_1 = zk.dot(delta);
   cout<<"u_bar_1:"<<u_bar_1<<"\n"; 
    // basic integral error 
    this->integratedAltitudeError += (posZCmd-posZ)*dt;
    cout<<"***(posZCmd-posZ)*dt:"<<(posZCmd-posZ)*dt<<"\n";
    cout<<"***integratedAltitudeError:"<<this->integratedAltitudeError<<"\n";
    u_bar_1 += this->integratedAltitudeError * KiPosZ;
    cout<<"u_bar_1:"<<u_bar_1<<"\n";
    
    thrust = -mass*(u_bar_1-CONST_GRAVITY)/R(2,2);
    
    cout<<"R(2,2):"<<R(2,2)<<"------thrust:"<<thrust<<"\n";
    cout<<"------------mass*CONST_GRAVITY:"<<mass*CONST_GRAVITY<<"\n";
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmdFF)
{
  // Calculate a desired horizontal acceleration based on 
  //  desired lateral position/velocity/acceleration and current pose
  // INPUTS: 
  //   posCmd: desired position, in NED [m]
  //   velCmd: desired velocity, in NED [m/s]
  //   pos: current position, NED [m]
  //   vel: current velocity, NED [m/s]
  //   accelCmdFF: feed-forward acceleration, NED [m/s2]
  // OUTPUT:
  //   return a V3F with desired horizontal accelerations. 
  //     the Z component should be 0
  // HINTS: 
  //  - use the gain parameters kpPosXY and kpVelXY
  //  - make sure you limit the maximum horizontal velocity and acceleration
  //    to maxSpeedXY and maxAccelXY

  // make sure we don't have any incoming z-component
  accelCmdFF.z = 0;
  velCmd.z = 0;
  posCmd.z = pos.z;

  // we initialize the returned desired acceleration to the feed-forward value.
  // Make sure to _add_, not simply replace, the result of your controller
  // to this variable
  V3F accelCmd = accelCmdFF;
 
    ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
    
    
    V3F pd = V3F(kpPosXY, kpVelXY, 1.f);
    posCmd  = -posCmd;
    pos  = -pos;
    velCmd  = -velCmd;
    vel = -vel;
     
    
    V3F deltaX = V3F((posCmd-pos)[0], CONSTRAIN((velCmd-vel)[0], -maxSpeedXY, maxSpeedXY), accelCmd[0]);
    cout<<"posCmd.x:"<<posCmd[0]<<" pos.x:"<< pos[0]<< " deltaX:"<<deltaX[0] <<"\n";
    cout<<"velCmd.x:"<<velCmd[0]<<" vel.x:"<< vel[0]<< " deltaV:"<<deltaX[1] <<"\n";
    V3F deltaY = V3F((posCmd-pos)[1], CONSTRAIN((velCmd-vel)[1], -maxSpeedXY, maxSpeedXY), accelCmd[1]);
    cout<<"posCmd.y:"<<posCmd[1]<<" pos.y:"<< pos[1]<< " deltaY:"<<deltaY[0] <<"\n";
    cout<<"velCmd.y:"<<velCmd[1]<<" vel.y:"<< vel[1]<<"\n";
    
    accelCmd = V3F(pd.dot(deltaX), pd.dot(deltaY), 0.f);
    cout<<"*******accelCmd.x:"<<accelCmd.x<<" accelCmd.y:"<<accelCmd.y<<"\n";
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return accelCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{
  // Calculate a desired yaw rate to control yaw to yawCmd
  // INPUTS: 
  //   yawCmd: commanded yaw [rad]
  //   yaw: current yaw [rad]
  // OUTPUT:
  //   return a desired yaw rate [rad/s]
  // HINTS: 
  //  - use fmodf(foo,b) to unwrap a radian angle measure float foo to range [0,b]. 
  //  - use the yaw control gain parameter kpYaw
     
    float yawRateCmd=0; 
    
    cout<<"++++kpYaw:"<< kpYaw<< " yawCmd:"<< yawCmd<<" yaw:"<<yaw<<"\n";
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

    
    yawRateCmd = kpYaw * (yawCmd - yaw);
    cout<<"yawRateCmd:"<<yawRateCmd<<"\n";
    
 
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return yawRateCmd;

}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
  curTrajPoint = GetNextTrajectoryPoint(simTime);

  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);
    
  // reserve some thrust margin for angle control
  float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
  collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);
  cout<<"miniMortorThrust:"<<(minMotorThrust+ thrustMargin)*4.f<<"\n";
   cout<<"collThrustCmd:"<<collThrustCmd<<"\n";
    
  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);
  
  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

  V3F desMoment = BodyRateControl(desOmega, estOmega);

  return GenerateMotorCommands(collThrustCmd, desMoment);
}
