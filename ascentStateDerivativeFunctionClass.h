/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      160518    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */


#ifndef ASCENTSTATEDERIVATIVEFUNCTIONCLASS_H
#define ASCENTSTATEDERIVATIVEFUNCTIONCLASS_H


#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <sstream>


#include <Eigen/Core>
#include <Eigen/Geometry>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <string>
#include <cmath>

//#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include <tudatApplications/thesisProject/linearAlgebraTypesUpdated.h>
#include <Tudat/Mathematics/BasicMathematics/coordinateConversions.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h>
//#include <tudatApplications/thesisProject/referenceFrameTransformationsUpdated.h>

/// Testing the celestial body class ///
#include <thesisProject/celestialBody.h>            // Final version

/// Testing the vehicle class ///
#include <thesisProject/MarsAscentVehicle.h>    // Final version

/// Testing the current state and time and its updater ///
#include <thesisProject/stateAndTime.h>             // Final version

/// Testing the auxiliary equations ///
#include <thesisProject/Auxiliary.h>                // Original test file

/// Testing the other required functions ///
#include <thesisProject/projectLibraries/otherRequiredFunctions.h>              // Original test file

class ascentStateDerivativeFunctionClass
{
public:
    ascentStateDerivativeFunctionClass(const celestialBody& planet_, const MarsAscentVehicle& MAV_){

        // Create variables to be used in this function

            Mars = planet_;                               // The celestial body class
            MAV = MAV_;                                 // The MAV class
//            StateAndTime stateAndTime = currentStateAndTime_;      // The current state and time class

            rotationalVelocityMars = Mars.rotationalVelocity();
            primeMeridianAngle = Mars.primeMeridianAngle();
            inertialFrameTime = Mars.inertialFrameTime();

    }         // NOTE TO SELF: WHEN INITIALIZING A NEW CLASS WITHOUT ANY PRESETS, MAKE SURE TO INCLUDE {} BEHIND (). IN THIS CASE THAT WOULD LOOK LIKE THIS: ascentStateDerivativeFunctionClass(){}

    /// Debug ///

    const double getTestValue(){

        const double testValue = 1;

        return testValue;
    }

    /// Debug ///

    tudat::basic_mathematics::Vector7d ascentStateDerivativeFunction(const double currentTime_, const tudat::basic_mathematics::Vector7d& currentState){



    /// Compute current spherical state ///

//        tudat::basic_mathematics::Vector7d currentState = stateAndTime.getCurrentState();

        const double xPosition = currentState(0);            // x position coordinate definition
        const double yPosition = currentState(1);            // y position coordinate definition
        const double zPosition = currentState(2);            // z position coordinate definition
        const double xVelocity = currentState(3);            // x velocity coordinate definition
        const double yVelocity = currentState(4);            // y velocity coordinate definition
        const double zVelocity = currentState(5);            // z velocity coordinate definition
        const double massMAV = currentState(6);              // MAV mass definition
        const double currentTime = currentTime_;   // current time definition

    // Computations

        const double Radius = sqrt(xPosition*xPosition+yPosition*yPosition+zPosition*zPosition);         // r [km]

        const double inertialVelocity = sqrt(xVelocity*xVelocity+yVelocity*yVelocity+zVelocity*zVelocity);       // V_I [km/s]

        const double inertialLongitude = atan2(yPosition,xPosition);         // lambda [rad]

        const double Latitude = asin(zPosition/Radius);              // delta [rad]

        const double rotationalLongitude = inertialLongitude-rotationalVelocityMars*(inertialFrameTime+currentTime)+primeMeridianAngle;  // tau [rad]

        double inertialLongitudeChange_;       // lambda_dot [rad/s] (placeholder)

        // Avoiding singularities
        if ((xPosition*xPosition+yPosition*yPosition) == 0){

            inertialLongitudeChange_ = 0;
        }
        else {
            inertialLongitudeChange_ = (xPosition*yVelocity-yPosition*xVelocity)/(xPosition*xPosition+yPosition*yPosition);
        };

        const double inertialLongitudeChange = inertialLongitudeChange_;        // lambda_dot [rad/s] (actual parameter)

        double rotationalLongitudeChange_ = inertialLongitudeChange-rotationalVelocityMars;     // tau_dot [rad/s] (placeholder)

        if (rotationalLongitudeChange_<=1e-15){  // Setting the accuracy to 1e-15 to avoid problems in the beginning with rounding errors...

            rotationalLongitudeChange_ = 0;

        };

        const double rotationalLongitudeChange = rotationalLongitudeChange_;    // tau_dot [rad/s] (actual parameter)

    /*    /// Debug ///

        std::cout<<"rotationalVelocityMars = "<<rotationalVelocityMars<<std::endl;
        std::cout<<"inertialLongitudeChange = "<<inertialLongitudeChange<<std::endl;
        std::cout<<"rotationalVelocityMars-7.088e-05 = "<<rotationalVelocityMars-7.088e-05<<std::endl;
        std::cout<<"inertialLongitudeChange-7.088e-05 = "<<inertialLongitudeChange-7.088e-05<<std::endl;
        std::cout<<"(xPosition*yVelocity-yPosition*xVelocity) = "<<(xPosition*yVelocity-yPosition*xVelocity)<<std::endl;
        std::cout<<"(xPosition*xPosition+yPosition*yPosition) = "<<(xPosition*xPosition+yPosition*yPosition)<<std::endl;
    //*/


        const double RadiusChange = (xPosition*xVelocity+yPosition*yVelocity+zPosition*zVelocity)/(Radius);      // radial velocity [km/s]

        double LatitudeChange_; // delta_dot [rad/s] (placeholder)

        if ((Radius*Radius*sqrt(1-(zPosition/Radius)*(zPosition/Radius))) == 0){
            LatitudeChange_ = 0;
        }
        else{
            LatitudeChange_ = (Radius*zVelocity-zPosition*RadiusChange)/(Radius*Radius*sqrt(1-(zPosition/Radius)*(zPosition/Radius)));
        };

        const double LatitudeChange = LatitudeChange_;  // delta_dot [rad/s] (actual parameter)

        const double localMarsRotationalVelocity = rotationalVelocityMars*Radius*cos(Latitude);  // V_M [km/s]

        const double inertialFlightPathAngle = asin(RadiusChange/inertialVelocity);          // gamma_I [rad]

        const double inertialAzimuth = atan2((inertialLongitudeChange*cos(Latitude)),LatitudeChange);    // chi_I [rad]

        const double rotationalVelocity = sqrt(localMarsRotationalVelocity*localMarsRotationalVelocity+inertialVelocity*inertialVelocity-2*localMarsRotationalVelocity*inertialVelocity*cos(inertialFlightPathAngle)*sin(inertialAzimuth));  // V_R [km/s]

        double rotationalFlightPathAngle_; // gamma_R [rad]  (placeholder)

        if (rotationalVelocity == 0){       // Setting the initial flight path angle in the rotational frame to 90 deg (or pi/s)

            rotationalFlightPathAngle_ = tudat::mathematical_constants::LONG_PI/2;
        }
        else {
            rotationalFlightPathAngle_ = asin(RadiusChange/rotationalVelocity);
        };

        const double rotationalFlightPathAngle = rotationalFlightPathAngle_;    // gamma_R [rad] (actual parameter)



        const double rotationalAzimuth = atan2((rotationalLongitudeChange*cos(Latitude)),LatitudeChange);    // chi_R [rad]

    /*    // Check output
        std::cout<<"Radius = "<<Radius<<std::endl;
        std::cout<<"inertialVelocity = "<<inertialVelocity<<std::endl;
        std::cout<<"inertialLongitude = "<<inertialLongitude<<std::endl;
        std::cout<<"Latitude = "<<Latitude<<std::endl;
        std::cout<<"rotationalLongitude = "<<rotationalLongitude<<std::endl;
        std::cout<<"inertialLongitudeChange = "<<inertialLongitudeChange<<std::endl;
        std::cout<<"rotationalLongitudeChange = "<<rotationalLongitudeChange<<std::endl;
        std::cout<<"RadiusChange = "<<RadiusChange<<std::endl;
        std::cout<<"LatitudeChange = "<<LatitudeChange<<std::endl;
        std::cout<<"localMarsRotationalVelocity = "<<localMarsRotationalVelocity<<std::endl;
        std::cout<<"inertialFlightPathAngle = "<<inertialFlightPathAngle<<std::endl;
        std::cout<<"inertialAzimuth = "<<inertialAzimuth<<std::endl;
        std::cout<<"rotationalVelocity = "<<rotationalVelocity<<std::endl;
        std::cout<<"rotationalFlightPathAngle = "<<rotationalFlightPathAngle<<std::endl;
        std::cout<<"rotationalAzimuth = "<<rotationalAzimuth<<std::endl;

    //*/


    /// Testing the local air temperature function ///

        const double currentAltitude = Radius-Mars.bodyReferenceRadius();

        const double currentTemperature = air_temperature::airTemperature(Mars.temperaturePolyCoefficients(), Mars.temperatureAltitudeRanges(),currentAltitude);

    /*   // Check output
        std::cout<<"currentAltitude = "<<currentAltitude<<std::endl;
        std::cout<<"currentTemperature = "<<currentTemperature<<std::endl;
        std::cout<<"Radius = "<<Radius<<std::endl;
        std::cout<<"R_MOLA = "<<Mars.bodyReferenceRadius()<<std::endl;
        std::cout<<"Radius-3395.4 = "<<Radius-3395.4<<std::endl;
        std::cout<<"R_MOLA-3396 = "<<Mars.bodyReferenceRadius()-3396<<std::endl;
    //*/

    /// Testing the local air density function ///

        const double currentDensity= air_density::airDensity(Mars.densityPolyCoefficients(),  currentAltitude);

    //    std::cout<<"The current air density = "<<currentDensity<<std::endl;

    /// Testing the ascentDragForce function ///

    /*    // Initially only with machNumber
        const double machNumber = Drag::ascentDragForce(rotationalVelocity,currentTemperature,Mars.adiabeticIndex(),Mars.specificGasConstant());

        std::cout<<"The current Mach number = "<<machNumber<<std::endl;


    /// Testing the dragCoefficient function ///

        const double currentDragCoefficient = Drag::dragCoefficient(machNumber,MAV.dragCoefficientPolyCoefficients(),MAV.dragCoefficientMachRanges());

        std::cout<<"The current drag coefficient = "<<currentDragCoefficient<<std::endl;

        // Drag coefficient function as part of the drag function
        const double currentDragCoefficient = Drag::ascentDragForce(rotationalVelocity,currentTemperature,Mars.adiabeticIndex(),Mars.specificGasConstant(),
                                                                    MAV.dragCoefficientPolyCoefficients(),MAV.dragCoefficientMachRanges());


        std::cout<<"The current drag coefficient = "<<currentDragCoefficient<<std::endl;
    //*/

        const double currentDrag = Drag::ascentDragForce(rotationalVelocity,currentTemperature,Mars.adiabeticIndex(),Mars.specificGasConstant(),
                                                         MAV.dragCoefficientPolyCoefficients(),MAV.dragCoefficientMachRanges(),
                                                         MAV.referenceArea(), currentDensity);

    //    std::cout<<"currentDrag = "<<currentDrag<<std::endl;



        /// Thrust acceleration in B-frame ///   thrustAccelerationsBframe

        const Eigen::Vector3d thrustAccelerationsPframe = Eigen::Vector3d((MAV.Thrust()/massMAV),0,0);            // THIS HAS TO BE CHANGED IN THE FUTURE TO INCLUDE A WIDE RANGE OF THRUST AZIMUTH AND ELEVATION ANGLES!!!

        const double thrustAzimuthTestDeg = 0;             // thrust azimuth gimbal angle [Deg] 10 for testing
        const double thrustElevationTestDeg = 0;            // thrust elevation gimbal angle [Deg] 5 for testing

        const double thrustAzimuthTest = deg2rad(thrustAzimuthTestDeg);     // thrust azimuth gimbal angle [rad]
        const double thrustElevationTest = deg2rad(thrustElevationTestDeg); // thrust elevation gimbal angle [rad]


        const Eigen::Vector3d thrustAccelerationsBframe = getPropulsionToBodyFrameTransformationMatrix(thrustAzimuthTest,thrustElevationTest)*thrustAccelerationsPframe;

//        std::cout<<"The thrust accelerations in the B-frame are "<<thrustAccelerationsBframe<<std::endl;

        /// Drag acceleration in B-frame ///

        const Eigen::Vector3d dragAccelerationsBframe = Eigen::Vector3d((-currentDrag/massMAV),0,0);

//        std::cout<<"The drag accelerations in the B-frame are "<<dragAccelerationsBframe<<std::endl;

    /// Transfer to the inertial frame ///
    ///
        /// Accelerations in the V-frame ///

        // Thrust accelerations from B-frame to V-frame
        const Eigen::Vector3d thrustAccelerationsVframe = tudat::reference_frames::getTrajectoryToLocalVerticalFrameTransformationMatrix(rotationalFlightPathAngle,rotationalAzimuth)*thrustAccelerationsBframe;

        // Drag accelerations from B-frame to V-frame
        const Eigen::Vector3d dragAccelerationsVframe = tudat::reference_frames::getTrajectoryToLocalVerticalFrameTransformationMatrix(rotationalFlightPathAngle,rotationalAzimuth)*dragAccelerationsBframe;

    //    std::cout<<"The thrust accelerations in the V-frame are "<<thrustAccelerationsVframe<<std::endl;
    //    std::cout<<"The drag accelerations in the V-frame are "<<dragAccelerationsVframe<<std::endl;


        /// Accelerations in the R-frame ///

        // Thrust accelerations from V-frame to R-frame
        const Eigen::Vector3d thrustAccelerationsRframe = tudat::reference_frames::getLocalVerticalToRotatingPlanetocentricFrameTransformationMatrix(rotationalLongitude,Latitude)*thrustAccelerationsVframe;

        // Drag accelerations from V-frame to R-frame
        const Eigen::Vector3d dragAccelerationsRframe = tudat::reference_frames::getLocalVerticalToRotatingPlanetocentricFrameTransformationMatrix(rotationalLongitude,Latitude)*dragAccelerationsVframe;

    //    std::cout<<"The thrust accelerations in the R-frame are "<<thrustAccelerationsRframe<<std::endl;
    //    std::cout<<"The drag accelerations in the R-frame are "<<dragAccelerationsRframe<<std::endl;

        /// Accelerations in the I-frame ///

        const double angleItoR = rotationalVelocityMars*(inertialFrameTime+currentTime)-primeMeridianAngle;

        // Thrust accelerations from R-frame to I-frame
        const Eigen::Vector3d thrustAccelerationsIframe = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(angleItoR)*thrustAccelerationsRframe;

        // Drag accelerations from R-frame to I-frame
        const Eigen::Vector3d dragAccelerationsIframe = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(angleItoR)*dragAccelerationsRframe;

    //    std::cout<<"The thrust accelerations in the I-frame are "<<thrustAccelerationsIframe<<std::endl;
    //    std::cout<<"The drag accelerations in the I-frame are "<<dragAccelerationsIframe<<std::endl;

    //    std::cout<<"The thrust+drag accelerations in the I-frame are "<<thrustAccelerationsIframe+dragAccelerationsIframe<<std::endl;


    /// Compute gravitational acceleration ///

        Eigen::Vector3d gravAccelerationsIframe_;        // Define the placeholder gravitational acceleration vector

        for (int i = 0; i<3;i++){

            gravAccelerationsIframe_(i)=-Mars.standardGravitationalParameter()*(currentState(i)/(pow(Radius,3)));

    //        std::cout<<"gravAccelerationsIframe_("<<i<<") = "<<gravAccelerationsIframe_(i)<<std::endl;
        };

        const Eigen::Vector3d gravAccelerationsIframe = gravAccelerationsIframe_;         // Actual gravitational acceleration vector


    //    std::cout<<"The gravitational accelerations in the I-frame are "<<gravAccelerationsIframe<<std::endl;

    /// Compute total acceleration ///

        const Eigen::Vector3d totalAccelerationsIframe = gravAccelerationsIframe+dragAccelerationsIframe+thrustAccelerationsIframe;



    /// Compute the mass flow rate ///

        const double massFlowRate = -MAV.Thrust()/(tudat::physical_constants::STANDARD_EARTH_GRAVITATIONAL_ACCELERATION*MAV.specificImpulse());



//        std::cout<<"The velocities in the I-frame are "<<stateAndTime.getCurrentVelocity()<<std::endl;
//        std::cout<<"The total accelerations in the I-frame are "<<totalAccelerationsIframe<<std::endl;
//        std::cout<<"The mass flow rate = "<<massFlowRate<<std::endl;


    /// Define the output vector and fill it ///


        tudat::basic_mathematics::Vector7d stateDerivativeVector;

        stateDerivativeVector(0) = xVelocity;
        stateDerivativeVector(1) = yVelocity;
        stateDerivativeVector(2) = zVelocity;
        stateDerivativeVector(3) = totalAccelerationsIframe(0);
        stateDerivativeVector(4) = totalAccelerationsIframe(1);
        stateDerivativeVector(5) = totalAccelerationsIframe(2);
        stateDerivativeVector(6) = massFlowRate;

        return stateDerivativeVector;
        } // end of the state derivative function


private:

    // Create variables to be used in this function

        celestialBody Mars;                               // The celestial body class       // This initilizes a new Mars, but later replaces it by the provided planet
        MarsAscentVehicle MAV;                                 // The MAV class             // This initializes a new MAV, but later replaces it by the provided MAV

        double rotationalVelocityMars;                  // The rotational velocity of Mars [rad/s]
        double primeMeridianAngle;                      // The angle of the prime Meridian of Mars at time of inertial frame set [rad]
        double inertialFrameTime;                       // The time at inertial frame set [s]


}; // end of class

#endif // ASCENTSTATEDERIVATIVEFUNCTIONCLASS_H