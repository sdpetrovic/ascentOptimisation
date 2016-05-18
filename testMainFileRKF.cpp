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
 *      160517    S.D. Petrovic     File created
 *      160518    S.D. Petrovic     Changed the state derivative function header/source file structure into a class. Now the input for the function itself is simply time and state.
 *                                  This way, it should be able to be used by the existing RKF integrator from the tudat library.
 *
 *    References
 *
 *    Notes
 *
 */

// This is a test main file to test the different class files, header/source files to see if any output is produced (verification)

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


/* Taylor Series:
 *
/// Testing the basic recurrence relations ///
#include <thesisProject/projectLibraries/basicRecurrenceRelations.h>               // Original test file

/// Testing all recurrence relations ///
#include <thesisProject/projectLibraries/allRecurrenceRelations.h>          // Original test file

/// Testing the stepSize class ///
#include <thesisProject/StepSize.h>             // Original test file

/// Testing the actual Taylor Series integration fucntion ///
#include <thesisProject/projectLibraries/TaylorSeriesIntegration.h>             // Original test file

//*/


// State Derivative Function

/// Testing the airTemperature functions ///
#include <thesisProject/projectLibraries/airTemperature.h>      // Original test file

/// Testing the airDensity function ///
#include <thesisProject/projectLibraries/airDensity.h>          // Original test file

/// Testing the ascentDragForce function ///
#include <thesisProject/projectLibraries/ascentDragForce.h>     // Original test file

/// Testing the dragCoefficient function ///
#include <thesisProject/projectLibraries/dragCoefficient.h>     // Original test file

/// Testing the ascentStateDerivativeFunction ///
#include <thesisProject/projectLibraries/ascentStateDerivativeFunction.h>   // Original test file

/// Testing the ascentStateDerivativeFunctionClass ///
#include <thesisProject/ascentStateDerivativeFunctionClass.h>       // Adapted file from thesisProject/projectLibraries/ascentStateDerivativeFunction.h

// testing




int main()

{

std::cout<<setprecision(15)<<"Setting output precision to 15"<<std::endl;

    /// Setting the Celestial Body class ///


//    // First test

//    celestialBody Mars;

    // Second test

//    const std::string planet = "Mars";
//    const std::string planet = "Venus";

    celestialBody Mars;

//    const double adiabeticIndex = Mars.adiabeticIndex();
//    const double specificGasConstant = Mars.specificGasConstant();
//    const double standardGravitationalParameter = Mars.standardGravitationalParameter();
    const double rotationalVelocityMars = Mars.rotationalVelocity();
    const double primeMeridianAngle = Mars.primeMeridianAngle();
    const double inertialFrameTime = Mars.inertialFrameTime();

   const double bodyReferenceRadius = Mars.bodyReferenceRadius();

//    celestialBody Mars(planet);
//    Mars.setPlanet(planet);  // Does not exist in the class anymore!


    /// Setting the vehicle class ///

    MarsAscentVehicle MAV;



  /// Initial conditions ///

    // Launch site characteristics

//    const double initialAltitude = -0.6e3;             // Starting altitude [m MOLA]
    const double initialAltitude = -0.6;                 // Starting altitude [km MOLA]
    const double initialLatitudeDeg = 21;               // Starting latitude [deg]
    const double initialLongitudeDeg = 74.5;            // Starting longitude [deg]

//    const double initialLatitude = initialLatitudeDeg*tudat::mathematical_constants::LONG_PI/180;       // Starting latitude [rad]
//    const double initialLongitude = initialLongitudeDeg*tudat::mathematical_constants::LONG_PI/180;     // Starting longitude [rad]

    const double initialLatitude = deg2rad(initialLatitudeDeg);       // Starting latitude [rad]
    const double initialLongitude = deg2rad(initialLongitudeDeg);     // Starting longitude [rad]

    const double initialRadius = bodyReferenceRadius+initialAltitude;               // Starting radius in km
//        const double initialRadius = bodyReferenceRadius+initialAltitude;               // Starting radius in m

        // Converting the initial spherical position to cartesian position using the standard convertSphericalToCartesian function of Tudat
        // Please note that this function requires the zenith angle as input which is pi/2-latitude!

    Eigen::Vector3d initialCartesianPositionRotationalFrame = Eigen::Vector3d::Zero(3);

      initialCartesianPositionRotationalFrame(0) = initialRadius*cos(initialLatitude)*cos(initialLongitude); // x_R
      initialCartesianPositionRotationalFrame(1) = initialRadius*cos(initialLatitude)*sin(initialLongitude); // y_R
      initialCartesianPositionRotationalFrame(2) = initialRadius*sin(initialLatitude); // z_R

    const Eigen::Vector3d initialCartesianPositionInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(rotationalVelocityMars*inertialFrameTime-primeMeridianAngle)*initialCartesianPositionRotationalFrame;


    // Compute initial velocity in y-direction as seen from the launch site in the inertial frame

    const Eigen::Vector3d initialVelocityLaunchSite = Eigen::Vector3d(0,(rotationalVelocityMars*initialRadius*cos(initialLatitude)),0);

    const Eigen::Vector3d initialVelocityInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(rotationalVelocityMars*inertialFrameTime-primeMeridianAngle+initialLongitude)*initialVelocityLaunchSite;

    /// Setting StateAndTime class using modified vector ///

    tudat::basic_mathematics::Vector7d aState;

    aState(0) = initialCartesianPositionInertialFrame(0);       // x,y,z position in [km]
    aState(1) = initialCartesianPositionInertialFrame(1);
    aState(2) = initialCartesianPositionInertialFrame(2);
    aState(3) = initialVelocityInertialFrame(0);                // x,y,z velocity in [km/s]
    aState(4) = initialVelocityInertialFrame(1);
    aState(5) = initialVelocityInertialFrame(2);
    aState(6) = 227;  // Mass [kg] from literature study

    StateAndTime stateAndTime(aState);        // Creating the current state class using the namespace and class directly





///////////////////////////////////////////////////////////////////////////////////
////////////////////// Testing the State Derivative Function //////////////////////
///////////////////////////////////////////////////////////////////////////////////


    /// Testing with the state derivative function header and source file
    const tudat::basic_mathematics::Vector7d stateDerivativeVector = ascentStateDerivativeFunction(Mars,MAV,stateAndTime);

    std::cout<<"The state derivative vector is "<<stateDerivativeVector<<std::endl;


    /// Testing with the state derivative function class

    // Initial test
//    ascentStateDerivativeFunctionClass stateDerivativeFunctionClass;    // Initialize the class

//    const tudat::basic_mathematics::Vector7d stateDerivativeClassVector = stateDerivativeFunctionClass.ascentStateDerivativeFunction(Mars,MAV,stateAndTime);

//    std::cout<<"The state derivative vector obtained from the class is "<<stateDerivativeClassVector<<std::endl;

    // Full complete test

    std::cout<<"It should not be re-initilizing the class, especially here..."<<std::endl;

    ascentStateDerivativeFunctionClass stateDerivativeFunctionClass(Mars,MAV);     // Initialize the class

    const tudat::basic_mathematics::Vector7d stateDerivativeClassVector = stateDerivativeFunctionClass.ascentStateDerivativeFunction(stateAndTime.getCurrentTime(),stateAndTime.getCurrentState());

    std::cout<<"The state derivative vector obtained from the class is "<<stateDerivativeClassVector<<std::endl;


///// Testing the implementation in the integrator ///

//    const double initialStepSize = 0.2;     // Using the same initial step-size as defined for TSI





    return 0;
}


