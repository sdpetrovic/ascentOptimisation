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
 *      160919    S.D. Petrovic     File created
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
#include <stdio.h>

#include <quadmath.h>   // Quad precision...

#include <time.h>   // To determine the current computer time
#include <sys/time.h> // To determine the current computer time

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <string>
#include <cmath>

// Used for the RKF integrator
#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
//#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
//#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
//#include <boost/make_shared.hpp>
//#include <boost/shared_ptr.hpp>
#include <Tudat/Mathematics/NumericalIntegrators/euler.h>
//#include <boost/bind.hpp>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h>

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


// TSI
/// Testing the auxiliary equations ///
//#include <thesisProject/Auxiliary.h>                // Original test file
#include <thesisProject/AuxiliarySpherical.h>                // Spherical test file


/// Testing the basic recurrence relations ///
#include <thesisProject/projectLibraries/basicRecurrenceRelations.h>               // Original test file

/// Testing all recurrence relations ///
//#include <thesisProject/projectLibraries/allRecurrenceRelations.h>          // Original test file
#include <thesisProject/projectLibraries/allRecurrenceRelationsSpherical.h>          // Spherical test file

/// Testing the stepSize class ///
#include <thesisProject/StepSize.h>             // Original test file

/// Testing the actual Taylor Series integration fucntion ///
//#include <thesisProject/projectLibraries/TaylorSeriesIntegration.h>             // Original test file
#include <thesisProject/projectLibraries/TaylorSeriesIntegrationSpherical.h>             // Original test file

/// Testing the other required functions ///
#include <thesisProject/projectLibraries/otherRequiredFunctions.h>              // Original test file

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
//#include <thesisProject/projectLibraries/ascentStateDerivativeFunction.h>   // Original test file

/// Testing the ascentStateDerivativeFunctionClass ///
#include <thesisProject/ascentStateDerivativeFunctionClass.h>       // Adapted file from thesisProject/projectLibraries/ascentStateDerivativeFunction.h

// Circularisation
#include <tudat/Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>


// testing
#include <thesisProject/projectLibraries/trajectoryIntegration.h>  // Test the trajectory integration file



int main()

{





    //// Old stuff ////



    /// Setting the desired end orbit ///

    double desiredOrbitalAltitude = 320.0; // Desired orbital altitude in kmMOLA (320 km is default)
//    const double desiredEccentricity = 0.0; // Desired orbital eccentricity
    double desiredInclinationDeg = 45.0; // Desired orbital inclination [deg] (45 is default)





    //////////////////////////////////////////////////////// Test Cases ////////////////////////////////////////////////////////

//    Test case from Woolley 2015 (case 10 SSTO)
    std::cout<<"Test case 1: Woolley 2015 SSTO"<<std::endl;
    desiredOrbitalAltitude = 390.0; // Desired orbital altitude in kmMOLA (320 km is default)
    desiredInclinationDeg = 45.0; // Desired orbital inclination [deg] (45 is default)

    const double massMAV = 267.4;                      // Set the MAV GLOM in [kg]
    const double thrust = 3.56;                        // Set the MAV thrust in [kN]
    const double specificImpulse = 256;                // Set the MAV specific impulse [s]
    const double initialBurnTime = 142.5;       // Set the burn time from launch till coast [s]
//    const double burnOutAngle = deg2rad(6.0);  // Set the burn out angle (flight-path angle at end of first burn) [rad]
//    const double finalBurnOutMass = 60.7;       // Set the final burn out mass (empty mass + OS mass + excess propellant mass) [kg]
    const double initialLongitudeDeg = 74.5;     // Set the launch latitude in [deg] (tau)
    const double initialLatitudeDeg = 0.0;        // Starting latitude [deg] initial condition in (delta)
    const double HeadingAngleDeg = 90.0;  // Set the launch azimuth [deg] (psi)
    const double initialAltitude = -0.6;          // Starting altitude [km MOLA] initial condition is -0.6 km MOLA
    const double initialGroundVelocity = 0.00001;          // Starting velocity in km/s (is suppose to be 0.0...) 0.00001 default at initial step-size of 0.01 sec




//// Test case from Joel (hybrid) case7_3_2016_v33
// std::cout<<"Test case 2: Joel (hybrid) case7_3_2016_v33"<<std::endl;
//    desiredOrbitalAltitude = 3875.19000000064-bodyReferenceRadius; // Desired orbital altitude in kmMOLA (320 km is default)
//    desiredInclination = deg2rad(92.6899999999988); // Desired orbital inclination (45 is default)

//    MAV.setMAVmass(288.95985303149);                      // Set the MAV GLOM in [kg]
//    MAV.setThrust(6.01868886452604);                        // Set the MAV thrust in [kN]
//    MAV.setThrustResetValue(MAV.Thrust());      // Set the reset value equal to the original given thrust
//    MAV.setSpecificImpulse(315.9);                // Set the MAV specific impulse [s]
//    const double initialBurnTime = 99.361911852794;       // Set the burn time from launch till coast [s]
////    const double burnOutAngle = deg2rad(6.0);  // Set the burn out angle (flight-path angle at end of first burn) [rad]
////    const double finalBurnOutMass = 60.7;       // Set the final burn out mass (empty mass + OS mass + excess propellant mass) [kg]
//    const double initialLongitudeDeg = 90.0;     // Set the launch latitude in [deg] (tau)
//    const double initialLatitudeDeg = 0.0;        // Starting latitude [deg] initial condition in (delta)
//    const double HeadingAngle = deg2rad(90.0);  // Set the launch azimuth [rad] (psi)
//    Mars.setUpperAltitudeBound(desiredOrbitalAltitude); // Set the upper bound of the altitude [km] for the temperature
//    MAV.setUpdatedFinalAltitude(desiredOrbitalAltitude); // Set the upper bound of the altitude [km] for the thrust angles
//    const double initialAltitude = -0.6;          // Starting altitude [km MOLA] initial condition is -0.6 km MOLA
//    const double initialGroundVelocity = 0.000001;          // Starting velocity in km/s (is suppose to be 0.0...) 0.00001 default at initial step-size of 0.01 sec


//// Test case old verification
//    desiredOrbitalAltitude = 320.0; // Desired orbital altitude in kmMOLA (320 km is default)
//    desiredInclination = deg2rad(45.0); // Desired orbital inclination (45 is default)

////    MAV.setMAVmass(267.4);                      // Set the MAV GLOM in [kg]
////    MAV.setThrust(3.56);                        // Set the MAV thrust in [kN]
////    MAV.setThrustResetValue(MAV.Thrust());      // Set the reset value equal to the original given thrust
////    MAV.setSpecificImpulse(256);                // Set the MAV specific impulse [s]
//    const double initialBurnTime = 142.5;       // Set the burn time from launch till coast [s]
////    const double burnOutAngle = deg2rad(6.0);  // Set the burn out angle (flight-path angle at end of first burn) [rad]
////    const double finalBurnOutMass = 60.7;       // Set the final burn out mass (empty mass + OS mass + excess propellant mass) [kg]
//    const double initialLongitudeDeg = 0.0;     // Set the launch latitude in [deg] (tau)
//    const double initialLatitudeDeg = 0.0;        // Starting latitude [deg] initial condition in (delta)
//    const double HeadingAngle = deg2rad(90.0);  // Set the launch azimuth [rad] (psi)
////    Mars.setUpperAltitudeBound(desiredOrbitalAltitude); // Set the upper bound of the altitude [km] for the temperature
////    MAV.setUpdatedFinalAltitude(desiredOrbitalAltitude); // Set the upper bound of the altitude [km] for the thrust angles
//    const double initialAltitude = -0.6;          // Starting altitude [km MOLA] initial condition is -0.6 km MOLA
//    const double initialGroundVelocity = 0.0;          // Starting velocity in km/s (is suppose to be 0.0...) 0.00001 default at initial step-size of 0.01 sec



    //////////////////////////////////////////////////////// Test Cases ////////////////////////////////////////////////////////

    /// Comparison?
//    const bool comparison = true;

    /// Set initial flight path angle and heading angle
    const double FlightPathAngleDeg = 89.0;     // Set flight-path angle in deg --> Default = 90.0 deg


  /// Initial conditions /// a.k.a. control centre

//    const double initialBurnTime = 68.63; // Burn time from launch till coast
    const double setEndTime = 2000.0;  // Integration end time  // 77 sec for a remainder mass of about 100 kg  // 200 sec for free fall // 2000 for test cases
    const double RKFinitiaterTime = 1.0;    // Time that the RKF integrator is used for TSI initially

    /// TSI settings ///
    const int maxOrder = 20; // Eventually want order 20 (testing is 8)
    /// TSI settings ///

    /// Integration settings ///
    const double chosenLocalErrorTolerance = 1e-8;      // The chosen local error tolerance for TSI
    const double chosenStepSize = 0.01; // The chosen initial step-size for TSI



    Eigen::MatrixXd outputMatrix = performIntegration(desiredOrbitalAltitude,desiredInclinationDeg,initialAltitude,initialLatitudeDeg,initialLongitudeDeg,
                                                      FlightPathAngleDeg,HeadingAngleDeg,initialGroundVelocity,massMAV,thrust,specificImpulse,initialBurnTime,maxOrder,
                                                      chosenLocalErrorTolerance,chosenStepSize,setEndTime,RKFinitiaterTime);

    std::cout<<"outputMatrix = "<<outputMatrix<<std::endl;


    return 0;
}


