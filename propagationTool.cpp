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

    /// Input file ///

//    std::string nameOfFile = "test.txt"; // Test file
    std::string nameOfFile = "test1Ryan.txt"; // Test file from Ryan
//    std::string nameOfFile = "test2Joel.txt"; // Test file from Joel



/// Read the input file ///

    // This requires the complete path in order to work!
    std::string pathToWorkingFile = "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/07.InputFiles/";

    std::string completePathToInput = pathToWorkingFile+nameOfFile;

//    std::cout<<"completePathToInput = "<<completePathToInput<<std::endl;

    // This requires the complete path in order to work!
//    std::ifstream inputFile("/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/07.InputFiles/test.txt");
    std::ifstream inputFile(completePathToInput);


    // Read the first line
    std::string firstLine, variable, equalSign;
    std::getline(inputFile,firstLine);
    std::cout<<firstLine<<std::endl;

    double number;

    Eigen::VectorXd inputVectorValues = Eigen::VectorXd::Zero(1); // Create the first vector of the input matrix with all the values


    if (inputFile.is_open())
    {

        int line = 0; // Line to print in vector
        while(inputFile >> variable >> equalSign >> number){
//            std::cout<<variable<<" "<<equalSign<<" "<<number<<std::endl;

            inputVectorValues.conservativeResize(line+1); // This resizes the vector and makes it bigger to include all the values in the input file

            inputVectorValues(line) = number;

            line++;
        }



//        std::cout<<"inputVectorValues"<<'\n'<<inputVectorValues<<std::endl;

        inputFile.close();
    }
    else {std::cout<<"It didn't open..."<<std::endl;}

/// Set input values ///

    const double desiredOrbitalAltitude =   inputVectorValues(0);
    const double desiredInclinationDeg =    inputVectorValues(1);
    const double initialAltitude =          inputVectorValues(2);
    const double initialLatitudeDeg =       inputVectorValues(3);
    const double initialLongitudeDeg =      inputVectorValues(4);

    const double FlightPathAngleDeg =       inputVectorValues(5);
//    const double FlightPathAngleDeg = 89;

    const double HeadingAngleDeg =          inputVectorValues(6);
    const double initialGroundVelocity =    inputVectorValues(7);
    const double massMAV =                  inputVectorValues(8);
    const double thrust =                   inputVectorValues(9);
    const double specificImpulse =          inputVectorValues(10);
    const double initialBurnTime =          inputVectorValues(11);
    const double constantThrustElevationAngle = inputVectorValues(12);
    const double constantThrustAzimuthAngle = inputVectorValues(13);

//    const int maxOrder =                    inputVectorValues(14);
    const int maxOrder = 20;

//    const double chosenLocalErrorTolerance = inputVectorValues(15);
    const double chosenLocalErrorTolerance = 1e-15;

    const double chosenStepSize =           inputVectorValues(16);


//    const double setEndTime =               inputVectorValues(17);
    const double setEndTime = 1796;

    const double RKFinitiaterTime =         inputVectorValues(18);
    const bool rotatingPlanet =             inputVectorValues(19);
    const bool GravityAcc =                 inputVectorValues(20);
    const bool ThrustAcc =                  inputVectorValues(21);
    const bool DragAcc =                    inputVectorValues(22);
    const bool comparison =                 inputVectorValues(23);

/////////////////////////////////////////////////////////////////////// Actual propagation ///////////////////////////////////////////////////////////////////////

/// Perform the integration ///




        Eigen::MatrixXd outputMatrix = performIntegration(desiredOrbitalAltitude,desiredInclinationDeg,initialAltitude,initialLatitudeDeg,initialLongitudeDeg,
                                                          FlightPathAngleDeg,HeadingAngleDeg,initialGroundVelocity,massMAV,thrust,specificImpulse,initialBurnTime,
                                                          constantThrustElevationAngle,constantThrustAzimuthAngle,maxOrder,
                                                          chosenLocalErrorTolerance,chosenStepSize,setEndTime,RKFinitiaterTime,rotatingPlanet,GravityAcc,ThrustAcc,DragAcc,comparison);

    std::cout<<"outputMatrix = "<<outputMatrix<<std::endl;


    return 0;
}


