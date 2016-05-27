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
 *      160519    S.D. Petrovic     File created
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

//#include <thesisProject/celestialbody.h>          // Original test file
//#include <thesisProject/celestialBodyTestOne.h>     // First test
//#include <thesisProject/celestialBodyTestTwo.h>     // Second test
//#include <thesisProject/celestialBodyTestThree.h>   // Third test
//#include <thesisProject/celestialBodyTestFour.h>   // Fourth test

#include <thesisProject/celestialBody.h>            // Final version

/// Testing the vehicle class ///

//#include <thesisProject/MarsAscentVehicleTestOne.h> // First Test
//#include <thesisProject/MarsAscentVehicleTestTwo.h> // Second Test

#include <thesisProject/MarsAscentVehicle.h>    // Final version

/// Testing the current state and time and its updater ///

//#include <thesisProject/StateAndTime.h>           // Original test file
//#include <thesisProject/stateAndTimeTestTwo.h>      // Second test file
//#include <thesisProject/stateAndTimeTestThree.h>    // Third test file

#include <thesisProject/stateAndTime.h>             // Final version


// TSI
/// Testing the auxiliary equations ///
#include <thesisProject/Auxiliary.h>                // Original test file


/// Testing the basic recurrence relations ///
#include <thesisProject/projectLibraries/basicRecurrenceRelations.h>               // Original test file

/// Testing all recurrence relations ///
#include <thesisProject/projectLibraries/allRecurrenceRelations.h>          // Original test file

/// Testing the stepSize class ///
#include <thesisProject/StepSize.h>             // Original test file

/// Testing the actual Taylor Series integration fucntion ///
#include <thesisProject/projectLibraries/TaylorSeriesIntegration.h>             // Original test file

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

    // No Gravity

    Mars.setStandardGravitationalParameter(0);

    if (Mars.standardGravitationalParameter() == 0){
        std::cout<<"NO GRAVITY"<<std::endl;
    }

//    const double adiabeticIndex = Mars.adiabeticIndex();
//    const double specificGasConstant = Mars.specificGasConstant();
//    const double standardGravitationalParameter = Mars.standardGravitationalParameter();
    const double rotationalVelocity = Mars.rotationalVelocity();
    const double primeMeridianAngle = Mars.primeMeridianAngle();
    const double inertialFrameTime = Mars.inertialFrameTime();

   const double bodyReferenceRadius = Mars.bodyReferenceRadius();




    /// Setting the vehicle class ///

    MarsAscentVehicle MAV;

    // No Drag

    MAV.setReferenceArea(0);

    if (MAV.referenceArea() == 0){
        std::cout<<"NO DRAG"<<std::endl;
    }

    // No Thrust

    MAV.setThrust(0);

    if (MAV.Thrust() == 0){
        std::cout<<"NO THRUST"<<std::endl;
    }

  /// Initial conditions ///


    /// TSI settings ///
    const int maxOrder = 20; // Eventually want order 20 (testing is 8)
    const double chosenLocalErrorTolerance = 1e-8;      // The chosen local error tolerance for TSI
    /// TSI settings ///


    // Launch site characteristics

//    const double initialAltitude = -0.6e3;             // Starting altitude [m MOLA]
    const double initialAltitude = -0.6;                 // Starting altitude [km MOLA] initial condition is -0.6 km MOLA
    std::cout<<"The initial altitude = "<<initialAltitude<<std::endl;
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

    const Eigen::Vector3d initialCartesianPositionInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(rotationalVelocity*inertialFrameTime-primeMeridianAngle)*initialCartesianPositionRotationalFrame;


    // Compute initial velocity in y-direction as seen from the launch site in the inertial frame

    const Eigen::Vector3d initialVelocityLaunchSite = Eigen::Vector3d(0,(rotationalVelocity*initialRadius*cos(initialLatitude)),0);

    const Eigen::Vector3d initialVelocityInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(rotationalVelocity*inertialFrameTime-primeMeridianAngle+initialLongitude)*initialVelocityLaunchSite;

    /// Setting StateAndTime class using modified vector ///

    tudat::basic_mathematics::Vector7d aState;

    aState(0) = initialCartesianPositionInertialFrame(0);
    aState(1) = initialCartesianPositionInertialFrame(1);
    aState(2) = initialCartesianPositionInertialFrame(2);
    aState(3) = initialVelocityInertialFrame(0);
    aState(4) = initialVelocityInertialFrame(1);
    aState(5) = initialVelocityInertialFrame(2);
    aState(6) = 227;  // Mass [kg] from literature study

//    /// Debug ///

//        aState(0) = 852.252774466749;
//        aState(1) = 3073.12422474535;
//        aState(2) = 1224.18491564675;
//        aState(3) = -0.21782304504995;
//        aState(4) = 0.0604076766542031;
//        aState(5) = 0;
//        aState(6) = 227;  // Mass [kg] from literature study

//    /// Debug //

    StateAndTime currentStateAndTime(aState);        // Creating the current state class using the namespace and class directly

    std::cout<<"aState = "<<aState<<std::endl;
//    std::cout<<"x1-852.252774466749 = "<<aState(0)-852.252774466749<<std::endl;
//    std::cout<<"x1-850 = "<<aState(0)-850<<std::endl;
//    std::cout<<"x2-3073.12422474535 = "<<aState(1)-3073.12422474535<<std::endl;
//    std::cout<<"x3-1224.18491564675 = "<<aState(2)-1224.18491564675<<std::endl;
//    std::cout<<"x4+0.21782304504995 = "<<aState(3)+0.21782304504995<<std::endl;
//    std::cout<<"x5-0.0604076766542031 = "<<aState(4)-0.0604076766542031<<std::endl;
//    std::cout<<"x6-0 = "<<aState(5)-0<<std::endl;
//    std::cout<<"x7-227 = "<<aState(6)-227<<std::endl;


/////////////////////////////////////////////////////////////////////
////////////////////// Testing the integrators //////////////////////
/////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////
/*////////////////////// Testing the Taylor series integrator //////////////////////
//////////////////////////////////////////////////////////////////////////////////

/// Setting the data collection file for TSI and inserting the first values ///

    // Set directory where output files will be stored. THIS REQUIRES THE COMPLETE PATH IN ORDER TO WORK!!
    const std::string outputDirectory = "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/01.integrationResults/TSI/";


    // Set output format for matrix output.
    Eigen::IOFormat csvFormat( 15, 0, ", ", "\n" );

    // Set absolute path to file containing the Taylor Series Coefficients.
    std::string dataAbsolutePath = outputDirectory + "test6FullIntegrationTSIstateAndTime.csv";

    // Create a row vector for the storing of the data
    Eigen::MatrixXd outputVector = Eigen::MatrixXd::Zero(1,8); // Create a row vector for the storing of the data

    // Getting the initial conditions for storage
    const tudat::basic_mathematics::Vector7d initialState = currentStateAndTime.getCurrentState();

    // Filling the output vector
    outputVector(0,0) = currentStateAndTime.getCurrentTime();   // Storing the initial time
    outputVector(0,1) = initialState(0);   // Storing the initial x position
    outputVector(0,2) = initialState(1);   // Storing the initial y position
    outputVector(0,3) = initialState(2);   // Storing the initial z position
    outputVector(0,4) = initialState(3);   // Storing the initial x velocity
    outputVector(0,5) = initialState(4);   // Storing the initial y velocity
    outputVector(0,6) = initialState(5);   // Storing the initial z velocity
    outputVector(0,7) = initialState(6);   // Storing the initial MAV mass

    // Storing the data

    std::ifstream ifile(dataAbsolutePath.c_str()); // Check it as an input file

    bool fexists = false;   // Set the default to "It does not exist"

    if (ifile){         // Attempt to open the file


       fexists = true;      // If the file can be opened it must exist

       ifile.close();   // Close the file

    }


    // If so: error and create temporary file, if not: create new file and put data in

    if (fexists == true){


        /// Get time ///

        time_t rawtime;
        struct tm * timeinfo;

        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
//        printf ( "Current local time and date: %s", asctime (timeinfo) );         // (from stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c and from stackoverflow.com/questions/1442116/how-to-get-date-and-time-value-in-c-program)
//        std::cout<<"time = "<<time ( &rawtime )<<std::endl;
//        std::cout<< (timeinfo->tm_year+1900) << "-"
//                 << (timeinfo->tm_mon + 1) << "-"
//                 << timeinfo->tm_mday << " "
//                 << timeinfo->tm_hour << ":"
//                 << timeinfo->tm_min << ":"
//                 << timeinfo->tm_sec <<std::endl;

        ostringstream ConvertYear;
        ConvertYear << (timeinfo->tm_year+1900);
        std::string currentYear = ConvertYear.str();

        ostringstream ConvertMonth;
        ConvertMonth << (timeinfo->tm_mon+1);
        std::string currentMonth = ConvertMonth.str();

        ostringstream ConvertDay;
        ConvertDay << timeinfo->tm_mday;
        std::string currentDay = ConvertDay.str();

        ostringstream ConvertHour;
        ConvertHour << timeinfo->tm_hour;
        std::string currentHour = ConvertHour.str();

//        std::cout<<"currentHour = "<<currentHour<<std::endl;

        ostringstream ConvertMin;
        ConvertMin << timeinfo->tm_min;
        std::string currentMin = ConvertMin.str();

//        std::cout<<"currentMin = "<<currentMin<<std::endl;

        ostringstream ConvertSec;
        ConvertSec << timeinfo->tm_sec;
        std::string currentSec = ConvertSec.str();

//        std::cout<<"currentSec = "<<currentSec<<std::endl;

        // Making sure each one of the representation has at two numbers
        if(currentMonth.size() == 1){
            currentMonth = "0" + currentMonth;
        }


        if(currentDay.size() == 1){
            currentDay = "0" + currentDay;
        }


        if(currentSec.size() == 1){
            currentSec = "0" + currentSec;
        }

        if (currentMin.size() == 1){
            currentMin = "0" + currentMin;
        }

        if (currentHour.size() == 1){
            currentHour = "0" + currentHour;
        }


//        std::cout<<"The length of currentSec = "<<currentSec.size()<<" and the value = "<<currentSec<<std::endl;

        std::string ComputerTimeString = currentYear + "-" + currentMonth + "-" + currentDay + "_" + currentHour + ":" + currentMin + ":" + currentSec;  // Convert to string and store

        std::string newFileName = "backupTSIFileAtDateAndTime_" + ComputerTimeString + ".csv";

    std::cerr<<"The file name that you have chosen already exists, a new file with name "<<newFileName<<" will be created to store the data for now"<<std::endl;

    // Set new absolute path to file containing the data.
    dataAbsolutePath = outputDirectory + newFileName;

    // Export the data.
    std::ofstream exportFile1( dataAbsolutePath.c_str( ) ); // Make the new file
    std::cout<<"New file called "<<dataAbsolutePath<<" has been created"<<std::endl;
    exportFile1 << outputVector.format( csvFormat );          // Store the new values
    exportFile1.close( );   // Close the file


}
        else{

        // Export the data.
        std::ofstream exportFile1( dataAbsolutePath.c_str( ) ); // Make the new file
        std::cout<<"New file called "<<dataAbsolutePath<<" has been created"<<std::endl;
        exportFile1 << outputVector.format( csvFormat );          // Store the new values
        exportFile1.close( );   // Close the file
    };



/// Defining the order and initializing the StepSize class ///


    std::cout<<"The order of TSI = "<<maxOrder<<std::endl;

        StepSize stepSize; // Initializing the stepSize class. THIS SHOULD BE DONE BEFORE THE START OF THE INTEGRATION!!!!!

        stepSize.setLocalErrorTolerance(chosenLocalErrorTolerance);  // Setting the local error tolerance to the lowest possible value in order to compare to RKF7(8) and the others

//        /// Debug ///

//        std::cout<<"StepSize has been initialized"<<std::endl;
//        stepSize.setCurrentStepSize(1);
//        std::cout<<"The current step-size = "<<stepSize.getCurrentStepSize()<<std::endl;
//        stepSize.setCurrentStepSize(5);
//        std::cout<<"The current step-size = "<<stepSize.getCurrentStepSize()<<std::endl;

//        /// Debug ///

//        const double currentStepSize = stepSize.getCurrentStepSize();

//        std::cout<<"The current stepSize = "<<currentStepSize<<std::endl;

/// Performing the actual TSI integration ///

        // Define storing matrix for the intermediate values
        Eigen::MatrixXd dataStoringMatrix(1,8); // The size of this matrix will change in the do-loop

        // Set the end time
        const double endTime = 10; // sec

//        std::cout<<"It works till here 1"<<std::endl;


/// The integeration do-loop ///

        // Set initial running time. This is updated after each step that the numerical integrator takes.
     double runningTime = 0.0;
     int count = 0;

//	do
//	{

//     std::cout<<"It works till here 2"<<std::endl;

     for (int i = 0; i<2; i++){
         /// Debug ///
    std::cout<<"The current step-size = "<<stepSize.getCurrentStepSize()<<std::endl;
    std::cout<<"The current runningTime = "<<runningTime<<std::endl;
    std::cout<<"std::fabs(endTime-runningTime) = "<<std::fabs(endTime-runningTime)<<std::endl;
    std::cout<<"std::fabs( stepSize.getCurrentStepSize() ) * ( 1.0 + std::numeric_limits< double >::epsilon( ) ) = "<<std::fabs( stepSize.getCurrentStepSize() ) * ( 1.0 + std::numeric_limits< double >::epsilon( ) )<<std::endl;
         /// Debug ///

	if ( std::fabs( endTime - runningTime )
                             <= std::fabs( stepSize.getCurrentStepSize() ) * ( 1.0 + std::numeric_limits< double >::epsilon( ) ) )
                        {
        std::cout<<"It is indeed smaller than the step-size"<<std::endl;
                            stepSize.setCurrentStepSize(endTime - runningTime);
                        }
        std::cout<<"The new step-size = "<<stepSize.getCurrentStepSize()<<std::endl;
        Eigen::VectorXd updatedStateAndTimeVector = performTaylorSeriesIntegrationStep(Mars, MAV, currentStateAndTime, stepSize, maxOrder); /// The actual integration step
        // This function has the output: updated position, updated velocity, updated mass and updated time

        std::cout<<"updatedStateAndTimeVector = "<<updatedStateAndTimeVector<<std::endl;


        // Check to see if the class has been updated from within the TSI function
//        const double nextStepSize = stepSize.getCurrentStepSize();

//        std::cout<<"The next stepSize = "<<nextStepSize<<std::endl;

/// Storing the values ///

        outputVector = Eigen::MatrixXd::Zero(1,8); // Setting the output vector to zero again to be sure.

        // Filling the output vector
        outputVector(0,0) = updatedStateAndTimeVector(7);   // Storing the updated time
        outputVector(0,1) = updatedStateAndTimeVector(0);   // Storing the updated x position
        outputVector(0,2) = updatedStateAndTimeVector(1);   // Storing the updated y position
        outputVector(0,3) = updatedStateAndTimeVector(2);   // Storing the updated z position
        outputVector(0,4) = updatedStateAndTimeVector(3);   // Storing the updated x velocity
        outputVector(0,5) = updatedStateAndTimeVector(4);   // Storing the updated y velocity
        outputVector(0,6) = updatedStateAndTimeVector(5);   // Storing the updated z velocity
        outputVector(0,7) = updatedStateAndTimeVector(6);   // Storing the updated MAV mass

        // Store the new values in the data storage matrix

        if (count == 0){

          dataStoringMatrix.row(count) = outputVector.row(0); // Filling the matrix
        }
        else{
            dataStoringMatrix.conservativeResize(count+1,8); // Making the matrix bigger in order to store more values

            dataStoringMatrix.row(count) = outputVector.row(0); // Filling the matrix
        }

        // Updating the current state and time class!!!
        tudat::basic_mathematics::Vector7d currentStateVector; // Create the current state and time vector

        // Fill the curent state and time vector
        currentStateVector(0) = updatedStateAndTimeVector(0);  // Updated x position
        currentStateVector(1) = updatedStateAndTimeVector(1);  // Updated y position
        currentStateVector(2) = updatedStateAndTimeVector(2);  // Updated z position
        currentStateVector(3) = updatedStateAndTimeVector(3);  // Updated x velocity
        currentStateVector(4) = updatedStateAndTimeVector(4);  // Updated y velocity
        currentStateVector(5) = updatedStateAndTimeVector(5);  // Updated z velocity
        currentStateVector(6) = updatedStateAndTimeVector(6);  // Updated MAV mass

        runningTime = updatedStateAndTimeVector(7);             // Updated time

        currentStateAndTime.setCurrentStateAndTime(currentStateVector,runningTime); // Update the current state and time class!



     count++;
     };

//    }while( !( endTime - runningTime <= std::numeric_limits< double >::epsilon( ) ) );

        /// Adding the values to the file ///

        // Check if the file already exists.


        std::ifstream ifile2(dataAbsolutePath.c_str()); // Check it as an input file

        fexists = false;   // Set the default to "It does not exist"

        if (ifile2){         // Attempt to open the file


           fexists = true;      // If the file can be opened it must exist

           ifile2.close();   // Close the file

        }


        // If so: append, if not: error

        if (fexists == true){

            // Export the Taylor Series Coefficients matrix.
            std::ofstream exportFile1;                          // Define the file as an output file


            exportFile1.open(dataAbsolutePath.c_str(),std::ios_base::app);      // Open the file in append mode

            exportFile1 << "\n";                                            // Make sure the new matrix start on a new line

            exportFile1 << dataStoringMatrix.format( csvFormat ); // Add the new values

            std::cout<<"The file called "<<dataAbsolutePath<<" has been appended"<<std::endl;


            exportFile1.close( );   // Close the file
}
            else{

            std::cerr<<"Error: values could not be stored because storage file does not exist"<<std::endl;
        };
//*/

//        std::cout<<"It works till here 3"<<std::endl;
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////// Testing the RKF and other higher order integrators //////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

                    /// Setting the data collection file for RKF and inserting the first values ///

                        // Set directory where output files will be stored. THIS REQUIRES THE COMPLETE PATH IN ORDER TO WORK!!
                        const std::string outputDirectory = "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/01.integrationResults/RKF/";


                        // Set output format for matrix output.
                        Eigen::IOFormat csvFormat( 15, 0, ", ", "\n" );

                        // Set absolute path to file containing the data.
                        std::string dataAbsolutePath = outputDirectory + "test4FullIntegrationRKF7(8)stateAndTime.csv";

                        // Create a row vector for the storing of the data
                        Eigen::MatrixXd outputVector = Eigen::MatrixXd::Zero(1,8); // Create a row vector for the storing of the data

                        // Getting the initial conditions for storage
                        const tudat::basic_mathematics::Vector7d initialState = currentStateAndTime.getCurrentState();

                        // Filling the output vector
                        outputVector(0,0) = currentStateAndTime.getCurrentTime();   // Storing the initial time
                        outputVector(0,1) = initialState(0);   // Storing the initial x position
                        outputVector(0,2) = initialState(1);   // Storing the initial y position
                        outputVector(0,3) = initialState(2);   // Storing the initial z position
                        outputVector(0,4) = initialState(3);   // Storing the initial x velocity
                        outputVector(0,5) = initialState(4);   // Storing the initial y velocity
                        outputVector(0,6) = initialState(5);   // Storing the initial z velocity
                        outputVector(0,7) = initialState(6);   // Storing the initial MAV mass


                        // Storing the data

                        std::ifstream ifile(dataAbsolutePath.c_str()); // Check it as an input file

                        bool fexists = false;   // Set the default to "It does not exist"

                        if (ifile){         // Attempt to open the file


                           fexists = true;      // If the file can be opened it must exist

                           ifile.close();   // Close the file

                        }


                        // If so: error and create temporary file, if not: create new file and put data in

                        if (fexists == true){


                            // Get time //

                            time_t rawtime;
                            struct tm * timeinfo;

                            time ( &rawtime );
                            timeinfo = localtime ( &rawtime );
                    //        printf ( "Current local time and date: %s", asctime (timeinfo) );         // (from stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c and from stackoverflow.com/questions/1442116/how-to-get-date-and-time-value-in-c-program)
                    //        std::cout<<"time = "<<time ( &rawtime )<<std::endl;
                    //        std::cout<< (timeinfo->tm_year+1900) << "-"
                    //                 << (timeinfo->tm_mon + 1) << "-"
                    //                 << timeinfo->tm_mday << " "
                    //                 << timeinfo->tm_hour << ":"
                    //                 << timeinfo->tm_min << ":"
                    //                 << timeinfo->tm_sec <<std::endl;

                            ostringstream ConvertYear;
                            ConvertYear << (timeinfo->tm_year+1900);
                            std::string currentYear = ConvertYear.str();

                            ostringstream ConvertMonth;
                            ConvertMonth << (timeinfo->tm_mon+1);
                            std::string currentMonth = ConvertMonth.str();

                            ostringstream ConvertDay;
                            ConvertDay << timeinfo->tm_mday;
                            std::string currentDay = ConvertDay.str();

                            ostringstream ConvertHour;
                            ConvertHour << timeinfo->tm_hour;
                            std::string currentHour = ConvertHour.str();

                    //        std::cout<<"currentHour = "<<currentHour<<std::endl;

                            ostringstream ConvertMin;
                            ConvertMin << timeinfo->tm_min;
                            std::string currentMin = ConvertMin.str();

                    //        std::cout<<"currentMin = "<<currentMin<<std::endl;

                            ostringstream ConvertSec;
                            ConvertSec << timeinfo->tm_sec;
                            std::string currentSec = ConvertSec.str();


                            // Making sure each one of the representation has at two numbers
                            if(currentMonth.size() == 1){
                                currentMonth = "0" + currentMonth;
                            }


                            if(currentDay.size() == 1){
                                currentDay = "0" + currentDay;
                            }

                            if(currentSec.size() == 1){
                                currentSec = "0" + currentSec;
                            }

                            if (currentMin.size() == 1){
                                currentMin = "0" + currentMin;
                            }

                            if (currentHour.size() == 1){
                                currentHour = "0" + currentHour;
                            }


                    //        std::cout<<"The length of currentSec = "<<currentSec.size()<<" and the value = "<<currentSec<<std::endl;

                            // Create new file name

                            std::string ComputerTimeString = currentYear + "-" + currentMonth + "-" + currentDay + "_" + currentHour + ":" + currentMin + ":" + currentSec;  // Convert to string and store

                            std::string newFileName = "backupRKFFileAtDateAndTime_" + ComputerTimeString + ".csv";

                        std::cerr<<"The file name that you have chosen already exists, a new file with name "<<newFileName<<" will be created to store the data for now"<<std::endl;

                        // Set new absolute path to file containing the data.
                        dataAbsolutePath = outputDirectory + newFileName;

                        // Export the data.
                        std::ofstream exportFile1( dataAbsolutePath.c_str( ) ); // Make the new file
                        std::cout<<"New file called "<<dataAbsolutePath<<" has been created"<<std::endl;
                        exportFile1 << outputVector.format( csvFormat );          // Store the new values
                        exportFile1.close( );   // Close the file


                    }
                            else{

                            // Export the data.
                            std::ofstream exportFile1( dataAbsolutePath.c_str( ) ); // Make the new file
                            std::cout<<"New file called "<<dataAbsolutePath<<" has been created"<<std::endl;
                            exportFile1 << outputVector.format( csvFormat );          // Store the new values
                            exportFile1.close( );   // Close the file
                        };



                    /// Testing with the state derivative function class


                    // Full complete test

                    ascentStateDerivativeFunctionClass stateDerivativeFunctionClass(Mars,MAV);     // Initialize the class

                    // Just using boost::bind

                    boost::function< tudat::basic_mathematics::Vector7d( const double, const tudat::basic_mathematics::Vector7d ) > stateDerivativeFunction // Using boost::function to create the passing function
                            =boost::bind( &ascentStateDerivativeFunctionClass::ascentStateDerivativeFunction, &stateDerivativeFunctionClass, _1, _2);       // Then using boost::bind to bind the function as per class stateDerivativeFunctionClass which requires two inputs:
                                                                                                                                                            // _1 which is the current time and _2 which is the current state


                 ///// Testing the implementation in the integrator ///

                    // Initial conditions.
                    const double initialTime = currentStateAndTime.getCurrentTime();                            // Time.
                //    Eigen::VectorXd initialState = currentStateAndTime.getCurrentState(); // State: start with zero velocity at the origin.

                    const double endTime = 1000;     // Using the same initial step-size as defined for TSI

                    // Step-size settings.
                    // The minimum and maximum step-size are set such that the input data is fully accepted by the
                    // integrator, to determine the steps to be taken.
                    const double zeroMinimumStepSize = std::numeric_limits< double >::epsilon( );
                    const double infiniteMaximumStepSize = std::numeric_limits< double >::infinity( );
                    double stepSize = 0.2;          // Using the same initial step-size as defined for TSI

                    // Tolerances.
                    const double relativeTolerance = 1e-8;     // 1e-14 is used by TSI, original setting was 1e-15
                    const double absoluteTolerance = 1e-8;     // 1e-14 is used by TSI, original setting was 1e-15

                    // For RKF7(8) a step-size of 0.2 is only used if the tolerances are 1e-3.... and is accepted till 1e-8
                    // For RKF4(5) a step-size of 0.2 is only used if the tolerances are 1e-7.... and is accepted till 1e-10
                    // For DP8(7) a step-size of 0.2 is only used if the tolerances are 1e-7.... and is accepted till 1e-9

//                    /// RungeKutta4 numerical integrator.

//                    tudat::numerical_integrators::RungeKutta4IntegratorXd RK4integrator(
//                        stateDerivativeFunction, initialTime, initialState );

//                    // Integrate to the specified end time.
//                    Eigen::VectorXd RK4endState = RK4integrator.integrateTo( endTime, stepSize );

//                    std::cout<<"The RK4 end state is "<<RK4endState<<std::endl;

                    /// Runge-Kutta-Fehlberg 7(8) integrator.
                       tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integrator(
                                   tudat::numerical_integrators::RungeKuttaCoefficients::get(
                                       tudat::numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg78),
                                   stateDerivativeFunction, initialTime, initialState, zeroMinimumStepSize,
                                   infiniteMaximumStepSize, relativeTolerance, absoluteTolerance );

                //       /// Runge-Kutta-Fehlberg 4(5) integrator.
                //          tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integrator(
                //                      tudat::numerical_integrators::RungeKuttaCoefficients::get(
                //                          tudat::numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg45),
                //                      stateDerivativeFunction, initialTime, initialState, zeroMinimumStepSize,
                //                      infiniteMaximumStepSize, relativeTolerance, absoluteTolerance );

                //          /// Dormand-Prince 8(7) integrator.
                //             tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integrator(
                //                         tudat::numerical_integrators::RungeKuttaCoefficients::get(
                //                             tudat::numerical_integrators::RungeKuttaCoefficients::rungeKutta87DormandPrince),
                //                         stateDerivativeFunction, initialTime, initialState, zeroMinimumStepSize,
                //                         infiniteMaximumStepSize, relativeTolerance, absoluteTolerance );



                       // Set initial running time. This is updated after each step that the numerical integrator takes.
                    double runningTime = 0.0;
                    int count = 0;


                    // Define storing matrix for the intermediate values
                    Eigen::MatrixXd dataStoringMatrix(1,8); // The size of this matrix will change in the do-loop

                    do
                    {
                        // Make sure the integrator does not integrate beyond the end time.



                        if ( std::fabs( endTime - runningTime )
                             <= std::fabs( stepSize ) * ( 1.0 + std::numeric_limits< double >::epsilon( ) ) )
                        {
                            stepSize = endTime - integrator.getCurrentIndependentVariable( );
                        }

//                        double prevStepSize = stepSize;

                        // Perform a single integration step. Then update the step-size and running time.
                        integrator.performIntegrationStep( stepSize );
                        stepSize = integrator.getNextStepSize( );
                        runningTime = integrator.getCurrentIndependentVariable( );

                        Eigen::VectorXd currentState = integrator.getCurrentState();

//                        std::cout<<"The current stepSize is "<<prevStepSize<<" s"<<std::endl;
                        std::cout<<"The current running time is "<<runningTime<<std::endl;

                        /// Debug ///

                        const double rotationalVelocityMars = Mars.rotationalVelocity();
                        const double primeMeridianAngle = Mars.primeMeridianAngle();
                        const double inertialFrameTime = Mars.inertialFrameTime();

                        const double xPosition = currentState(0);            // x position coordinate definition
                        const double yPosition = currentState(1);            // y position coordinate definition
                        const double zPosition = currentState(2);            // z position coordinate definition
                        const double xVelocity = currentState(3);            // x velocity coordinate definition
                        const double yVelocity = currentState(4);            // y velocity coordinate definition
                        const double zVelocity = currentState(5);            // z velocity coordinate definition
                        const double massMAV = currentState(6);              // MAV mass definition
                        const double currentTime = runningTime;   // current time definition

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

                        // Check output
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

                        /// Debug ///

                        /// Storing the values ///

                                outputVector = Eigen::MatrixXd::Zero(1,8); // Setting the output vector to zero again to be sure.


                                // Filling the output vector
                                outputVector(0,0) = integrator.getCurrentIndependentVariable();   // Storing the updated time
                                outputVector(0,1) = currentState(0);   // Storing the updated x position
                                outputVector(0,2) = currentState(1);   // Storing the updated y position
                                outputVector(0,3) = currentState(2);   // Storing the updated z position
                                outputVector(0,4) = currentState(3);   // Storing the updated x velocity
                                outputVector(0,5) = currentState(4);   // Storing the updated y velocity
                                outputVector(0,6) = currentState(5);   // Storing the updated z velocity
                                outputVector(0,7) = currentState(6);   // Storing the updated MAV mass

                        // Store the new values in the data storage matrix

                        if (count == 0){

                          dataStoringMatrix.row(count) = outputVector.row(0); // Filling the matrix
                        }
                        else{
                            dataStoringMatrix.conservativeResize(count+1,8); // Making the matrix bigger in order to store more values

                            dataStoringMatrix.row(count) = outputVector.row(0); // Filling the matrix
                        }


                        count++;

                    }while( !( endTime - runningTime <= std::numeric_limits< double >::epsilon( ) ) );

                    /// Storing the values to the file ///


                    // Check if the file already exists.


                    std::ifstream ifile2(dataAbsolutePath.c_str()); // Check it as an input file

                    fexists = false;   // Set the default to "It does not exist"

                    if (ifile2){         // Attempt to open the file


                       fexists = true;      // If the file can be opened it must exist

                       ifile2.close();   // Close the file

                    }


                    // If so: append, if not: create new file and put data in

                    if (fexists == true){

                        // Export the Taylor Series Coefficients matrix.
                        std::ofstream exportFile1;                          // Define the file as an output file


                        exportFile1.open(dataAbsolutePath.c_str(),std::ios_base::app);      // Open the file in append mode

                        exportFile1 << "\n";                                            // Make sure the new matrix start on a new line

                        exportFile1 << dataStoringMatrix.format( csvFormat ); // Add the new values

                        std::cout<<"The file called "<<dataAbsolutePath<<" has been appended"<<std::endl;


                        exportFile1.close( );   // Close the file
            }
                        else{

                        std::cerr<<"Error: values could not be stored because storage file does not exist"<<std::endl;
                    };

                    // The result of the integration.
                    Eigen::VectorXd endState = integrator.getCurrentState( );

                    std::cout<<"Final number of integration steps is "<<count<<std::endl;
                    std::cout<<"The end state is "<<endState<<std::endl;

                    //*/

    return 0;
}

