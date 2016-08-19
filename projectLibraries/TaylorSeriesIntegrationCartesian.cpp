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
 *      160505    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */

#include "TaylorSeriesIntegrationCartesian.h"



/* ///// Some required functions ///// (Moved to separate header and source file)

/// deg2rad and rad2deg ///

const double deg2rad(const double deg){

    const double rad = deg*tudat::mathematical_constants::LONG_PI/180;

    return rad;
}

const double rad2deg(const double rad){

    const double deg = rad*180/tudat::mathematical_constants::LONG_PI;

    return deg;
}


/// B-P frame transformations ///



//! Get transformation quaternion from the Body (B) to the Propulsion (P) frame.
Eigen::Quaterniond getBodyToPropulsionFrameTransformationQuaternion(
        const double thrustAzimuth, const double thrustElevation )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd RotationAroundZaxis = Eigen::AngleAxisd(
                -1.0 * thrustAzimuth, Eigen::Vector3d::UnitZ( ) );
    Eigen::AngleAxisd RotationAroundYaxis = Eigen::AngleAxisd(
                -1.0 * thrustElevation,
                Eigen::Vector3d::UnitY( ) );
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
                ( RotationAroundYaxis * RotationAroundZaxis ) );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}


//! Get transformation matrix from the Body (B) to the Propulsion (P) frame.
Eigen::Matrix3d getBodyToPropulsionFrameTransformationMatrix(
    const double thrustAzimuth, const double thrustElevation )
{
    return getBodyToPropulsionFrameTransformationQuaternion(
            thrustAzimuth, thrustElevation ).toRotationMatrix( );
}

//! Get transformation matrix from the Propulsion (P) to the Body (B) frame.
Eigen::Matrix3d getPropulsionToBodyFrameTransformationMatrix(
    const double thrustAzimuth, const double thrustElevation )
{
    return getBodyToPropulsionFrameTransformationMatrix(
            thrustAzimuth, thrustElevation ).transpose( );
}

//! Get transformation quaternion from the Propulsion (P) to the Body (B) frame.
Eigen::Quaterniond getPropulsionToBodyFrameTransformationQuaternion(
        const double thrustAzimuth, const double thrustElevation )
{
    return getBodyToPropulsionFrameTransformationQuaternion(
            thrustAzimuth, thrustElevation ).inverse( );
}
//*/

/// Taylor series integration step ///


Eigen::VectorXd performCartesianTaylorSeriesIntegrationStep(const celestialBody& planet_, const MarsAscentVehicle& MAV_, const StateAndTime& currentStateAndTime_, StepSize& stepSize, const double maxOrder_ ,
                                                   const double FlightPathAngle_,
                                                   const double HeadingAngle_){

    // The stepSize class is the only class that will be updated directly from this function!

        // Create variables to be used in this function

    celestialBody planet = planet_;                               // The celestial body class
    MarsAscentVehicle MAV = MAV_;                                 // The MAV class
    StateAndTime currentStateAndTime = currentStateAndTime_;      // The current state and time class
    const double maxOrder = maxOrder_;                                  // The maximum order K (default value is 20)

    // Getting the required information from the classes

    // celestialBody
    const double adiabeticIndex = planet.adiabeticIndex();                  // The adiabetic index gamma_a
    const double specificGasConstant = planet.specificGasConstant();        // The specific gas constant Rstar
    const double standardGravitationalParameter = planet.standardGravitationalParameter();  // The standard gravitational parameter mu
    const double rotationalVelocity = planet.rotationalVelocity();          // The rotational velocity of the planet Omega
    const double primeMeridianAngle = planet.primeMeridianAngle();          // The angle of the prime meridian when the inertial frame was set omega_p
    const double inertialFrameTime = planet.inertialFrameTime();            // Time between start of simulation and when inertial frame was set

    const double bodyReferenceRadius = planet.bodyReferenceRadius();         // Reference radius of the planet R


    const Eigen::MatrixXd temperaturePolyCoefficients = planet.temperaturePolyCoefficients();   // Temperature polynomial coefficients for temperature curve P_Tn
    const Eigen::MatrixXd temperatureAltitudeRanges = planet.temperatureAltitudeRanges();       // Temperature curve altitude ranges h
    const Eigen::VectorXd densityPolyCoefficients = planet.densityPolyCoefficients();           // Density polynomial coefficients for density curve P_rho n

    // MarsAscentVehicle
    const double Thrust = MAV.Thrust();                                                          // T     engine nominal thrust
    const double specificImpulse = MAV.specificImpulse();                                        // Isp     engine nominal specific impulse
    const double referenceArea = MAV.referenceArea();                                             // S   vehicle reference area
    const Eigen::MatrixXd dragCoefficientPolyCoefficients = MAV.dragCoefficientPolyCoefficients();            // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
    const Eigen::MatrixXd dragCoefficientMachRanges = MAV.dragCoefficientMachRanges();                       // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient
    const Eigen::MatrixXd thrustAzimuth = MAV.thrustAzimuth();                                // psiT   these are the thrust azimuth-gimbal angles as a function of time (including the time ranges)
    const Eigen::MatrixXd thrustElevation = MAV.thrustElevation();                                // epsilonT   these are the thrust elevation-gimbal angles as a function of time (including the time ranges)

    // StateAndTime

    const tudat::basic_mathematics::Vector7d currentState = currentStateAndTime.getCurrentState(); // The complete state including position, velocity and mass
    const double currentMass = currentStateAndTime.getCurrentMass();                                // The mass seperately
    const double currentTime = currentStateAndTime.getCurrentTime();                                // The current time

    double currentStepSize = stepSize.getCurrentStepSize();                      // The current step-size

    // Specified initial conditions

    const double FlightPathAngle = FlightPathAngle_;            // Set flight-path angle
    const double HeadingAngle = HeadingAngle_;                  // Set heading angle

//    std::cout<<"Flight-path angle in integration.cpp = "<<FlightPathAngle<<std::endl;
//    std::cout<<"Heading angle in integration.cpp = "<<HeadingAngle<<std::endl;


    //std::cout<<"This works right 1?"<<std::endl;

//////// Computations //////////

    /// Thrust acceleration in B-frame ///   thrustAccelerationsBframe // initial (for u)

    const Eigen::Vector3d thrustAccelerationsPframe = Eigen::Vector3d((Thrust/currentMass),0,0);            // THIS HAS TO BE CHANGED IN THE FUTURE TO INCLUDE A WIDE RANGE OF THRUST AZIMUTH AND ELEVATION ANGLES!!!

    Eigen::MatrixXd thrustAzimuthMatrix = MAV.thrustAzimuth();
    Eigen::MatrixXd thrustElevationMatrix = MAV.thrustElevation();

//    const double thrustAzimuthTestDeg = 0;             // thrust azimuth gimbal angle [Deg] 10 for testing
//    const double thrustElevationTestDeg = 0;            // thrust elevation gimbal angle [Deg] 5 for testing

//    const double thrustAzimuthTest = deg2rad(thrustAzimuthTestDeg);     // thrust azimuth gimbal angle [rad]
//    const double thrustElevationTest = deg2rad(thrustElevationTestDeg); // thrust elevation gimbal angle [rad]

      const double thrustAzimuthTest = thrustAzimuthMatrix(0,2);             // thrust azimuth gimbal angle [rad]
      const double thrustElevationTest = thrustElevationMatrix(0,2);            // thrust elevation gimbal angle [rad]


    const Eigen::Vector3d thrustAccelerationsBframe = getPropulsionToBodyFrameTransformationMatrix(thrustAzimuthTest,thrustElevationTest)*thrustAccelerationsPframe;

//    std::cout<<"thrustAccelerationsBframe = "<<thrustAccelerationsBframe<<std::endl;

//std::cout<<"This works right 2?"<<std::endl;

    /// Computing the auxiliary equations, derivatives and functions ///

    AuxiliaryCartesian Aux(adiabeticIndex, specificGasConstant,standardGravitationalParameter, rotationalVelocity, primeMeridianAngle,
              inertialFrameTime, bodyReferenceRadius, temperaturePolyCoefficients, temperatureAltitudeRanges,
              densityPolyCoefficients, Thrust, thrustAzimuthMatrix, thrustElevationMatrix, specificImpulse,
              referenceArea, dragCoefficientPolyCoefficients, dragCoefficientMachRanges);

    //std::cout<<"This works right 3?"<<std::endl;

    // Set the initial launch angles
    Aux.setFlightPathAngleAndHeadingAngle(FlightPathAngle,HeadingAngle);    // Used to specify a certain flight-path angle other than 90 degrees and heading angle other than 0 degrees at t = 0 sec

    // Compute the auxiliary equations

    Eigen::VectorXd auxiliaryEquations =  Aux.getCartesianAuxiliaryEquations(currentState,currentTime,thrustAccelerationsBframe);

//    std::cout<<"The auxiliaryEquations are "<<auxiliaryEquations<<std::endl;

    /// Debug ///
//    if (currentTime == 0.2){
//        auxiliaryEquations(23) = 1;
//    }
    /// Debug ///

//std::cout<<"This works right 4?"<<std::endl;

    // Compute the auxiliary derivatives

    Eigen::VectorXd auxiliaryDerivatives = Aux.getCartesianAuxiliaryDerivatives(currentState,currentTime,thrustAccelerationsBframe,auxiliaryEquations);

//    std::cout<<"The auxiliaryDerivatives are "<<auxiliaryDerivatives<<std::endl;

    // Compute the auxiliary functions

    Eigen::MatrixXd auxiliaryFunctions = Aux.getCartesianAuxiliaryFunctions(currentState,currentTime,thrustAccelerationsBframe,auxiliaryEquations,auxiliaryDerivatives);

//    std::cout<<"The auxiliaryFunctions are "<<auxiliaryFunctions<<std::endl;

    /// Computing the Taylor Coefficients ///

    Eigen::MatrixXd TaylorCoefficients = getCartesianTaylorCoefficients(adiabeticIndex, specificGasConstant, standardGravitationalParameter, rotationalVelocity, primeMeridianAngle,
                          inertialFrameTime, bodyReferenceRadius,temperaturePolyCoefficients, temperatureAltitudeRanges,
                          densityPolyCoefficients, Thrust, specificImpulse,
                          referenceArea, dragCoefficientPolyCoefficients, dragCoefficientMachRanges,
            thrustAccelerationsBframe,
            thrustAzimuthMatrix,
            thrustElevationMatrix,
            auxiliaryEquations,
            auxiliaryDerivatives,
            auxiliaryFunctions,
            currentTime,
            maxOrder);

    /// Storing the Taylor Coefficients to a file ///

    //std::cout<<"Does this even work1?"<<std::endl;

//        Eigen::MatrixXd TaylorCoefficientsOutputMatrix = Eigen::MatrixXd::Zero(7,maxOrder+1);       // Create an output matrix for the file without the first empty row

//        TaylorCoefficientsOutputMatrix.row(0) = TaylorCoefficients.row(1);                  // The first line entries are the maxOrder+1 Taylor Series Coefficients for     the position in the x-direction
//        TaylorCoefficientsOutputMatrix.row(1) = TaylorCoefficients.row(2);                  // The second line entries are the maxOrder+1 Taylor Series Coefficients for    the position in the y-direction
//        TaylorCoefficientsOutputMatrix.row(2) = TaylorCoefficients.row(3);                  // The third line entries are the maxOrder+1 Taylor Series Coefficients for     the position in the z-direction
//        TaylorCoefficientsOutputMatrix.row(3) = TaylorCoefficients.row(4);                  // The fourth line entries are the maxOrder+1 Taylor Series Coefficients for    the velocity in the x-direction
//        TaylorCoefficientsOutputMatrix.row(4) = TaylorCoefficients.row(5);                  // The fifth line entries are the maxOrder+1 Taylor Series Coefficients for     the velocity in the y-direction
//        TaylorCoefficientsOutputMatrix.row(5) = TaylorCoefficients.row(6);                  // The sixth line entries are the maxOrder+1 Taylor Series Coefficients for     the velocity in the z-direction
//        TaylorCoefficientsOutputMatrix.row(6) = TaylorCoefficients.row(7);                  // The seventh line entries are the maxOrder+1 Taylor Series Coefficients for   the mass

        Eigen::MatrixXd TaylorCoefficientsOutputMatrix = Eigen::MatrixXd::Zero(8,maxOrder+1);       // Create an output matrix for the file without the first empty row

        TaylorCoefficientsOutputMatrix(0,0) = currentTime;                                  // The first line entries give the current time as primary entry and zero as the rest
        TaylorCoefficientsOutputMatrix.row(1) = TaylorCoefficients.row(1);                  // The second line entries are the maxOrder+1 Taylor Series Coefficients for     the position in the x-direction
        TaylorCoefficientsOutputMatrix.row(2) = TaylorCoefficients.row(2);                  // The third line entries are the maxOrder+1 Taylor Series Coefficients for      the position in the y-direction
        TaylorCoefficientsOutputMatrix.row(3) = TaylorCoefficients.row(3);                  // The fourth line entries are the maxOrder+1 Taylor Series Coefficients for     the position in the z-direction
        TaylorCoefficientsOutputMatrix.row(4) = TaylorCoefficients.row(4);                  // The fifth line entries are the maxOrder+1 Taylor Series Coefficients for      the velocity in the x-direction
        TaylorCoefficientsOutputMatrix.row(5) = TaylorCoefficients.row(5);                  // The sixth line entries are the maxOrder+1 Taylor Series Coefficients for      the velocity in the y-direction
        TaylorCoefficientsOutputMatrix.row(6) = TaylorCoefficients.row(6);                  // The seventh line entries are the maxOrder+1 Taylor Series Coefficients for    the velocity in the z-direction
        TaylorCoefficientsOutputMatrix.row(7) = TaylorCoefficients.row(7);                  // The eighth line entries are the maxOrder+1 Taylor Series Coefficients for     the mass


        /// Start debug ///

        //std::cout<<"Does this even work2? :S"<<std::endl;

/*        std::cout<<"x-position coefficients are: "<<TaylorCoefficientsOutputMatrix.row(0)<<std::endl;
        std::cout<<"y-position coefficients are: "<<TaylorCoefficientsOutputMatrix.row(1)<<std::endl;
        std::cout<<"z-position coefficients are: "<<TaylorCoefficientsOutputMatrix.row(2)<<std::endl;
        std::cout<<"x-velocity coefficients are: "<<TaylorCoefficientsOutputMatrix.row(3)<<std::endl;
        std::cout<<"y-velocity coefficients are: "<<TaylorCoefficientsOutputMatrix.row(4)<<std::endl;
        std::cout<<"z-velocity coefficients are: "<<TaylorCoefficientsOutputMatrix.row(5)<<std::endl;
        std::cout<<"mass coefficients are: "<<TaylorCoefficientsOutputMatrix.row(6)<<std::endl;

 //*/       /// End debug ///



        // Set directory where output files will be stored. THIS REQUIRES THE COMPLETE PATH IN ORDER TO WORK!!
        const std::string outputDirectory = "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/02.taylorSeriesCoefficientsOutputFolder/";


        // Set output format for matrix output.
        Eigen::IOFormat csvFormat( 15, 0, ", ", "\n" );


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

        std::string newFileName = "CartesianTaylorSeriesCoefficientsFileAtDateAndTime_" + ComputerTimeString + ".csv";

//    std::cerr<<"The file name that you have chosen already exists, a new file with name "<<newFileName<<" will be created to store the data for now"<<std::endl;

//    // Set new absolute path to file containing the data.
//    dataAbsolutePathTSI = outputDirectoryTSI + newFileName;

    // Set absolute path to file containing the Taylor Series Coefficients.
    const std::string taylorSeriesCoefficientsAbsolutePath = outputDirectory + newFileName;


        // Set absolute path to file containing the Taylor Series Coefficients.
//        const std::string taylorSeriesCoefficientsAbsolutePath = outputDirectory + "test1TaylorSeriesCoefficients(bugSearch31-05-2016).csv";



        // Check if the file already exists.


        std::ifstream ifile(taylorSeriesCoefficientsAbsolutePath.c_str()); // Check it as an input file

        bool fexists = false;   // Set the default to "It does not exist"

        if (ifile){         // Attempt to open the file


           fexists = true;      // If the file can be opened it must exist

           ifile.close();   // Close the file

        }


        // If so: append, if not: create new file and put data in

        if (fexists == true){

            // Export the Taylor Series Coefficients matrix.
            std::ofstream exportFile1;                          // Define the file as an output file


            exportFile1.open(taylorSeriesCoefficientsAbsolutePath.c_str(),std::ios_base::app);      // Open the file in append mode

            exportFile1 << "\n";                                            // Make sure the new matrix start on a new line

            exportFile1 << TaylorCoefficientsOutputMatrix.format( csvFormat ); // Add the new values

//            std::cout<<"The file called "<<taylorSeriesCoefficientsAbsolutePath<<" has been appended"<<std::endl;


            exportFile1.close( );   // Close the file
}
            else{

            // Export the Taylor Series Coefficients matrix.
            std::ofstream exportFile1( taylorSeriesCoefficientsAbsolutePath.c_str( ) ); // Make the new file
            std::cout<<"New file called "<<taylorSeriesCoefficientsAbsolutePath<<" has been created"<<std::endl;
            exportFile1 << TaylorCoefficientsOutputMatrix.format( csvFormat );          // Store the new values
            exportFile1.close( );   // Close the file
        };

//std::cout<<"Does this even work3?"<<std::endl;

        /// Performing the actual Taylor Series expansion for every state variable ///



        tudat::basic_mathematics::Vector7d updatedState = tudat::basic_mathematics::Vector7d::Zero();        // Create a vector for the updatedState and setting it to zero

        for (int n = 0; n<updatedState.size();n++){                 // All variables

        for (int k = 0; k<maxOrder+1;k++){                      // Taylor series summation

            updatedState(n) += TaylorCoefficients((n+1),k)*pow(currentStepSize,k);      // Perform one step of the taylor series expansion and then add it to the previous step

//            if (n==3){
//           std::cout<<"updatedState(3) = "<<updatedState(3)<<std::endl;
//           std::cout<<"currentStepSize = "<<currentStepSize<<std::endl;
//           std::cout<<"TaylorCoefficients((n+1),k) = "<<TaylorCoefficients((n+1),k)<<std::endl;

//            }

        } // Taylor series summation

}   // All variables







        /// Determine if a new section has been reached for the temperature and correct timestep accordingly ///

    //        int currentSectionT = 0; // Default value
//            double originalAltitude = currentState(0)-bodyReferenceRadius; // h
            double originalAltitude = sqrt(currentState(0)*currentState(0)+currentState(1)*currentState(1)+currentState(2)*currentState(2))-bodyReferenceRadius; //h
            double limitAltitude = temperatureAltitudeRanges(0,1); // Default value

                    for (int i=0; i < 3+1; i++){

            if (temperatureAltitudeRanges(i,0) <= originalAltitude && originalAltitude < temperatureAltitudeRanges(i,1)){

    //            currentSectionT = i;
    //            std::cout<<"It actually does go inside the limitAltitude section"<<std::endl;
                 limitAltitude = temperatureAltitudeRanges(i,1);
    //             std::cout<<"limitAltitude = "<<limitAltitude<<std::endl;


            }
            else {
    //            std::cout<<"The originalAltitude = "<<originalAltitude<<", which is lower than the lowest altitude."<<std::endl;
            }

            /// Debug ///

    //        std::cout<<"originalAltitude = "<<originalAltitude<<std::endl;
    //        std::cout<<"temperatureAltitudeRanges(i,0) = "<<temperatureAltitudeRanges(i,0)<<std::endl;
    //        std::cout<<"temperatureAltitudeRanges(i,1) = "<<temperatureAltitudeRanges(i,1)<<std::endl;


            /// Debug ///

    };



                    // Set the initial values for the time-step domain

                    double tFrom = currentTime;    // t0,0
                    double fFrom = originalAltitude-limitAltitude; // f(t0,0)
                    double fFromDot;        // fdot(t0,0)
                    double fFromDotX;       // fdotx (t0,0)
                    double fFromDotY;       // fdoty (t0,0)
                    double fFromDotZ;       // fdotz (t0,0)
                    double tTo = currentTime+currentStepSize;  // t1,0
//                    double newAltitude = updatedState(0)-bodyReferenceRadius; // h_new
                    double newAltitude = sqrt(updatedState(0)*updatedState(0)+updatedState(1)*updatedState(1)+updatedState(2)*updatedState(2))-bodyReferenceRadius; // h_new


                    double fTo = newAltitude-limitAltitude; // f(t1,0)

                    /// Debug ///

                    std::cout<<"limitAltitude = "<<limitAltitude<<std::endl;
                    std::cout<<"originalAltitude = "<<originalAltitude<<std::endl;
                    std::cout<<"newAltitude = "<<newAltitude<<std::endl;
                    std::cout<<"fFrom = "<<fFrom<<std::endl;
                    std::cout<<"fTo = "<<fTo<<std::endl;
                    std::cout<<"currentStepSize before = "<<currentStepSize<<std::endl;
                    std::cout<<"initialNewTime = "<<currentTime+currentStepSize<<std::endl;

                    /// Debug ///


    bool signIsPositive = false; // Check in which direction, up or down, the curve is heading

            if (fFrom<0){
                signIsPositive = true;
            }


    bool altitudeAccept = false; // Default value is that the altitude is not accepted and has therefore gone beyond the altitude limit

    if (fTo <0){    // If the newAltitude is below the limitAltitude (so fTO<0) then the new altitude is still within the same section and can be accepted

        altitudeAccept = true;
    }
    else{

            std::cout<<"//////////////////////////////////////////////////////////////////////// Do loop has been started ////////////////////////////////////////////////////////////////////////"<<std::endl;
            do{ // If the altitude has gone beyond the limit altitude, this do loop will determine a new "currentStepSize"


            if (newAltitude - limitAltitude <= 1e-6 && newAltitude - limitAltitude >= 0.0){   // Checking if the convergence condition has been met and an answer has been found

                altitudeAccept = true;
    }
            else{
                // Compute an updated time step //

                fFromDot = 0.0; // Reset
                fFromDotX = 0.0; // Reset
                fFromDotY = 0.0; // Reset
                fFromDotZ = 0.0; // Reset

                for (int k = 1; k < maxOrder+1; k++){
                fFromDotX += k*TaylorCoefficients(1,k)*pow(currentStepSize,(k-1)) ;      // Compute tFromDotX
                fFromDotY += k*TaylorCoefficients(2,k)*pow(currentStepSize,(k-1)) ;      // Compute tFromDotY
                fFromDotZ += k*TaylorCoefficients(3,k)*pow(currentStepSize,(k-1)) ;      // Compute tFromDotZ

    }
                fFromDot = (currentState(0)*fFromDotX+currentState(1)*fFromDotY+currentState(2)*fFromDotZ)/(sqrt(updatedState(0)*updatedState(0)+updatedState(1)*updatedState(1)+updatedState(2)*updatedState(2))); // Compute tFromDot which is the derivative of the radius vector

                double alpha = (fTo-fFrom)/((tTo-tFrom)*(tTo-tFrom))-fFromDot/(tTo-tFrom);  // Basically computing fFromDoubleDot/2 (or coefficient)

                double tFromNew;
                if (signIsPositive == true){    // Compute the new "from" time
            tFromNew = tFrom + (-fFromDot+sqrt(fFromDot*fFromDot+4.0*alpha*fFrom))/(2.0*alpha);

            /// Debug ///

            std::cout<<"Debug"<<std::endl;
            std::cout<<"tTo = "<<tTo<<std::endl;
            std::cout<<"tFrom = "<<tFrom<<std::endl;
            std::cout<<"(tTo-tFrom) = "<<(tTo-tFrom)<<std::endl;

            std::cout<<"(fTo-fFrom) = "<<(fTo-fFrom)<<std::endl;
            std::cout<<"fFromDot = "<<fFromDot<<std::endl;

            std::cout<<"(fTo-fFrom)/((tTo-tFrom)*(tTo-tFrom)) = "<<(fTo-fFrom)/((tTo-tFrom)*(tTo-tFrom))<<std::endl;
            std::cout<<"fFromDot/(tTo-tFrom) = "<<fFromDot/(tTo-tFrom)<<std::endl;
            std::cout<<"(tTo-tFrom) = "<<(tTo-tFrom)<<std::endl;
            std::cout<<"((tTo-tFrom)*(tTo-tFrom)) = "<<((tTo-tFrom)*(tTo-tFrom))<<std::endl;


            std::cout<<"(-fFromDot+sqrt(fFromDot*fFromDot+4.0*alpha*fFrom))/(2.0*alpha) = "<<(-fFromDot+sqrt(fFromDot*fFromDot+4.0*alpha*fFrom))/(2.0*alpha)<<std::endl;
            std::cout<<"alpha = "<<alpha<<std::endl;
            std::cout<<"sqrt(fFromDot*fFromDot+4.0*alpha*fFrom) = "<<sqrt(fFromDot*fFromDot+4.0*alpha*fFrom)<<std::endl;
            std::cout<<"fFromDot*fFromDot+4.0*alpha*fFrom = "<<fFromDot*fFromDot+4.0*alpha*fFrom<<std::endl;

//            std::cout<<"Debug"<<std::endl;
            std::cout<<" "<<std::endl;

            /// Debug ///

    } // +
                else {
                   tFromNew = tFrom + (-fFromDot-sqrt(fFromDot*fFromDot+4.0*alpha*fFrom))/(2.0*alpha);
                } // -

                currentStepSize = currentTime-tFromNew; // Update the new step-size to compute the new value for the altitude using the Taylor Coefficients

                /// Debug ///
//                std::cout<<"signIsPositive = "<<signIsPositive<<std::endl;
//                std::cout<<"fFromDot = "<<fFromDot<<std::endl;
//                std::cout<<"tFromNew = "<<tFromNew<<std::endl;
//                std::cout<<"tFrom = "<<tFrom<<std::endl;
//                std::cout<<"currentStepSize = currentTime-tFromNew = "<<currentStepSize<<std::endl;
////                std::cout<<" "<<std::endl;
                /// Debug ///


            tudat::basic_mathematics::Vector7d updatedStateNew = tudat::basic_mathematics::Vector7d::Zero();        // Create a vector for the updatedState and setting it to zero

                for (int n = 0; n<updatedStateNew.size();n++){                 // All variables

                for (int k = 0; k<maxOrder+1;k++){                      // Taylor series summation

                    updatedStateNew(n) += TaylorCoefficients((n+1),k)*pow(currentStepSize,k);      // Perform one step of the taylor series expansion and then add it to the previous step


                } // Taylor series summation

                }   // All variables

                updatedState = updatedStateNew; // update the updated state

//                newAltitude = updatedState(0)-bodyReferenceRadius;  // Determine the new final altitude at the end of the integration step (if this step-size would be the actual step-size for this integration step)
                newAltitude = sqrt(updatedState(0)*updatedState(0)+updatedState(1)*updatedState(1)+updatedState(2)*updatedState(2))-bodyReferenceRadius;  // Determine the new final altitude at the end of the integration step (if this step-size would be the actual step-size for this integration step)


                double fFromNew = newAltitude-limitAltitude; // Determine the new "from" function value

                if ((fFrom/(fabs(fFrom))) != (fFromNew/(fabs(fFromNew)))){  // If the function value of the new "from" value is on the other side of the root, then the old "from" values are now the new "to" values (evaluate back to the beginning). Otherwise the "to" values don't change.
                    fTo = fFrom;
                    tTo = tFrom;

                    /// Debug ///
//                    std::cout<<"fTo and tTo changed!"<<std::endl;
//                    std::cout<<"tTo = "<<tTo<<std::endl;
//                    std::cout<<"tFromNew = "<<tFromNew<<std::endl;
//                    std::cout<<"tFrom = "<<tFrom<<std::endl;
                    /// Debug ///

                }

                fFrom = fFromNew;   // Update the new "from" values
                tFrom = tFromNew;

                /// Debug ///
//                std::cout<<"Debug"<<std::endl;
//                std::cout<<" "<<std::endl;
                /// Debug ///






            } // update!

    }while(altitudeAccept == false);

            std::cout<<"currentStepSize after = "<<currentStepSize<<std::endl;
            std::cout<<"//////////////////////////////////////////////////////////////////////// Do loop has ended ////////////////////////////////////////////////////////////////////////"<<std::endl;

    } // end of else



//*/









//std::cout<<"Does this even work4?"<<std::endl;
        double updatedTime = currentTime+currentStepSize;           // Create the updated time variable


        /// Updating the step-size ///

        tudat::basic_mathematics::Vector7d penultimateCoefficients;         // Define the Xn(K-1) vector
        tudat::basic_mathematics::Vector7d lastCoefficients;                // Define the Xn(K) vector

        for (int i = 0; i<7; i++){                                      // Fill the vectors

            penultimateCoefficients(i)= TaylorCoefficients((i+1),maxOrder-1);       // Xn(K-1)
            lastCoefficients(i)= TaylorCoefficients((i+1),(maxOrder));              // Xn(K)
        }

//        std::cout<<"TaylorSeriesIntegration works till here 1"<<std::endl;

        stepSize.updateStepSizeUsingIteration(penultimateCoefficients, lastCoefficients, maxOrder); // Determining the new step-size and updating the current step-size to that

//        std::cout<<"TaylorSeriesIntegration works till here 2"<<std::endl;

        // Please note that instead of updateStepSizeUsingInteration, you can also use updateStepSizeUsingPreviousStepSize, however, this method does not work if the maximum truncation error estimate is larger
        // than the local error tolerance and requires the previous step-size

        //std::cout<<"Does this even work5?"<<std::endl;

        /// Setting the ouput vector ///

        Eigen::VectorXd updatedStateAndTime(8);

        updatedStateAndTime(0) = updatedState(0);   // Updated position in the x-direction
        updatedStateAndTime(1) = updatedState(1);   // Updated position in the y-direction
        updatedStateAndTime(2) = updatedState(2);   // Updated position in the z-direction
        updatedStateAndTime(3) = updatedState(3);   // Updated velocity in the x-direction
        updatedStateAndTime(4) = updatedState(4);   // Updated velocity in the y-direction
        updatedStateAndTime(5) = updatedState(5);   // Updated velocity in the z-direction
        updatedStateAndTime(6) = updatedState(6);   // Updated mass
        updatedStateAndTime(7) = updatedTime;   // Updated time

//        std::cout<<"updatedStateAndTime = "<<updatedStateAndTime<<std::endl;

//        std::cout<<"TaylorSeriesIntegration works till here 3"<<std::endl;

    return updatedStateAndTime;

} // End of the function
