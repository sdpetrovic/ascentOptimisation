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
 *      160413    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */





#ifndef AUXILIARY_H
#define AUXILIARY_H


#include <iostream>
#include <cmath>

#include <Eigen/Core>

#include <tudat/Tudat/Mathematics/BasicMathematics/linearAlgebra.h>
#include <tudatApplications/thesisProject/linearAlgebraTypesUpdated.h>

using namespace std;

class Auxiliary

        /* This class will describe the different auxiliary equations, derivatives and functions
         * These are represented by:
         *
         *  - xn                    Auxiliary equation with n is 1, ..., 34
         *  - un                    Auxiliary derivatives n is 1, ..., 34
         *  - wn,m                  Auxiliary functions n is 1, ..., 34 and m is 1, ..., 24
         *
         */





{
public:

    /* In this case, the constructor does only takes celestial body and vehicle constant input. The class function will contain the variable parameters
     *
     *    // The diferent celestial body constant parameters and polynomial coefficient parameter matrices
     *
     * const double adiabeticIndex                                   // gamma_a      adiabetic index
     * const double specificGasConstant                        // Rstar    [m^2/(s^2*K)]    specific gas constant
     * const double standardGravitationalParameter  // mu_M     [m^3/s^2]    standard gravitational parameter
     * const double rotationalVelocity                          // rotational velocity of Mars  [rad/s]
     * const double primeMeridianAngle                           // OmegaP   [rad]   relative angle between the prime meridian and the x-axis
     * const double inertialFrameTime                             // t0       [s]    time between the start time and the time that the inertial frame was set
     * const Eigen::MatrixXd temperaturePolyCoefficients // PTn    temperature polynomial coefficients
     * const Eigen::MatrixXd temperatureAltitudeRanges     // altitude range per section for the temperature-altitude curve [km MOLA]
     * const Eigen::VectorXd densityPolyCoefficients         // Prho n density polynomial coefficients
     *
     *     // The differnt vehicle constant parameters and polynomial coefficients
     *
     * const double Thrust                                                         // T   [N]  engine nominal thrust
     * const double specificImpulse                                         // Isp [s]    engine nominal specific impulse
     * const double referenceArea                                            // S [m^2]  vehicle reference area
     * const Eigen::MatrixXd dragCoefficientPolyCoefficients            // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
     * const Eigen::MatrixXd dragCoefficientMachRanges                        // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient
     *
     *
     *
     */


    Auxiliary(const double adiabeticIndex_, const double specificGasConstant_, const double standardGravitationalParameter_, const double rotationalVelocity_, const double primeMeridianAngle_,
              const double inertialFrameTime_, const double bodyReferenceRadius_, const Eigen::MatrixXd temperaturePolyCoefficients_, const Eigen::MatrixXd temperatureAltitudeRanges_,
              const Eigen::VectorXd densityPolyCoefficients_, const double Thrust_, const double specificImpulse_,
              const double referenceArea_, const Eigen::MatrixXd dragCoefficientPolyCoefficients_, const Eigen::MatrixXd dragCoefficientMachRanges_){

            // Set the diferent celestial body constant parameters and polynomial coefficient parameter matrices

         adiabeticIndex = adiabeticIndex_;                                   // gamma_a      adiabetic index
         specificGasConstant = specificGasConstant_;                        // Rstar    [m^2/(s^2*K)]    specific gas constant
         standardGravitationalParameter = standardGravitationalParameter_;  // mu_M     [m^3/s^2]    standard gravitational parameter
         rotationalVelocity = rotationalVelocity_;                         // rotational velocity of Mars  [rad/s]
         primeMeridianAngle = primeMeridianAngle_;                          // OmegaP   [rad]   relative angle between the prime meridian and the x-axis
         inertialFrameTime = inertialFrameTime_;                            // t0       [s]    time between the start time and the time that the inertial frame was set
         bodyReferenceRadius = bodyReferenceRadius_;                                  // Rm       [m]     MOLA radius of Mars

         temperaturePolyCoefficients = temperaturePolyCoefficients_; // PTn    temperature polynomial coefficients
         temperatureAltitudeRanges = temperatureAltitudeRanges_;    // altitude range per section for the temperature-altitude curve [km MOLA]
         densityPolyCoefficients = densityPolyCoefficients_;         // Prho n density polynomial coefficients

             // Set the differnt vehicle constant parameters and polynomial coefficients

         Thrust = Thrust_;                                                         // T   [N]  engine nominal thrust
         specificImpulse = specificImpulse_;                                        // Isp [s]    engine nominal specific impulse
         referenceArea = referenceArea_;                                           // S [m^2]  vehicle reference area
         dragCoefficientPolyCoefficients = dragCoefficientPolyCoefficients_;           // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
         dragCoefficientMachRanges = dragCoefficientMachRanges_;                       // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient







    }




    /// Compute the Auxiliary Equations ///

    /* In c++ the convention is that the first entry of a vector is 0, however in this case this first entry will be filled up by w4,2 such that the entry corresponds to the associated auxiliary equation number.
     *
     * The getAuxiliaryEquations function takes three inputs:
     *
     *  - const tudat::basic_mathematics::Vector7d& aState      Which is the current state
     *  - const double time                                     Which is the current time
     *  - const Eigen::Vector3d& thrustAccelerationsBframe      Which are the thrust accelerations in the Bframe
     *
     */

    Eigen::VectorXd getAuxiliaryEquations( const tudat::basic_mathematics::Vector7d& aState, const double time, const Eigen::Vector3d& thrustAccelerationsBframe){

        auxiliaryEquationsVector = Eigen::VectorXd::Zero(45);       // Setting the complete vector and filling it with zeros for now

        // The following expressions are described in the order in which the equations have to be computed corresponding to the respective vector entry

        auxiliaryEquationsVector(1) = aState(0);              // x1
        auxiliaryEquationsVector(2) = aState(1);              // x2
        auxiliaryEquationsVector(3) = aState(2);              // x3
        auxiliaryEquationsVector(4) = aState(3);              // x4
        auxiliaryEquationsVector(5) = aState(4);              // x5
        auxiliaryEquationsVector(6) = aState(5);              // x6
        auxiliaryEquationsVector(7) = aState(6);              // x7

        auxiliaryEquationsVector(8) = auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2)+
                auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3) ;              // x8

        auxiliaryEquationsVector(9) = pow(auxiliaryEquationsVector(8), 1.5);              // x9

        auxiliaryEquationsVector(10) = rotationalVelocity*(inertialFrameTime+time)-primeMeridianAngle ;              // x10

        auxiliaryEquationsVector(19) = auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2);              // x19

        auxiliaryEquationsVector(21) = auxiliaryEquationsVector(4)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(5)*auxiliaryEquationsVector(5)+
                auxiliaryEquationsVector(6)*auxiliaryEquationsVector(6) ;              // x21

        auxiliaryEquationsVector(26) = 2*(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5)+
                                          auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6));              // x26

        auxiliaryEquationsVector(11) = atan2(auxiliaryEquationsVector(2),auxiliaryEquationsVector(1))-auxiliaryEquationsVector(10);              // x11

        auxiliaryEquationsVector(20) = pow(auxiliaryEquationsVector(8), 0.5);              // x20

        auxiliaryEquationsVector(36) = pow(auxiliaryEquationsVector(21), 0.5) ;                // x36

        auxiliaryEquationsVector(12) = asin(auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)) ;              // x12

        auxiliaryEquationsVector(25) = auxiliaryEquationsVector(25)/(2*auxiliaryEquationsVector(20));              // x25

        // Please note that the altitude h (x31) is expressed in km MOLA (which is also the input for the density and temperature curves!)
        auxiliaryEquationsVector(31) = (auxiliaryEquationsVector(20)-bodyReferenceRadius)/1000;              // x31 [km]!!!

        auxiliaryEquationsVector(24) = (auxiliaryEquationsVector(20)*auxiliaryEquationsVector(6)-auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25))/
                (auxiliaryEquationsVector(8)*sqrt(1-pow((auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)),2)));              // x24

        // Computing the polynomial fit using the altitude and fit parameters for density
        for (int i = 0; i < 10;i++) {

        auxiliaryEquationsVector(30) += pow(auxiliaryEquationsVector(31),i)*densityPolyCoefficients(i);              // x30
};

        // Determine which section of the temperature curve needs to be used and what the corresponding order is

        if (temperatureAltitudeRanges(0,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(0,1)){

        sectionT = 0;
        powerT = 1;

        }
        else if (temperatureAltitudeRanges(1,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(1,1)){

        sectionT = 1;
        powerT = 3;

        }
        else if (temperatureAltitudeRanges(2,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(2,1)){

        sectionT = 2;
        powerT = 6;

        }
        else if (temperatureAltitudeRanges(3,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(3,1)){

            sectionT = 3;
            powerT = 8;
        }
        else if (temperatureAltitudeRanges(4,0) <= auxiliaryEquationsVector(31)){

            sectionT = 4;
            powerT = 0;
        }
        else {


            std::cerr<<"The current altitude: "<<auxiliaryEquationsVector(31)<<" [km MOLA] is not a valid altitude (lower than the lowest reference altitude)"<<std::endl;

                       sectionT = 0;
                        powerT = 1;

        };

        // Computing the polynomial fit using the altitude and fit parameters for temperature
        for (int i=0; i < powerT;i++){

        auxiliaryEquationsVector(34) += pow(auxiliaryEquationsVector(31),i)*temperaturePolyCoefficients(sectionT,i);              // x34

};
        auxiliaryEquationsVector(35) = rotationalVelocity*auxiliaryEquationsVector(20)*cos(auxiliaryEquationsVector(12));                // x35

        auxiliaryEquationsVector(37) = auxiliaryEquationsVector(25)/auxiliaryEquationsVector(36);                // x37

        auxiliaryEquationsVector(18) = auxiliaryEquationsVector(20)*auxiliaryEquationsVector(24);              // x18

        auxiliaryEquationsVector(28) = exp(auxiliaryEquationsVector(30));              // x28

        auxiliaryEquationsVector(33) = sqrt(adiabeticIndex*specificGasConstant*auxiliaryEquationsVector(34));              // x33

        auxiliaryEquationsVector(38) = asin(auxiliaryEquationsVector(37));                // x38

        auxiliaryEquationsVector(41) = auxiliaryEquationsVector(35)*auxiliaryEquationsVector(36);                // x41

        auxiliaryEquationsVector(44) = auxiliaryEquationsVector(36)*cos(auxiliaryEquationsVector(38));                // x44

        auxiliaryEquationsVector(39) = auxiliaryEquationsVector(18)/auxiliaryEquationsVector(44);                // x39

        auxiliaryEquationsVector(40) = acos(auxiliaryEquationsVector(39));                // x40

        auxiliaryEquationsVector(42) = cos(auxiliaryEquationsVector(38)*sin(auxiliaryEquationsVector(40)));                // x42

        auxiliaryEquationsVector(43) = auxiliaryEquationsVector(41)*auxiliaryEquationsVector(42);                // x43

        auxiliaryEquationsVector(15) = sqrt(auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)+auxiliaryEquationsVector(21)-2*auxiliaryEquationsVector(43));              // x15

        auxiliaryEquationsVector(23) = auxiliaryEquationsVector(25)/auxiliaryEquationsVector(15);              // x23

        auxiliaryEquationsVector(14) = asin(auxiliaryEquationsVector(23));              // x14

        auxiliaryEquationsVector(16) = cos(auxiliaryEquationsVector(14));              // x16

        auxiliaryEquationsVector(32) = auxiliaryEquationsVector(15)/auxiliaryEquationsVector(33);              // x32

        auxiliaryEquationsVector(17) = auxiliaryEquationsVector(16)*auxiliaryEquationsVector(15);              // x17

        // Determine which section of the drag coefficient curve needs to be used

        for (int i=0; i < 5; i++){

            if (dragCoefficientMachRanges(i,0) <= auxiliaryEquationsVector(32) && auxiliaryEquationsVector(32) < dragCoefficientMachRanges(i,1)){

                sectionCD = i;
            }

        };

        auxiliaryEquationsVector(29) = dragCoefficientPolyCoefficients(sectionCD,1)*auxiliaryEquationsVector(32)+dragCoefficientPolyCoefficients(sectionCD,0);              // x29

        auxiliaryEquationsVector(22) = auxiliaryEquationsVector(18)/auxiliaryEquationsVector(17);              // x22

        auxiliaryEquationsVector(27) = 0.5*referenceArea*auxiliaryEquationsVector(28)*auxiliaryEquationsVector(15)*auxiliaryEquationsVector(15)*auxiliaryEquationsVector(29);              // x27

        auxiliaryEquationsVector(13) = acos(auxiliaryEquationsVector(22));              // x13

        auxiliaryEquationsVector(0) = thrustAccelerationsBframe(0)-(auxiliaryEquationsVector(27)/auxiliaryEquationsVector(7));              // w4,2





// auxiliaryEquationsVector()





       return auxiliaryEquationsVector;
    }







private:

    // The diferent celestial body constant parameters and polynomial coefficient parameter matrices

 double adiabeticIndex;                                   // gamma_a      adiabetic index
 double specificGasConstant;                        // Rstar    [m^2/(s^2*K)]    specific gas constant
 double standardGravitationalParameter;  // mu_M     [m^3/s^2]    standard gravitational parameter
 double rotationalVelocity;                         // rotational velocity of Mars  [rad/s]
 double primeMeridianAngle;                          // OmegaP   [rad]   relative angle between the prime meridian and the x-axis
 double inertialFrameTime;                            // t0       [s]    time between the start time and the time that the inertial frame was set
 double bodyReferenceRadius;                         // Rm    [m]       MOLA radius of Mars

 Eigen::MatrixXd temperaturePolyCoefficients; // PTn    temperature polynomial coefficients
 Eigen::MatrixXd temperatureAltitudeRanges;    // altitude range per section for the temperature-altitude curve [km MOLA]
 Eigen::VectorXd densityPolyCoefficients;         // Prho n density polynomial coefficients

     // The different vehicle constant parameters and polynomial coefficients

 double Thrust;                                                         // T   [N]  engine nominal thrust
 double specificImpulse;                                        // Isp [s]    engine nominal specific impulse
 double referenceArea;                                           // S [m^2]  vehicle reference area
 Eigen::MatrixXd dragCoefficientPolyCoefficients;           // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
 Eigen::MatrixXd dragCoefficientMachRanges;                       // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient



    Eigen::VectorXd auxiliaryEquationsVector;               // The vector containing the auxiliary equations xn
    Eigen::VectorXd auxiliaryDerivativesVector;             // The vector containing the auxiliary derivatives un
    Eigen::VectorXd auxiliaryFunctionsVector;               // The vector containgin the auxiliary functions wn,m



    // Additional in-class used variables


    int sectionT;                                   // This variable holds the "(section number -1)" for the temperature curve fit
    int powerT;                                     // This variable holds the section corresponding order for the temperature curve fit

    int sectionCD;                                  // This variable holds the "(section number -1)" for the drag coefficient curve fit

};

#endif // AUXILIARY_H
