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
 *      160518    S.D. Petrovic     Fixed the mistake I made with the transformation matrix T_IB which resulted in a mistake in u6
 *      160520    S.D. Petrovic     Fixed mistake in x25 where it said x25 = x25/(2*x20) which should be x25 = x26/(2*x20). Also, u41 had a - instead of a +.
 *                                  Also, at u45 and u21 the tolerance had to include an abs function!
 *      160526    S.D. Petrovic     Added W4,0 to be able to properly evaluate the recurrence relation of W4,2
 *      160527    S.D. Petrovic     Corrected mistake in x42
 *      160531    S.D. Petrovic     Added more ways to deal with rounding errors
 *      160602    S.D. Petrovic     Added the thrust auxiliary functions
 *      160603    S.D. Petrovic     Updated u10 to be OmegaM instead of 0!
 *      160618    S.D. Petrovic     Found mistake in u24 and updated it and the corresponding auxiliary functions
 *      160622    S.D. Petrovic     Found a huge mistake in u15 (extra acc were not taken into account, so it was giving zero values != possible), implemented new equations...
 *
 *    References
 *
 *    Notes
 *
 */





#ifndef AUXILIARY_H
#define AUXILIARY_H


#include <iostream>
#include <iomanip>
#include <cmath>

#include <Eigen/Core>

#include <tudat/Tudat/Mathematics/BasicMathematics/linearAlgebra.h>
#include <tudatApplications/thesisProject/linearAlgebraTypesUpdated.h>

#include <tudatApplications/thesisProject/physicalConstantsUpdated.h>

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

    /* In this case, the constructor only takes celestial body and vehicle constant input. The class function will contain the variable parameters
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
              const Eigen::VectorXd densityPolyCoefficients_, const double Thrust_, const Eigen::MatrixXd thrustAzimuthMatrix_, const Eigen::MatrixXd thrustElevationMatrix_, const double specificImpulse_,
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
         thrustAzimuthMatrix = thrustAzimuthMatrix_;                                // psi_T [rad] thrust azimuth angles
         thrustElevationMatrix = thrustElevationMatrix_;                            // epsilon_T [rad] thrust elevation angles
         specificImpulse = specificImpulse_;                                        // Isp [s]    engine nominal specific impulse
         referenceArea = referenceArea_;                                           // S [m^2]  vehicle reference area
         dragCoefficientPolyCoefficients = dragCoefficientPolyCoefficients_;           // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
         dragCoefficientMachRanges = dragCoefficientMachRanges_;                       // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient


std::cout<<"verticalInertialFlightPathAngleSet original = "<<verticalInertialFlightPathAngleSet<<std::endl;
std::cout<<"verticalRotationalFlightPathAngleSet original = "<<verticalRotationalFlightPathAngleSet<<std::endl;



        // Set the booleans to false in case of faulty memory assignment and mistakes in the deletion of the previous class
       rotationalFlightPathAngleSet = false;         // All of these are used to let the program know that a predefined angle was set and that that angle should be used (initially)
       inertialFlightPathAngleSet = false;
       rotationalHeadingAngleSet = false;
       inertialHeadingAngleSet = false;


       verticalRotationalFlightPathAngleSet = false;       // All of these are used for the vertical ascent case
       verticalInertialFlightPathAngleSet = false;
       verticalRotationalHeadingAngleSet = false;
       verticalInertialHeadingAngleSet = false;





std::cout<<"verticalInertialFlightPathAngleSet original 2 = "<<verticalInertialFlightPathAngleSet<<std::endl;




    }

/// Set functions ///

    void setRotationalFlightPathAngle(const double rotationalFlightPathAngle_){
     rotationalFlightPathAngle = rotationalFlightPathAngle_;
     rotationalFlightPathAngleSet = true;
     if (rotationalFlightPathAngle_ ==tudat::mathematical_constants::LONG_PI/2){
         verticalRotationalFlightPathAngleSet = true;
     }

    }         // Rotational flight path angle in rad

    void setInertialFlightPathAngle(const double inertialFlightPathAngle_){
        inertialFlightPathAngle = inertialFlightPathAngle_;
        inertialFlightPathAngleSet = true;
        if (inertialFlightPathAngle_ ==tudat::mathematical_constants::LONG_PI/2){
            verticalInertialFlightPathAngleSet = true;
        }
    }           // Inertial flight path angle in rad

    void setRotationalHeadingAngle(const double rotationalHeadingAngle_){
        rotationalHeadingAngle = rotationalHeadingAngle_;
        rotationalHeadingAngleSet = true;
    }            // Rotational heading angle in rad

    void setInertialHeadingAngle(const double inertialHeadingAngle_){
        inertialHeadingAngle = inertialHeadingAngle_;
        inertialHeadingAngleSet = true;
    }              // Inertial heading angle in rad

    /// Test functions to revert changes if tolerance is reached



//////////////////////////////////////////////// Auxiliary Equations //////////////////////////////////////////////////////////////////////

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
        std::cout<<"verticalInertialFlightPathAngleSet eq 1 = "<<verticalInertialFlightPathAngleSet<<std::endl;
        std::cout<<"verticalRotationalFlightPathAngleSet eq 1 = "<<verticalRotationalFlightPathAngleSet<<std::endl;

        auxiliaryEquationsVector = Eigen::VectorXd::Zero(56);       // Setting the complete vector and filling it with zeros for now

        // The following expressions are described in the order in which the equations have to be computed corresponding to the respective vector entry

        auxiliaryEquationsVector(1) = aState(0);              // x1
        auxiliaryEquationsVector(2) = aState(1);              // x2
        auxiliaryEquationsVector(3) = aState(2);              // x3
        auxiliaryEquationsVector(4) = aState(3);              // x4
        auxiliaryEquationsVector(5) = aState(4);              // x5
        auxiliaryEquationsVector(6) = aState(5);              // x6
        auxiliaryEquationsVector(7) = aState(6);              // x7

       //std::cout<<"Surely this works 1..."<<std::endl;

        auxiliaryEquationsVector(8) = auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2)+
                auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3) ;              // x8

        auxiliaryEquationsVector(9) = pow(auxiliaryEquationsVector(8), 1.5);              // x9

        auxiliaryEquationsVector(10) = rotationalVelocity*(inertialFrameTime+time)-primeMeridianAngle ;              // x10

        auxiliaryEquationsVector(19) = auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2);              // x19

        auxiliaryEquationsVector(21) = auxiliaryEquationsVector(4)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(5)*auxiliaryEquationsVector(5)+
                auxiliaryEquationsVector(6)*auxiliaryEquationsVector(6) ;              // x21

//        std::cout<<"x21 = "<<auxiliaryEquationsVector(21)<<std::endl;

        auxiliaryEquationsVector(26) = 2*(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5)+
                                          auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6));              // x26

//        std::cout<<"auxiliaryEquationsVector(26) = "<<auxiliaryEquationsVector(26)<<std::endl;

        // Setting accuracy for x26 to 1 mm^2/s
        if (abs(auxiliaryEquationsVector(26))<=1e-12){
            auxiliaryEquationsVector(26) = 0;
        };

        ///Debug///
        //std::cout<<"Surely this works 2..."<<std::endl;
//        auxiliaryEquationsVector(26) = 2*((auxiliaryEquationsVector(1)/1e6)*(auxiliaryEquationsVector(4))+(auxiliaryEquationsVector(2)/1e6)*(auxiliaryEquationsVector(5))+
//                                          (auxiliaryEquationsVector(3)/1e6)*(auxiliaryEquationsVector(6)))*1e6;              // x26

//      std::cout<<setprecision(15)<<"x26 = "<<auxiliaryEquationsVector(26)<<std::endl;

//        std::cout<<"x1 = "<<auxiliaryEquationsVector(1)<<std::endl;
//        std::cout<<"x1-852.252774466749 = "<<auxiliaryEquationsVector(1)-852.252774466749<<std::endl;
//        std::cout<<setprecision(15)<<"x1*x4 = "<<auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4)<<std::endl;
//        std::cout<<setprecision(15)<<"x1*x4+185.640294486616 = "<<auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4)+185.640294486616<<std::endl;
//        std::cout<<setprecision(15)<<"x2*x5 = "<<auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5)<<std::endl;
//        std::cout<<setprecision(15)<<"x2*x5-185.640294486616 = "<<auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5)-185.640294486616<<std::endl;
//        std::cout<<setprecision(15)<<"x3*x6 = "<<auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6)<<std::endl;
//        std::cout<<setprecision(15)<<"x1*x4+x2*x5 = "<<(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4))+(auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5))<<std::endl;
//        std::cout<<setprecision(15)<<"x1*x4+x2*x5+x3*x6 = "<<auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5)+
//                   auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6)<<std::endl;
//        ///Debug///

        auxiliaryEquationsVector(49) = atan2(auxiliaryEquationsVector(2),auxiliaryEquationsVector(1));      // x 49

//        auxiliaryEquationsVector(11) = atan2(auxiliaryEquationsVector(2),auxiliaryEquationsVector(1))-auxiliaryEquationsVector(10);              // x11
        auxiliaryEquationsVector(11) = auxiliaryEquationsVector(49)-auxiliaryEquationsVector(10);



//        auxiliaryEquationsVector(20) = pow(auxiliaryEquationsVector(8), 0.5);              // x20

        auxiliaryEquationsVector(20) = sqrt(auxiliaryEquationsVector(8));               // x20

//        auxiliaryEquationsVector(36) = pow(auxiliaryEquationsVector(21), 0.5) ;                // x36

        auxiliaryEquationsVector(36) = sqrt(auxiliaryEquationsVector(21));              // x36

        auxiliaryEquationsVector(52) = (auxiliaryEquationsVector(4)+rotationalVelocity*auxiliaryEquationsVector(2))*cos(auxiliaryEquationsVector(10))+
                sin(auxiliaryEquationsVector(10))*(auxiliaryEquationsVector(5)-rotationalVelocity*auxiliaryEquationsVector(1));     // x52

        auxiliaryEquationsVector(53) = -(auxiliaryEquationsVector(4)+rotationalVelocity*auxiliaryEquationsVector(2))*sin(auxiliaryEquationsVector(10))+
                cos(auxiliaryEquationsVector(10))*(auxiliaryEquationsVector(5)-rotationalVelocity*auxiliaryEquationsVector(1));     // x53

        auxiliaryEquationsVector(54) = auxiliaryEquationsVector(4)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(5)*auxiliaryEquationsVector(5)+
                rotationalVelocity*rotationalVelocity*auxiliaryEquationsVector(19)+2.0*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)
                                                                                                         -auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1));     // x54

        auxiliaryEquationsVector(55) = auxiliaryEquationsVector(53)*auxiliaryEquationsVector(50)-auxiliaryEquationsVector(52)*auxiliaryEquationsVector(51);     // x55



/*
        std::cout<<"x21-x21 = "<<auxiliaryEquationsVector(21)-auxiliaryEquationsVector(21)<<std::endl;
        std::cout<<"x21-sqrt(x21*x21) = "<<auxiliaryEquationsVector(21)-sqrt(auxiliaryEquationsVector(21)*auxiliaryEquationsVector(21))<<std::endl;
        std::cout<<"x21-sqrt(x21)*sqrt(x21) = "<<auxiliaryEquationsVector(21)-sqrt(auxiliaryEquationsVector(21))*sqrt(auxiliaryEquationsVector(21))<<std::endl;

        std::cout<<"4-2*2 = "<<4-2*2<<std::endl;
        std::cout<<"4.365416516-sqrt(4.365416516)*sqrt(4.365416516) = "<<4.365416516-sqrt(4.365416516)*sqrt(4.365416516)<<std::endl;

                std::cout<<"x21-x36^2 = "<<auxiliaryEquationsVector(21)-auxiliaryEquationsVector(36)*auxiliaryEquationsVector(36)<<std::endl;
*/




/*        /// Debug ///
        //std::cout<<"Surely this works 3..."<<std::endl;

        std::cout<<"x45 = "<<auxiliaryEquationsVector(45)<<std::endl;
        std::cout<<"x45-7.088e-05 = "<<auxiliaryEquationsVector(45)-7.088e-05<<std::endl;
        std::cout<<"x1*x5-x2*x4 = "<<auxiliaryEquationsVector(1)*auxiliaryEquationsVector(5)-auxiliaryEquationsVector(2)*auxiliaryEquationsVector(4)<<std::endl;
        std::cout<<"x19 = "<<auxiliaryEquationsVector(19)<<std::endl;

        std::cout<<"x1*x5-x2*x4-712.211649225101 = "<<auxiliaryEquationsVector(1)*auxiliaryEquationsVector(5)-auxiliaryEquationsVector(2)*auxiliaryEquationsVector(4)-712.211649225101<<std::endl;
        std::cout<<"x19-10048132.7486611 = "<<auxiliaryEquationsVector(19)-10048132.7486611<<std::endl;
        std::cout<<"(x1*x5-x2*x4-712.211649225101)/(x19-10048132.7486611) = "<<(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(5)-auxiliaryEquationsVector(2)*auxiliaryEquationsVector(4)-712.211649225101)/
                   (auxiliaryEquationsVector(19)-10048132.7486611)<<std::endl;
        std::cout<<"712.211649225101/10048132.7486611 = "<<712.211649225101/10048132.7486611<<std::endl;

        std::cout<<"x45 in Mm = "<<(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(5)*1E-6-auxiliaryEquationsVector(2)*auxiliaryEquationsVector(4)*1E-6)/(auxiliaryEquationsVector(19)*1E-6)<<std::endl;
        std::cout<<"x45 in Mm-7.088e-05 = "<<(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(5)*1E-6-auxiliaryEquationsVector(2)*auxiliaryEquationsVector(4)*1E-6)/(auxiliaryEquationsVector(19)*1E-6)-7.088e-05<<std::endl;
//*/

        auxiliaryEquationsVector(12) = asin(auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)) ;              // x12


        auxiliaryEquationsVector(50) = auxiliaryEquationsVector(20)*cos(auxiliaryEquationsVector(12)*cos(auxiliaryEquationsVector(11)));  // x50

        auxiliaryEquationsVector(51) = auxiliaryEquationsVector(20)*cos(auxiliaryEquationsVector(12)*sin(auxiliaryEquationsVector(11)));  // x51





//        auxiliaryEquationsVector(25) = (auxiliaryEquationsVector(26)*1e-6)/(2*auxiliaryEquationsVector(20)*1e-6);              // x25
        // If the inertial flight path angle is set to 90 degrees, then x25 = x36
        if (verticalInertialFlightPathAngleSet == true){
            auxiliaryEquationsVector(25) = auxiliaryEquationsVector(36);
        }
        else{
        auxiliaryEquationsVector(25) = (auxiliaryEquationsVector(26))/(2.0*auxiliaryEquationsVector(20));              // x25
        }

        // Acount for when the angle was not predefined but is still 90 degrees
        if (verticalInertialFlightPathAngleSet == false && auxiliaryEquationsVector(25) == auxiliaryEquationsVector(36)){
            std::cout<<"verticalInertialFlightPathAngleSet eq 2 = "<<verticalInertialFlightPathAngleSet<<std::endl;
            verticalInertialFlightPathAngleSet = true;
        }

        std::cout<<"verticalInertialFlightPathAngleSet eq 3 = "<<verticalInertialFlightPathAngleSet<<std::endl;

        // Please note that the altitude h (x31) is expressed in km MOLA (which is also the input for the density and temperature curves!)
//        auxiliaryEquationsVector(31) = (auxiliaryEquationsVector(20)-bodyReferenceRadius)/1000;              // x31 [km]!!!
//        auxiliaryEquationsVector(31) = (auxiliaryEquationsVector(20)*1e6-bodyReferenceRadius*1e6)/1e6;              // x31 [km]!!!
        auxiliaryEquationsVector(31) = (auxiliaryEquationsVector(20)-bodyReferenceRadius);              // x31 [km]!!!
/*
        std::cout<<"bodyReferenceRadius = "<<bodyReferenceRadius<<std::endl;
        std::cout<<"Computed initial radius = "<<auxiliaryEquationsVector(20)<<std::endl;
        std::cout<<"Computed altitude in m is "<<auxiliaryEquationsVector(20)-bodyReferenceRadius<<std::endl;
        std::cout<<"Computed altitude in km is "<<auxiliaryEquationsVector(31)<<std::endl;
*/

        // Avoid singularities and account for the case of vertical flight (no change in longitude)
        if (auxiliaryEquationsVector(19) == 0 || verticalInertialFlightPathAngleSet == true){

            auxiliaryEquationsVector(45) = 0;
        }
        else {
        auxiliaryEquationsVector(45) = (auxiliaryEquationsVector(1)*auxiliaryEquationsVector(5)-auxiliaryEquationsVector(2)*auxiliaryEquationsVector(4))/auxiliaryEquationsVector(19);               // x45
//        auxiliaryEquationsVector(45) = ((auxiliaryEquationsVector(1)/1e6)*auxiliaryEquationsVector(5)-(auxiliaryEquationsVector(2)/1e6)*auxiliaryEquationsVector(4))/(auxiliaryEquationsVector(19)/1e6); // x45 in 1e6 km/1e6 km
};


        // If the spacecraft flies over a pole i.e. if +-z=r (or x3=x20) then the change in latitude is undefined and has to be set equal to zero here. This is done to avoid singularities.
        // Also, if the MAV is in vertical flight, there is no change in latitude either

        if (auxiliaryEquationsVector(3)==auxiliaryEquationsVector(20) || -auxiliaryEquationsVector(3)==auxiliaryEquationsVector(20) || verticalInertialFlightPathAngleSet == true){

            auxiliaryEquationsVector(24) = 0;
        }
        else {

        auxiliaryEquationsVector(24) = (auxiliaryEquationsVector(20)*auxiliaryEquationsVector(6)-auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25))/
                (auxiliaryEquationsVector(8)*sqrt(1-(auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20))*(auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20))));              // x24

        };

//        std::cout<<"Altitude = "<<auxiliaryEquationsVector(31)<<std::endl;

//        std::cout<<"Lowest altitude range = "<<temperatureAltitudeRanges(0,0)<<std::endl;

        // Computing the polynomial fit using the altitude and fit parameters for density
        for (int i = 0; i < 10+1;i++) {

        auxiliaryEquationsVector(30) += pow(auxiliaryEquationsVector(31),i)*densityPolyCoefficients(i);              // x30
};

        // Determine which section of the temperature curve needs to be used and what the corresponding order is
        // Also, because a computer is less than perfect, a small correction is made to the lower bound of the first section to make sure that the initial altitude is still valid

        if ((temperatureAltitudeRanges(0,0)-0.000000000001) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(0,1)){

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
        //std::cout<<"Surely this works 4..."<<std::endl;

        // Computing the polynomial fit using the altitude and fit parameters for temperature
        for (int i=0; i < powerT+1;i++){

        auxiliaryEquationsVector(34) += pow(auxiliaryEquationsVector(31),i)*temperaturePolyCoefficients(sectionT,i);              // x34

//        std::cout<<"x34 interval "<<i<<" = "<<auxiliaryEquationsVector(34)<<std::endl;

};
        // Avoid cosine round-off errors and acounting for non-rotating Mars
        if (abs(cos(auxiliaryEquationsVector(12)))<6.2e-17 || rotationalVelocity == 0.0){
            auxiliaryEquationsVector(35) = 0;
        }
        else {
        auxiliaryEquationsVector(35) = rotationalVelocity*auxiliaryEquationsVector(20)*cos(auxiliaryEquationsVector(12));                // x35
}
        // Avoid singularities
        if (auxiliaryEquationsVector(36) == 0){
            auxiliaryEquationsVector(37) = 0;
        }
        // And dealing with round-off errors
        else if ((auxiliaryEquationsVector(25)/auxiliaryEquationsVector(36)) < -1.0 && ((auxiliaryEquationsVector(25)/auxiliaryEquationsVector(36))+1.0) > -1e-15){
            auxiliaryEquationsVector(37) = -1.0;
            std::cout<<"Rounding x37 to -1"<<std::endl;
        }
        else if ((auxiliaryEquationsVector(25)/auxiliaryEquationsVector(36)) > 1.0 && ((auxiliaryEquationsVector(25)/auxiliaryEquationsVector(36))-1.0) < 1e-15){
            auxiliaryEquationsVector(37) = 1.0;
            std::cout<<"Rounding x37 to 1"<<std::endl;
        }
        // And dealing with vertical flight where x37 is 1
        else if (verticalInertialFlightPathAngleSet == true){
            auxiliaryEquationsVector(37) = 1.0;
        }
        else {
        auxiliaryEquationsVector(37) = auxiliaryEquationsVector(25)/auxiliaryEquationsVector(36);                // x37
}
//        /// Debug ///

//        std::cout<<"(auxiliaryEquationsVector(25)/auxiliaryEquationsVector(36)) = "<<(auxiliaryEquationsVector(25)/auxiliaryEquationsVector(36))<<std::endl;
//        std::cout<<"x25-x36 = "<<auxiliaryEquationsVector(25)-auxiliaryEquationsVector(36)<<std::endl;

//        std::cout<<"x37 = "<<auxiliaryEquationsVector(37)<<std::endl;

//        std::cout<<"x36 = "<<auxiliaryEquationsVector(36)<<std::endl;
//        std::cout<<"x21 = "<<auxiliaryEquationsVector(21)<<std::endl;
//        std::cout<<"sqrt(x21) = "<<sqrt(auxiliaryEquationsVector(21))<<std::endl;
//        std::cout<<"x4 = "<<auxiliaryEquationsVector(4)<<std::endl;
//        std::cout<<"x4-x36 = "<<auxiliaryEquationsVector(4)-auxiliaryEquationsVector(36)<<std::endl;
//        std::cout<<"x25 = "<<auxiliaryEquationsVector(25)<<std::endl;
//        std::cout<<"x20 = "<<auxiliaryEquationsVector(20)<<std::endl;
//        std::cout<<"x1 = "<<auxiliaryEquationsVector(1)<<std::endl;
//        std::cout<<"x1-x20 = "<<auxiliaryEquationsVector(1)-auxiliaryEquationsVector(20)<<std::endl;
//        std::cout<<"x26 = "<<auxiliaryEquationsVector(26)<<std::endl;
//        std::cout<<"2*x1*x4 = "<<2*auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4)<<std::endl;
//        std::cout<<"x25-x4 = "<<auxiliaryEquationsVector(25)-auxiliaryEquationsVector(4)<<std::endl;
//        std::cout<<"2*x1*x4/2*x1 = "<<(2*auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4))/(2*auxiliaryEquationsVector(1))<<std::endl;
//        std::cout<<"(2*x1*x4/2*x1)-x4 = "<<((2*auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4))/(2*auxiliaryEquationsVector(1)))-auxiliaryEquationsVector(4)<<std::endl;
//        std::cout<<"2*x4/2 = "<<(2*auxiliaryEquationsVector(4))/(2)<<std::endl;
//        std::cout<<"(2*x4/2)-x4 = "<<((2*auxiliaryEquationsVector(4))/(2))-auxiliaryEquationsVector(4)<<std::endl;
//        std::cout<<"x26/(2*x20)-x4 = "<<(auxiliaryEquationsVector(26)/(2*auxiliaryEquationsVector(20)))-auxiliaryEquationsVector(4)<<std::endl;
//        std::cout<<"((x1*x4)/(x1))-x4 = "<<(((auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4)))/(auxiliaryEquationsVector(1)))-auxiliaryEquationsVector(4)<<std::endl;
//        std::cout<<"((x1*x4*1e-6)/(x1*1e-6))-x4 = "<<(((auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4))*1e-6)/(auxiliaryEquationsVector(1)*1e-6))-auxiliaryEquationsVector(4)<<std::endl;
//        std::cout<<"(x1/x1)*x4-x4 = "<<(auxiliaryEquationsVector(1)/auxiliaryEquationsVector(1))*auxiliaryEquationsVector(4)-auxiliaryEquationsVector(4)<<std::endl;

//        /// Debug ///

        if (rotationalVelocity == 0.0){ // For non-rotating Mars
            auxiliaryEquationsVector(46) = auxiliaryEquationsVector(45);
        }
        else {
        auxiliaryEquationsVector(46) = auxiliaryEquationsVector(45)-rotationalVelocity;             // x46
        }

//        // Set tolerance for the velocity in case of rounding errors... It is set such that the the accuracy is 10 micro-metres/sec of surface movement
//        if (abs(auxiliaryEquationsVector(46))<=1e-11){
//            auxiliaryEquationsVector(46) = 0;
//        }

//std::cout<<"Surely this works 5..."<<std::endl;
        // Avoid cosine rounding errors and account for vertical flight where x45 is zero
                if (abs(cos(auxiliaryEquationsVector(12)))<6.2e-17 || verticalInertialFlightPathAngleSet == true){
                  auxiliaryEquationsVector(47) = 0;
                }
                else {
        auxiliaryEquationsVector(47) = cos(auxiliaryEquationsVector(12))*auxiliaryEquationsVector(45);               // x47
                }

//        auxiliaryEquationsVector(18) = auxiliaryEquationsVector(20)*auxiliaryEquationsVector(24);              // x18


        auxiliaryEquationsVector(28) = exp(auxiliaryEquationsVector(30));              // x28

        auxiliaryEquationsVector(33) = sqrt(adiabeticIndex*specificGasConstant*auxiliaryEquationsVector(34));              // x33


//        std::cout<<"x33 (or speed of sound) = "<<auxiliaryEquationsVector(33)<<std::endl;

        if (verticalInertialFlightPathAngleSet == true ){ // Set equal to 90 degrees if vertical launch
            auxiliaryEquationsVector(38) = tudat::mathematical_constants::LONG_PI/2.0;
        }
        else if (inertialFlightPathAngleSet == true && time == 0.0){ // If the inertial flight path angle was set initially then x38 is that angle for t=0
            auxiliaryEquationsVector(38) = inertialFlightPathAngle;
        }
        else {
        auxiliaryEquationsVector(38) = asin(auxiliaryEquationsVector(37));                // x38
        }

        std::cout<<"x38 = "<<auxiliaryEquationsVector(38)<<std::endl;

        if (verticalInertialFlightPathAngleSet == true){ // If vertical flight then set the heading angle to zero
            auxiliaryEquationsVector(40) = 0.0;
        }
        else if (time == 0.0 && inertialHeadingAngleSet == true){ // if not vertical flight but heading is still defined, then set heading angle for t = 0.0
            auxiliaryEquationsVector(40) = inertialHeadingAngle;
        }
        else {
        auxiliaryEquationsVector(40) = atan2(auxiliaryEquationsVector(47),auxiliaryEquationsVector(24));               // x40
        }

        std::cout<<"x40 = "<<auxiliaryEquationsVector(40)<<std::endl;

        if (rotationalVelocity == 0.0 || auxiliaryEquationsVector(35) == 0.0){
            auxiliaryEquationsVector(41) = 0.0;
        }
        else {
        auxiliaryEquationsVector(41) = auxiliaryEquationsVector(35)*auxiliaryEquationsVector(36);                // x41
        }

        // Avoid cosine rounding errors
                if (abs(cos(auxiliaryEquationsVector(12)))<6.2e-17){
                  auxiliaryEquationsVector(48) = 0;
                }
                else {
        auxiliaryEquationsVector(48) = cos(auxiliaryEquationsVector(12))*auxiliaryEquationsVector(46);               // x48
                }
//std::cout<<"Surely this works 6..."<<std::endl;
        /*
//        auxiliaryEquationsVector(44) = auxiliaryEquationsVector(36)*cos(auxiliaryEquationsVector(38));                // x44

//        auxiliaryEquationsVector(39) = auxiliaryEquationsVector(18)/auxiliaryEquationsVector(44);                // x39

//        auxiliaryEquationsVector(40) = acos(auxiliaryEquationsVector(39));                // x40
*/


        // Avoid cosine rounding errors
                if (abs(cos(auxiliaryEquationsVector(38)))<6.2e-17){
                  auxiliaryEquationsVector(42) = 0;
                }
                else if (verticalInertialFlightPathAngleSet == true || auxiliaryEquationsVector(40) == 0){// Vertical ascent
                    auxiliaryEquationsVector(42) = 0;
                    }
                else {
        auxiliaryEquationsVector(42) = cos(auxiliaryEquationsVector(38))*sin(auxiliaryEquationsVector(40));                // x42
}


//        /// Debug ///

             //std::cout<<"Surely this works 6.1..."<<std::endl;

//        std::cout<<"cos(x38) = "<<cos(auxiliaryEquationsVector(38))<<std::endl;
//                std::cout<<"x40 = "<<auxiliaryEquationsVector(40)<<std::endl;
//        std::cout<<"sin(x40) = "<<sin(auxiliaryEquationsVector(40))<<std::endl;
//        std::cout<<"cos(x38)*sin(x40) = "<<cos(auxiliaryEquationsVector(38))*sin(auxiliaryEquationsVector(40))<<std::endl;
//        std::cout<<"x42 = "<<auxiliaryEquationsVector(42)<<std::endl;
//        std::cout<<"x42-cos(x38)*sin(x40) = "<<auxiliaryEquationsVector(42)-cos(auxiliaryEquationsVector(38))*sin(auxiliaryEquationsVector(40))<<std::endl;
//        std::cout<<"cos(x38)*sin(x40)-0.99999999989952-4.44089209850063e-16 = "<<cos(auxiliaryEquationsVector(38))*sin(auxiliaryEquationsVector(40))-0.99999999989952-4.44089209850063e-16<<std::endl;
//        std::cout<<"x42-0.99999999989952-4.44089209850063e-16 = "<<auxiliaryEquationsVector(42)-0.99999999989952-4.44089209850063e-16<<std::endl;

//        /// Debug ///

//        auxiliaryEquationsVector(43) = auxiliaryEquationsVector(41)*auxiliaryEquationsVector(42);                // x43

                // Avoid cosine rounding errors
                        if (abs(cos(auxiliaryEquationsVector(38)))<6.2e-17){
                          auxiliaryEquationsVector(43) = 0;
                        }
                        else if (verticalInertialFlightPathAngleSet == true || auxiliaryEquationsVector(41) == 0.0 || auxiliaryEquationsVector(42) == 0.0){ // Vertical ascent
                                auxiliaryEquationsVector(43) = 0.0;
                        }
                        else if (auxiliaryEquationsVector(41) == 0.0 || rotationalVelocity == 0.0){ // For non-rotating Mars
                            auxiliaryEquationsVector(43) = 0.0;
                        }
                        else {
        auxiliaryEquationsVector(43) = auxiliaryEquationsVector(41)*cos(auxiliaryEquationsVector(38))*sin(auxiliaryEquationsVector(40));                // x43
}

                        //std::cout<<"Surely this works 6.2..."<<std::endl;
                        std::cout<<"x43 = "<<auxiliaryEquationsVector(43)<<std::endl;
//                        std::cout<<"u43 = "<<auxiliaryDerivativesVector(43)<<std::endl;
                        std::cout<<"verticalInertialFlightPathAngleSet eq 4 = "<<verticalInertialFlightPathAngleSet<<std::endl;
                        std::cout<<"verticalRotationalFlightPathAngleSet eq 2 = "<<verticalRotationalFlightPathAngleSet<<std::endl;
                        std::cout<<"x35 = "<<auxiliaryEquationsVector(35)<<std::endl;
                        std::cout<<"x21+Omega_M^2*x19+2*Omega_M*(x4*x2-x5*x1) = "<<auxiliaryEquationsVector(21)+rotationalVelocity*rotationalVelocity*auxiliaryEquationsVector(19)+
                                   2.0*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)-auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1))<<std::endl;

                        std::cout<<"x15 = "<<sqrt(auxiliaryEquationsVector(21)+rotationalVelocity*rotationalVelocity*auxiliaryEquationsVector(19)+
                                                  2.0*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)-auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1)))<<std::endl;
                        std::cout<<"x25 = "<<auxiliaryEquationsVector(25)<<std::endl;

                        // Account for rounding errors
                        if (auxiliaryEquationsVector(21)+rotationalVelocity*rotationalVelocity*auxiliaryEquationsVector(19)+
                                2.0*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)-auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1))<0){
                            auxiliaryEquationsVector(15) = 0.0;
                        }
                        else{
                        auxiliaryEquationsVector(15) = sqrt(auxiliaryEquationsVector(21)+rotationalVelocity*rotationalVelocity*auxiliaryEquationsVector(19)+
                                                            2.0*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)-auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1)));        // x15
                        }

                        std::cout<<"x25-x15 = "<<auxiliaryEquationsVector(25)-auxiliaryEquationsVector(15)<<std::endl;
                        std::cout<<"x25/x15 = "<<auxiliaryEquationsVector(25)/auxiliaryEquationsVector(15)<<std::endl;

//        // If vertical ascent
//                        if (verticalInertialFlightPathAngleSet == true || auxiliaryEquationsVector(43) == 0.0){ //std::cout<<"It goes here ..."<<std::endl;
//                            if (auxiliaryEquationsVector(35) == 0.0 || rotationalVelocity == 0.0){ // And non-rotating Mars
////                                std::cout<<"It goes here 1"<<std::endl;
//                                auxiliaryEquationsVector(15) = auxiliaryEquationsVector(36);

//                            }
//                            else{
////                                std::cout<<"It goes here 2"<<std::endl;
//                            auxiliaryEquationsVector(15) = sqrt(auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)+auxiliaryEquationsVector(21));

//                        }}
//                        else if (verticalRotationalFlightPathAngleSet == true){ // Vertical flight
////                            std::cout<<"It goes here 3"<<std::endl;
//                                auxiliaryEquationsVector(15) = auxiliaryEquationsVector(25);


//                        }

//                        else{
////                            std::cout<<"It goes here 4"<<std::endl;
//        auxiliaryEquationsVector(15) = sqrt(auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)+auxiliaryEquationsVector(21)-2.0*auxiliaryEquationsVector(43));              // x15

//                        }

                        //std::cout<<"Surely this works 6.3..."<<std::endl;

         /*               /// Debug ///

                        std::cout<<"x25 = "<<auxiliaryEquationsVector(25)<<std::endl;
                        std::cout<<"x15 = "<<auxiliaryEquationsVector(15)<<std::endl;
                        std::cout<<"x15 new = "<<sqrt(auxiliaryEquationsVector(21)+rotationalVelocity*rotationalVelocity*(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+
                                                                                                                          auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2))+
                                                      2.0*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)-auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1)))<<std::endl;
                        std::cout<<"x15 new ini sqr = "<<(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(5)*auxiliaryEquationsVector(5)+
                                                          auxiliaryEquationsVector(6)*auxiliaryEquationsVector(6))+rotationalVelocity*rotationalVelocity*(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+
                                                                                                                         auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2))+
                                                     2*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)-auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1))<<std::endl;
                        std::cout<<"x15 new ini = "<<sqrt((auxiliaryEquationsVector(4)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(5)*auxiliaryEquationsVector(5)+
                                                           auxiliaryEquationsVector(6)*auxiliaryEquationsVector(6))+rotationalVelocity*rotationalVelocity*(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+
                                                                                                                          auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2))+
                                                      2.0*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)-auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1)))<<std::endl;
                        std::cout<<"x15 new ini (VI = VM) sqr = "<<(2.0*rotationalVelocity*rotationalVelocity*(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+
                                                                                                                          auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2))+
                                                      2.0*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)-auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1)))<<std::endl;
                        std::cout<<"x25 ini = "<<(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5)+
                                                  auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6))/(sqrt(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2)+
                                                                                                                 auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3)))<<std::endl;
                        std::cout<<"x25 ini / x15 new ini = "<<(((auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5)+
                                                 auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6))/(sqrt(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2)+
                                                                                                                auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3))))/
                                   (sqrt((auxiliaryEquationsVector(4)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(5)*auxiliaryEquationsVector(5)+
                                                                                                                                                                                  auxiliaryEquationsVector(6)*auxiliaryEquationsVector(6))+rotationalVelocity*rotationalVelocity*(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+
                                                                                                                                                                                                                                                 auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2))+
                                                                                                                                                                             2.0*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)-auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1)))))<<std::endl;
                        std::cout<<"FPA_R = "<<asin(((auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5)+
                                                 auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6))/(sqrt(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2)+
                                                                                                                auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3))))/
                                   (sqrt((auxiliaryEquationsVector(4)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(5)*auxiliaryEquationsVector(5)+
                                                                                                                                                                                  auxiliaryEquationsVector(6)*auxiliaryEquationsVector(6))+rotationalVelocity*rotationalVelocity*(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+
                                                                                                                                                                                                                                                 auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2))+
                                                                                                                                                                             2.0*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)-auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1)))))<<std::endl;
                        std::cout<<"x25-x15 = "<<auxiliaryEquationsVector(25)-auxiliaryEquationsVector(15)<<std::endl;

          //*/              /// Debug ///


                        // Acount for when the angle was not predefined but is still 90 degrees
                        if (verticalRotationalFlightPathAngleSet == false && auxiliaryEquationsVector(25) == auxiliaryEquationsVector(15)){
                            std::cout<<"verticalRotationalFlightPathAngleSet eq 3 = "<<verticalRotationalFlightPathAngleSet<<std::endl;
                            verticalRotationalFlightPathAngleSet = true;
                            std::cout<<"verticalRotationalFlightPathAngleSet eq 4 = "<<verticalRotationalFlightPathAngleSet<<std::endl;
                        }

//        // Set tolerance for the velocity in case of rounding errors... It is set such that the the accuracy is 10 micro-metres/sec
//        if (sqrt(abs(auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)+auxiliaryEquationsVector(21)-2*auxiliaryEquationsVector(43)))<=1E-8){

//            auxiliaryEquationsVector(15) = 0;


//        }

        /// Debug ///

//        std::cout<<"r_dot_I-V_I*sin(gamma_I) = "<<auxiliaryEquationsVector(25)-auxiliaryEquationsVector(36)*sin(auxiliaryEquationsVector(38))<<std::endl;
//        std::cout<<"x15 = "<<auxiliaryEquationsVector(15)<<std::endl;
//        std::cout<<"x25/x15 = "<<auxiliaryEquationsVector(25)/auxiliaryEquationsVector(15)<<std::endl;
//        std::cout<<"r_dot_R/x15 = "<<(auxiliaryEquationsVector(25)-auxiliaryEquationsVector(36)*sin(auxiliaryEquationsVector(38)))/auxiliaryEquationsVector(15)<<std::endl;


        //std::cout<<"Surely this works 7..."<<std::endl;

//        std::cout<<"x15 = "<<auxiliaryEquationsVector(15)<<std::endl;
//        else {
//        auxiliaryEquationsVector(15) = sqrt(auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)+auxiliaryEquationsVector(21)-2*auxiliaryEquationsVector(43));              // x15
//        auxiliaryEquationsVector(15) = sqrt(((auxiliaryEquationsVector(35)/1e6)*(auxiliaryEquationsVector(35))+(auxiliaryEquationsVector(21)/1e6)-2*(auxiliaryEquationsVector(43)/1e6))*1e6);              // x15
//};



/* // What the actual f*ck?!
        std::cout<<"((auxiliaryEquationsVector(35)/1e6)*(auxiliaryEquationsVector(35)/1e6)+(auxiliaryEquationsVector(21)/1e6)-2*(auxiliaryEquationsVector(43)/1e6))*1e6 = "<<((auxiliaryEquationsVector(35)/1e6)*(auxiliaryEquationsVector(35)/1e6)+(auxiliaryEquationsVector(21)/1e6)-2*(auxiliaryEquationsVector(43)/1e6))*1e6<<std::endl;
        std::cout<<"(auxiliaryEquationsVector(35)/1e6)*(auxiliaryEquationsVector(35)) = "<<(auxiliaryEquationsVector(35)/1e6)*(auxiliaryEquationsVector(35))<<std::endl;
        std::cout<<"(auxiliaryEquationsVector(21)/1e6) = "<<(auxiliaryEquationsVector(21)/1e6)<<std::endl;
        std::cout<<"2*(auxiliaryEquationsVector(43)/1e6) = "<<2*(auxiliaryEquationsVector(43)/1e6)<<std::endl;
        std::cout<<"(auxiliaryEquationsVector(35)/1e6) = "<<(auxiliaryEquationsVector(35)/1e6)<<std::endl;

        std::cout<<"x25 = "<<auxiliaryEquationsVector(25)<<std::endl;
        std::cout<<"x26 = "<<auxiliaryEquationsVector(26)<<std::endl;
        std::cout<<"x35 = "<<auxiliaryEquationsVector(35)<<std::endl;
        std::cout<<"x36 = "<<auxiliaryEquationsVector(36)<<std::endl;
        std::cout<<"x35-x36 = "<<auxiliaryEquationsVector(35)-auxiliaryEquationsVector(36)<<std::endl;
        std::cout<<"x42 = "<<auxiliaryEquationsVector(42)<<std::endl;
        std::cout<<"x41 = "<<auxiliaryEquationsVector(41)<<std::endl;
        std::cout<<"x35*x36 = "<<auxiliaryEquationsVector(35)*auxiliaryEquationsVector(36)<<std::endl;
        std::cout<<"x35^2-x43 = "<<auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)-auxiliaryEquationsVector(43)<<std::endl;
        std::cout<<"x21-x43 = "<<auxiliaryEquationsVector(21)-auxiliaryEquationsVector(43)<<std::endl;
        std::cout<<"x35*x35-x43+x21-x43 = "<<auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)-auxiliaryEquationsVector(43)+auxiliaryEquationsVector(21)-auxiliaryEquationsVector(43)<<std::endl;
        std::cout<<"x35*x35+x21-x43-x43 = "<<auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)+auxiliaryEquationsVector(21)-auxiliaryEquationsVector(43)-auxiliaryEquationsVector(43)<<std::endl;
        std::cout<<"x35*x35+x21-2*x43 = "<<auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)+auxiliaryEquationsVector(21)-2*auxiliaryEquationsVector(43)<<std::endl;
        std::cout<<"x35*x35-x43-x43+x21 = "<<auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)-auxiliaryEquationsVector(43)-auxiliaryEquationsVector(43)+auxiliaryEquationsVector(21)<<std::endl;
        std::cout<<"x35*x35+x21-2*x43 = "<<auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)+auxiliaryEquationsVector(21)-2*auxiliaryEquationsVector(43)<<std::endl;
        std::cout<<"x35*x35+x21-2*41*42 = "<<auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)+auxiliaryEquationsVector(21)-2*auxiliaryEquationsVector(41)*auxiliaryEquationsVector(42)<<std::endl;
        std::cout<<"x35*x35+x21-2*41*cos(x38)*sin(x40) = "<<auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)+auxiliaryEquationsVector(21)-2*auxiliaryEquationsVector(41)*cos(auxiliaryEquationsVector(38))*sin(auxiliaryEquationsVector(40))<<std::endl;
        std::cout<<"x35*x35+x21-2*41*sin(x40)*cos(x38) = "<<auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)+auxiliaryEquationsVector(21)-2*auxiliaryEquationsVector(41)*sin(auxiliaryEquationsVector(40))*cos(auxiliaryEquationsVector(38))<<std::endl;
        std::cout<<"-2*41*42 = "<<-2*auxiliaryEquationsVector(41)*auxiliaryEquationsVector(42)<<std::endl;
        std::cout<<"-2*41*cos(x38)*sin(x40) = "<<-2*auxiliaryEquationsVector(41)*cos(auxiliaryEquationsVector(38))*sin(auxiliaryEquationsVector(40))<<std::endl;
        std::cout<<"x35*x35+x21 = "<<auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)+auxiliaryEquationsVector(21)<<std::endl;
        std::cout<<"x35*x35+x21-0.10096312339415 = "<<auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)+auxiliaryEquationsVector(21)-0.10096312339415<<std::endl;

        std::cout<<"-2*41*42+0.10096312339415 = "<<-2*auxiliaryEquationsVector(41)*auxiliaryEquationsVector(42)+0.10096312339415<<std::endl;
        std::cout<<"-2*41*cos(x38)*sin(x40)+0.10096312339415 = "<<-2*auxiliaryEquationsVector(41)*cos(auxiliaryEquationsVector(38))*sin(auxiliaryEquationsVector(40))+0.10096312339415<<std::endl;


        std::cout<<"x35 = "<<auxiliaryEquationsVector(35)<<std::endl;
        std::cout<<"x35^2 = "<<auxiliaryEquationsVector(35)*auxiliaryEquationsVector(35)<<std::endl;


        std::cout<<"x21 = "<<auxiliaryEquationsVector(21)<<std::endl;
        std::cout<<"x43 = "<<auxiliaryEquationsVector(43)<<std::endl;

        std::cout<<"x15 (or V_R) = "<<auxiliaryEquationsVector(15)<<std::endl;
//*/


        // If the velocity of the MAV in the rotating frame V_R (or x15) = 0 m/s then the flight path angle in that frame (gamma_R or x14) is not defined.
        // Therefore to avoid singularities, this flight path angle is set equal to 90 degrees: x23 = 1.

        if (auxiliaryEquationsVector(15)==0){

            auxiliaryEquationsVector(23) = 1.0;
        }
        // And dealing with round-off errors
        else if ((auxiliaryEquationsVector(25)/auxiliaryEquationsVector(15)) < -1.0 && ((auxiliaryEquationsVector(25)/auxiliaryEquationsVector(15))+1.0) > -1e-15){
            auxiliaryEquationsVector(23) = -1.0;
            std::cout<<"Rounding x23 to -1"<<std::endl;
        }
        else if ((auxiliaryEquationsVector(25)/auxiliaryEquationsVector(15)) > 1.0 && ((auxiliaryEquationsVector(25)/auxiliaryEquationsVector(15))-1.0) < 1e-15){
            auxiliaryEquationsVector(23) = 1.0;
            std::cout<<"Rounding x23 to 1"<<std::endl;
        }
        // And dealing with vertical flight where x23 is 1
        else if (verticalRotationalFlightPathAngleSet == true){
            auxiliaryEquationsVector(23) = 1.0;
        }
        else {
        auxiliaryEquationsVector(23) = auxiliaryEquationsVector(25)/auxiliaryEquationsVector(15);              // x23
};
//        auxiliaryEquationsVector(23) = ;  // Debug purposes!

/// Debug ///
//std::cout<<"Surely this works 8..."<<std::endl;
//        std::cout<<"(auxiliaryEquationsVector(25)/auxiliaryEquationsVector(15)) = "<<(auxiliaryEquationsVector(25)/auxiliaryEquationsVector(15))<<std::endl;

//std::cout<<"x23 = "<<auxiliaryEquationsVector(23)<<std::endl;
//std::cout<<"x23^2 = "<<auxiliaryEquationsVector(23)*auxiliaryEquationsVector(23)<<std::endl;
//std::cout<<"x23^2-1 = "<<auxiliaryEquationsVector(23)*auxiliaryEquationsVector(23)-1<<std::endl;
//std::cout<<"x23 - 0.999999999899838 = "<<auxiliaryEquationsVector(23)-0.999999999899838<<std::endl;
//std::cout<<"x15 = "<<auxiliaryEquationsVector(15)<<std::endl;
//std::cout<<"x25 = "<<auxiliaryEquationsVector(25)<<std::endl;
//std::cout<<"The whole vector = "<<auxiliaryEquationsVector<<std::endl;

/// Debug ///

        if (verticalRotationalFlightPathAngleSet == true){ // Set equal to 90 degrees if vertical launch
            auxiliaryEquationsVector(14) = tudat::mathematical_constants::LONG_PI/2.0;
        }
        else if (rotationalFlightPathAngleSet == true && time == 0.0){ // If the rotational flight path angle was set initially then x14 is that angle for t=0
            auxiliaryEquationsVector(14) = rotationalFlightPathAngle;
        }
        else if (rotationalVelocity == 0.0){ // For non-rotating Mars
            auxiliaryEquationsVector(14) = auxiliaryEquationsVector(38);
        }
        else{
        auxiliaryEquationsVector(14) = asin(auxiliaryEquationsVector(23));              // x14
        }

        if (verticalRotationalFlightPathAngleSet == true){ // Heading angle defined 0 if vertical flight
            auxiliaryEquationsVector(13) = 0.0;
        }
        else if (rotationalHeadingAngleSet == true && time == 0.0){ // if not vertical flight but heading angle is still defined then set equal to heading at t = 0
            auxiliaryEquationsVector(13) = rotationalHeadingAngle;
        }
        else if (rotationalVelocity == 0.0){ // For non-rotating Mars
            auxiliaryEquationsVector(13) = auxiliaryEquationsVector(40);
        }
        else{
        auxiliaryEquationsVector(13) =  atan2(auxiliaryEquationsVector(48),auxiliaryEquationsVector(24));                   // x13
        }

//        // Dealing with the inaccuracy in pi
//        if (auxiliaryEquationsVector(23) == 1 || auxiliaryEquationsVector(23) == -1){
//            auxiliaryEquationsVector(16) = 0;
//        }
//        else{
//        auxiliaryEquationsVector(16) = cos(auxiliaryEquationsVector(14));              // x16
//        };
//        // Avoid cosine rounding errors
//                if (abs(auxiliaryEquationsVector(16))<6.2e-17){
//                  auxiliaryEquationsVector(16) = 0;
//                }

        auxiliaryEquationsVector(32) = auxiliaryEquationsVector(15)/auxiliaryEquationsVector(33);              // x32

//        std::cout<<"x32 (or Mach number) = "<<auxiliaryEquationsVector(32)<<std::endl;

//        auxiliaryEquationsVector(17) = auxiliaryEquationsVector(16)*auxiliaryEquationsVector(15);              // x17

        // Determine which section of the drag coefficient curve needs to be used

        for (int i=0; i < 5+1; i++){

            if (dragCoefficientMachRanges(i,0) <= auxiliaryEquationsVector(32) && auxiliaryEquationsVector(32) < dragCoefficientMachRanges(i,1)){

                sectionCD = i;

//                std::cout<<"Required section for CD is section "<<sectionCD<<std::endl;
            }


        };

//        std::cout<<"As I said sectionCD is "<<sectionCD<<std::endl;

//        std::cout<<"Yes?"<<std::endl;

        auxiliaryEquationsVector(29) = dragCoefficientPolyCoefficients(sectionCD,1)*auxiliaryEquationsVector(32)+dragCoefficientPolyCoefficients(sectionCD,0);              // x29
//        std::cout<<"No?"<<std::endl;

//std::cout<<"Surely this works 9..."<<std::endl;
        /*
        // Similar to the flight path angle, if V_R (or x15) = 0 m/s then the azimuth angle is not defined. It is also not defined if the flight path angle is +-90 degrees (i.e. if x23 = +-1)
        // To avoid singularities x22 is set to be equal to 1, such that the azimuth angle is 0.

//        if (auxiliaryEquationsVector(15) == 0 || auxiliaryEquationsVector(23) == 1 || auxiliaryEquationsVector(23) == -1){

//            auxiliaryEquationsVector(22) = 1;
//        }
//        else {

//        auxiliaryEquationsVector(22) = auxiliaryEquationsVector(18)/auxiliaryEquationsVector(17);              // x22
//}; */


//std::cout<<"Surely this works 10..."<<std::endl;
        auxiliaryEquationsVector(27) = 0.5*referenceArea*auxiliaryEquationsVector(28)*auxiliaryEquationsVector(15)*auxiliaryEquationsVector(15)*auxiliaryEquationsVector(29);              // x27

        /// Debug ///

//        std::cout<<"x27 = "<<auxiliaryEquationsVector(27)<<std::endl;
//        std::cout<<"referenceArea = "<<referenceArea<<std::endl;

        /// Debug ///
//std::cout<<"Surely this works 11..."<<std::endl;
//        auxiliaryEquationsVector(13) = acos(auxiliaryEquationsVector(22));              // x13

        auxiliaryEquationsVector(0) = thrustAccelerationsBframe(0)-(auxiliaryEquationsVector(27)/auxiliaryEquationsVector(7));              // w4,2



// Set vertical ascent to false again
//        verticalInertialFlightPathAngleSet = false;
        verticalInertialFlightPathAngleSet = NULL;
        std::cout<<"verticalInertialFlightPathAngleSet eq 5 = "<<verticalInertialFlightPathAngleSet<<std::endl;
//        verticalRotationalFlightPathAngleSet = false;
        verticalRotationalFlightPathAngleSet = NULL;
        std::cout<<"verticalRotationalFlightPathAngleSet eq 5 = "<<verticalRotationalFlightPathAngleSet<<std::endl;


// auxiliaryEquationsVector() = ;               // x


//std::cout<<"Surely this works end..."<<std::endl;


       return auxiliaryEquationsVector;
    }

//////////////////////////////////////////////// Auxiliary Derivatives //////////////////////////////////////////////////////////////////////

    Eigen::VectorXd getAuxiliaryDerivatives( const tudat::basic_mathematics::Vector7d& aState, const double time, const Eigen::Vector3d& thrustAccelerationsBframe, const Eigen::VectorXd& auxiliaryEquationsVector){

        std::cout<<"verticalInertialFlightPathAngleSet der 1 = "<<verticalInertialFlightPathAngleSet<<std::endl;
        std::cout<<"verticalRotationalFlightPathAngleSet der 1 = "<<verticalRotationalFlightPathAngleSet<<std::endl;

    auxiliaryDerivativesVector = Eigen::VectorXd::Zero(56);       // Setting the complete vector and filling it with zeros for now

    // The following expressions are described in the order in which the equations have to be computed corresponding to the respective vector entry
    // Which in this case means that the first entry of the vector is 0 and is not used.

//    auxiliaryDerivativesVector(10) = 0;
    auxiliaryDerivativesVector(10) = rotationalVelocity;                // u10

    // Avoid cosine rounding errors
            if (abs(cos(auxiliaryEquationsVector(12)))<6.2e-17){
              auxiliaryDerivativesVector(35) = rotationalVelocity*(-auxiliaryEquationsVector(20)*auxiliaryEquationsVector(24)*sin(auxiliaryEquationsVector(12)));
            }
            else if (rotationalVelocity == 0.0){ // For non-rotating Mars
                auxiliaryDerivativesVector(35) = 0.0;
            }
            else {
    auxiliaryDerivativesVector(35) = rotationalVelocity*(cos(auxiliaryEquationsVector(12))*auxiliaryEquationsVector(25)-auxiliaryEquationsVector(20)*auxiliaryEquationsVector(24)*sin(auxiliaryEquationsVector(12)));                // u35
}
    auxiliaryDerivativesVector(1) = auxiliaryEquationsVector(4);                // u1

    auxiliaryDerivativesVector(2) = auxiliaryEquationsVector(5);                // u2

    auxiliaryDerivativesVector(3) = auxiliaryEquationsVector(6);                // u3



    // Avoid cosine rounding errors
    double cx10x11;

//    if (abs(cos(auxiliaryEquationsVector(10)+auxiliaryEquationsVector(11)))<6.2e-17){
//        cx10x11 = 0;
//    }
//    else {
//        cx10x11 = cos(auxiliaryEquationsVector(10)+auxiliaryEquationsVector(11));
//    }

    if (abs(cos(auxiliaryEquationsVector(49)))<6.2e-17){  // Using x49
        cx10x11 = 0;
    }
    else {
        cx10x11 = cos(auxiliaryEquationsVector(49));
    }

    double cx12;

    if (abs(cos(auxiliaryEquationsVector(12)))<6.2e-17){
        cx12 = 0;
    }
    else {
        cx12 = cos(auxiliaryEquationsVector(12));
    }

    double cx13;

    if (abs(cos(auxiliaryEquationsVector(13)))<6.2e-17){
        cx13 = 0;
    }
    else {
        cx13 = cos(auxiliaryEquationsVector(13));
    }

    double cx14;

    if (abs(cos(auxiliaryEquationsVector(14)))<6.2e-17){
        cx14 = 0;
//        std::cout<<"cx14 has been set to 0"<<std::endl;
    }
    else {
        cx14 = cos(auxiliaryEquationsVector(14));
    }

    // Same for sin

   double sx10x11;

    if (abs(sin(auxiliaryEquationsVector(49)))<1.2e-16){ // Using x49

        sx10x11 = 0;

    }
    else {
        sx10x11 = sin(auxiliaryEquationsVector(49));
    }



    auxiliaryDerivativesVector(4) = -standardGravitationalParameter*(auxiliaryEquationsVector(1)/auxiliaryEquationsVector(9))+auxiliaryEquationsVector(0)*
            (cx10x11*(-sin(auxiliaryEquationsVector(12))*cx13*cx14+cx12*sin(auxiliaryEquationsVector(14)))-
                 sx10x11*sin(auxiliaryEquationsVector(13))*cx14)+
            thrustAccelerationsBframe(1)*(cx10x11*sin(auxiliaryEquationsVector(12))*sin(auxiliaryEquationsVector(13))-sx10x11*cx13)+
            thrustAccelerationsBframe(2)*(cx10x11*(-sin(auxiliaryEquationsVector(12))*cx13*sin(auxiliaryEquationsVector(14))-
                                                                                                          cx12*cx14)-
                                          sx10x11*sin(auxiliaryEquationsVector(13))*sin(auxiliaryEquationsVector(14)));                // u4




    auxiliaryDerivativesVector(5) = -standardGravitationalParameter*(auxiliaryEquationsVector(2)/auxiliaryEquationsVector(9))+auxiliaryEquationsVector(0)*
            (sx10x11*(-sin(auxiliaryEquationsVector(12))*cx13*cx14+cx12*sin(auxiliaryEquationsVector(14)))+
                 cx10x11*sin(auxiliaryEquationsVector(13))*cx14)+
            thrustAccelerationsBframe(1)*(sx10x11*sin(auxiliaryEquationsVector(12))*sin(auxiliaryEquationsVector(13))+cx10x11*cx13)+
            thrustAccelerationsBframe(2)*(sx10x11*(-sin(auxiliaryEquationsVector(12))*cx13*sin(auxiliaryEquationsVector(14))-
                                                                                                          cx12*cx14)+
                                          cx10x11*sin(auxiliaryEquationsVector(13))*sin(auxiliaryEquationsVector(14)));                // u5

    auxiliaryDerivativesVector(6) = -standardGravitationalParameter*(auxiliaryEquationsVector(3)/auxiliaryEquationsVector(9))+auxiliaryEquationsVector(0)*
            (cx12*cx13*cx14+sin(auxiliaryEquationsVector(12))*sin(auxiliaryEquationsVector(14)))-
            thrustAccelerationsBframe(1)*cx12*sin(auxiliaryEquationsVector(13))+
            thrustAccelerationsBframe(2)*(cx12*cx13*sin(auxiliaryEquationsVector(14))-sin(auxiliaryEquationsVector(12))*cx14);                // u6
//}
    // u6 changed becuase of the mistake found in the complete transformation matrix:
    // (cos(auxiliaryEquationsVector(12))*cos(auxiliaryEquationsVector(13))*cos(auxiliaryEquationsVector(14))-cos(auxiliaryEquationsVector(12))*sin(auxiliaryEquationsVector(14))) [wrong] =>
    // (cos(auxiliaryEquationsVector(12))*cos(auxiliaryEquationsVector(13))*cos(auxiliaryEquationsVector(14))+sin(auxiliaryEquationsVector(12))*sin(auxiliaryEquationsVector(14))) [correct]

//    std::cout<<"The initial acceleration in the x-direction = "<<auxiliaryDerivativesVector(4)<<std::endl;
//    std::cout<<"The initial acceleration in the y-direction = "<<auxiliaryDerivativesVector(5)<<std::endl;
//    std::cout<<"The initial acceleration in the z-direction = "<<auxiliaryDerivativesVector(6)<<std::endl;

    /// Debug ///

/*    std::cout<<"Thrust in I-frame x-direction = "<<auxiliaryDerivativesVector(4)<<std::endl;
    std::cout<<"Thrust in I-frame y-direction = "<<auxiliaryDerivativesVector(5)<<std::endl;
    std::cout<<"Thrust in I-frame z-direction = "<<auxiliaryDerivativesVector(6)<<std::endl;

    std::cout<<"u4 = "<<auxiliaryDerivativesVector(4)<<std::endl;
    std::cout<<"u5 = "<<auxiliaryDerivativesVector(5)<<std::endl;
    std::cout<<"u6 = "<<auxiliaryDerivativesVector(6)<<std::endl;
//    std::cout<<"u4-4.89819672896927 = "<<auxiliaryDerivativesVector(4)-4.89819672896927<<std::endl;
//    std::cout<<"u5-17.6623268076567 = "<<auxiliaryDerivativesVector(5)-17.6623268076567<<std::endl;
//    std::cout<<"u6+23.1286117670555 = "<<auxiliaryDerivativesVector(6)+23.1286117670555<<std::endl;
    std::cout<<"First part of u5 = "<<auxiliaryEquationsVector(0)*
               (sx10x11*(-sin(auxiliaryEquationsVector(12))*cx13*cx14+cx12*sin(auxiliaryEquationsVector(14)))+
                    cx10x11*sin(auxiliaryEquationsVector(13))*cx14)<<std::endl;
        std::cout<<"Second u5 = "<<thrustAccelerationsBframe(1)*(sx10x11*sin(auxiliaryEquationsVector(12))*sin(auxiliaryEquationsVector(13))+cx10x11*cx13)<<std::endl;
    std::cout<<"Third u5 = "<<thrustAccelerationsBframe(2)*(sx10x11*(-sin(auxiliaryEquationsVector(12))*cx13*sin(auxiliaryEquationsVector(14))-
                                                                     cx12*cx14)+
     cx10x11*sin(auxiliaryEquationsVector(13))*sin(auxiliaryEquationsVector(14)))<<std::endl;
    std::cout<<"First part of u6 = "<<auxiliaryEquationsVector(0)*
               (cx12*cx13*cx14+sin(auxiliaryEquationsVector(12))*sin(auxiliaryEquationsVector(14)))<<std::endl;
    std::cout<<"Second u6 = "<<thrustAccelerationsBframe(1)*cx12*sin(auxiliaryEquationsVector(13))<<std::endl;
    std::cout<<"Third u6 = "<<thrustAccelerationsBframe(2)*(cx12*cx13*sin(auxiliaryEquationsVector(14))-sin(auxiliaryEquationsVector(12))*cx14)<<std::endl;
    std::cout<<"w4,0 = "<<auxiliaryEquationsVector(0)<<std::endl;
    std::cout<<"cx12 = "<<cx12<<std::endl;
    std::cout<<"cx13 = "<<cx13<<std::endl;
    std::cout<<"cx14 = "<<cx14<<std::endl;
    std::cout<<"sx12 = "<<sin(auxiliaryEquationsVector(12))<<std::endl;
    std::cout<<"sx14 = "<<sin(auxiliaryEquationsVector(14))<<std::endl;
    std::cout<<"x12 = "<<auxiliaryEquationsVector(12)<<std::endl;
    std::cout<<"x13 = "<<auxiliaryEquationsVector(13)<<std::endl;
    std::cout<<"x14 = "<<auxiliaryEquationsVector(14)<<std::endl;


//*/

    auxiliaryDerivativesVector(50) = auxiliaryEquationsVector(52);      // u50

    auxiliaryDerivativesVector(51) = auxiliaryEquationsVector(53);      // u51

    auxiliaryDerivativesVector(52) = cos(auxiliaryEquationsVector(10))*(auxiliaryDerivativesVector(4)+rotationalVelocity*auxiliaryEquationsVector(5))-
            auxiliaryDerivativesVector(10)*sin(auxiliaryEquationsVector(10))*(auxiliaryEquationsVector(4)+rotationalVelocity*auxiliaryEquationsVector(2))+
            sin(auxiliaryEquationsVector(10))*(auxiliaryDerivativesVector(5)-rotationalVelocity*auxiliaryEquationsVector(4))+
            auxiliaryDerivativesVector(10)*cos(auxiliaryEquationsVector(10))*(auxiliaryEquationsVector(5)-rotationalVelocity*auxiliaryEquationsVector(1));  // u52

    auxiliaryDerivativesVector(53) = -sin(auxiliaryEquationsVector(10))*(auxiliaryDerivativesVector(4)+rotationalVelocity*auxiliaryEquationsVector(5))-
            auxiliaryDerivativesVector(10)*cos(auxiliaryEquationsVector(10))*(auxiliaryEquationsVector(4)+rotationalVelocity*auxiliaryEquationsVector(2))+
            cos(auxiliaryEquationsVector(10))*(auxiliaryDerivativesVector(5)-rotationalVelocity*auxiliaryEquationsVector(4))+
            auxiliaryDerivativesVector(10)*sin(auxiliaryEquationsVector(10))*(auxiliaryEquationsVector(5)-rotationalVelocity*auxiliaryEquationsVector(1));  // u53

    auxiliaryDerivativesVector(54) = 2.0*auxiliaryEquationsVector(4)*auxiliaryDerivativesVector(4)+2*auxiliaryEquationsVector(5)*auxiliaryDerivativesVector(5)+
            rotationalVelocity*rotationalVelocity*auxiliaryDerivativesVector(19)+2.0*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(5)+
                                                                                                         auxiliaryDerivativesVector(4)*auxiliaryEquationsVector(2)-
                                                                                                         auxiliaryEquationsVector(5)*auxiliaryEquationsVector(4)-
                                                                                                         auxiliaryDerivativesVector(5)*auxiliaryEquationsVector(1));       // u54

    auxiliaryDerivativesVector(55) = auxiliaryDerivativesVector(53)*auxiliaryEquationsVector(50)-auxiliaryDerivativesVector(52)*auxiliaryEquationsVector(51);       // u55


    auxiliaryDerivativesVector(7) = -(Thrust/(tudat::physical_constants::STANDARD_EARTH_GRAVITATIONAL_ACCELERATION*specificImpulse));                // u7

    auxiliaryDerivativesVector(8) = 2.0*(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5)+auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6));                // u8

    auxiliaryDerivativesVector(19) = 2.0*(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5));                // u19

    auxiliaryDerivativesVector(9) = 1.5*((auxiliaryEquationsVector(9)*auxiliaryDerivativesVector(8))/auxiliaryEquationsVector(8));                // u9

/*
 * // If the spacecraft flies over a pole i.e. if +-z=r (or x3=x20) then the change in latitude is undefined and has to be set equal to zero here. This is done to avoid singularities.

    if (auxiliaryEquationsVector(3)==auxiliaryEquationsVector(20) || -auxiliaryEquationsVector(3)==auxiliaryEquationsVector(20)){

        auxiliaryDerivativesVector(11) = 0;
    }
    else {
    auxiliaryDerivativesVector(11) = (auxiliaryEquationsVector(1)*auxiliaryEquationsVector(5)-auxiliaryEquationsVector(2)*auxiliaryEquationsVector(4))/auxiliaryEquationsVector(19);                // u11

    };
*/

    auxiliaryDerivativesVector(49) = auxiliaryEquationsVector(45);                      // u49

//    auxiliaryDerivativesVector(11) = auxiliaryEquationsVector(45)-rotationalVelocity;                // u11
    auxiliaryDerivativesVector(11) = auxiliaryEquationsVector(46);                // u11

//    auxiliaryDerivativesVector(20) = auxiliaryEquationsVector(26)/(2.0*auxiliaryEquationsVector(20));                // u20 = \dot(r)
    auxiliaryDerivativesVector(20) = auxiliaryEquationsVector(25);                // u20 = \dot(r)


    // If the spacecraft flies over a pole i.e. if +-z=r (or x3=x20) then the change in latitude is undefined and has to be set equal to zero here. This is done to avoid singularities.

    if (auxiliaryEquationsVector(3)==auxiliaryEquationsVector(20) || -auxiliaryEquationsVector(3)==auxiliaryEquationsVector(20)){

        auxiliaryDerivativesVector(12) = 0;
    }
    else {
//    auxiliaryDerivativesVector(12) = (auxiliaryEquationsVector(20)*auxiliaryEquationsVector(6)-auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25))/
//            (auxiliaryEquationsVector(8)*sqrt(1.0-(auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20))*(auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20))));                // u12
        auxiliaryDerivativesVector(12) = auxiliaryEquationsVector(24);                // u12
    };

    auxiliaryDerivativesVector(31) = auxiliaryDerivativesVector(20);                // u31


    // Avoiding singularities
    if (auxiliaryEquationsVector(19) == 0){
//std::cout<<"The equation goes here 1"<<std::endl;
        auxiliaryDerivativesVector(45) = 0;
    }
//    else if(abs(auxiliaryEquationsVector(1)*auxiliaryDerivativesVector(5)-auxiliaryEquationsVector(2)*auxiliaryDerivativesVector(4))<=1E-8){                        // Set tolerance for the differences to 1E-8

////        std::cout<<"The equation goes here 2"<<std::endl;
////        std::cout<<"(x1*u5-x2*u4) = "<<(auxiliaryEquationsVector(1)*auxiliaryDerivativesVector(5)-auxiliaryEquationsVector(2)*auxiliaryDerivativesVector(4))<<std::endl;

//        auxiliaryDerivativesVector(45) = 0;
////        std::cout<<"Yep"<<std::endl;
//    }
    else if (verticalInertialFlightPathAngleSet == true){ // For vertical flight there is not acceleration for longitude
        auxiliaryDerivativesVector(45) = 0.0;
    }
    else {
//std::cout<<"The equation goes here 3"<<std::endl;
    auxiliaryDerivativesVector(45) = (auxiliaryEquationsVector(19)*(auxiliaryEquationsVector(1)*auxiliaryDerivativesVector(5)-auxiliaryEquationsVector(2)*auxiliaryDerivativesVector(4))
                                      -auxiliaryDerivativesVector(19)*(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(5)-auxiliaryEquationsVector(2)*auxiliaryEquationsVector(4)))/
            (auxiliaryEquationsVector(19)*auxiliaryEquationsVector(19));                // u45
    };

    /// Debug ///
/*
    std::cout<<"x1*u5 = "<<auxiliaryEquationsVector(1)*auxiliaryDerivativesVector(5)<<std::endl;
    std::cout<<"x2*u4 = "<<auxiliaryEquationsVector(2)*auxiliaryDerivativesVector(4)<<std::endl;
    std::cout<<"x1*u5-x2*u4 = "<<auxiliaryEquationsVector(1)*auxiliaryDerivativesVector(5)-auxiliaryEquationsVector(2)*auxiliaryDerivativesVector(4)<<std::endl;
    std::cout<<"x1*u5-20.4150429110416 = "<<auxiliaryEquationsVector(1)*auxiliaryDerivativesVector(5)-20.4150429110416<<std::endl;
    std::cout<<"x2*u4-20.4150429110416 = "<<auxiliaryEquationsVector(2)*auxiliaryDerivativesVector(4)-20.4150429110416<<std::endl;
    std::cout<<"x1 = "<<auxiliaryEquationsVector(1)<<std::endl;
    std::cout<<"u5 = "<<auxiliaryDerivativesVector(5)<<std::endl;
    std::cout<<"x2 = "<<auxiliaryEquationsVector(2)<<std::endl;
    std::cout<<"u4 = "<<auxiliaryDerivativesVector(4)<<std::endl;
    std::cout<<"x1-847113.311014168 = "<<auxiliaryEquationsVector(1)-847113.311014168<<std::endl;
    std::cout<<"u5-17.6623268076567 = "<<auxiliaryDerivativesVector(5)-17.6623268076567<<std::endl;
    std::cout<<"x2-3054591.91823781 = "<<auxiliaryEquationsVector(2)-3054591.91823781<<std::endl;
    std::cout<<"u4-4.89819672896927 = "<<auxiliaryDerivativesVector(4)-4.89819672896927<<std::endl;
    std::cout<<"x2/x9 = "<<auxiliaryEquationsVector(2)/auxiliaryEquationsVector(9)<<std::endl;
    std::cout<<"auxiliaryDerivativesVector(45) = (auxiliaryEquationsVector(19)*(auxiliaryEquationsVector(1)*auxiliaryDerivativesVector(5)-auxiliaryEquationsVector(2)*auxiliaryDerivativesVector(4))-auxiliaryDerivativesVector(19)*(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(5)-auxiliaryEquationsVector(2)*auxiliaryEquationsVector(4)))/(auxiliaryEquationsVector(19)*auxiliaryEquationsVector(19)) = "<<
               (auxiliaryEquationsVector(19)*(auxiliaryEquationsVector(1)*auxiliaryDerivativesVector(5)-auxiliaryEquationsVector(2)*auxiliaryDerivativesVector(4))
                                                 -auxiliaryDerivativesVector(19)*(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(5)-auxiliaryEquationsVector(2)*auxiliaryEquationsVector(4)))/
                       (auxiliaryEquationsVector(19)*auxiliaryEquationsVector(19))<<std::endl;
//*/

    auxiliaryDerivativesVector(21) = 2.0*(auxiliaryEquationsVector(4)*auxiliaryDerivativesVector(4)+auxiliaryEquationsVector(5)*auxiliaryDerivativesVector(5)+auxiliaryEquationsVector(6)*auxiliaryDerivativesVector(6));                // u21

//    // Set tolerance for the derivative of the inertial velocity squared to 1E-10
//    if (abs(auxiliaryDerivativesVector(21))<=1E-10){

//        auxiliaryDerivativesVector(21) = 0;
//    }


/*
    /// Debug ///

    std::cout<<"x4*u4  = "<<auxiliaryEquationsVector(4)*auxiliaryDerivativesVector(4)<<std::endl;
    std::cout<<"x5*u5  = "<<auxiliaryEquationsVector(5)*auxiliaryDerivativesVector(5)<<std::endl;
    std::cout<<"x6*u6  = "<<auxiliaryEquationsVector(6)*auxiliaryDerivativesVector(6)<<std::endl;

    std::cout<<"x4*u4+x5*u5  = "<<auxiliaryEquationsVector(4)*auxiliaryDerivativesVector(4)+auxiliaryEquationsVector(5)*auxiliaryDerivativesVector(5)<<std::endl;
    std::cout<<"(x4*u4)/(x5*u5) = "<<(auxiliaryEquationsVector(4)*auxiliaryDerivativesVector(4))/(auxiliaryEquationsVector(5)*auxiliaryDerivativesVector(5))<<std::endl;

    std::cout<<"x4u4+1060.50600304257 = "<<auxiliaryEquationsVector(4)*auxiliaryDerivativesVector(4)+1060.50600304257<<std::endl;
    std::cout<<"x5u5-1060.50600304257 = "<<auxiliaryEquationsVector(5)*auxiliaryDerivativesVector(5)-1060.50600304257<<std::endl;
*/


    auxiliaryDerivativesVector(26) = 2.0*auxiliaryEquationsVector(21)+2.0*(auxiliaryEquationsVector(1)*auxiliaryDerivativesVector(4)+auxiliaryEquationsVector(2)*auxiliaryDerivativesVector(5)+auxiliaryEquationsVector(3)*auxiliaryDerivativesVector(6));                // u26


    // Computing the polynomial fit derivative using the altitude and fit parameters for density
    for (int i = 1; i < 10+1;i++) {

    auxiliaryDerivativesVector(30) += auxiliaryDerivativesVector(31)*i*pow(auxiliaryEquationsVector(31),i-1)*densityPolyCoefficients(i);              // u30
};


    // Determine which section of the temperature curve needs to be used and what the corresponding order is
    // Also, because a computer is less than perfect, a small correction is made to the lower bound of the first section to make sure that the initial altitude is still valid

    if (temperatureAltitudeRanges(0,0)-0.000000000001 <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(0,1)){

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
    for (int i=1; i < powerT+1;i++){

    auxiliaryDerivativesVector(34) += auxiliaryDerivativesVector(31)*i*pow(auxiliaryEquationsVector(31),i-1)*temperaturePolyCoefficients(sectionT,i);              // u34

};


    auxiliaryDerivativesVector(28) = auxiliaryDerivativesVector(30)*exp(auxiliaryEquationsVector(30)) ;                // u28

    auxiliaryDerivativesVector(33) = ((adiabeticIndex*specificGasConstant)/(2.0*auxiliaryEquationsVector(33)))*auxiliaryDerivativesVector(34);                // u33

    auxiliaryDerivativesVector(36) = sqrt(auxiliaryDerivativesVector(4)*auxiliaryDerivativesVector(4)+auxiliaryDerivativesVector(5)*auxiliaryDerivativesVector(5)+
                                          auxiliaryDerivativesVector(6)*auxiliaryDerivativesVector(6));          // u36

    if (verticalInertialFlightPathAngleSet == true){ // Vertical flight, then the radius changes with the same pase as the inertial velocity
        auxiliaryDerivativesVector(25) = auxiliaryDerivativesVector(36);
    }
    else {
    auxiliaryDerivativesVector(25) = (2.0*auxiliaryEquationsVector(8)*auxiliaryDerivativesVector(26)-auxiliaryEquationsVector(26)*auxiliaryEquationsVector(26))/(4.0*auxiliaryEquationsVector(9));                // u25
    }

    /// Debug ///
/*
    std::cout<<"x28*u26 = "<<auxiliaryEquationsVector(8)*auxiliaryDerivativesVector(26)<<std::endl;
    std::cout<<"x26^2 = "<<auxiliaryEquationsVector(26)*auxiliaryEquationsVector(26)<<std::endl;
    std::cout<<"x9 = "<<auxiliaryEquationsVector(9)<<std::endl;
*/

//    // Avoid singularities
//    if (auxiliaryEquationsVector(36) ==0){
//        auxiliaryDerivativesVector(36) = 0;
//    }
//    else{
//    auxiliaryDerivativesVector(36) = auxiliaryDerivativesVector(21)/(2.0*auxiliaryEquationsVector(36));                // u36
//}







    auxiliaryDerivativesVector(46) = auxiliaryDerivativesVector(45);                // u46

    // Avoid cosine rounding errors
    if (auxiliaryDerivativesVector(45) == 0.0 || verticalInertialFlightPathAngleSet == true){ // In vertical flight, u47 is zero
        auxiliaryDerivativesVector(47) = 0.0;
    }
    else if (abs(cos(auxiliaryEquationsVector(12)))<6.2e-17){
               auxiliaryDerivativesVector(47) = -auxiliaryEquationsVector(45)*auxiliaryDerivativesVector(12)*sin(auxiliaryEquationsVector(12));
            }

            else {
    auxiliaryDerivativesVector(47) = auxiliaryDerivativesVector(45)*cos(auxiliaryEquationsVector(12))-auxiliaryEquationsVector(45)*auxiliaryDerivativesVector(12)*sin(auxiliaryEquationsVector(12));                // u47
}

    // If the spacecraft flies over a pole i.e. if +-z=r (or x3=x20) then the change in latitude is undefined and has to be set equal to zero here. This is done to avoid singularities.

    if (auxiliaryEquationsVector(3)==auxiliaryEquationsVector(20) || -auxiliaryEquationsVector(3)==auxiliaryEquationsVector(20)){

        auxiliaryDerivativesVector(24) = 0;
    }
    else if (verticalInertialFlightPathAngleSet == true){ // Vertical flight has no acceleration in the latitude
        auxiliaryDerivativesVector(24) = 0;
    }
    else {
//    auxiliaryDerivativesVector(24) = (auxiliaryEquationsVector(6)*auxiliaryDerivativesVector(20)+auxiliaryEquationsVector(20)*auxiliaryDerivativesVector(6)-auxiliaryEquationsVector(6)*auxiliaryEquationsVector(25)-auxiliaryEquationsVector(3)*auxiliaryDerivativesVector(25))/
//            (auxiliaryEquationsVector(8)*sqrt(1-pow((auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)),2.0)))-(((2.0*auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3)*auxiliaryDerivativesVector(25))/(pow(auxiliaryEquationsVector(20),3.0))-
//                                                                                                                      (2.0*auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6))/auxiliaryEquationsVector(8))*(auxiliaryEquationsVector(20)*auxiliaryEquationsVector(6)-
//                                                                                                                                                                                                                auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25)))/
//            (2.0*auxiliaryEquationsVector(8)*pow((1.0-pow((auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)),2.0)),1.5))-
//            (auxiliaryDerivativesVector(8)*(auxiliaryEquationsVector(20)*auxiliaryEquationsVector(6)-auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25)))/
//            (auxiliaryEquationsVector(8)*auxiliaryEquationsVector(8)*sqrt(1.0-pow((auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)),2.0)));                // u24

    auxiliaryDerivativesVector(24) = (auxiliaryEquationsVector(20)*auxiliaryDerivativesVector(6)-auxiliaryEquationsVector(3)*auxiliaryDerivativesVector(25))/
            (auxiliaryEquationsVector(8)*sqrt(1-pow((auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)),2.0)))-(((2.0*auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25))/(pow(auxiliaryEquationsVector(20),3.0))-
                                                                                                                      (2.0*auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6))/auxiliaryEquationsVector(8))*(auxiliaryEquationsVector(20)*auxiliaryEquationsVector(6)-
                                                                                                                                                                                                                auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25)))/
            (2.0*auxiliaryEquationsVector(8)*pow((1.0-pow((auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)),2.0)),1.5))-
            (auxiliaryDerivativesVector(8)*(auxiliaryEquationsVector(20)*auxiliaryEquationsVector(6)-auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25)))/
            (auxiliaryEquationsVector(8)*auxiliaryEquationsVector(8)*sqrt(1.0-pow((auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)),2.0)));                // u24
    };



    /// Debug ///

/*    std::cout<<"u24 old = "<<(auxiliaryEquationsVector(6)*auxiliaryDerivativesVector(20)+auxiliaryEquationsVector(20)*auxiliaryDerivativesVector(6)-auxiliaryEquationsVector(6)*auxiliaryEquationsVector(25)-auxiliaryEquationsVector(3)*auxiliaryDerivativesVector(25))/
            (auxiliaryEquationsVector(8)*sqrt(1-pow((auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)),2)))-(((2*auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3)*auxiliaryDerivativesVector(25))/(pow(auxiliaryEquationsVector(20),3))-
                                                                                                                      (2*auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6))/auxiliaryEquationsVector(8))*(auxiliaryEquationsVector(20)*auxiliaryEquationsVector(6)-
                                                                                                                                                                                                                auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25)))/
            (2*auxiliaryEquationsVector(8)*pow((1-pow((auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)),2)),1.5))-
            (auxiliaryDerivativesVector(8)*(auxiliaryEquationsVector(20)*auxiliaryEquationsVector(6)-auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25)))/
            (auxiliaryEquationsVector(8)*auxiliaryEquationsVector(8)*sqrt(1-pow((auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)),2)))<<std::endl;
    std::cout<<"u24 new = "<<(((auxiliaryDerivativesVector(6)/auxiliaryEquationsVector(20))-((2*auxiliaryEquationsVector(6)*auxiliaryEquationsVector(25))/auxiliaryEquationsVector(8))-
                              ((auxiliaryEquationsVector(3)*auxiliaryDerivativesVector(25))/auxiliaryEquationsVector(8))+((2*auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25)*auxiliaryEquationsVector(25))/
                                                                                                                          (auxiliaryEquationsVector(8)*auxiliaryEquationsVector(20))))/
            sqrt(1-((auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3))/auxiliaryEquationsVector(8))))-
            ((((2*auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25))/(auxiliaryEquationsVector(8)*auxiliaryEquationsVector(20)))-
              ((2*auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6))/auxiliaryEquationsVector(8)))*((auxiliaryEquationsVector(6)/auxiliaryEquationsVector(20))-((auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25))/auxiliaryEquationsVector(8)))/
             (2*pow((1-((auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3))/auxiliaryEquationsVector(8))),(3/2))))<<std::endl;

    double x = 2.0;
    double fx = x*x;
    double fxp = 2.0*x;
    double fxdp = 2.0;
    double gx = x*x*x;
    double gxp = 3.0*x*x;
    double gxdp = 6.0*x;
    double gxs = x*x*x*x*x*x;
    double gxsp = 6.0*x*x*x*x*x;
    double hx = fx/gx;
    double hxp = (gx*fxp-fx*gxp)/gxs;
    double hxdp = (gxs*fxdp-gx*(2.0*fxp*gxp+fx*gxdp)+2.0*fx*gxp*gxp)/(gxs*gx);

    std::cout<<"x = "<<x<<std::endl;
    std::cout<<"fx = "<<fx<<std::endl;
    std::cout<<"fxp = "<<fxp<<std::endl;
    std::cout<<"fxdp = "<<fxdp<<std::endl;
    std::cout<<"gx = "<<gx<<std::endl;
    std::cout<<"gxp = "<<gxp<<std::endl;
    std::cout<<"gxdp = "<<gxdp<<std::endl;
    std::cout<<"gxs = "<<gxs<<std::endl;
    std::cout<<"hx = "<<hx<<std::endl;
    std::cout<<"hxp = "<<hxp<<std::endl;
    std::cout<<"hxdp = "<<hxdp<<std::endl;
    std::cout<<"asin(h(x))' = "<<hxp/sqrt(1-hx*hx)<<std::endl;
    std::cout<<"asin(f(x)/g(x))' = "<<(gx*fxp-fx*gxp)/(gxs*sqrt(1-((fx*fx)/gxs)))<<std::endl;

    std::cout<<"u24 example old = "<<(fxp*gx+gx*fxdp-fxp*gxp-fx*gxdp)/
            (gxs*sqrt(1-pow((fx/gx),2)))-(((2*fx*fx*gxdp)/(pow(gx,3))-
                                                                                                                      (2*fx*fxp)/gxs)*(gx*fxp-
                                                                                                                                                                                                                fx*gxp))/
            (2*gxs*pow((1-pow((fx/gx),2)),1.5))-
            (gxs*(gx*fxp-fx*gxp))/
            (gxs*gxs*sqrt(1-pow((fx/gx),2)))<<std::endl;

    std::cout<<"u24 example old adjusted = "<<(gx*fxdp+fxp*gxp-fxp*gxp-fx*gxdp)/(gxs*sqrt(1-pow((fx/gx),2)))-(((2*fx*fx*gxp)/(pow(gx,3))-(2*fx*fxp)/gxs)*(gx*fxp-fx*gxp))/
            (2*gxs*pow((1-pow((fx/gx),2)),1.5))-
            (gxsp*(gx*fxp-fx*gxp))/
            (gxs*gxs*sqrt(1-pow((fx/gx),2)))<<std::endl;

    std::cout<<"u24 example original newly written = "<<((gx*fxdp+fxp*gxp-fx*gxdp-gxp*fxp)/(gxs*sqrt(1.0-((fx*fx)/gxs))))-((gxsp*(fxp*gx-gxp*fx))/(gxs*gxs*sqrt(1.0-((fx*fx)/gxs))))-
            (((((2.0*fx*fx*gxp)/(gx*gx*gx))-((2.0*fx*fxp)/gxs))*(fxp*gx-gxp*fx))/(2.0*gxs*pow((1.0-((fx*fx)/gxs)),1.5)))<<std::endl;

    std::cout<<"u24 example original newly written opt = "<<((gx*fxdp-fx*gxdp)/(gxs*sqrt(1.0-((fx*fx)/gxs))))-((gxsp*(fxp*gx-gxp*fx))/(gxs*gxs*sqrt(1.0-((fx*fx)/gxs))))-
            (((((2.0*fx*fx*gxp)/(gx*gx*gx))-((2.0*fx*fxp)/gxs))*(fxp*gx-gxp*fx))/(2.0*gxs*pow((1.0-((fx*fx)/gxs)),1.5)))<<std::endl;

    std::cout<<"u24 example old newly written = "<<((gx*fxdp-fx*gxdp)/(gx*gx*sqrt(1.0-((fx*fx)/(gx*gx)))))-
            ((2.0*gxp*(gx*fxp-fx*gxp))/(gx*gx*gx*sqrt(1.0-((fx*fx)/(gx*gx)))))-
            (((((2.0*fx*fx*gxp)/(gx*gx*gx))-((2.0*fx*fxp)/(gx*gx)))*(gx*fxp-fx*gxp))/(2.0*gx*gx*pow((1.0-((fx*fx)/(gx*gx))),1.5)))<<std::endl;

    std::cout<<"u24 example new = "<<(((fxdp/gx)-((2.0*fxp*gxp)/gxs)-
                              ((fx*gxdp)/gxs)+((2.0*fx*gxp*gxp)/
                                                                                                                          (gxs*gx)))/
            sqrt(1.0-((fx*fx)/gxs)))-
            ((((2.0*fx*fx*gxp)/(gxs*gx))-
              ((2.0*fx*fxp)/gxs))*((fxp/gx)-((fx*gxp)/gxs))/
             (2.0*pow((1.0-((fx*fx)/gxs)),1.5)))<<std::endl;

    std::cout<<"u24 example simple = "<<(2.0*fx-1.0)/(sqrt(1.0-(1.0/fx))*gx*(fx-1.0))<<std::endl;

    std::cout<<"u24 example h(x) pow = "<<(hx*hxp*hxp-(hx*hx-1)*hxdp)/pow((1-hx*hx),1.5)<<std::endl;
    std::cout<<"u24 example h(x) = "<<(hx*hxp*hxp-(hx*hx-1.0)*hxdp)/((1.0-hx*hx)*sqrt(1.0-hx*hx))<<std::endl;
    std::cout<<"pow((1-hx*hx),(3/2)) = "<<pow((1-hx*hx),(3/2))<<std::endl;
    std::cout<<"((1-hx*hx)*sqrt(1-hx*hx)) = "<<((1-hx*hx)*sqrt(1-hx*hx))<<std::endl;
    std::cout<<"pow(0.5,2) = "<<pow(0.5,2)<<std::endl;
    std::cout<<"pow(0.5,(3/2)) = "<<pow(0.5,(3/2))<<std::endl;
    double pfunc = 3/2;
    double three = 3;
    double two = 2;
    double threeOverTwo = three/two;
    std::cout<<"pfunc = 3/2 = "<<pfunc<<std::endl;
    std::cout<<"threeOverTwo = "<<threeOverTwo<<std::endl;
    std::cout<<"pow(0.5,threeOverTwo) = "<<pow(0.5,threeOverTwo)<<std::endl;
    std::cout<<"pow(0.5,(pfunc)) = "<<pow(0.5,(pfunc))<<std::endl;
    std::cout<<"pow(0.5,1.5) = "<<pow(0.5,1.5)<<std::endl;
    std::cout<<"0.5*sqrt(0.5) = "<<0.5*sqrt(0.5)<<std::endl;

    std::cout<<"Part 1 of x24 = "<<(auxiliaryEquationsVector(6)*auxiliaryDerivativesVector(20)+auxiliaryEquationsVector(20)*auxiliaryDerivativesVector(6)-auxiliaryEquationsVector(6)*auxiliaryEquationsVector(25)-auxiliaryEquationsVector(3)*auxiliaryDerivativesVector(25))/
               (auxiliaryEquationsVector(8)*sqrt(1-pow((auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)),2)))<<std::endl;

    std::cout<<"Part 1.1 of x24 = "<<(auxiliaryEquationsVector(6)*auxiliaryDerivativesVector(20)+auxiliaryEquationsVector(20)*auxiliaryDerivativesVector(6)-auxiliaryEquationsVector(6)*auxiliaryEquationsVector(25)-auxiliaryEquationsVector(3)*auxiliaryDerivativesVector(25))<<std::endl;
    std::cout<<"Part 1.2 of x24 = "<<(auxiliaryEquationsVector(8)*sqrt(1-pow((auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)),2)))<<std::endl;

    std::cout<<"Part 2 of x24 = "<<(((2*auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3)*auxiliaryDerivativesVector(25))/(pow(auxiliaryEquationsVector(20),3))-
                                     (2*auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6))/auxiliaryEquationsVector(8))*(auxiliaryEquationsVector(20)*auxiliaryEquationsVector(6)-
                                                                                                                               auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25)))/
(2*auxiliaryEquationsVector(8)*pow((1-pow((auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)),2)),1.5))<<std::endl;

    std::cout<<"Part 3 of x24 = "<<(auxiliaryDerivativesVector(8)*(auxiliaryEquationsVector(20)*auxiliaryEquationsVector(6)-auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25)))/
               (auxiliaryEquationsVector(8)*sqrt(1-pow((auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20)),2)))<<std::endl;
*/





    // Avoid singularities and account for vertical flight (there is not going to be a difference between r_dot and V_I in vertical flight)
    if (auxiliaryEquationsVector(36)==0 || verticalInertialFlightPathAngleSet == true){
        auxiliaryDerivativesVector(37) = 0;
    }
    else {
    auxiliaryDerivativesVector(37) = (auxiliaryEquationsVector(36)*auxiliaryDerivativesVector(25)-auxiliaryEquationsVector(25)*auxiliaryDerivativesVector(36))/(auxiliaryEquationsVector(36)*auxiliaryEquationsVector(36));                // u37
}
    if (rotationalVelocity == 0.0){ // Non-rotating Mars
        auxiliaryDerivativesVector(41) = 0.0;
    }
    else {
    auxiliaryDerivativesVector(41) = auxiliaryEquationsVector(36)*auxiliaryDerivativesVector(35)+auxiliaryEquationsVector(35)*auxiliaryDerivativesVector(36);                // u41
    }

    // Avoid cosine rounding errors
                if (abs(cos(auxiliaryEquationsVector(12)))<6.2e-17){
                    auxiliaryDerivativesVector(48) = -auxiliaryEquationsVector(46)*auxiliaryDerivativesVector(12)*sin(auxiliaryEquationsVector(12));
                }
                else {
    auxiliaryDerivativesVector(48) = auxiliaryDerivativesVector(46)*cos(auxiliaryEquationsVector(12))-auxiliaryEquationsVector(46)*auxiliaryDerivativesVector(12)*sin(auxiliaryEquationsVector(12));                // u47
}

//    auxiliaryDerivativesVector(18) = auxiliaryEquationsVector(20)*auxiliaryDerivativesVector(24)+auxiliaryEquationsVector(24)*auxiliaryEquationsVector(25);                // u18

    // Avoiding singularities and accounting for vertical flight
    if ((auxiliaryEquationsVector(47)*auxiliaryEquationsVector(47)+auxiliaryEquationsVector(24)*auxiliaryEquationsVector(24))==0 || verticalInertialFlightPathAngleSet == true){

        auxiliaryDerivativesVector(40) = 0;


    }
    else {
    auxiliaryDerivativesVector(40) = (auxiliaryEquationsVector(24)*auxiliaryDerivativesVector(47)-auxiliaryEquationsVector(47)*auxiliaryDerivativesVector(24))/
            (auxiliaryEquationsVector(47)*auxiliaryEquationsVector(47)+auxiliaryEquationsVector(24)*auxiliaryEquationsVector(24));                // u40

    };

    /// Debug ///
/*
    std::cout<<"x24 = "<<auxiliaryEquationsVector(24)<<std::endl;
    std::cout<<"u47 = "<<auxiliaryDerivativesVector(47)<<std::endl;
    std::cout<<"x47 = "<<auxiliaryEquationsVector(47)<<std::endl;
    std::cout<<"u24 = "<<auxiliaryDerivativesVector(24)<<std::endl;
*/








    // If x37 = +-1 then the derivative is undefined (should not happen)
    if (auxiliaryEquationsVector(37)==1){

        auxiliaryDerivativesVector(38) = 0;
        verticalInertialFlightPathAngleSet = true; // Declare vertical flight in inertial frame
        std::cout<<"verticalInertialFlightPathAngleSet der 2 = "<<verticalInertialFlightPathAngleSet<<std::endl;
    }
    else if (auxiliaryEquationsVector(37)==-1){
        auxiliaryDerivativesVector(38) = 0;
    }
    else {
    auxiliaryDerivativesVector(38) = auxiliaryDerivativesVector(37)/sqrt(1.0-auxiliaryEquationsVector(37)*auxiliaryEquationsVector(37));                // u38
};
/*
    auxiliaryDerivativesVector(44) = auxiliaryDerivativesVector(36)*cos(auxiliaryEquationsVector(38))-auxiliaryEquationsVector(36)*auxiliaryDerivativesVector(38)*sin(auxiliaryEquationsVector(38));                // u44

    auxiliaryDerivativesVector(39) = (auxiliaryEquationsVector(44)*auxiliaryDerivativesVector(18)-auxiliaryEquationsVector(18)*auxiliaryDerivativesVector(44))/(auxiliaryEquationsVector(44)*auxiliaryEquationsVector(44));                // u39

    // If x39 = +-1 then the derivative is undefined
    if (auxiliaryEquationsVector(39)==1 || auxiliaryEquationsVector(39)==-1){

        auxiliaryDerivativesVector(40) = 0;
    }
    else {
    auxiliaryDerivativesVector(40) = -auxiliaryDerivativesVector(39)/sqrt(1-auxiliaryEquationsVector(39)*auxiliaryEquationsVector(39));                // u40
};
*/
    // Account for vertical flight
    if (verticalInertialFlightPathAngleSet == true){
        auxiliaryDerivativesVector(42) = 0.0;
    }
    // Avoid cosine rounding errors
    else if (abs(cos(auxiliaryEquationsVector(38)))<6.2e-17 || abs(cos(auxiliaryEquationsVector(40)))<6.2e-17){
                     auxiliaryDerivativesVector(42) = -auxiliaryDerivativesVector(38)*sin(auxiliaryEquationsVector(38))*sin(auxiliaryEquationsVector(40));
                }
                else {
    auxiliaryDerivativesVector(42) = cos(auxiliaryEquationsVector(38))*cos(auxiliaryEquationsVector(40))*auxiliaryDerivativesVector(40)-auxiliaryDerivativesVector(38)*sin(auxiliaryEquationsVector(38))*sin(auxiliaryEquationsVector(40));                // u42
}
//    /// Debug ///
//    std::cout<<"u42 = "<<auxiliaryDerivativesVector(42)<<std::endl;
//    std::cout<<"cos(x38)*cos(x40)*u40-u38*sin(x38)*sin(x40) = "<<cos(auxiliaryEquationsVector(38))*cos(auxiliaryEquationsVector(40))*auxiliaryDerivativesVector(40)-auxiliaryDerivativesVector(38)*sin(auxiliaryEquationsVector(38))*sin(auxiliaryEquationsVector(40))<<std::endl;
//    std::cout<<"cos(x38)*cos(x40)*u40 = "<<cos(auxiliaryEquationsVector(38))*cos(auxiliaryEquationsVector(40))*auxiliaryDerivativesVector(40)<<std::endl;
//    std::cout<<"u38*sin(x38)*sin(x40) = "<<auxiliaryDerivativesVector(38)*sin(auxiliaryEquationsVector(38))*sin(auxiliaryEquationsVector(40))<<std::endl;
//    std::cout<<"cos(x38) = "<<cos(auxiliaryEquationsVector(38))<<std::endl;
//    std::cout<<"cos(x40) = "<<cos(auxiliaryEquationsVector(40))<<std::endl;
//    std::cout<<"u40 = "<<auxiliaryDerivativesVector(40)<<std::endl;
//    std::cout<<"x38 = "<<auxiliaryEquationsVector(38)<<std::endl;
//    std::cout<<"x40 = "<<auxiliaryEquationsVector(40)<<std::endl;
//    std::cout<<"x40 - 1.5707963267949 = "<<auxiliaryEquationsVector(40)-1.5707963267949<<std::endl;
//    /// Debug ///

    if (verticalInertialFlightPathAngleSet == true || rotationalVelocity == 0.0){ // For vertical flight x43 does not change (it also does not change when dealing with a non-rotating Mars)
        auxiliaryDerivativesVector(43) = 0.0;
    }
    else {
    auxiliaryDerivativesVector(43) = auxiliaryEquationsVector(41)*auxiliaryDerivativesVector(42)+auxiliaryEquationsVector(42)*auxiliaryDerivativesVector(41);                // u43
    }


//    // If V_R = 0 m/s then the derivative is undefined
//    if (auxiliaryEquationsVector(15)==0){

//        auxiliaryDerivativesVector(15) = 0;
//    }
//    else if (auxiliaryEquationsVector(23) == 1.0){ // If vertical flight then r_dot and V_R dot change at the same rate (are equal)
//        auxiliaryDerivativesVector(15) = auxiliaryDerivativesVector(25);
//    }
//    else {
//    auxiliaryDerivativesVector(15) = (2.0*auxiliaryEquationsVector(35)*auxiliaryDerivativesVector(35)+auxiliaryDerivativesVector(21)-2.0*auxiliaryDerivativesVector(43))/(2.0*auxiliaryEquationsVector(15));                // u15
//};


    auxiliaryDerivativesVector(15) = rotationalVelocity*sqrt(rotationalVelocity*rotationalVelocity*auxiliaryEquationsVector(19)+4.0*auxiliaryEquationsVector(54)+
                                                             4.0*rotationalVelocity*auxiliaryEquationsVector(55))-(standardGravitationalParameter/auxiliaryEquationsVector(8))+
            (Thrust/auxiliaryEquationsVector(7))+(auxiliaryEquationsVector(27)/auxiliaryEquationsVector(7));   // u15






    // If V_R = 0 m/s then the derivative is undefined
    if (auxiliaryEquationsVector(15)==0){

        auxiliaryDerivativesVector(23) = 0;
    }
    else {
    auxiliaryDerivativesVector(23) = (auxiliaryEquationsVector(15)*auxiliaryDerivativesVector(25)-auxiliaryEquationsVector(25)*auxiliaryDerivativesVector(15))/(auxiliaryEquationsVector(15)*auxiliaryEquationsVector(15));                // u23
};

    auxiliaryDerivativesVector(32) = (auxiliaryEquationsVector(33)*auxiliaryDerivativesVector(15)-auxiliaryEquationsVector(15)*auxiliaryDerivativesVector(33))/(auxiliaryEquationsVector(33)*auxiliaryEquationsVector(33));                // u32


    // If x23 = +-1 then the derivative is undefined
    if (auxiliaryEquationsVector(23)==1){

        auxiliaryDerivativesVector(14) = 0;
        verticalRotationalFlightPathAngleSet = true;    // Declare vertical flight in rotating frame
        std::cout<<"verticalRotationalFlightPathAngleSet der 2 = "<<verticalRotationalFlightPathAngleSet<<std::endl;
    }
    else if (auxiliaryEquationsVector(23)==-1){
        auxiliaryDerivativesVector(14) = 0;
    }
    else {
    auxiliaryDerivativesVector(14) = auxiliaryDerivativesVector(23)/sqrt(1.0-auxiliaryEquationsVector(23)*auxiliaryEquationsVector(23));                // u14
};

    // Avoiding singularities

    if ((auxiliaryEquationsVector(48)*auxiliaryEquationsVector(48)+auxiliaryEquationsVector(24)*auxiliaryEquationsVector(24))==0){

        auxiliaryDerivativesVector(13) = 0;
    }
    else if (verticalRotationalFlightPathAngleSet == true){ // For vertical flight
        auxiliaryDerivativesVector(13) = 0.0;
    }
    else{
    auxiliaryDerivativesVector(13) = (auxiliaryEquationsVector(24)*auxiliaryDerivativesVector(48)-auxiliaryEquationsVector(48)*auxiliaryDerivativesVector(24))/
            (auxiliaryEquationsVector(48)*auxiliaryEquationsVector(48)+auxiliaryEquationsVector(24)*auxiliaryEquationsVector(24));                // u13

    };

    // Determine which section of the drag coefficient curve needs to be used

    for (int i=0; i < 5+1; i++){

        if (dragCoefficientMachRanges(i,0) <= auxiliaryEquationsVector(32) && auxiliaryEquationsVector(32) < dragCoefficientMachRanges(i,1)){

            sectionCD = i;
        }

    };

    auxiliaryDerivativesVector(29) = dragCoefficientPolyCoefficients(sectionCD,1)*auxiliaryDerivativesVector(32);              // u29



//    auxiliaryDerivativesVector(16) = -sin(auxiliaryEquationsVector(14))*auxiliaryDerivativesVector(14);                // u16

    auxiliaryDerivativesVector(27) = 0.5*referenceArea*auxiliaryEquationsVector(15)*(auxiliaryEquationsVector(15)*(auxiliaryEquationsVector(29)*auxiliaryDerivativesVector(28)+auxiliaryEquationsVector(28)*auxiliaryDerivativesVector(29))+
                                                                                     2.0*auxiliaryEquationsVector(28)*auxiliaryEquationsVector(29)*auxiliaryDerivativesVector(15));                // u27
/*
    auxiliaryDerivativesVector(17) = auxiliaryEquationsVector(16)*auxiliaryDerivativesVector(15)-auxiliaryEquationsVector(15)*auxiliaryDerivativesVector(16);                // u17

    // If x17 = 0 then the derivative is undefined
    if (auxiliaryEquationsVector(17)==0){

        auxiliaryDerivativesVector(22) = 0;
    }
    else {
    auxiliaryDerivativesVector(22) = (auxiliaryEquationsVector(17)*auxiliaryDerivativesVector(18)-auxiliaryEquationsVector(18)*auxiliaryDerivativesVector(17))/(auxiliaryEquationsVector(17)*auxiliaryEquationsVector(17));                // u22
};

    // If x22 = +-1 then the derivative is undefined
    if (auxiliaryEquationsVector(22)==1 || auxiliaryEquationsVector(22)==-1){

        auxiliaryDerivativesVector(13) = 0;
    }
    else {
    auxiliaryDerivativesVector(13) = -auxiliaryDerivativesVector(22)/sqrt(1-auxiliaryEquationsVector(22)*auxiliaryEquationsVector(22));                // u13
};
*/

// auxiliaryDerivativesVector() = ;                // u

    // Set vertical ascent to false again
    //        verticalInertialFlightPathAngleSet = false;
            verticalInertialFlightPathAngleSet = NULL;
            std::cout<<"verticalInertialFlightPathAngleSet der 3 = "<<verticalInertialFlightPathAngleSet<<std::endl;
    //        verticalRotationalFlightPathAngleSet = false;
            verticalRotationalFlightPathAngleSet = NULL;
            std::cout<<"verticalRotationalFlightPathAngleSet der 3 = "<<verticalRotationalFlightPathAngleSet<<std::endl;


    return auxiliaryDerivativesVector;
}

//////////////////////////////////////////////// Auxiliary Functions //////////////////////////////////////////////////////////////////////

Eigen::MatrixXd getAuxiliaryFunctions( const tudat::basic_mathematics::Vector7d& aState, const double time, const Eigen::Vector3d& thrustAccelerationsBframe, const Eigen::VectorXd& auxiliaryEquationsVector,
                                         const Eigen::VectorXd& auxiliaryDerivativesVector){

    auxiliaryFunctionsMatrix = Eigen::MatrixXd::Zero(56,39);       // Setting the complete matrix and filling it with zeros for now

//    Eigen::MatrixXd thrustAzimuthMatrix = MAV.thrustAzimuth();      // Setting the thrust angle matrices
//    Eigen::MatrixXd thrustElevationMatrix = MAV.thrustElevation();

    // The following expressions are described in the order in which the equations have to be computed corresponding to the respective vector entry
    // Which in this case means that the not all positions in the matrix will be used. The other values will simply be 0.


    // w4
    auxiliaryFunctionsMatrix(4,0) = auxiliaryEquationsVector(27)/auxiliaryEquationsVector(7);   // Added because of the mistake found in the recurrence relation of W4,2
    auxiliaryFunctionsMatrix(4,1) = auxiliaryEquationsVector(1)/auxiliaryEquationsVector(9);
    auxiliaryFunctionsMatrix(4,2) = auxiliaryEquationsVector(0);
//    auxiliaryFunctionsMatrix(4,3) = cos(auxiliaryEquationsVector(10)+auxiliaryEquationsVector(11));
    auxiliaryFunctionsMatrix(4,3) = cos(auxiliaryEquationsVector(49));
    // Avoid cosine rounding errors
            if (abs(auxiliaryFunctionsMatrix(4,3))<6.2e-17){
              auxiliaryFunctionsMatrix(4,3) = 0;
            }
    auxiliaryFunctionsMatrix(4,4) = sin(auxiliaryEquationsVector(12));
//    std::cout<<"w4,4 (sin(x12)) = "<<auxiliaryFunctionsMatrix(4,4)<<std::endl;
    auxiliaryFunctionsMatrix(4,5) = cos(auxiliaryEquationsVector(13));
    // Avoid cosine rounding errors
            if (abs(auxiliaryFunctionsMatrix(4,5))<6.2e-17){
              auxiliaryFunctionsMatrix(4,5) = 0;
            }
    auxiliaryFunctionsMatrix(4,6) = cos(auxiliaryEquationsVector(12));
    // Avoid cosine rounding errors
            if (abs(auxiliaryFunctionsMatrix(4,6))<6.2e-17){
              auxiliaryFunctionsMatrix(4,6) = 0;
            }
//    std::cout<<"w4,6 (cos(x12)) = "<<auxiliaryFunctionsMatrix(4,6)<<std::endl;
    auxiliaryFunctionsMatrix(4,7) = sin(auxiliaryEquationsVector(14));

    auxiliaryFunctionsMatrix(4,38) = cos(auxiliaryEquationsVector(14));
    // Avoid cosine rounding errors
            if (abs(auxiliaryFunctionsMatrix(4,38))<6.2e-17){
              auxiliaryFunctionsMatrix(4,38) = 0;
            }
//    auxiliaryFunctionsMatrix(4,8) = sin(auxiliaryEquationsVector(10)+auxiliaryEquationsVector(11));
    auxiliaryFunctionsMatrix(4,8) = sin(auxiliaryEquationsVector(49));
    // Avoid sine rounding errors
    if (abs(auxiliaryFunctionsMatrix(4,8))<1.22e-16){
        auxiliaryFunctionsMatrix(4,8) = 0;
    }
    auxiliaryFunctionsMatrix(4,9) = sin(auxiliaryEquationsVector(13));
    auxiliaryFunctionsMatrix(4,10) = auxiliaryFunctionsMatrix(4,4)*auxiliaryFunctionsMatrix(4,5);
    auxiliaryFunctionsMatrix(4,11) = auxiliaryFunctionsMatrix(4,6)*auxiliaryFunctionsMatrix(4,7);
//    auxiliaryFunctionsMatrix(4,12) = auxiliaryFunctionsMatrix(4,9)*auxiliaryEquationsVector(16);
    auxiliaryFunctionsMatrix(4,12) = auxiliaryFunctionsMatrix(4,9)*auxiliaryFunctionsMatrix(4,38);
    auxiliaryFunctionsMatrix(4,13) = auxiliaryFunctionsMatrix(4,4)*auxiliaryFunctionsMatrix(4,9);
    auxiliaryFunctionsMatrix(4,14) = auxiliaryFunctionsMatrix(4,8)*auxiliaryFunctionsMatrix(4,5);
//    auxiliaryFunctionsMatrix(4,15) = auxiliaryFunctionsMatrix(4,6)*auxiliaryEquationsVector(16);
    auxiliaryFunctionsMatrix(4,15) = auxiliaryFunctionsMatrix(4,6)*auxiliaryFunctionsMatrix(4,38);
    auxiliaryFunctionsMatrix(4,16) = auxiliaryFunctionsMatrix(4,9)*auxiliaryFunctionsMatrix(4,7);
//    auxiliaryFunctionsMatrix(4,17) = auxiliaryFunctionsMatrix(4,10)*auxiliaryEquationsVector(16);
    auxiliaryFunctionsMatrix(4,17) = auxiliaryFunctionsMatrix(4,10)*auxiliaryFunctionsMatrix(4,38);
    auxiliaryFunctionsMatrix(4,18) = auxiliaryFunctionsMatrix(4,8)*auxiliaryFunctionsMatrix(4,12);
    auxiliaryFunctionsMatrix(4,19) = auxiliaryFunctionsMatrix(4,3)*auxiliaryFunctionsMatrix(4,13);
    auxiliaryFunctionsMatrix(4,20) = auxiliaryFunctionsMatrix(4,10)*auxiliaryFunctionsMatrix(4,7);
    auxiliaryFunctionsMatrix(4,21) = auxiliaryFunctionsMatrix(4,8)*auxiliaryFunctionsMatrix(4,16);
    auxiliaryFunctionsMatrix(4,22) = auxiliaryFunctionsMatrix(4,3)*(auxiliaryFunctionsMatrix(4,11)-auxiliaryFunctionsMatrix(4,17));
    auxiliaryFunctionsMatrix(4,23) = auxiliaryFunctionsMatrix(4,3)*(-auxiliaryFunctionsMatrix(4,20)-auxiliaryFunctionsMatrix(4,15));
    auxiliaryFunctionsMatrix(4,24) = auxiliaryFunctionsMatrix(4,2)*(auxiliaryFunctionsMatrix(4,22)-auxiliaryFunctionsMatrix(4,18));

    auxiliaryFunctionsMatrix(4,25) = Thrust/auxiliaryEquationsVector(7);
    auxiliaryFunctionsMatrix(4,26) = cos(thrustAzimuthMatrix(0,2));
    // Avoid cosine rounding errors
            if (abs(auxiliaryFunctionsMatrix(4,26))<6.2e-17){
              auxiliaryFunctionsMatrix(4,26) = 0;
            }

    auxiliaryFunctionsMatrix(4,27) = cos(thrustElevationMatrix(0,2));
    // Avoid cosine rounding errors
            if (abs(auxiliaryFunctionsMatrix(4,27))<6.2e-17){
              auxiliaryFunctionsMatrix(4,27) = 0;
            }

    auxiliaryFunctionsMatrix(4,28) = sin(thrustAzimuthMatrix(0,2));
    auxiliaryFunctionsMatrix(4,29) = sin(thrustElevationMatrix(0,2));
    auxiliaryFunctionsMatrix(4,30) = auxiliaryFunctionsMatrix(4,26)*auxiliaryFunctionsMatrix(4,27);
    auxiliaryFunctionsMatrix(4,31) = auxiliaryFunctionsMatrix(4,28)*auxiliaryFunctionsMatrix(4,27);
    auxiliaryFunctionsMatrix(4,32) = auxiliaryFunctionsMatrix(4,25)*auxiliaryFunctionsMatrix(4,30);
    auxiliaryFunctionsMatrix(4,33) = auxiliaryFunctionsMatrix(4,25)*auxiliaryFunctionsMatrix(4,31);
    auxiliaryFunctionsMatrix(4,34) = auxiliaryFunctionsMatrix(4,25)*auxiliaryFunctionsMatrix(4,29);
    auxiliaryFunctionsMatrix(4,35) = auxiliaryFunctionsMatrix(4,33)*(auxiliaryFunctionsMatrix(4,19)-auxiliaryFunctionsMatrix(4,14));
    auxiliaryFunctionsMatrix(4,36) = auxiliaryFunctionsMatrix(4,34)*(auxiliaryFunctionsMatrix(4,23)-auxiliaryFunctionsMatrix(4,21));
    auxiliaryFunctionsMatrix(4,37) = 1.0/auxiliaryEquationsVector(7);



    // w5
    auxiliaryFunctionsMatrix(5,1) = auxiliaryEquationsVector(2)/auxiliaryEquationsVector(9);
    auxiliaryFunctionsMatrix(5,2) = auxiliaryFunctionsMatrix(4,8)*(auxiliaryFunctionsMatrix(4,11)-auxiliaryFunctionsMatrix(4,17));
    auxiliaryFunctionsMatrix(5,3) = auxiliaryFunctionsMatrix(4,3)*auxiliaryFunctionsMatrix(4,12);
    auxiliaryFunctionsMatrix(5,4) = auxiliaryFunctionsMatrix(4,8)*auxiliaryFunctionsMatrix(4,13);
    auxiliaryFunctionsMatrix(5,5) = auxiliaryFunctionsMatrix(4,3)*auxiliaryFunctionsMatrix(4,5);
    auxiliaryFunctionsMatrix(5,6) = auxiliaryFunctionsMatrix(4,8)*(-auxiliaryFunctionsMatrix(4,20)-auxiliaryFunctionsMatrix(4,11));
    auxiliaryFunctionsMatrix(5,7) = auxiliaryFunctionsMatrix(4,3)*auxiliaryFunctionsMatrix(4,16);
    auxiliaryFunctionsMatrix(5,8) = auxiliaryFunctionsMatrix(4,2)*(auxiliaryFunctionsMatrix(5,2)+auxiliaryFunctionsMatrix(5,3));
    auxiliaryFunctionsMatrix(5,9) = auxiliaryFunctionsMatrix(4,33)*(auxiliaryFunctionsMatrix(5,4)+auxiliaryFunctionsMatrix(5,5));
    auxiliaryFunctionsMatrix(5,10) = auxiliaryFunctionsMatrix(4,34)*(auxiliaryFunctionsMatrix(5,6)+auxiliaryFunctionsMatrix(5,7));


    // w6
    auxiliaryFunctionsMatrix(6,0) = auxiliaryFunctionsMatrix(4,4)*auxiliaryFunctionsMatrix(4,7);  // Added because of the mistake found in the complete transformation matrix
    auxiliaryFunctionsMatrix(6,1) = auxiliaryEquationsVector(3)/auxiliaryEquationsVector(9);
    auxiliaryFunctionsMatrix(6,2) = auxiliaryFunctionsMatrix(4,5)*auxiliaryFunctionsMatrix(4,15);
    auxiliaryFunctionsMatrix(6,3) = auxiliaryFunctionsMatrix(4,6)*auxiliaryFunctionsMatrix(4,9);
    auxiliaryFunctionsMatrix(6,4) = auxiliaryFunctionsMatrix(4,5)*auxiliaryFunctionsMatrix(4,11);
//    auxiliaryFunctionsMatrix(6,5) = auxiliaryFunctionsMatrix(4,4)*auxiliaryEquationsVector(16);
    auxiliaryFunctionsMatrix(6,5) = auxiliaryFunctionsMatrix(4,4)*auxiliaryFunctionsMatrix(4,38);
    auxiliaryFunctionsMatrix(6,6) = auxiliaryFunctionsMatrix(4,2)*(auxiliaryFunctionsMatrix(6,2)+auxiliaryFunctionsMatrix(6,0)); // Changed becuase of the mistake found in the complete transformation matrix
    auxiliaryFunctionsMatrix(6,7) = auxiliaryFunctionsMatrix(4,33)*auxiliaryFunctionsMatrix(6,3);
    auxiliaryFunctionsMatrix(6,8) = auxiliaryFunctionsMatrix(4,34)*(auxiliaryFunctionsMatrix(6,4)-auxiliaryFunctionsMatrix(6,5));



    // w8
    auxiliaryFunctionsMatrix(8,1) = auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4);
    auxiliaryFunctionsMatrix(8,2) = auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5);
    auxiliaryFunctionsMatrix(8,3) = auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6);





    // w9
    auxiliaryFunctionsMatrix(9,1) = (auxiliaryEquationsVector(9)*auxiliaryDerivativesVector(8))/auxiliaryEquationsVector(8);




/*    // w11
    auxiliaryFunctionsMatrix(11,1) = auxiliaryEquationsVector(1)*auxiliaryEquationsVector(5);
    auxiliaryFunctionsMatrix(11,2) = auxiliaryEquationsVector(2)*auxiliaryEquationsVector(4);


    // Avoiding singularities
    if (auxiliaryEquationsVector(19) == 0){

        auxiliaryFunctionsMatrix(11,3) = 0;
    }
    else {
    auxiliaryFunctionsMatrix(11,3) = (auxiliaryFunctionsMatrix(11,1)-auxiliaryFunctionsMatrix(11,2))/auxiliaryEquationsVector(19);

};
*/


    // w12
    auxiliaryFunctionsMatrix(12,1) = auxiliaryEquationsVector(20)*auxiliaryEquationsVector(6);
    auxiliaryFunctionsMatrix(12,2) = auxiliaryEquationsVector(3)*auxiliaryEquationsVector(25);
    auxiliaryFunctionsMatrix(12,3) = auxiliaryEquationsVector(3)/auxiliaryEquationsVector(20);
    auxiliaryFunctionsMatrix(12,4) = auxiliaryFunctionsMatrix(12,3)*auxiliaryFunctionsMatrix(12,3);
//    std::cout<<"w12,3 = "<<auxiliaryFunctionsMatrix(12,3)<<std::endl;
//    std::cout<<"w12,4 = "<<auxiliaryFunctionsMatrix(12,4)<<std::endl;

    auxiliaryFunctionsMatrix(12,5) = sqrt(1-auxiliaryFunctionsMatrix(12,4));
    auxiliaryFunctionsMatrix(12,6) = auxiliaryEquationsVector(8)*auxiliaryFunctionsMatrix(12,5);


 /*
    std::cout<<"x8 = "<<auxiliaryEquationsVector(8)<<std::endl;
    std::cout<<"w12,5 = "<<auxiliaryFunctionsMatrix(12,5)<<std::endl;
    std::cout<<"x8-11528741.16 = "<<auxiliaryEquationsVector(8)-11528741.16<<std::endl;
    std::cout<<"w12,5-0.933580426497202 = "<<auxiliaryFunctionsMatrix(12,5)-0.933580426497202<<std::endl;

    std::cout<<"w12,6 = "<<auxiliaryFunctionsMatrix(12,6)<<std::endl;
    std::cout<<"w12,6-10763007.0891286 = "<<auxiliaryFunctionsMatrix(12,6)-10763007.0891286<<std::endl;
//*/

    // Avoiding singularities
    if (auxiliaryFunctionsMatrix(12,6) == 0){

        auxiliaryFunctionsMatrix(12,7) = 0;
    }
    else {

    auxiliaryFunctionsMatrix(12,7) = (auxiliaryFunctionsMatrix(12,1)-auxiliaryFunctionsMatrix(12,2))/(auxiliaryFunctionsMatrix(12,6));

};


/*
    // w13
    auxiliaryFunctionsMatrix(13,1) = auxiliaryEquationsVector(22)*auxiliaryEquationsVector(22);
    auxiliaryFunctionsMatrix(13,2) = sqrt(1-auxiliaryFunctionsMatrix(13,1));


    // Avoiding singularities
    if (auxiliaryFunctionsMatrix(13,2) == 0){

        auxiliaryFunctionsMatrix(13,3) = 0;
    }
    else {
    auxiliaryFunctionsMatrix(13,3) = (auxiliaryDerivativesVector(22))/auxiliaryFunctionsMatrix(13,2);
};
*/

    // w13
    auxiliaryFunctionsMatrix(13,1) = auxiliaryEquationsVector(24)*auxiliaryDerivativesVector(48);
    auxiliaryFunctionsMatrix(13,2) = auxiliaryEquationsVector(48)*auxiliaryDerivativesVector(24);
    auxiliaryFunctionsMatrix(13,3) = auxiliaryEquationsVector(48)*auxiliaryEquationsVector(48);
    auxiliaryFunctionsMatrix(13,4) = auxiliaryEquationsVector(24)*auxiliaryEquationsVector(24);

    // Avoiding singularities
    if (auxiliaryFunctionsMatrix(13,3)+auxiliaryFunctionsMatrix(13,4) == 0){

        auxiliaryFunctionsMatrix(13,5) = 0;
    }
    else {
    auxiliaryFunctionsMatrix(13,5) = (auxiliaryFunctionsMatrix(13,1)-auxiliaryFunctionsMatrix(13,2))/(auxiliaryFunctionsMatrix(13,3)+auxiliaryFunctionsMatrix(13,4));
};



    // w14
    auxiliaryFunctionsMatrix(14,1) = auxiliaryEquationsVector(23)*auxiliaryEquationsVector(23);
    auxiliaryFunctionsMatrix(14,2) = sqrt(1.0-auxiliaryFunctionsMatrix(14,1));
//    auxiliaryFunctionsMatrix(14,2) = 1/sqrt(1-auxiliaryFunctionsMatrix(14,1));

//    std::cout<<"1-w14,1 = "<<1-auxiliaryFunctionsMatrix(14,1)<<std::endl;

    // Avoiding singularities
    if (auxiliaryFunctionsMatrix(14,2) == 0){

        auxiliaryFunctionsMatrix(14,3) = 0;
    }
    else {
    auxiliaryFunctionsMatrix(14,3) = (auxiliaryDerivativesVector(23))/auxiliaryFunctionsMatrix(14,2);
};

//        // Avoiding singularities
//        if (1-auxiliaryFunctionsMatrix(14,1) == 0){

//            auxiliaryFunctionsMatrix(14,2) = 0;
//            auxiliaryFunctionsMatrix(14,3) = 0;
//        }
//        else {
//        auxiliaryFunctionsMatrix(14,3) = (auxiliaryDerivativesVector(23))*auxiliaryFunctionsMatrix(14,2);
//    };



    // w15
    auxiliaryFunctionsMatrix(15,1) = sqrt(rotationalVelocity*rotationalVelocity*auxiliaryEquationsVector(19)+4.0*auxiliaryEquationsVector(54)+4.0*rotationalVelocity*auxiliaryEquationsVector(55));
    auxiliaryFunctionsMatrix(15,2) = 1/auxiliaryEquationsVector(8);

//    auxiliaryFunctionsMatrix(15,1) = auxiliaryEquationsVector(35)*auxiliaryDerivativesVector(35);

//    // Avoiding singularities
//    if (auxiliaryEquationsVector(15) == 0){

//        auxiliaryFunctionsMatrix(15,2) = 0;
//    }
//    else {
//    auxiliaryFunctionsMatrix(15,2) = (2.0*auxiliaryFunctionsMatrix(15,1)+auxiliaryDerivativesVector(21)-2.0*auxiliaryDerivativesVector(43))/(auxiliaryEquationsVector(15));
//};

/*    /// Debug ///

    std::cout<<"w15,2 = "<<auxiliaryFunctionsMatrix(15,2)<<std::endl;
         std::cout<<"x15 = "<<auxiliaryEquationsVector(15)<<std::endl;
    std::cout<<"2*auxiliaryFunctionsMatrix(15,1)+auxiliaryDerivativesVector(21)-2*auxiliaryDerivativesVector(43) = "<<2*auxiliaryFunctionsMatrix(15,1)+auxiliaryDerivativesVector(21)-2*auxiliaryDerivativesVector(43)<<std::endl;
    std::cout<<"calc w15,2 = "<<(2*auxiliaryFunctionsMatrix(15,1)+auxiliaryDerivativesVector(21)-2*auxiliaryDerivativesVector(43))/(0)<<std::endl;
*/


    // w16
    auxiliaryFunctionsMatrix(16,1) = auxiliaryFunctionsMatrix(4,7)*auxiliaryDerivativesVector(14);





/*    // w17
    auxiliaryFunctionsMatrix(17,1) = auxiliaryEquationsVector(16)*auxiliaryDerivativesVector(15);
    auxiliaryFunctionsMatrix(17,2) = auxiliaryEquationsVector(15)*auxiliaryDerivativesVector(16);
*/




/*    // w18
    auxiliaryFunctionsMatrix(18,1) = auxiliaryEquationsVector(20)*auxiliaryDerivativesVector(24);
    auxiliaryFunctionsMatrix(18,2) = auxiliaryEquationsVector(24)*auxiliaryEquationsVector(25);
*/




    // w20
    auxiliaryFunctionsMatrix(20,1) = auxiliaryEquationsVector(26)/auxiliaryEquationsVector(20);




    // w21
    auxiliaryFunctionsMatrix(21,1) = auxiliaryEquationsVector(4)*auxiliaryDerivativesVector(4);
    auxiliaryFunctionsMatrix(21,2) = auxiliaryEquationsVector(5)*auxiliaryDerivativesVector(5);
    auxiliaryFunctionsMatrix(21,3) = auxiliaryEquationsVector(6)*auxiliaryDerivativesVector(6);





/*    // w22
    auxiliaryFunctionsMatrix(22,1) = auxiliaryEquationsVector(17)*auxiliaryDerivativesVector(18);
    auxiliaryFunctionsMatrix(22,2) = auxiliaryEquationsVector(18)*auxiliaryDerivativesVector(17);
    auxiliaryFunctionsMatrix(22,3) = auxiliaryEquationsVector(17)*auxiliaryEquationsVector(17);

    // Avoiding singularities
    if (auxiliaryFunctionsMatrix(22,3) == 0){

        auxiliaryFunctionsMatrix(22,4) = 0;
    }
    else {
    auxiliaryFunctionsMatrix(22,4) = (auxiliaryFunctionsMatrix(22,1)-auxiliaryFunctionsMatrix(22,2))/(auxiliaryFunctionsMatrix(22,3));
};

*/


    // w23
    auxiliaryFunctionsMatrix(23,1) = auxiliaryEquationsVector(15)*auxiliaryDerivativesVector(25);
    auxiliaryFunctionsMatrix(23,2) = auxiliaryEquationsVector(25)*auxiliaryDerivativesVector(15);
    auxiliaryFunctionsMatrix(23,3) = auxiliaryEquationsVector(15)*auxiliaryEquationsVector(15);

    // Avoiding singularities
    if (auxiliaryFunctionsMatrix(23,3) == 0){

        auxiliaryFunctionsMatrix(23,4) = 0;
    }
    else {
    auxiliaryFunctionsMatrix(23,4) = (auxiliaryFunctionsMatrix(23,1)-auxiliaryFunctionsMatrix(23,2))/(auxiliaryFunctionsMatrix(23,3));
};




    // w24
//    auxiliaryFunctionsMatrix(24,1) = auxiliaryEquationsVector(6)*auxiliaryDerivativesVector(20);      // Not needed anymore
    auxiliaryFunctionsMatrix(24,2) = auxiliaryEquationsVector(20)*auxiliaryDerivativesVector(6);
//    auxiliaryFunctionsMatrix(24,3) = auxiliaryEquationsVector(25)*auxiliaryEquationsVector(6);        // Not needed anymore
    auxiliaryFunctionsMatrix(24,4) = auxiliaryEquationsVector(3)*auxiliaryDerivativesVector(25);

    // Avoiding singularities
    if (auxiliaryFunctionsMatrix(12,6) == 0){

        auxiliaryFunctionsMatrix(24,5) = 0;
    }
    else {
//    auxiliaryFunctionsMatrix(24,5) = (auxiliaryFunctionsMatrix(24,1)+auxiliaryFunctionsMatrix(24,2)-auxiliaryFunctionsMatrix(24,3)-auxiliaryFunctionsMatrix(24,4))/(auxiliaryFunctionsMatrix(12,6));  // Wrong
        auxiliaryFunctionsMatrix(24,5) = (auxiliaryFunctionsMatrix(24,2)-auxiliaryFunctionsMatrix(24,4))/(auxiliaryFunctionsMatrix(12,6));

    };

//    auxiliaryFunctionsMatrix(24,6) = auxiliaryEquationsVector(3)*auxiliaryFunctionsMatrix(24,4);  // Wrong
    auxiliaryFunctionsMatrix(24,6) = auxiliaryEquationsVector(3)*auxiliaryFunctionsMatrix(12,2);
    auxiliaryFunctionsMatrix(24,7) = pow(auxiliaryEquationsVector(20),3.0);
    auxiliaryFunctionsMatrix(24,8) = auxiliaryFunctionsMatrix(8,3)/auxiliaryEquationsVector(8);
    auxiliaryFunctionsMatrix(24,9) = auxiliaryFunctionsMatrix(24,6)/auxiliaryFunctionsMatrix(24,7);
    auxiliaryFunctionsMatrix(24,10) = (2.0*auxiliaryFunctionsMatrix(24,9)-2.0*auxiliaryFunctionsMatrix(24,8))*(auxiliaryFunctionsMatrix(12,1)-auxiliaryFunctionsMatrix(12,2));
    auxiliaryFunctionsMatrix(24,11) = auxiliaryDerivativesVector(8)*(auxiliaryFunctionsMatrix(12,1)-auxiliaryFunctionsMatrix(12,2));
    auxiliaryFunctionsMatrix(24,12) = pow((1.0-auxiliaryFunctionsMatrix(12,4)),1.5);
    auxiliaryFunctionsMatrix(24,13) = auxiliaryEquationsVector(8)*auxiliaryEquationsVector(8);
    auxiliaryFunctionsMatrix(24,14) = auxiliaryEquationsVector(8)*auxiliaryFunctionsMatrix(24,12);
    auxiliaryFunctionsMatrix(24,15) = auxiliaryFunctionsMatrix(24,13)*auxiliaryFunctionsMatrix(12,5);

    // Avoiding singularities
    if (auxiliaryFunctionsMatrix(24,14) == 0){

        auxiliaryFunctionsMatrix(24,16) = 0;
    }
    else {
    auxiliaryFunctionsMatrix(24,16) = auxiliaryFunctionsMatrix(24,10)/auxiliaryFunctionsMatrix(24,14);
    };

    // Avoiding singularities
    if (auxiliaryFunctionsMatrix(24,15) == 0){

        auxiliaryFunctionsMatrix(24,17) = 0;
    }
    else {
    auxiliaryFunctionsMatrix(24,17) = auxiliaryFunctionsMatrix(24,11)/auxiliaryFunctionsMatrix(24,15);
};





    // w25
    auxiliaryFunctionsMatrix(25,1) = auxiliaryEquationsVector(8)*auxiliaryDerivativesVector(26);
    auxiliaryFunctionsMatrix(25,2) = auxiliaryEquationsVector(26)*auxiliaryEquationsVector(26);
    auxiliaryFunctionsMatrix(25,3) = (2.0*auxiliaryFunctionsMatrix(25,1)-auxiliaryFunctionsMatrix(25,2))/(auxiliaryEquationsVector(9));





    // w26
    auxiliaryFunctionsMatrix(26,1) = auxiliaryEquationsVector(1)*auxiliaryDerivativesVector(4);
    auxiliaryFunctionsMatrix(26,2) = auxiliaryEquationsVector(2)*auxiliaryDerivativesVector(5);
    auxiliaryFunctionsMatrix(26,3) = auxiliaryEquationsVector(3)*auxiliaryDerivativesVector(6);





    // w27
    auxiliaryFunctionsMatrix(27,1) = auxiliaryEquationsVector(29)*auxiliaryDerivativesVector(28);
    auxiliaryFunctionsMatrix(27,2) = auxiliaryEquationsVector(28)*auxiliaryDerivativesVector(29);
    auxiliaryFunctionsMatrix(27,3) = auxiliaryEquationsVector(28)*auxiliaryEquationsVector(29);
    auxiliaryFunctionsMatrix(27,4) = auxiliaryEquationsVector(15)*(auxiliaryFunctionsMatrix(27,1)+auxiliaryFunctionsMatrix(27,2));
    auxiliaryFunctionsMatrix(27,5) = auxiliaryFunctionsMatrix(27,3)*auxiliaryDerivativesVector(15);
    auxiliaryFunctionsMatrix(27,6) = auxiliaryEquationsVector(15)*(auxiliaryFunctionsMatrix(27,4)+auxiliaryFunctionsMatrix(27,5));





    // w28
    auxiliaryFunctionsMatrix(28,1) = auxiliaryDerivativesVector(30)*auxiliaryEquationsVector(28);






    // w30
    auxiliaryFunctionsMatrix(30,1) = pow(auxiliaryEquationsVector(31),9.0);
    auxiliaryFunctionsMatrix(30,2) = pow(auxiliaryEquationsVector(31),8.0);
    auxiliaryFunctionsMatrix(30,3) = pow(auxiliaryEquationsVector(31),7.0);
    auxiliaryFunctionsMatrix(30,4) = pow(auxiliaryEquationsVector(31),6.0);
    auxiliaryFunctionsMatrix(30,5) = pow(auxiliaryEquationsVector(31),5.0);
    auxiliaryFunctionsMatrix(30,6) = pow(auxiliaryEquationsVector(31),4.0);
    auxiliaryFunctionsMatrix(30,7) = pow(auxiliaryEquationsVector(31),3.0);
    auxiliaryFunctionsMatrix(30,8) = pow(auxiliaryEquationsVector(31),2.0);

    for (int i=2; i<10+1; i++){

        if (i==2){

            auxiliaryFunctionsMatrix(30,9) = auxiliaryDerivativesVector(31)*(2.0*densityPolyCoefficients(2)*auxiliaryEquationsVector(31)+densityPolyCoefficients(1));
        }
        else {

           auxiliaryFunctionsMatrix(30,9) += auxiliaryDerivativesVector(31)*i*densityPolyCoefficients(i)*auxiliaryFunctionsMatrix(30,(11-i));

        };
    };






    // w32
    auxiliaryFunctionsMatrix(32,1) = auxiliaryEquationsVector(33)*auxiliaryDerivativesVector(15);
    auxiliaryFunctionsMatrix(32,2) = auxiliaryEquationsVector(15)*auxiliaryDerivativesVector(33);
    auxiliaryFunctionsMatrix(32,3) = auxiliaryEquationsVector(33)*auxiliaryEquationsVector(33);
    auxiliaryFunctionsMatrix(32,4) = (auxiliaryFunctionsMatrix(32,1)-auxiliaryFunctionsMatrix(32,2))/(auxiliaryFunctionsMatrix(32,3));





    // w33
    auxiliaryFunctionsMatrix(33,1) = auxiliaryDerivativesVector(34)/auxiliaryEquationsVector(33);




    // w34
    // First it has to be determined whether or not these functions have to be computed. If not, they remain zero. If so they only one of them is computed.

    if (temperatureAltitudeRanges(1,0)<=auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31)<temperatureAltitudeRanges(1,1)){

      auxiliaryFunctionsMatrix(34,2) = auxiliaryDerivativesVector(31)*(3.0*temperaturePolyCoefficients(3,2)*auxiliaryFunctionsMatrix(30,8)+2.0*temperaturePolyCoefficients(2,2)*auxiliaryEquationsVector(31)+temperaturePolyCoefficients(1,2));
    }
    else if (temperatureAltitudeRanges(2,0)<=auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31)<temperatureAltitudeRanges(2,1)){

        for (int i=2; i<6+1;i++){

            if (i==2){

                auxiliaryFunctionsMatrix(34,3) = auxiliaryDerivativesVector(31)*(2*temperaturePolyCoefficients(2,3)*auxiliaryEquationsVector(31)+temperaturePolyCoefficients(1,3));
            }


            else {

                auxiliaryFunctionsMatrix(34,3) += auxiliaryDerivativesVector(31)*i*temperaturePolyCoefficients(i,3)*auxiliaryFunctionsMatrix(30,(11-i));
            };

        };



    }
    else if (temperatureAltitudeRanges(3,0)<=auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31)<temperatureAltitudeRanges(3,1)){

        for (int i=2; i<8+1;i++){

            if (i==2){

                auxiliaryFunctionsMatrix(34,4) = auxiliaryDerivativesVector(31)*(2.0*temperaturePolyCoefficients(2,4)*auxiliaryEquationsVector(31)+temperaturePolyCoefficients(1,4));
            }


            else {

                auxiliaryFunctionsMatrix(34,4) += auxiliaryDerivativesVector(31)*i*temperaturePolyCoefficients(i,4)*auxiliaryFunctionsMatrix(30,(11-i));
            };

        };
};



    // w35
    auxiliaryFunctionsMatrix(35,1) = auxiliaryFunctionsMatrix(4,6)*auxiliaryEquationsVector(25);
    auxiliaryFunctionsMatrix(35,2) = auxiliaryFunctionsMatrix(4,4)*auxiliaryEquationsVector(24);
    auxiliaryFunctionsMatrix(35,3) = auxiliaryFunctionsMatrix(35,2)*auxiliaryEquationsVector(12);





    // w36
    auxiliaryFunctionsMatrix(36,1) = auxiliaryDerivativesVector(21)/auxiliaryEquationsVector(36);




    // w37
    auxiliaryFunctionsMatrix(37,1) = auxiliaryEquationsVector(36)*auxiliaryDerivativesVector(25);
    auxiliaryFunctionsMatrix(37,2) = auxiliaryEquationsVector(25)*auxiliaryDerivativesVector(36);
    auxiliaryFunctionsMatrix(37,3) = auxiliaryEquationsVector(36)*auxiliaryEquationsVector(36);
    auxiliaryFunctionsMatrix(37,4) = (auxiliaryFunctionsMatrix(37,1)-auxiliaryFunctionsMatrix(37,2))/(auxiliaryFunctionsMatrix(37,3));





    // w38
    auxiliaryFunctionsMatrix(38,1) = auxiliaryEquationsVector(37)*auxiliaryEquationsVector(37);
    auxiliaryFunctionsMatrix(38,2) = sqrt(1.0-auxiliaryFunctionsMatrix(38,1));

    // Avoiding singularities
    if (auxiliaryFunctionsMatrix(38,2) == 0){

        auxiliaryFunctionsMatrix(38,3) = 0;
    }
    else {
    auxiliaryFunctionsMatrix(38,3) = auxiliaryDerivativesVector(37)/auxiliaryFunctionsMatrix(38,2);
};




/*    // w39
    auxiliaryFunctionsMatrix(39,1) = auxiliaryEquationsVector(44)*auxiliaryDerivativesVector(18);
    auxiliaryFunctionsMatrix(39,2) = auxiliaryEquationsVector(18)*auxiliaryDerivativesVector(44);
    auxiliaryFunctionsMatrix(39,3) = auxiliaryEquationsVector(44)*auxiliaryEquationsVector(44);
    auxiliaryFunctionsMatrix(39,4) = (auxiliaryFunctionsMatrix(39,1)-auxiliaryFunctionsMatrix(39,2))/(auxiliaryFunctionsMatrix(39,3));

*/


/*
    // w40
    auxiliaryFunctionsMatrix(40,1) = auxiliaryEquationsVector(39)*auxiliaryEquationsVector(39);
    auxiliaryFunctionsMatrix(40,2) = sqrt(1-auxiliaryFunctionsMatrix(40,1));

    // Avoiding singularities
    if (auxiliaryFunctionsMatrix(40,2) == 0){

        auxiliaryFunctionsMatrix(40,3) = 0;
    }
    else {
    auxiliaryFunctionsMatrix(40,3) = auxiliaryDerivativesVector(39)/auxiliaryFunctionsMatrix(40,2);
};

*/

    // w40
    auxiliaryFunctionsMatrix(40,1) = auxiliaryEquationsVector(24)*auxiliaryDerivativesVector(47);
    auxiliaryFunctionsMatrix(40,2) = auxiliaryEquationsVector(47)*auxiliaryDerivativesVector(24);
    auxiliaryFunctionsMatrix(40,3) = auxiliaryEquationsVector(47)*auxiliaryEquationsVector(47);


    // Avoiding singularities
    if (auxiliaryFunctionsMatrix(40,3)+auxiliaryFunctionsMatrix(40,4) == 0){

        auxiliaryFunctionsMatrix(40,4) = 0;
    }
    else {
    auxiliaryFunctionsMatrix(40,4) = (auxiliaryFunctionsMatrix(40,1)-auxiliaryFunctionsMatrix(40,2))/(auxiliaryFunctionsMatrix(40,3)+auxiliaryFunctionsMatrix(13,4));
};


    // w41
    auxiliaryFunctionsMatrix(41,1) = auxiliaryEquationsVector(36)*auxiliaryDerivativesVector(35);
    auxiliaryFunctionsMatrix(41,2) = auxiliaryEquationsVector(35)*auxiliaryDerivativesVector(36);





    // w42
    auxiliaryFunctionsMatrix(42,1) = cos(auxiliaryEquationsVector(38));
    // Avoid cosine rounding errors
            if (abs(auxiliaryFunctionsMatrix(42,1))<6.2e-17){
              auxiliaryFunctionsMatrix(42,1) = 0;
            }
    auxiliaryFunctionsMatrix(42,2) = cos(auxiliaryEquationsVector(40));
    // Avoid cosine rounding errors
            if (abs(auxiliaryFunctionsMatrix(42,2))<6.2e-17){
              auxiliaryFunctionsMatrix(42,2) = 0;
            }
    auxiliaryFunctionsMatrix(42,3) = sin(auxiliaryEquationsVector(38));
    auxiliaryFunctionsMatrix(42,4) = sin(auxiliaryEquationsVector(40));
    auxiliaryFunctionsMatrix(42,5) = auxiliaryFunctionsMatrix(42,2)*auxiliaryDerivativesVector(40);
    auxiliaryFunctionsMatrix(42,6) = auxiliaryFunctionsMatrix(42,3)*auxiliaryDerivativesVector(38);
    auxiliaryFunctionsMatrix(42,7) = auxiliaryFunctionsMatrix(42,1)*auxiliaryFunctionsMatrix(42,5);
    auxiliaryFunctionsMatrix(42,8) = auxiliaryFunctionsMatrix(42,6)*auxiliaryFunctionsMatrix(42,4);





    // w43
    auxiliaryFunctionsMatrix(43,1) = auxiliaryEquationsVector(41)*auxiliaryDerivativesVector(42);
    auxiliaryFunctionsMatrix(43,2) = auxiliaryEquationsVector(42)*auxiliaryDerivativesVector(41);





/*    // w44
    auxiliaryFunctionsMatrix(44,1) = auxiliaryEquationsVector(36)*auxiliaryDerivativesVector(38);
    auxiliaryFunctionsMatrix(44,2) = auxiliaryDerivativesVector(36)*auxiliaryFunctionsMatrix(42,1);
    auxiliaryFunctionsMatrix(44,3) = auxiliaryFunctionsMatrix(44,1)*auxiliaryFunctionsMatrix(42,3);
*/

    // w45
    auxiliaryFunctionsMatrix(45,1) = auxiliaryEquationsVector(1)*auxiliaryDerivativesVector(5);
    auxiliaryFunctionsMatrix(45,2) = auxiliaryEquationsVector(2)*auxiliaryDerivativesVector(4);
    auxiliaryFunctionsMatrix(45,3) = auxiliaryEquationsVector(1)*auxiliaryEquationsVector(5);
    auxiliaryFunctionsMatrix(45,4) = auxiliaryEquationsVector(2)*auxiliaryEquationsVector(4);
    auxiliaryFunctionsMatrix(45,5) = auxiliaryEquationsVector(19)*auxiliaryEquationsVector(19);
    auxiliaryFunctionsMatrix(45,6) = auxiliaryEquationsVector(19)*(auxiliaryFunctionsMatrix(45,1)-auxiliaryFunctionsMatrix(45,2));
    auxiliaryFunctionsMatrix(45,7) = auxiliaryDerivativesVector(19)*(auxiliaryFunctionsMatrix(45,3)-auxiliaryFunctionsMatrix(45,4));

    // Avoiding singularities
    if (auxiliaryFunctionsMatrix(45,5) == 0){

        auxiliaryFunctionsMatrix(45,8) = 0;
    }
    else {
    auxiliaryFunctionsMatrix(45,8) = (auxiliaryFunctionsMatrix(45,6)-auxiliaryFunctionsMatrix(45,7))/auxiliaryFunctionsMatrix(45,5);
};

    // w47
    auxiliaryFunctionsMatrix(47,1) =  auxiliaryDerivativesVector(45)*auxiliaryFunctionsMatrix(4,6);
    auxiliaryFunctionsMatrix(47,2) =  auxiliaryDerivativesVector(12)*auxiliaryFunctionsMatrix(4,4);
    auxiliaryFunctionsMatrix(47,3) =  auxiliaryEquationsVector(45)*auxiliaryFunctionsMatrix(47,2);

    // w48
    auxiliaryFunctionsMatrix(48,1) =  auxiliaryDerivativesVector(46)*auxiliaryFunctionsMatrix(4,6);
    auxiliaryFunctionsMatrix(48,2) =  auxiliaryEquationsVector(46)*auxiliaryFunctionsMatrix(47,2);


    // w52
    auxiliaryFunctionsMatrix(52,1) = cos(auxiliaryEquationsVector(10));
    auxiliaryFunctionsMatrix(52,2) = sin(auxiliaryEquationsVector(10));
    auxiliaryFunctionsMatrix(52,3) = auxiliaryDerivativesVector(10)*auxiliaryFunctionsMatrix(52,2);
    auxiliaryFunctionsMatrix(52,4) = auxiliaryDerivativesVector(10)*auxiliaryFunctionsMatrix(52,1);
    auxiliaryFunctionsMatrix(52,5) = auxiliaryFunctionsMatrix(52,1)*(auxiliaryDerivativesVector(4)+rotationalVelocity*auxiliaryEquationsVector(5));
    auxiliaryFunctionsMatrix(52,6) = auxiliaryFunctionsMatrix(52,3)*(auxiliaryEquationsVector(4)+rotationalVelocity*auxiliaryEquationsVector(2));
    auxiliaryFunctionsMatrix(52,7) = auxiliaryFunctionsMatrix(52,2)*(auxiliaryDerivativesVector(5)-rotationalVelocity*auxiliaryEquationsVector(4));
    auxiliaryFunctionsMatrix(52,8) = auxiliaryFunctionsMatrix(52,4)*(auxiliaryEquationsVector(5)-rotationalVelocity*auxiliaryEquationsVector(1));



    // w53
    auxiliaryFunctionsMatrix(53,1) = auxiliaryFunctionsMatrix(52,2)*(auxiliaryDerivativesVector(4)+rotationalVelocity*auxiliaryEquationsVector(5));
    auxiliaryFunctionsMatrix(53,2) = auxiliaryFunctionsMatrix(52,4)*(auxiliaryEquationsVector(4)+rotationalVelocity*auxiliaryEquationsVector(2));
    auxiliaryFunctionsMatrix(53,3) = auxiliaryFunctionsMatrix(52,1)*(auxiliaryDerivativesVector(5)-rotationalVelocity*auxiliaryEquationsVector(4));
    auxiliaryFunctionsMatrix(53,4) = auxiliaryFunctionsMatrix(52,3)*(auxiliaryEquationsVector(5)-rotationalVelocity*auxiliaryEquationsVector(1));


    // w55
    auxiliaryFunctionsMatrix(55,1) = auxiliaryDerivativesVector(53)*auxiliaryEquationsVector(50);
    auxiliaryFunctionsMatrix(55,2) = auxiliaryDerivativesVector(52)*auxiliaryEquationsVector(51);



// auxiliaryFunctionsMatrix(,)


    return auxiliaryFunctionsMatrix;

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
 Eigen::MatrixXd thrustAzimuthMatrix;                           // psi_T [rad] thrust azimuth angles
 Eigen::MatrixXd thrustElevationMatrix;                         // epsilon_T [rad] thrust elevation angles
 double specificImpulse;                                        // Isp [s]    engine nominal specific impulse
 double referenceArea;                                           // S [m^2]  vehicle reference area
 Eigen::MatrixXd dragCoefficientPolyCoefficients;           // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
 Eigen::MatrixXd dragCoefficientMachRanges;                       // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient

    // Set functions

 double rotationalFlightPathAngle;         // Rotational flight path angle in rad
 double inertialFlightPathAngle;           // Inertial flight path angle in rad
 double rotationalHeadingAngle;            // Rotational heading angle in rad
 double inertialHeadingAngle;              // Inertial heading angle in rad

 bool rotationalFlightPathAngleSet;         // All of these are used to let the program know that a predefined angle was set and that that angle should be used (initially)
 bool inertialFlightPathAngleSet;
 bool rotationalHeadingAngleSet;
 bool inertialHeadingAngleSet;


 bool verticalRotationalFlightPathAngleSet;       // All of these are used for the vertical ascent case
 bool verticalInertialFlightPathAngleSet;
 bool verticalRotationalHeadingAngleSet;
 bool verticalInertialHeadingAngleSet;





    // Additional in-class used variables


    int sectionT;                                   // This variable holds the "(section number -1)" for the temperature curve fit
    int powerT;                                     // This variable holds the section corresponding order for the temperature curve fit

    int sectionCD;                                  // This variable holds the "(section number -1)" for the drag coefficient curve fit

    Eigen::VectorXd auxiliaryEquationsVector;               // The vector containing the auxiliary equations xn
    Eigen::VectorXd auxiliaryDerivativesVector;             // The vector containing the auxiliary derivatives un
    Eigen::MatrixXd auxiliaryFunctionsMatrix;               // The matrix containing the auxiliary functions wn,m

};

#endif // AUXILIARY_H
