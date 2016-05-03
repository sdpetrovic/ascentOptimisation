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
 *      160428    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */



#include "allRecurrenceRelations.h"



Eigen::MatrixXd getTaylorCoefficients(const double adiabeticIndex, const double specificGasConstant, const double standardGravitationalParameter, const double rotationalVelocity, const double primeMeridianAngle,
                      const double inertialFrameTime, const double bodyReferenceRadius, const Eigen::MatrixXd temperaturePolyCoefficients, const Eigen::MatrixXd temperatureAltitudeRanges,
                      const Eigen::VectorXd densityPolyCoefficients, const double Thrust, const double specificImpulse,
                      const double referenceArea, const Eigen::MatrixXd dragCoefficientPolyCoefficients, const Eigen::MatrixXd dragCoefficientMachRanges,
        const Eigen::VectorXd& thrustAccelerationsBframe,
        const Eigen::VectorXd& initialEquationsVector,
        const Eigen::VectorXd& initialDerivativesVector,
        const Eigen::MatrixXd& initialFunctionsMatrix,
        const double currentTime,
        const int maxOrder){


    /// Set initial values ///


    Eigen::VectorXd W9IntermediateVector = Eigen::VectorXd::Zero(maxOrder);         // Intermediate vector created to be able to use the basic recurrence relations (should have done this from the start...)


    XMatrix = Eigen::MatrixXd::Zero(initialEquationsVector.size(),maxOrder);
    UMatrix = Eigen::MatrixXd::Zero(initialDerivativesVector.size(),maxOrder);

    for (int i = 1; i < initialEquationsVector.size()+1; i++){

      XMatrix(i,0) = initialEquationsVector(i);
      XMatrix(i,1) = initialDerivativesVector(i);
      UMatrix(i,0) = initialDerivativesVector(i);

    };




    WVector4_1 = Eigen::VectorXd::Zero(maxOrder);     // W4,1
    WVector4_2 = Eigen::VectorXd::Zero(maxOrder);     // W4,2
    WVector4_3 = Eigen::VectorXd::Zero(maxOrder);     // W4,3
    WVector4_4 = Eigen::VectorXd::Zero(maxOrder);     // W4,4
    WVector4_5 = Eigen::VectorXd::Zero(maxOrder);     // W4,5
    WVector4_6 = Eigen::VectorXd::Zero(maxOrder);     // W4,6
    WVector4_7 = Eigen::VectorXd::Zero(maxOrder);     // W4,7
    WVector4_8 = Eigen::VectorXd::Zero(maxOrder);     // W4,8
    WVector4_9 = Eigen::VectorXd::Zero(maxOrder);     // W4,9
    WVector4_10 = Eigen::VectorXd::Zero(maxOrder);     // W4,10

    WVector4_11 = Eigen::VectorXd::Zero(maxOrder);     // W4,11
    WVector4_12 = Eigen::VectorXd::Zero(maxOrder);     // W4,12
    WVector4_13 = Eigen::VectorXd::Zero(maxOrder);     // W4,13
    WVector4_14 = Eigen::VectorXd::Zero(maxOrder);     // W4,14
    WVector4_15 = Eigen::VectorXd::Zero(maxOrder);     // W4,15
    WVector4_16 = Eigen::VectorXd::Zero(maxOrder);     // W4,16
    WVector4_17 = Eigen::VectorXd::Zero(maxOrder);     // W4,17
    WVector4_18 = Eigen::VectorXd::Zero(maxOrder);     // W4,18
    WVector4_19 = Eigen::VectorXd::Zero(maxOrder);     // W4,19
    WVector4_20 = Eigen::VectorXd::Zero(maxOrder);     // W4,20

    WVector4_21 = Eigen::VectorXd::Zero(maxOrder);     // W4,21
    WVector4_22 = Eigen::VectorXd::Zero(maxOrder);     // W4,22
    WVector4_23 = Eigen::VectorXd::Zero(maxOrder);     // W4,23
    WVector4_24 = Eigen::VectorXd::Zero(maxOrder);     // W4,24


    WVector5_1 = Eigen::VectorXd::Zero(maxOrder);     // W5,1
    WVector5_2 = Eigen::VectorXd::Zero(maxOrder);     // W5,2
    WVector5_3 = Eigen::VectorXd::Zero(maxOrder);     // W5,3
    WVector5_4 = Eigen::VectorXd::Zero(maxOrder);     // W5,4
    WVector5_5 = Eigen::VectorXd::Zero(maxOrder);     // W5,5
    WVector5_6 = Eigen::VectorXd::Zero(maxOrder);     // W5,6
    WVector5_7 = Eigen::VectorXd::Zero(maxOrder);     // W5,7
    WVector5_8 = Eigen::VectorXd::Zero(maxOrder);     // W5,8


    WVector6_1 = Eigen::VectorXd::Zero(maxOrder);     // W6,1
    WVector6_2 = Eigen::VectorXd::Zero(maxOrder);     // W6,2
    WVector6_3 = Eigen::VectorXd::Zero(maxOrder);     // W6,3
    WVector6_4 = Eigen::VectorXd::Zero(maxOrder);     // W6,4
    WVector6_5 = Eigen::VectorXd::Zero(maxOrder);     // W6,5
    WVector6_6 = Eigen::VectorXd::Zero(maxOrder);     // W6,6


    WVector8_1 = Eigen::VectorXd::Zero(maxOrder);     // W8,1
    WVector8_2 = Eigen::VectorXd::Zero(maxOrder);     // W8,2
    WVector8_3 = Eigen::VectorXd::Zero(maxOrder);     // W8,3


    WVector9 = Eigen::VectorXd::Zero(maxOrder);     // W9


    WVector12_1 = Eigen::VectorXd::Zero(maxOrder);     // W12,1
    WVector12_2 = Eigen::VectorXd::Zero(maxOrder);     // W12,2
    WVector12_3 = Eigen::VectorXd::Zero(maxOrder);     // W12,3
    WVector12_4 = Eigen::VectorXd::Zero(maxOrder);     // W12,4
    WVector12_5 = Eigen::VectorXd::Zero(maxOrder);     // W12,5
    WVector12_6 = Eigen::VectorXd::Zero(maxOrder);     // W12,6
    WVector12_7 = Eigen::VectorXd::Zero(maxOrder);     // W12,7


    WVector13_1 = Eigen::VectorXd::Zero(maxOrder);     // W13,1
    WVector13_2 = Eigen::VectorXd::Zero(maxOrder);     // W13,2
    WVector13_3 = Eigen::VectorXd::Zero(maxOrder);     // W13,3
    WVector13_4 = Eigen::VectorXd::Zero(maxOrder);     // W13,4
    WVector13_5 = Eigen::VectorXd::Zero(maxOrder);     // W13,5


    WVector14_1 = Eigen::VectorXd::Zero(maxOrder);     // W14,1
    WVector14_2 = Eigen::VectorXd::Zero(maxOrder);     // W14,2
    WVector14_3 = Eigen::VectorXd::Zero(maxOrder);     // W14,3


    WVector15_1 = Eigen::VectorXd::Zero(maxOrder);     // W15,1
    WVector15_2 = Eigen::VectorXd::Zero(maxOrder);     // W15,2


    WVector16 = Eigen::VectorXd::Zero(maxOrder);     // W16


    WVector20 = Eigen::VectorXd::Zero(maxOrder);     // W20


    WVector21_1 = Eigen::VectorXd::Zero(maxOrder);     // W21,1
    WVector21_2 = Eigen::VectorXd::Zero(maxOrder);     // W21,2
    WVector21_3 = Eigen::VectorXd::Zero(maxOrder);     // W21,3


    WVector23_1 = Eigen::VectorXd::Zero(maxOrder);     // W23,1
    WVector23_2 = Eigen::VectorXd::Zero(maxOrder);     // W23,2
    WVector23_3 = Eigen::VectorXd::Zero(maxOrder);     // W23,3
    WVector23_4 = Eigen::VectorXd::Zero(maxOrder);     // W23,4


    WVector24_1 = Eigen::VectorXd::Zero(maxOrder);     // W24,1
    WVector24_2 = Eigen::VectorXd::Zero(maxOrder);     // W24,2
    WVector24_3 = Eigen::VectorXd::Zero(maxOrder);     // W24,3
    WVector24_4 = Eigen::VectorXd::Zero(maxOrder);     // W24,4
    WVector24_5 = Eigen::VectorXd::Zero(maxOrder);     // W24,5
    WVector24_6 = Eigen::VectorXd::Zero(maxOrder);     // W24,6
    WVector24_7 = Eigen::VectorXd::Zero(maxOrder);     // W24,7
    WVector24_8 = Eigen::VectorXd::Zero(maxOrder);     // W24,8
    WVector24_9 = Eigen::VectorXd::Zero(maxOrder);     // W24,9
    WVector24_10 = Eigen::VectorXd::Zero(maxOrder);     // W24,10

    WVector24_11 = Eigen::VectorXd::Zero(maxOrder);     // W24,11
    WVector24_12 = Eigen::VectorXd::Zero(maxOrder);     // W24,12
    WVector24_13 = Eigen::VectorXd::Zero(maxOrder);     // W24,13
    WVector24_14 = Eigen::VectorXd::Zero(maxOrder);     // W24,14
    WVector24_15 = Eigen::VectorXd::Zero(maxOrder);     // W24,15
    WVector24_16 = Eigen::VectorXd::Zero(maxOrder);     // W24,16
    WVector24_17 = Eigen::VectorXd::Zero(maxOrder);     // W24,17


    WVector25_1 = Eigen::VectorXd::Zero(maxOrder);     // W25,1
    WVector25_2 = Eigen::VectorXd::Zero(maxOrder);     // W25,2
    WVector25_3 = Eigen::VectorXd::Zero(maxOrder);     // W25,3


    WVector26_1 = Eigen::VectorXd::Zero(maxOrder);     // W26,1
    WVector26_2 = Eigen::VectorXd::Zero(maxOrder);     // W26,2
    WVector26_3 = Eigen::VectorXd::Zero(maxOrder);     // W26,3


    WVector27_1 = Eigen::VectorXd::Zero(maxOrder);     // W27,1
    WVector27_2 = Eigen::VectorXd::Zero(maxOrder);     // W27,2
    WVector27_3 = Eigen::VectorXd::Zero(maxOrder);     // W27,3
    WVector27_4 = Eigen::VectorXd::Zero(maxOrder);     // W27,4
    WVector27_5 = Eigen::VectorXd::Zero(maxOrder);     // W27,5
    WVector27_6 = Eigen::VectorXd::Zero(maxOrder);     // W27,6


    WVector28 = Eigen::VectorXd::Zero(maxOrder);     // W28


    WVector30_1 = Eigen::VectorXd::Zero(maxOrder);     // W30,1
    WVector30_2 = Eigen::VectorXd::Zero(maxOrder);     // W30,2
    WVector30_3 = Eigen::VectorXd::Zero(maxOrder);     // W30,3
    WVector30_4 = Eigen::VectorXd::Zero(maxOrder);     // W30,4
    WVector30_5 = Eigen::VectorXd::Zero(maxOrder);     // W30,5
    WVector30_6 = Eigen::VectorXd::Zero(maxOrder);     // W30,6
    WVector30_7 = Eigen::VectorXd::Zero(maxOrder);     // W30,7
    WVector30_8 = Eigen::VectorXd::Zero(maxOrder);     // W30,8
    WVector30_9 = Eigen::VectorXd::Zero(maxOrder);     // W30,9


    WVector32_1 = Eigen::VectorXd::Zero(maxOrder);     // W32,1
    WVector32_2 = Eigen::VectorXd::Zero(maxOrder);     // W32,2
    WVector32_3 = Eigen::VectorXd::Zero(maxOrder);     // W32,3
    WVector32_4 = Eigen::VectorXd::Zero(maxOrder);     // W32,4


    WVector33 = Eigen::VectorXd::Zero(maxOrder);     // W33


    WVector34_2 = Eigen::VectorXd::Zero(maxOrder);     // W34,2
    WVector34_3 = Eigen::VectorXd::Zero(maxOrder);     // W34,3
    WVector34_4 = Eigen::VectorXd::Zero(maxOrder);     // W34,4


    WVector35_1 = Eigen::VectorXd::Zero(maxOrder);     // W35,1
    WVector35_2 = Eigen::VectorXd::Zero(maxOrder);     // W35,2
    WVector35_3 = Eigen::VectorXd::Zero(maxOrder);     // W35,3


    WVector36 = Eigen::VectorXd::Zero(maxOrder);     // W36


    WVector37_1 = Eigen::VectorXd::Zero(maxOrder);     // W37,1
    WVector37_2 = Eigen::VectorXd::Zero(maxOrder);     // W37,2
    WVector37_3 = Eigen::VectorXd::Zero(maxOrder);     // W37,3
    WVector37_4 = Eigen::VectorXd::Zero(maxOrder);     // W37,4


    WVector38_1 = Eigen::VectorXd::Zero(maxOrder);     // W38,1
    WVector38_2 = Eigen::VectorXd::Zero(maxOrder);     // W38,2
    WVector38_3 = Eigen::VectorXd::Zero(maxOrder);     // W38,3


    WVector40_1 = Eigen::VectorXd::Zero(maxOrder);     // W40,1
    WVector40_2 = Eigen::VectorXd::Zero(maxOrder);     // W40,2
    WVector40_3 = Eigen::VectorXd::Zero(maxOrder);     // W40,3
    WVector40_4 = Eigen::VectorXd::Zero(maxOrder);     // W40,4


    WVector41_1 = Eigen::VectorXd::Zero(maxOrder);     // W41,1
    WVector41_2 = Eigen::VectorXd::Zero(maxOrder);     // W41,2


    WVector42_1 = Eigen::VectorXd::Zero(maxOrder);     // W42,1
    WVector42_2 = Eigen::VectorXd::Zero(maxOrder);     // W42,2
    WVector42_3 = Eigen::VectorXd::Zero(maxOrder);     // W42,3
    WVector42_4 = Eigen::VectorXd::Zero(maxOrder);     // W42,4
    WVector42_5 = Eigen::VectorXd::Zero(maxOrder);     // W42,5
    WVector42_6 = Eigen::VectorXd::Zero(maxOrder);     // W42,6
    WVector42_7 = Eigen::VectorXd::Zero(maxOrder);     // W42,7
    WVector42_8 = Eigen::VectorXd::Zero(maxOrder);     // W42,8



    WVector43_1 = Eigen::VectorXd::Zero(maxOrder);     // W43,1
    WVector43_2 = Eigen::VectorXd::Zero(maxOrder);     // W43,2


    WVector45_1 = Eigen::VectorXd::Zero(maxOrder);     // W45,1
    WVector45_2 = Eigen::VectorXd::Zero(maxOrder);     // W45,2
    WVector45_3 = Eigen::VectorXd::Zero(maxOrder);     // W45,3
    WVector45_4 = Eigen::VectorXd::Zero(maxOrder);     // W45,4
    WVector45_5 = Eigen::VectorXd::Zero(maxOrder);     // W45,5
    WVector45_6 = Eigen::VectorXd::Zero(maxOrder);     // W45,6
    WVector45_7 = Eigen::VectorXd::Zero(maxOrder);     // W45,7
    WVector45_8 = Eigen::VectorXd::Zero(maxOrder);     // W45,8


    WVector47_1 = Eigen::VectorXd::Zero(maxOrder);     // W47,1
    WVector47_2 = Eigen::VectorXd::Zero(maxOrder);     // W47,2
    WVector47_3 = Eigen::VectorXd::Zero(maxOrder);     // W47,3


    WVector48_1 = Eigen::VectorXd::Zero(maxOrder);     // W48,1
    WVector48_2 = Eigen::VectorXd::Zero(maxOrder);     // W48,2


/// Recurrence computations ///

    // getMultiplicationRecurrenceRelation                  Multiplication
    // getDivisionRecurrenceRelation                        Division
    // getPowerRecurrenceRelation                           Power
    // getExponentialRecurrenceRelation                     Exponential
    // getCosineRecurrenceRelation                          Cosine
    // getSineRecurrenceRelation                            Sine


    for (int k = 1; k< maxOrder;k++){


        // 1
        UMatrix(1,k) = XMatrix(4,k);

        // 2
        UMatrix(2,k) = XMatrix(5,k);

        // 3
        UMatrix(3,k) = XMatrix(6,k);

        // 4

        WVector4_1(k) = getDivisionRecurrenceRelation(XMatrix.row(1),XMatrix.row(9),k); //  XMatrix(1,k)/XMatrix(9,k);
        WVector4_2(k) = thrustAccelerationsBframe(0)-getDivisionRecurrenceRelation(XMatrix.row(27),XMatrix.row(7),k);
        WVector4_3(k) = getCosineRecurrenceRelation((XMatrix.row(10)+XMatrix.row(11)),WVector4_8,k);
        WVector4_4(k) = getSineRecurrenceRelation(XMatrix.row(12),WVector4_6,k);
        WVector4_5(k) = getCosineRecurrenceRelation(XMatrix.row(13),WVector4_9,k);
        WVector4_6(k) = getCosineRecurrenceRelation(XMatrix.row(12),WVector4_4,k);
        WVector4_7(k) = getSineRecurrenceRelation(XMatrix.row(14),XMatrix.row(16),k);
        WVector4_8(k) = getSineRecurrenceRelation((XMatrix.row(10)+XMatrix.row(11)),WVector4_3,k);
        WVector4_9(k) = getSineRecurrenceRelation(XMatrix.row(13),WVector4_5,k);
        WVector4_10(k) = getMultiplicationRecurrenceRelation(WVector4_4,WVector4_5,k);
        WVector4_11(k) = getMultiplicationRecurrenceRelation(WVector4_6,WVector4_7,k);
        WVector4_12(k) = getMultiplicationRecurrenceRelation(WVector4_9,XMatrix.row(16),k);
        WVector4_13(k) = getMultiplicationRecurrenceRelation(WVector4_4,WVector4_9,k);
        WVector4_14(k) = getMultiplicationRecurrenceRelation(WVector4_8,WVector4_5,k);
        WVector4_15(k) = getMultiplicationRecurrenceRelation(WVector4_6,XMatrix.row(16),k);
        WVector4_16(k) = getMultiplicationRecurrenceRelation(WVector4_9,WVector4_7,k);
        WVector4_17(k) = getMultiplicationRecurrenceRelation(WVector4_10,XMatrix.row(16),k);
        WVector4_18(k) = getMultiplicationRecurrenceRelation(WVector4_8,WVector4_12,k);
        WVector4_19(k) = getMultiplicationRecurrenceRelation(WVector4_3,WVector4_13,k);
        WVector4_20(k) = getMultiplicationRecurrenceRelation(WVector4_10,WVector4_7,k);
        WVector4_21(k) = getMultiplicationRecurrenceRelation(WVector4_8,WVector4_16,k);
        WVector4_22(k) = getMultiplicationRecurrenceRelation(WVector4_3,(WVector4_11-WVector4_17),k);
        WVector4_23(k) = getMultiplicationRecurrenceRelation(WVector4_3,(-WVector4_20-WVector4_15),k);
        WVector4_24(k) = getMultiplicationRecurrenceRelation(WVector4_2,(WVector4_22-WVector4_18),k);

        UMatrix(4,k) = -standardGravitationalParameter*WVector4_1(k)+WVector4_24(k)+thrustAccelerationsBframe(1)*(WVector4_19(k)-WVector4_14(k))+thrustAccelerationsBframe(2)*(WVector4_23(k)-WVector4_21(k));

        // 5

        WVector5_1(k) = getDivisionRecurrenceRelation(XMatrix.row(2),XMatrix.row(9),k);
        WVector5_2(k) = getMultiplicationRecurrenceRelation(WVector4_8,(WVector4_11-WVector4_17),k);
        WVector5_3(k) = getMultiplicationRecurrenceRelation(WVector4_3,WVector4_12,k);
        WVector5_4(k) = getMultiplicationRecurrenceRelation(WVector4_8,WVector4_13,k);
        WVector5_5(k) = getMultiplicationRecurrenceRelation(WVector4_3,WVector4_5,k);
        WVector5_6(k) = getMultiplicationRecurrenceRelation(WVector4_8,(-WVector4_20-WVector4_11),k);
        WVector5_7(k) = getMultiplicationRecurrenceRelation(WVector4_3,WVector4_16,k);
        WVector5_8(k) = getMultiplicationRecurrenceRelation(WVector4_2,(WVector5_2+WVector5_3),k);


        UMatrix(5,k) = -standardGravitationalParameter*WVector5_1(k)+WVector5_8(k)+thrustAccelerationsBframe(1)*(WVector5_4(k)+WVector5_5(k))+thrustAccelerationsBframe(2)*(WVector5_6(k)+WVector5_7(k));

        // 6

        WVector6_1(k) = getDivisionRecurrenceRelation(XMatrix.row(3),XMatrix.row(9),k);
        WVector6_2(k) = getMultiplicationRecurrenceRelation(WVector4_5,WVector4_15,k);
        WVector6_3(k) = getMultiplicationRecurrenceRelation(WVector4_6,WVector4_9,k);
        WVector6_4(k) = getMultiplicationRecurrenceRelation(WVector4_5,WVector4_11,k);
        WVector6_5(k) = getMultiplicationRecurrenceRelation(WVector4_4,XMatrix.row(16),k);
        WVector6_6(k) = getMultiplicationRecurrenceRelation(WVector4_2,(WVector6_2-WVector4_11),k);

        UMatrix(6,k) = -standardGravitationalParameter*WVector6_1(k)+WVector6_6(k)-thrustAccelerationsBframe(1)*WVector6_3(k)+thrustAccelerationsBframe(2)*(WVector6_4(k)-WVector6_5(k));

        // 7

//        UMatrix(7,k) =-thrust/(tudat::physical_constants::STANDARD_EARTH_GRAVITATIONAL_ACCELERATION*specificImpulse);
        UMatrix(7,k) = 0;

        // 8

        WVector8_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(1),XMatrix.row(4),k);
        WVector8_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(2),XMatrix.row(5),k);
        WVector8_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(3),XMatrix.row(6),k);

        UMatrix(8,k) = 2*(WVector8_1+WVector8_2+WVector8_3);

        // 10

        UMatrix(10,k) = 0;

        // 11

        UMatrix(11,k) = XMatrix(45,k)-rotationalVelocity;


        // 12

        WVector12_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(20),XMatrix.row(6),k);
        WVector12_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(3),XMatrix.row(25),k);
        WVector12_3(k) = getDivisionRecurrenceRelation(XMatrix.row(3),XMatrix.row(20),k);
        WVector12_4(k) = getMultiplicationRecurrenceRelation(WVector12_3,WVector12_3,k);
        WVector12_5(k) = getPowerRecurrenceRelation((1-WVector12_4),0.5,k);
        WVector12_6(k) = getMultiplicationRecurrenceRelation(XMatrix.row(8),WVector12_5,k);
        WVector12_7(k) = getDivisionRecurrenceRelation((WVector12_1-WVector12_2),WVector12_6,k);


        UMatrix(12,k) = WVector12_7(k);

        // 19

        UMatrix(19,k) = 2*(WVector8_1(k)+WVector8_2);

        // 20

        WVector20_1(k) = getDivisionRecurrenceRelation(XMatrix.row(26),XMatrix.row(20),k);


        UMatrix(20,k) = 0.5*WVector20(k);

        // 35

        UMatrix(35,k) = rotationalVelocity*(WVector35_1(k)-WVector35_3(k));

        // 9

        W9IntermediateVector(k) = getMultiplicationRecurrenceRelation(XMatrix.row(9),UMatrix.row(8),k);

        WVector9(k) = getDivisionRecurrenceRelation(W9IntermediateVector,XMatrix.row(8),k);

        UMatrix(9,k) = 1.5*WVector9(k);

        // 21

        UMatrix(21,k) = 2*(WVector21_1+WVector21_2+WVector21_3);

        // 26

        UMatrix(26,k) = 2*XMatrix(21,k)+2*(WVector26_1+WVector26_2+WVector26_3);

        // 31

        UMatrix(31,k) = UMatrix(20,k);

        // 45

        UMatrix(45,k) = WVector45_8(k);

        // 25

        UMatrix(25,k) = 0.25*WVector25_3(k);

        // 30

        UMatrix(30,k) = WVector30_9(k);

        // 34
        // First it has to be determined which equation has to be used and which values have to be computed

        if (temperatureAltitudeRanges(0,0)<=XMatrix(31,0) && XMatrix(31,0) < temperatureAltitudeRanges(0,1)){

            UMatrix(34,k) = temperaturePolyCoefficients(0,1)*UMatrix(31,k);
        }

        else if (temperatureAltitudeRanges(1,0)<=XMatrix(31,0) && XMatrix(31,0)<temperatureAltitudeRanges(1,1)){

            UMatrix(34,k) = WVector34_2(k);
        }

        else if (temperatureAltitudeRanges(2,0)<=XMatrix(31,0) && XMatrix(31,0)<temperatureAltitudeRanges(2,1)){

            UMatrix(34,k) = WVector34_3(k);

            }

        else if (temperatureAltitudeRanges(3,0)<=XMatrix(31,0) && XMatrix(31,0)<temperatureAltitudeRanges(3,1)){



            UMatrix(34,k) = WVector34,4(k);


        }
        else if (temperatureAltitudeRanges(4,0)<=XMatrix(31,0) && XMatrix(31,0)<temperatureAltitudeRanges(4,1)){

            UMatrix(34,k) = 0;


    };



        // 36

        UMatrix(36,k) = 0.5*WVector36(k);

        // 46

        UMatrix(46,k) = UMatrix(45,k);

        // 47

        UMatrix(47,k) = WVector47_1(k)-WVector47_3(k);

        // 24

        UMatrix(24,k) = WVector24_5(k)-0.5*WVector24_16(k)-WVector24_17(k);

        // 37

        UMatrix(37,k) = WVector37_4(k);

        // 41

        UMatrix(41,k) = WVector41_1(k)+WVector41_2(k);

        // 48

        UMatrix(48,k) = WVector48_1(k)-WVector48_2(k);

        // 13

        // w13
        WVector13_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(24),UMatrix.row(48),k);
        WVector13_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(48),UMatrix.row(24),k);
        WVector13_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(48),XMatrix.row(48),k);
        WVector13_4(k) = getMultiplicationRecurrenceRelation(XMatrix.row(24),XMatrix.row(24),k);
        WVector13_5(k) = getDivisionRecurrenceRelation((WVector13_1-WVector13_2),(WVector13_3+WVector13_4),k);


        UMatrix(13,k) = WVector13_5(k);

        // 38

        UMatrix(38,k) = WVector38_3(k);

        // 40

        UMatrix(40,k) = WVector40_4(k);

        // 42

        UMatrix(42,k) = WVector42_7(k) - WVector42_8(k);

        // 43

        UMatrix(43,k) = WVector43_1(k) + WVector43_2(k);

        // 15

        WVector15_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(35),UMatrix.row(35),k);
        WVector15_2(k) = getDivisionRecurrenceRelation((2*WVector15_1+UMatrix.row(21)-2*UMatrix.row(43)),XMatrix.row(15),k);

        UMatrix(15,k) = 0.5*WVector15_2(k);

        // 23

        UMatrix(23,k) = WVector23_4(k);

        // 32

        UMatrix(32,k) = WVector32_4(k);

        // 14

        WVector14_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(23),XMatrix.row(23),k);
        WVector14_2(k) = getPowerRecurrenceRelation((1-WVector14_1),0.5,k);
        WVector14_3(k) = getDivisionRecurrenceRelation(UMatrix.row(23),WVector14_2,k);


        UMatrix(14,k) = WVector14_3(k);

        // 29

        for (int i=0; i < 5+1; i++){

            if (dragCoefficientMachRanges(i,0) <= XMatrix(32,0) && XMatrix(32,0) < dragCoefficientMachRanges(i,1)){

                int sectionCD = i;
            }

        };

        UMatrix(29,k) = dragCoefficientPolyCoefficients(sectionCD,1)*UMatrix(32,k);


        // 16

        WVector16_1(k) = getMultiplicationRecurrenceRelation(WVector4_7,UMatrix.row(14),k);

        UMatrix(16,k) = -WVector16(k);

        // 27

        UMatrix(27,k) = 0.5*referenceArea*WVector27_6(k);



        /// Compute all auxiliary Equation coefficients from the Derivative coefficients ///

        for (int i = 1; i<i < initialEquationsVector.size()+1; i++){

            XMatrix(i,k+1) = UMatrix(i,k)/(k+1);
        }


//////////////////////





        // w21
        auxiliaryFunctionsMatrix(21,1) = auxiliaryEquationsVector(4)*auxiliaryDerivativesVector(4);
        auxiliaryFunctionsMatrix(21,2) = auxiliaryEquationsVector(5)*auxiliaryDerivativesVector(5);
        auxiliaryFunctionsMatrix(21,3) = auxiliaryEquationsVector(6)*auxiliaryDerivativesVector(6);








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
        auxiliaryFunctionsMatrix(24,1) = auxiliaryEquationsVector(6)*auxiliaryDerivativesVector(20);
        auxiliaryFunctionsMatrix(24,2) = auxiliaryEquationsVector(20)*auxiliaryDerivativesVector(6);
        auxiliaryFunctionsMatrix(24,3) = auxiliaryEquationsVector(25)*auxiliaryEquationsVector(6);
        auxiliaryFunctionsMatrix(24,4) = auxiliaryEquationsVector(3)*auxiliaryDerivativesVector(25);

        // Avoiding singularities
        if (auxiliaryFunctionsMatrix(12,6) == 0){

            auxiliaryFunctionsMatrix(24,5) = 0;
        }
        else {
        auxiliaryFunctionsMatrix(24,5) = (auxiliaryFunctionsMatrix(24,1)+auxiliaryFunctionsMatrix(24,2)-auxiliaryFunctionsMatrix(24,3)-auxiliaryFunctionsMatrix(24,4))/(auxiliaryFunctionsMatrix(12,6));

        };

        auxiliaryFunctionsMatrix(24,6) = auxiliaryEquationsVector(3)*auxiliaryFunctionsMatrix(24,4);
        auxiliaryFunctionsMatrix(24,7) = pow(auxiliaryEquationsVector(20),3);
        auxiliaryFunctionsMatrix(24,8) = auxiliaryFunctionsMatrix(8,3)/auxiliaryEquationsVector(8);
        auxiliaryFunctionsMatrix(24,9) = auxiliaryFunctionsMatrix(24,6)/auxiliaryFunctionsMatrix(24,7);
        auxiliaryFunctionsMatrix(24,10) = (2*auxiliaryFunctionsMatrix(24,9)-2*auxiliaryFunctionsMatrix(24,8))*(auxiliaryFunctionsMatrix(12,1)-auxiliaryFunctionsMatrix(12,2));
        auxiliaryFunctionsMatrix(24,11) = auxiliaryDerivativesVector(8)*(auxiliaryFunctionsMatrix(12,1)*auxiliaryFunctionsMatrix(12,2));
        auxiliaryFunctionsMatrix(24,12) = pow((1-auxiliaryFunctionsMatrix(12,4)),1.5);
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
        auxiliaryFunctionsMatrix(25,3) = (2*auxiliaryFunctionsMatrix(25,1)-auxiliaryFunctionsMatrix(25,2))/(auxiliaryEquationsVector(9));





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
        auxiliaryFunctionsMatrix(30,1) = pow(auxiliaryEquationsVector(31),9);
        auxiliaryFunctionsMatrix(30,2) = pow(auxiliaryEquationsVector(31),8);
        auxiliaryFunctionsMatrix(30,3) = pow(auxiliaryEquationsVector(31),7);
        auxiliaryFunctionsMatrix(30,4) = pow(auxiliaryEquationsVector(31),6);
        auxiliaryFunctionsMatrix(30,5) = pow(auxiliaryEquationsVector(31),5);
        auxiliaryFunctionsMatrix(30,6) = pow(auxiliaryEquationsVector(31),4);
        auxiliaryFunctionsMatrix(30,7) = pow(auxiliaryEquationsVector(31),3);
        auxiliaryFunctionsMatrix(30,8) = pow(auxiliaryEquationsVector(31),2);

        for (int i=2; i<10+1; i++){

            if (i==2){

                auxiliaryFunctionsMatrix(30,9) = auxiliaryDerivativesVector(31)*(2*densityPolyCoefficients(2)*auxiliaryEquationsVector(31)+densityPolyCoefficients(1));
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

          auxiliaryFunctionsMatrix(34,2) = auxiliaryDerivativesVector(31)*(3*temperaturePolyCoefficients(3,2)*auxiliaryFunctionsMatrix(30,8)+2*temperaturePolyCoefficients(2,2)*auxiliaryEquationsVector(31)+temperaturePolyCoefficients(1,2));
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

                    auxiliaryFunctionsMatrix(34,4) = auxiliaryDerivativesVector(31)*(2*temperaturePolyCoefficients(2,4)*auxiliaryEquationsVector(31)+temperaturePolyCoefficients(1,4));
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
        auxiliaryFunctionsMatrix(38,2) = sqrt(1-auxiliaryFunctionsMatrix(38,1));

        // Avoiding singularities
        if (auxiliaryFunctionsMatrix(38,2) == 0){

            auxiliaryFunctionsMatrix(38,3) = 0;
        }
        else {
        auxiliaryFunctionsMatrix(38,3) = auxiliaryDerivativesVector(37)/auxiliaryFunctionsMatrix(38,2);
    };





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
        auxiliaryFunctionsMatrix(42,2) = cos(auxiliaryEquationsVector(40));
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
    }




    /// Set return matrix ///

    stateTaylorCoefficients = Eigen::MatrixXd::Zero(8,maxOrder);


    return stateTaylorCoefficients;
}




