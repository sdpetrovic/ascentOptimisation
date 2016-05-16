#ifndef OTHERREQUIREDFUNCTIONS_H
#define OTHERREQUIREDFUNCTIONS_H


#include <iostream>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <string>
#include <cmath>

#include <tudatApplications/thesisProject/linearAlgebraTypesUpdated.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>



/// Some required functions ///


/// deg2rad ///
/// \brief deg2rad is a function to convert degrees to radians
/// \param deg
/// \return
///

const double deg2rad(const double deg);
//{

//    const double rad = deg*tudat::mathematical_constants::LONG_PI/180;

//    return rad;
//}

/// rad2deg ///
/// \brief rad2deg is a function to convert radians to degrees
/// \param rad
/// \return
///

const double rad2deg(const double rad);
//{

//    const double deg = rad*180/tudat::mathematical_constants::LONG_PI;

//    return deg;
//}


/// B-P frame transformations ///



//! Get transformation quaternion from the Body (B) to the Propulsion (P) frame.
Eigen::Quaterniond getBodyToPropulsionFrameTransformationQuaternion(
        const double thrustAzimuth, const double thrustElevation );
//{
//    // Compute transformation quaternion.
//    // Note the sign change, because how angleAxisd is defined.
//    Eigen::AngleAxisd RotationAroundZaxis = Eigen::AngleAxisd(
//                -1.0 * thrustAzimuth, Eigen::Vector3d::UnitZ( ) );
//    Eigen::AngleAxisd RotationAroundYaxis = Eigen::AngleAxisd(
//                -1.0 * thrustElevation,
//                Eigen::Vector3d::UnitY( ) );
//    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
//                ( RotationAroundYaxis * RotationAroundZaxis ) );

//    // Return transformation quaternion.
//    return frameTransformationQuaternion;
//}


//! Get transformation matrix from the Body (B) to the Propulsion (P) frame.
Eigen::Matrix3d getBodyToPropulsionFrameTransformationMatrix(
    const double thrustAzimuth, const double thrustElevation );
//{
//    return getBodyToPropulsionFrameTransformationQuaternion(
//            thrustAzimuth, thrustElevation ).toRotationMatrix( );
//}

//! Get transformation matrix from the Propulsion (P) to the Body (B) frame.
Eigen::Matrix3d getPropulsionToBodyFrameTransformationMatrix(
    const double thrustAzimuth, const double thrustElevation );
//{
//    return getBodyToPropulsionFrameTransformationMatrix(
//            thrustAzimuth, thrustElevation ).transpose( );
//}

//! Get transformation quaternion from the Propulsion (P) to the Body (B) frame.
Eigen::Quaterniond getPropulsionToBodyFrameTransformationQuaternion(
        const double thrustAzimuth, const double thrustElevation );
//{
//    return getBodyToPropulsionFrameTransformationQuaternion(
//            thrustAzimuth, thrustElevation ).inverse( );
//}


#endif // OTHERREQUIREDFUNCTIONS_H
