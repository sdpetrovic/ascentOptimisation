#ifndef AIRTEMPERATURE_H
#define AIRTEMPERATURE_H


#include <iostream>
#include <Eigen/Core>
#include <cmath>


namespace air_temperature
{

/// Air Temperature function ///
/// \brief airTemperature   Computes the current air temperature in [K]
/// \param temperaturePolyCoefficients  The polynomial coefficients for the temperature curve
/// \param temperatureAltitudeRanges    The altitudes defining each section of the temperature curve [km]
/// \param altitude r-R_MOLA [km]
/// \return
///
const double airTemperature(const Eigen::MatrixXd temperaturePolyCoefficients, const Eigen::MatrixXd temperatureAltitudeRanges, const double altitude);




} // end namespace air_temperature

#endif // AIRTEMPERATURE_H
