#ifndef AIRDENSITY_H
#define AIRDENSITY_H


#include <iostream>
#include <Eigen/Core>
#include <cmath>

namespace air_density
{

/// Air density function ///
/// \brief airDensity       Computes the current air density from the fitted density curve
/// \param densityPolyCoefficients  The polynomial coefficients for the log fit
/// \param altitude     The current altitude
/// \return
///

const double airDensity(const Eigen::MatrixXd densityPolyCoefficients, const double altitude);

} // end air_density namespace

#endif // AIRDENSITY_H
