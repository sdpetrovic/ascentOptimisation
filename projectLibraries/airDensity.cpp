#include "airDensity.h"

namespace air_density
{

/// Air density function ///

const double airDensity(const Eigen::MatrixXd densityPolyCoefficients, const double altitude){

double currentDensity_ = 0;  // Define the current log density and setting it to 0

for (int i = 0; i<11; i++){ // Compute the complete polynomial

    currentDensity_ += densityPolyCoefficients(i)*pow(altitude,i);

}

const double currentDensity = exp(currentDensity_); // Taking the exponential of the polynomial

return currentDensity;

} // end function airDensity

} // end air_density namespace

