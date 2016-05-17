#include "airTemperature.h"

namespace air_temperature
{

const double airTemperature(const Eigen::MatrixXd temperaturePolyCoefficients, const Eigen::MatrixXd temperatureAltitudeRanges, const double altitude){

int section = 0;    // Define the section and setting the default to 0

for (int i = 0; i<5; i++){      // Determine in which section of the curve the MAV is

    if (temperatureAltitudeRanges(i,0)<=altitude && altitude < temperatureAltitudeRanges(i,1)){
        section = i;
    };
};

int order = 0; // Define the order and setting the default to 0

// Determine the order
if (section == 0){
    order = 1;

}
else if (section == 1){
    order = 2;
}
else if (section == 2){
    order = 6;
}
else if(section == 3){
    order = 8;
}
else if (section == 4){
    order = 0;
};


// Compute the current temperature

double currentTemperature = 0; // Set current temperature to 0;

for (int j = 0; j<order+1; j++){

    currentTemperature += temperaturePolyCoefficients(section,j)*pow(altitude,j);
};

return currentTemperature;

} // end of function airTemperature



} // end namespace air_temperature
