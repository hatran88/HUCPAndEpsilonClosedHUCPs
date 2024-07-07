
#include <cmath>

#include "SetNumberPointFloat.h"


float Precision(float f, int places)
{
    float n = std::pow(10.0f, places);
    return std::round(f * n) / n;
}

