#ifndef CONST_H
#define CONST_H

#include <numbers>

constexpr double MIN_PS = 20000.;
constexpr double MAX_PS = 120000.;

constexpr double MIN_FIS = -100000.;
constexpr double MAX_FIS = 100000.;

constexpr double MIN_T = 150.;
constexpr double MAX_T = 400.;

constexpr double MIN_Q = 0.0;
constexpr double MAX_Q = 0.1;

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288 /* pi */
#endif

#ifndef RAD_CONVERT
#define RAD_CONVERT
constexpr double RAD2DEG = 180.0 / std::numbers::pi;  // conversion from radians to degrees
constexpr double DEG2RAD = std::numbers::pi / 180.0;  // conversion from degrees to radians
#endif

#endif /* CONST_H */
