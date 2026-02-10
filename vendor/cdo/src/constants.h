#ifndef _CONSTANTS_H
#define _CONSTANTS_H

// Thermodynamical constants adopted from ECMWF IFS-Code

constexpr double C_RKBOL = 1.380658e-23;   // Boltzmann constant in J/K
constexpr double C_RNAVO = 6.0221367e+23;  // Avogadro constant in 1/mol
constexpr double C_RMD = 28.9644;          // molecular weight of dry air
constexpr double C_RMV = 18.0153;          // molecular weight of water vapor
constexpr double C_R = C_RKBOL * C_RNAVO;
constexpr double C_RV = 1000. * C_R / C_RMV;

constexpr double C_EARTH_RD = 1000. * C_R / C_RMD;
constexpr double C_EARTH_RADIUS = 6371000.0;  // radius of the Earth in m
constexpr double C_EARTH_GRAV = 9.80665;

#define C_RG (1.0 / PlanetGrav)

constexpr double C_RCPV = 4.0 * C_RV;
#define C_RETV (C_RV / PlanetRD - 1.0)
constexpr double C_RCW = 4218.;        // specific water heat capacity ??
constexpr double C_RCS = 2106.;        // specific ice heat capacity ??
constexpr double C_RTT = 273.16;       // melting temperature of ice/snow
constexpr double C_RLVTT = 2.5008e+6;  // latent heat for vaporisation in J/kg
constexpr double C_RLSTT = 2.8345e+6;  // latent heat for sublimation in J/kg
constexpr double C_RESTT = 611.14;
#define C_RCPD (3.5 * PlanetRD);

constexpr double C_TIMES_RHOH2O = -333700000.0;

constexpr double PlanetRDDefault = C_EARTH_RD;
constexpr double PlanetRadiusDefault = C_EARTH_RADIUS;
constexpr double PlanetGravDefault = C_EARTH_GRAV;

extern double PlanetRD;
extern double PlanetRadius;
extern double PlanetGrav;

#endif /* _CONSTANTS_H */
