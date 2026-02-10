/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cmath>
#include "grid_convert.h"

double
lamrot_to_lam(double phirot, double lamrot, double polphi, double pollam, double polgam)
{
  /*
    Name of the original Fortran function: PHTOPHS

    This function converts lambda from one rotated system to lambda in another
    system. If the optional argument polgam is present, the other system can
    also be a rotated one, where polgam is the angle between the two north
    poles. If polgam is not present, the other system is the real geographical
    system.

    phirot : latitude in the rotated system
    lamrot : longitude in the rotated system (E>0)
    polphi : latitude of the rotated north pole
    pollam : longitude of the rotated north pole

    result : longitude in the geographical system
  */
  double zarg1, zarg2;

  double zsinpol = std::sin(DEG2RAD * polphi);
  double zcospol = std::cos(DEG2RAD * polphi);

  double zlampol = DEG2RAD * pollam;
  double zphirot = DEG2RAD * phirot;
  if (lamrot > 180.0) lamrot -= 360.0;
  double zlamrot = DEG2RAD * lamrot;

  if (polgam < 0) polgam += 360.0;
  if (polgam > 0)
  {
    double zgam = DEG2RAD * polgam;
    zarg1 = std::sin(zlampol)
                * (-zsinpol * std::cos(zphirot) * (std::cos(zlamrot) * std::cos(zgam) - std::sin(zlamrot) * std::sin(zgam))
                   + zcospol * std::sin(zphirot))
            - std::cos(zlampol) * std::cos(zphirot) * (std::sin(zlamrot) * std::cos(zgam) + std::cos(zlamrot) * std::sin(zgam));

    zarg2 = std::cos(zlampol)
                * (-zsinpol * std::cos(zphirot) * (std::cos(zlamrot) * std::cos(zgam) - std::sin(zlamrot) * std::sin(zgam))
                   + zcospol * std::sin(zphirot))
            + std::sin(zlampol) * std::cos(zphirot) * (std::sin(zlamrot) * std::cos(zgam) + std::cos(zlamrot) * std::sin(zgam));
  }
  else
  {
    zarg1 = std::sin(zlampol) * (-zsinpol * std::cos(zlamrot) * std::cos(zphirot) + zcospol * std::sin(zphirot))
            - std::cos(zlampol) * std::sin(zlamrot) * std::cos(zphirot);
    zarg2 = std::cos(zlampol) * (-zsinpol * std::cos(zlamrot) * std::cos(zphirot) + zcospol * std::sin(zphirot))
            + std::sin(zlampol) * std::sin(zlamrot) * std::cos(zphirot);
  }

  double result = 0;
  if (std::fabs(zarg2) > 0) result = RAD2DEG * std::atan2(zarg1, zarg2);
  if (std::fabs(result) < 9.e-14) result = 0;

  return result;
}

double
phirot_to_phi(double phirot, double lamrot, double polphi, double polgam)
{
  /*
    Name of the original Fortran function: PHSTOPH

    This function converts phi from one rotated system to phi in another
    system. If the optional argument polgam is present, the other system
    can also be a rotated one, where polgam is the angle between the two
    north poles.
    If polgam is not present, the other system is the real geographical
    system.

    phirot : latitude in the rotated system
    lamrot : longitude in the rotated system (E>0)
    polphi : latitude of the rotated north pole
    polgam : angle between the north poles of the systems

    result : latitude in the geographical system
  */
  double zarg;

  double zsinpol = std::sin(DEG2RAD * polphi);
  double zcospol = std::cos(DEG2RAD * polphi);

  double zphirot = DEG2RAD * phirot;
  if (lamrot > 180.0) lamrot -= 360.0;
  double zlamrot = DEG2RAD * lamrot;

  if (polgam < 0) polgam += 360.0;
  if (polgam > 0)
  {
    double zgam = DEG2RAD * polgam;
    zarg = zsinpol * std::sin(zphirot)
           + zcospol * std::cos(zphirot) * (std::cos(zlamrot) * std::cos(zgam) - std::sin(zgam) * std::sin(zlamrot));
  }
  else
    zarg = zcospol * std::cos(zphirot) * std::cos(zlamrot) + zsinpol * std::sin(zphirot);

  return RAD2DEG * std::asin(zarg);
}

static double
lam_to_lamrot(double phi, double rla, double polphi, double pollam)
{
  /*
    Name of the original Fortran function: RLSTORL

    Umrechnung von rla (geo. System) auf rlas (rot. System)

    phi    : Breite im geographischen System (N>0)
    rla    : Laenge im geographischen System (E>0)
    polphi : Geographische Breite des Nordpols des rot. Systems
    pollam : Geographische Laenge des Nordpols des rot. Systems

    result : Rotierte Laenge
  */
  double zsinpol = std::sin(DEG2RAD * polphi);
  double zcospol = std::cos(DEG2RAD * polphi);
  double zlampol = DEG2RAD * pollam;

  if (rla > 180.0) rla -= 360.0;

  double zrla = DEG2RAD * rla;
  double zphi = DEG2RAD * phi;

  double zarg1 = -std::sin(zrla - zlampol) * std::cos(zphi);
  double zarg2 = -zsinpol * std::cos(zphi) * std::cos(zrla - zlampol) + zcospol * std::sin(zphi);

  if (std::fabs(zarg2) < 1.0e-20) zarg2 = 1.0e-20;

  return RAD2DEG * std::atan2(zarg1, zarg2);
}

#ifdef TEST_GRID_ROT
static double
phi_to_phirot(double phi, double lam, double polphi, double pollam)
{
  /*
    Name of the original Fortran function: PHTOPHS

    Umrechnung von phi (geo. System) auf phis (rot. System)

    phi    : Breite im geographischen System (N>0)
    lam    : Laenge im geographischen System (E>0)
    polphi : Geographische Breite des Nordpols des rot. Systems
    pollam : Geographische Laenge des Nordpols des rot. Systems

    result : Rotierte Breite
  */
  double zsinpol = std::sin(DEG2RAD * polphi);
  double zcospol = std::cos(DEG2RAD * polphi);
  double zlampol = DEG2RAD * pollam;

  double zphi = DEG2RAD * phi;
  if (lam > 180.0) lam -= 360.0;
  double zlam = DEG2RAD * lam;

  double zarg = zcospol * std::cos(zphi) * std::cos(zlam - zlampol) + zsinpol * std::sin(zphi);

  return RAD2DEG * std::asin(zarg);
}
#endif

void
usvs_to_uv(double us, double vs, double phi, double rla, double polphi, double pollam, double *u, double *v)
{
  /*
    Umrechnen der windkomponenten us, vs im rotierten sphaerischen
    system in die windkomponenten u, v, im geographischen system

    us     : 'zonaler wind im rotierten system
    vs     : 'merid. wind im rotierten  system
    phi    : Breite im geographischen system (N>0)
    rla    : Laenge im geographischen system (E>0)
    polphi : Geographische breite des Nordpols des rot. Systems
    pollam : Geographische laenge des Nordpols des rot. Systems

    u      : zonaler wind im geographischen system
    v      : merid. wind im geographischen system
  */
  /* umrechnung von grad in bogenmass */
  double zpolphi = polphi * DEG2RAD;
  double zpollam = pollam * DEG2RAD;
  // Added by Uwe Schulzweida (17/11/2017)
  if (pollam < 0 && rla < pollam) rla += 360.0;
  // if ( pollam < 0 && rla < 0 ) rla += 360.0;
  double zrla = rla * DEG2RAD;
  double pollamd = pollam;
  if (pollamd < 0.0) pollamd += 360.0;

  // laenge im rotierten system berechnen
  double zrlas = lam_to_lamrot(phi, rla, polphi, pollam) * DEG2RAD;

  // winkel zbeta berechen (schnittwinkel der breitenkreise)
  double zarg = -std::sin(zpolphi) * std::sin(zrla - zpollam) * std::sin(zrlas) - std::cos(zrla - zpollam) * std::cos(zrlas);
  if (zarg > 1.0) zarg = 1.0;
  if (zarg < -1.0) zarg = -1.0;
  /*
  zbeta = std::acos(zarg);
  zbeta = sign(zbeta, -(rla - (pollamd-180.0)));
  */
  double zbeta = std::fabs(std::acos(zarg));
  // if ( -(rla - (pollamd-180.0)) < 0 ) zbeta = -zbeta;
  if ((-(rla - (pollamd - 180.0)) < 0) && (-(rla - (pollamd - 180.0)) >= -180)) zbeta = -zbeta;

  // us - wind transformieren
  *u = us * std::cos(zbeta) - vs * std::sin(zbeta);

  // vs - wind transformieren
  *v = us * std::sin(zbeta) + vs * std::cos(zbeta);
}

#ifdef TEST_GRID_ROT
int
main(void)
{
  double angle = 0.0;
  double polphi = 32.5;
  double pollam = -170.0;

  double y0 = 0.0;

  for (int i = 0; i < 10; ++i)
  {
    double x0 = i * 20.0;
    printf("rot in: %g %g\n", x0, y0);

    double x1 = lamrot_to_lam(y0, x0, polphi, pollam, angle);
    double y1 = phirot_to_phi(y0, x0, polphi, angle);
    printf("geo: %g %g\n", x1, y1);

    double x2 = lam_to_lamrot(y1, x1, polphi, pollam);
    double double y2 = phi_to_phirot(y1, x1, polphi, pollam);
    printf("rot out:%g %g\n", x2, y2);
  }

  double x1, x2;
  usvs_to_uv(30.0, 20.0, 30.0, 0.0, polphi, pollam, &x1, &x2);
  printf("usvs_to_uv: %g %g %g %g\n", polphi, pollam, x1, x2);
  printf("usvs_to_uv: 32.5 -170 26.3124 24.6507 <-- reference\n");

  return 0;
}
#endif
