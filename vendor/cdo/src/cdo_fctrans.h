#ifndef CDO_FCTRANS_H
#define CDO_FCTRANS_H

void fc2gp(const double *fc, double *gp, long nlat, long nlon, long nlev, long nfc);
void gp2fc(const double *gp, double *fc, long nlat, long nlon, long nlev, long nfc);

#endif
