#ifndef TRANSFORM_H
#define TRANSFORM_H

void after_legini_full(long ntr, long nlat, double *poli, double *pold, double *pdev, double *pol2, double *pol3, double *coslat);

// FFT
int fft_set(double *trigs, long *ifax, long n);
void fc2gp(const double *trig, const long *ifax, const double *fc, double *gp, long nlat, long nlon, long nlev, long nfc);
void gp2fc(const double *trig, const long *ifax, const double *gp, double *fc, long nlat, long nlon, long nlev, long nfc);

// Convert Spectral Array to new resolution
void sp2sp(const double *arrayIn, long truncIn, double *arrayOut, long truncOut);
void sp2fc(const double *sa, double *fa, const double *poli, long nlev, long nlat, long nfc, long nt);
void fc2sp(const double *fa, double *sa, const double *poli, long klev, long nlat, long nfc, long nt);

// Physc
void dv2ps(const double *div, double *pot, long nlev, long ntr);
void dv2uv(const double *d, const double *o, double *u, double *v, const double *f, const double *g, long nt, long nsp, long nlev);
void scaluv(double *fu, const double *rclat, long nlat, long lot);
void uv2dv(const double *fu, const double *fv, double *sd, double *sv, const double *pol2, const double *pol3, long klev, long nlat,
           long nt);
void geninx(long ntr, double *f, double *g);

#endif
