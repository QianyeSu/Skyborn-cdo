/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_output.h"
#include "cdo_fctrans.h"
#include "specspace.h"
#include <mpim_grid.h>

void
grid2spec(const SP_Transformation &spTrans, int gridIDin, Varray<double> const &arrayIn, int gridIDout, Varray<double> &arrayOut)
{
  auto const &fcTrans = spTrans.fcTrans;
  long nlev = 1;
  long ntr = gridInqTrunc(gridIDout);
  long nlon = gridInqXsize(gridIDin);
  long nlat = gridInqYsize(gridIDin);
  long waves = ntr + 1;
  long nfc = waves * 2;

  Varray<double> fpwork(nlat * nfc * nlev);

  if (fcTrans.use_fftw)
    gp2fc(arrayIn.data(), fpwork.data(), nlat, nlon, nlev, nfc);
  else
    gp2fc(fcTrans.vtrig.data(), fcTrans.ifax, arrayIn.data(), fpwork.data(), nlat, nlon, nlev, nfc);

  fc2sp(fpwork.data(), arrayOut.data(), spTrans.pold.data(), nlev, nlat, nfc, ntr);
}

void
spec2grid(const SP_Transformation &spTrans, int gridIDin, Varray<double> const &arrayIn, int gridIDout, Varray<double> &arrayOut)
{
  auto const &fcTrans = spTrans.fcTrans;
  long nlev = 1;
  long ntr = gridInqTrunc(gridIDin);
  long nlon = gridInqXsize(gridIDout);
  long nlat = gridInqYsize(gridIDout);
  long waves = ntr + 1;
  long nfc = waves * 2;

  Varray<double> fpwork(nlat * nfc * nlev);

  sp2fc(arrayIn.data(), fpwork.data(), spTrans.poli.data(), nlev, nlat, nfc, ntr);

  if (fcTrans.use_fftw)
    fc2gp(fpwork.data(), arrayOut.data(), nlat, nlon, nlev, nfc);
  else
    fc2gp(fcTrans.vtrig.data(), fcTrans.ifax, fpwork.data(), arrayOut.data(), nlat, nlon, nlev, nfc);
}

void
four2spec(const SP_Transformation &spTrans, int gridIDin, Varray<double> const &arrayIn, int gridIDout, Varray<double> &arrayOut)
{
  (void) gridIDin;
  long nlev = 1;
  long ntr = gridInqTrunc(gridIDout);
  long nlat = spTrans.nlat;
  long waves = ntr + 1;
  long nfc = waves * 2;

  fc2sp(arrayIn.data(), arrayOut.data(), spTrans.pold.data(), nlev, nlat, nfc, ntr);
}

void
spec2four(const SP_Transformation &spTrans, int gridIDin, Varray<double> const &arrayIn, int gridIDout, Varray<double> &arrayOut)
{
  long nlev = 1;
  long ntr = gridInqTrunc(gridIDin);
  long nfc = gridInqSize(gridIDout);
  long nlat = nfc_to_nlat(nfc, ntr);
  long waves = ntr + 1;
  nfc = waves * 2;

  sp2fc(arrayIn.data(), arrayOut.data(), spTrans.poli.data(), nlev, nlat, nfc, ntr);
}

void
four2grid(const FC_Transformation &fcTrans, int gridIDin, Varray<double> const &arrayIn, int gridIDout, Varray<double> &arrayOut)
{
  long nlev = 1;
  long ntr = gridInqTrunc(gridIDin);
  long nlon = gridInqXsize(gridIDout);
  long nlat = gridInqYsize(gridIDout);
  long waves = ntr + 1;
  long nfc = waves * 2;

  if (fcTrans.use_fftw)
    fc2gp(arrayIn.data(), arrayOut.data(), nlat, nlon, nlev, nfc);
  else
    fc2gp(fcTrans.vtrig.data(), fcTrans.ifax, arrayIn.data(), arrayOut.data(), nlat, nlon, nlev, nfc);
}

void
grid2four(const FC_Transformation &fcTrans, int gridIDin, Varray<double> const &arrayIn, int gridIDout, Varray<double> &arrayOut)
{
  long nlev = 1;
  long ntr = gridInqTrunc(gridIDout);
  long nlon = gridInqXsize(gridIDin);
  long nlat = gridInqYsize(gridIDin);
  long waves = ntr + 1;
  long nfc = waves * 2;

  if (fcTrans.use_fftw)
    gp2fc(arrayIn.data(), arrayOut.data(), nlat, nlon, nlev, nfc);
  else
    gp2fc(fcTrans.vtrig.data(), fcTrans.ifax, arrayIn.data(), arrayOut.data(), nlat, nlon, nlev, nfc);
}

void
spec2spec(int gridIDin, Varray<double> const &arrayIn, int gridIDout, Varray<double> &arrayOut)
{
  long ntrIn = gridInqTrunc(gridIDin);
  long ntrOut = gridInqTrunc(gridIDout);

  sp2sp(arrayIn.data(), ntrIn, arrayOut.data(), ntrOut);
}

void
speccut(int gridIDin, Varray<double> const &arrayIn, Varray<double> &arrayOut, Varray<int> const &waves)
{
  long ntr = gridInqTrunc(gridIDin);

  spcut(arrayIn.data(), arrayOut.data(), ntr, waves.data());
}

void
trans_uv2dv(const SP_Transformation &spTrans, long nlev, int gridID1, Varray<double> const &gu, Varray<double> const &gv,
            int gridID2, Varray<double> &sd, Varray<double> &svo)
{
  if (gridInqType(gridID1) != GRID_GAUSSIAN)
    cdo_abort("unexpected grid1 type: %s instead of Gaussian", gridNamePtr(gridInqType(gridID1)));

  if (gridInqType(gridID2) != GRID_SPECTRAL)
    cdo_abort("unexpected grid2 type: %s instead of spectral", gridNamePtr(gridInqType(gridID2)));

  auto const &fcTrans = spTrans.fcTrans;
  long ntr = gridInqTrunc(gridID2);
  long nlon = gridInqXsize(gridID1);
  long nlat = gridInqYsize(gridID1);
  long waves = ntr + 1;
  long nfc = waves * 2;

  Varray<double> fpwork1(nlat * nfc * nlev);
  Varray<double> fpwork2(nlat * nfc * nlev);

  if (fcTrans.use_fftw)
  {
    gp2fc(gu.data(), fpwork1.data(), nlat, nlon, nlev, nfc);
    gp2fc(gv.data(), fpwork2.data(), nlat, nlon, nlev, nfc);
  }
  else
  {
    gp2fc(fcTrans.vtrig.data(), fcTrans.ifax, gu.data(), fpwork1.data(), nlat, nlon, nlev, nfc);
    gp2fc(fcTrans.vtrig.data(), fcTrans.ifax, gv.data(), fpwork2.data(), nlat, nlon, nlev, nfc);
  }

  scaluv(fpwork1.data(), spTrans.coslat.data(), nlat, nfc * nlev);
  scaluv(fpwork2.data(), spTrans.coslat.data(), nlat, nfc * nlev);

  uv2dv(fpwork1.data(), fpwork2.data(), sd.data(), svo.data(), spTrans.pol2.data(), spTrans.pol3.data(), nlev, nlat, ntr);
}

void
trans_dv2uv(const SP_Transformation &spTrans, const DV_Transformation &dvTrans, long nlev, int gridID1, Varray<double> const &sd,
            Varray<double> const &svo, int gridID2, Varray<double> &gu, Varray<double> &gv)
{
  if (gridInqType(gridID1) != GRID_SPECTRAL)
    cdo_warning("unexpected grid1 type: %s instead of spectral", gridNamePtr(gridInqType(gridID1)));
  if (gridInqType(gridID2) != GRID_GAUSSIAN)
    cdo_warning("unexpected grid2 type: %s instead of Gaussian", gridNamePtr(gridInqType(gridID2)));

  auto const &fcTrans = spTrans.fcTrans;
  long ntr = gridInqTrunc(gridID1);
  long nlon = gridInqXsize(gridID2);
  long nlat = gridInqYsize(gridID2);
  long waves = ntr + 1;
  long nfc = waves * 2;
  long dimsp = (ntr + 1) * (ntr + 2);

  double *su = gu.data();
  double *sv = gv.data();

  dv2uv(sd.data(), svo.data(), su, sv, dvTrans.f1.data(), dvTrans.f2.data(), ntr, dimsp, nlev);

  Varray<double> fpwork(nlat * nfc * nlev);

  sp2fc(su, fpwork.data(), spTrans.poli.data(), nlev, nlat, nfc, ntr);
  scaluv(fpwork.data(), spTrans.rcoslat.data(), nlat, nfc * nlev);

  if (fcTrans.use_fftw)
    fc2gp(fpwork.data(), gu.data(), nlat, nlon, nlev, nfc);
  else
    fc2gp(fcTrans.vtrig.data(), fcTrans.ifax, fpwork.data(), gu.data(), nlat, nlon, nlev, nfc);

  sp2fc(sv, fpwork.data(), spTrans.poli.data(), nlev, nlat, nfc, ntr);
  scaluv(fpwork.data(), spTrans.rcoslat.data(), nlat, nfc * nlev);

  if (fcTrans.use_fftw)
    fc2gp(fpwork.data(), gv.data(), nlat, nlon, nlev, nfc);
  else
    fc2gp(fcTrans.vtrig.data(), fcTrans.ifax, fpwork.data(), gv.data(), nlat, nlon, nlev, nfc);
}
