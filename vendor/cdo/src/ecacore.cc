/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast
          Uwe Schulzweida

*/

#include <cstring>

#include <cdi.h>

#include <mpim_grid.h>
#include "process_int.h"
#include "ecacore.h"
#include "ecautil.h"
#include "util_date.h"
#include "datetime.h"
#include "field_functions.h"

constexpr int FIRST_VAR_ID = 0;

#define IS_NOT_SET(x) (x == nullptr)
#define IS_SET(x) (x != nullptr)

static void
init_field(Field &field, int gridID, double missval)
{
  field.grid = gridID;
  field.missval = missval;
}

static void
init_field(Field &field, int gridID, double missval, size_t gridsize)
{
  field.grid = gridID;
  field.missval = missval;
  field.resize(gridsize);
}

static void
init_field(Field &field, int gridID, double missval, size_t gridsize, double value)
{
  field.grid = gridID;
  field.missval = missval;
  field.resize(gridsize, value);
}

void
eca1(const ECA_REQUEST_1 &request)
{
  CdiDateTime ovDateTime{};
  CdiDateTime inDateTime21{};
  int itsID;
  int otsID;

  auto compareDate = request.compare_type;

  auto istreamID = cdo_open_read(0);

  auto ivlistID = cdo_stream_inq_vlist(istreamID);
  auto ovlistID = vlistCreate();

  VarList varList1(ivlistID);
  auto const &varFirst = varList1.vars[FIRST_VAR_ID];
  auto gridID = varFirst.gridID;
  auto zaxisID = varFirst.zaxisID;
  auto missval = varFirst.missval;

  auto ovarID1 = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);
  vlistDefVarMissval(ovlistID, ovarID1, missval);
  if (IS_SET(request.var1.name)) cdiDefKeyString(ovlistID, ovarID1, CDI_KEY_NAME, request.var1.name);
  if (IS_SET(request.var1.longname)) cdiDefKeyString(ovlistID, ovarID1, CDI_KEY_LONGNAME, request.var1.longname);
  if (IS_SET(request.var1.units)) cdiDefKeyString(ovlistID, ovarID1, CDI_KEY_UNITS, request.var1.units);

  if (IS_SET(request.var2.h2) || IS_SET(request.var2.h3))
  {
    auto ovarID2 = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);
    vlistDefVarMissval(ovlistID, ovarID2, missval);
    if (IS_SET(request.var2.name)) cdiDefKeyString(ovlistID, ovarID2, CDI_KEY_NAME, request.var2.name);
    if (IS_SET(request.var2.longname)) cdiDefKeyString(ovlistID, ovarID2, CDI_KEY_LONGNAME, request.var2.longname);
    if (IS_SET(request.var2.units)) cdiDefKeyString(ovlistID, ovarID2, CDI_KEY_UNITS, request.var2.units);
  }

  if (request.compare_type == 16) vlistDefNtsteps(ovlistID, 1);

  auto itaxisID = vlistInqTaxis(ivlistID);
  auto otaxisID = cdo_taxis_create(TAXIS_RELATIVE);
  taxisDefTunit(otaxisID, TUNIT_DAY);
  //  taxisDefTunit(otaxisID, TUNIT_MINUTE);
  //  taxisDefCalendar(otaxisID, CALENDAR_PROLEPTIC);
  taxisDefCalendar(otaxisID, taxisInqCalendar(itaxisID));
  taxisDefRdatetime(otaxisID, cdiDateTime_set(request.var1.refdate, 0));
  vlistDefTaxis(ovlistID, otaxisID);

  auto ostreamID = cdo_open_write(1);
  cdo_def_vlist(ostreamID, ovlistID);

  auto gridsize = gridInqSize(gridID);

  Field field1, field2, field3;
  field1.resize(gridsize);
  field3.resize(gridsize);
  if (IS_SET(request.var2.h2) || IS_SET(request.var2.h3)) field2.resize(gridsize);

  auto nlevels = zaxisInqSize(zaxisID);

  FieldVector var12(nlevels), samp1(nlevels), samp2(nlevels);
  FieldVector var13, var21, var23;
  if (IS_SET(request.var1.f3)) var13.resize(nlevels);
  if (IS_SET(request.var2.h2)) var21.resize(nlevels);
  if (IS_SET(request.var2.h3)) var23.resize(nlevels);

  for (int levelID = 0; levelID < nlevels; ++levelID)
  {
    init_field(var12[levelID], gridID, missval, gridsize);
    init_field(samp1[levelID], gridID, missval, gridsize);
    init_field(samp2[levelID], gridID, missval);

    if (IS_SET(request.var1.f3)) init_field(var13[levelID], gridID, missval, gridsize);
    if (IS_SET(request.var2.h2)) init_field(var21[levelID], gridID, missval, gridsize);
    if (IS_SET(request.var2.h3)) init_field(var23[levelID], gridID, missval, gridsize);
  }

  itsID = 0;
  otsID = 0;
  while (true)
  {
    int numFields = 0;
    long numSets = 0;
    while (true)
    {
      numFields = cdo_stream_inq_timestep(istreamID, itsID);
      if (numFields == 0) break;

      auto ivDateTime = taxisInqVdatetime(itaxisID);

      if (numSets == 0) inDateTime21 = ivDateTime;

      if (date_is_neq(ivDateTime, inDateTime21, compareDate))
      {
        cdo_add_steps(-1);
        break;
      }

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(istreamID);
        if (varID != FIRST_VAR_ID) continue;

        if (numSets == 0)
        {
          if (request.var1.f2 != &vfarnum2)
          {
            field_fill(var12[levelID], missval);
            var12[levelID].numMissVals = gridsize;
          }
          field_fill(samp1[levelID], missval);
          if (!samp2[levelID].empty()) field_fill(samp2[levelID], 0.0);
          if (IS_SET(request.var1.f3)) field_fill(var13[levelID], missval);
          if (IS_SET(request.var2.h2)) field_fill(var21[levelID], missval);
          if (IS_SET(request.var2.h3)) field_fill(var23[levelID], missval);

          samp1[levelID].numMissVals = gridsize;
          if (IS_SET(request.var1.f3)) var13[levelID].numMissVals = gridsize;
          if (IS_SET(request.var2.h2)) var21[levelID].numMissVals = gridsize;
          if (IS_SET(request.var2.h3)) var23[levelID].numMissVals = gridsize;
        }

        cdo_read_field(istreamID, field1);
        field1.grid = var12[levelID].grid;
        field1.missval = var12[levelID].missval;

        vfarnum(samp1[levelID], field1);

        if (IS_SET(request.var2.h2))
        {
          field2.vec_d = field1.vec_d;
          field2.numMissVals = field1.numMissVals;
          field2.grid = field1.grid;
          field2.missval = field1.missval;
        }

        if (IS_SET(request.var1.f1)) request.var1.f1(field1, request.var1.f1arg);

        if (field1.numMissVals || !samp2[levelID].empty())
        {
          if (samp2[levelID].empty())
          {
            samp2[levelID].resize(gridsize);
            field_fill(samp2[levelID], numSets);
          }
          for (size_t i = 0; i < gridsize; ++i)
          {
            if (fp_is_equal(field1.vec_d[i], field1.missval)) continue;
            samp2[levelID].vec_d[i]++;
          }
        }

        if (is_not_equal(request.var1.mulc, 0.0)) fieldc_mul(field1, request.var1.mulc);
        if (is_not_equal(request.var1.addc, 0.0)) fieldc_add(field1, request.var1.addc);

        if (IS_SET(request.var1.f3) && request.var1.f2 == &vfarnum2)
        {
          varray_copy(gridsize, var12[levelID].vec_d, field3.vec_d);
          field3.numMissVals = var12[levelID].numMissVals;
          field3.grid = var12[levelID].grid;
          field3.missval = var12[levelID].missval;
        }

        request.var1.f2(var12[levelID], field1);

        if (IS_SET(request.var2.h2) || IS_SET(request.var2.h3))
        {
          // if h2 is null, use the output of f2 as input for h1
          if (IS_NOT_SET(request.var2.h2))
          {
            varray_copy(gridsize, var12[levelID].vec_d, field2.vec_d);
            field2.numMissVals = var12[levelID].numMissVals;
            field2.grid = var12[levelID].grid;
            field2.missval = var12[levelID].missval;
          }

          if (IS_SET(request.var2.h1)) request.var2.h1(field2, request.var2.h1arg);

          if (IS_NOT_SET(request.var2.h2))
            request.var2.h3(var23[levelID], field2);
          else
          {
            request.var2.h2(var21[levelID], field2);
            if (IS_SET(request.var2.h3)) request.var2.h3(var23[levelID], var21[levelID]);
          }
        }

        if (IS_SET(request.var1.f3))
        {
          if (request.var1.f2 == &vfarnum2)
          {
            auto &array1 = field3.vec_d;
            auto const &array2 = var12[levelID].vec_d;
            auto missval2 = field2.missval;

            auto len = field1.size;
            if (len != field2.size) cdo_abort("Fields have different size (%s)", __func__);

            if (field2.numMissVals)
            {
              for (size_t i = 0; i < len; ++i)
              {
                if (fp_is_not_equal(array2[i], missval2) && fp_is_not_equal(array2[i], 0.0)) array1[i] = 0.0;
              }
            }
            else
            {
              for (size_t i = 0; i < len; ++i)
              {
                if (fp_is_not_equal(array2[i], 0.0)) array1[i] = 0.0;
              }
            }
            request.var1.f3(var13[levelID], field3);
          }
          else
            request.var1.f3(var13[levelID], var12[levelID]);
        }
      }

      ovDateTime = ivDateTime;
      numSets++;
      itsID++;
    }

    if (request.var1.f2 == &vfarnum2)
    {
      for (int levelID = 0; levelID < nlevels; ++levelID)
      {
        if (numFields == 0) { request.var1.f3(var13[levelID], var12[levelID]); }
        else
        {
          auto &array1 = var13[levelID].vec_d;
          auto const &array2 = var12[levelID].vec_d;
          auto missval2 = field2.missval;
          auto len = field1.size;

          for (size_t i = 0; i < len; ++i)
          {
            if (fp_is_equal(array2[i], numSets) || array2[i] > numSets) array1[i] = missval2;
          }
        }
      }
    }

    if (numFields == 0 && numSets == 0) break;

    if (request.var1.epilog == MEAN || request.var1.epilog == PERCENT_OF_TIME)
      for (int levelID = 0; levelID < nlevels; ++levelID)
      {
        auto &rvar = IS_SET(request.var1.f3) ? var13[levelID] : var12[levelID];

        if (samp2[levelID].empty())
          fieldc_div(rvar, numSets);
        else
          field2_div(rvar, samp2[levelID]);

        if (request.var1.epilog == PERCENT_OF_TIME) fieldc_mul(rvar, 100.0);
      }

    if (request.var1.refdate == 19550101)
    {
      taxisDefVdatetime(otaxisID, ovDateTime);
      cdo_def_timestep(ostreamID, otsID);
    }
    else
    {
      int year, month, day;
      cdiDate_decode(inDateTime21.date, &year, &month, &day);
      define_mid_of_time(request.compare_type, otaxisID, year, month, 12);
      cdo_def_timestep(ostreamID, otsID);
    }

    if (otsID && varFirst.timeType == TIME_CONSTANT) continue;

    for (int levelID = 0, varID = 0; levelID < nlevels; ++levelID)
    {
      auto &rvar = IS_SET(request.var1.f3) ? var13[levelID] : var12[levelID];

      vfarsel(rvar, samp1[levelID]);

      cdo_def_field(ostreamID, varID, levelID);
      cdo_write_field(ostreamID, rvar);
    }

    if (IS_SET(request.var2.h2) || IS_SET(request.var2.h3))
    {
      for (int levelID = 0, varID = 1; levelID < nlevels; ++levelID)
      {
        auto &rvar = IS_SET(request.var2.h3) ? var23[levelID] : var21[levelID];

        vfarsel(rvar, samp1[levelID]);

        cdo_def_field(ostreamID, varID, levelID);
        cdo_write_field(ostreamID, rvar);
      }
    }

    if (numFields == 0) break;
    otsID++;
  }

  cdo_stream_close(ostreamID);
  cdo_stream_close(istreamID);
}

void
eca2(const ECA_REQUEST_2 &request)
{
  CdiDateTime ovDateTime{};
  CdiDateTime inDateTime21{};
  int itsID;
  int otsID;

  auto compareDate = request.compare_type;

  auto istreamID1 = cdo_open_read(0);
  auto istreamID2 = cdo_open_read(1);

  auto ivlistID1 = cdo_stream_inq_vlist(istreamID1);
  auto ivlistID2 = cdo_stream_inq_vlist(istreamID2);
  auto ovlistID = vlistCreate();

  VarList varList1(ivlistID1);
  VarList varList2(ivlistID2);
  auto const &var1First = varList1.vars[FIRST_VAR_ID];
  auto const &var2First = varList2.vars[FIRST_VAR_ID];

  varList_compare(varList1, varList2);

  auto gridID = var1First.gridID;
  auto zaxisID = var1First.zaxisID;
  auto missval1 = var1First.missval;
  auto missval2 = var2First.missval;

  auto ovarID1 = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);
  vlistDefVarMissval(ovlistID, ovarID1, missval1);
  if (IS_SET(request.var1.name)) cdiDefKeyString(ovlistID, ovarID1, CDI_KEY_NAME, request.var1.name);
  if (IS_SET(request.var1.longname)) cdiDefKeyString(ovlistID, ovarID1, CDI_KEY_LONGNAME, request.var1.longname);
  if (IS_SET(request.var1.units)) cdiDefKeyString(ovlistID, ovarID1, CDI_KEY_UNITS, request.var1.units);

  if (IS_SET(request.var2.h2))
  {
    auto ovarID2 = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);
    vlistDefVarMissval(ovlistID, ovarID2, missval2);
    if (IS_SET(request.var2.name)) cdiDefKeyString(ovlistID, ovarID2, CDI_KEY_NAME, request.var2.name);
    if (IS_SET(request.var2.longname)) cdiDefKeyString(ovlistID, ovarID2, CDI_KEY_LONGNAME, request.var2.longname);
    if (IS_SET(request.var2.units)) cdiDefKeyString(ovlistID, ovarID2, CDI_KEY_UNITS, request.var2.units);
  }

  if (request.compare_type == 16) vlistDefNtsteps(ovlistID, 1);

  auto itaxisID1 = vlistInqTaxis(ivlistID1);
  auto otaxisID = cdo_taxis_create(TAXIS_RELATIVE);
  taxisDefTunit(otaxisID, TUNIT_DAY);
  //  taxisDefTunit(otaxisID, TUNIT_MINUTE);
  //  taxisDefCalendar(otaxisID, CALENDAR_PROLEPTIC);
  taxisDefCalendar(otaxisID, taxisInqCalendar(itaxisID1));
  taxisDefRdatetime(otaxisID, cdiDateTime_set(request.var1.refdate, 0));
  vlistDefTaxis(ovlistID, otaxisID);

  auto ostreamID = cdo_open_write(2);
  cdo_def_vlist(ostreamID, ovlistID);

  auto gridsize = gridInqSize(gridID);

  Field field1, field2, field3;
  field1.resize(gridsize);
  field2.resize(gridsize);
  field3.resize(gridsize);
  constexpr int MaxDays = 373;
  FieldVector2D varsData2[MaxDays];

  auto nlevels = zaxisInqSize(zaxisID);

  FieldVector var14(nlevels), samp1(nlevels), samp2(nlevels), samp3(nlevels);
  FieldVector total, var15, var22;
  if (request.var1.epilog == PERCENT_OF_TOTAL_AMOUNT) total.resize(nlevels);
  if (IS_SET(request.var1.f5)) var15.resize(nlevels);
  if (IS_SET(request.var2.h2)) var22.resize(nlevels);

  for (int levelID = 0; levelID < nlevels; ++levelID)
  {
    init_field(var14[levelID], gridID, missval1, gridsize);
    init_field(samp1[levelID], gridID, missval1, gridsize);
    init_field(samp2[levelID], gridID, missval1, gridsize);
    init_field(samp3[levelID], gridID, missval1);

    if (request.var1.epilog == PERCENT_OF_TOTAL_AMOUNT) init_field(total[levelID], gridID, missval1, gridsize);
    if (IS_SET(request.var1.f5)) init_field(var15[levelID], gridID, missval1, gridsize);
    if (IS_SET(request.var2.h2)) init_field(var22[levelID], gridID, missval1, gridsize);
  }

  itsID = 0;
  while (true)
  {
    auto numFields = cdo_stream_inq_timestep(istreamID2, itsID);
    if (numFields == 0) break;

    auto ivDateTime = taxisInqVdatetime(vlistInqTaxis(ivlistID2));

    auto dayOfYear = decode_day_of_year(ivDateTime.date);
    if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);

    if (!varsData2[dayOfYear].size()) field2D_init(varsData2[dayOfYear], varList2, FIELD_VEC);

    for (int fieldID = 0; fieldID < numFields; ++fieldID)
    {
      auto [varID, levelID] = cdo_inq_field(istreamID2);
      if (varID != FIRST_VAR_ID) continue;
      auto &rvar = varsData2[dayOfYear][0][levelID];
      cdo_read_field(istreamID2, rvar);
    }

    itsID++;
  }

  itsID = 0;
  otsID = 0;

  while (true)
  {
    int numFields = 0;
    long numSets = 0;
    while (true)
    {
      numFields = cdo_stream_inq_timestep(istreamID1, itsID);
      if (numFields == 0) break;

      auto ivDateTime = taxisInqVdatetime(itaxisID1);

      auto dayOfYear = decode_day_of_year(ivDateTime.date);
      if (!varsData2[dayOfYear].size()) cdo_abort("Input streams have different time values!");

      if (numSets == 0) inDateTime21 = ivDateTime;

      if (date_is_neq(ivDateTime, inDateTime21, compareDate))
      {
        cdo_add_steps(-1);
        break;
      }

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(istreamID1);
        if (varID != FIRST_VAR_ID) continue;

        if (numSets == 0)
        {
          if (request.var1.f4 != &vfarnum2)
          {
            field_fill(var14[levelID], missval1);
            var14[levelID].numMissVals = gridsize;
          }
          field_fill(samp1[levelID], missval1);
          field_fill(samp2[levelID], missval1);
          if (!samp3[levelID].empty()) field_fill(samp3[levelID], 0.0);
          if (request.var1.epilog == PERCENT_OF_TOTAL_AMOUNT) field_fill(total[levelID], 0.0);
          if (IS_SET(request.var1.f5)) field_fill(var15[levelID], missval1);
          if (IS_SET(request.var2.h2)) field_fill(var22[levelID], missval1);

          samp1[levelID].numMissVals = gridsize;
          samp2[levelID].numMissVals = gridsize;
          if (request.var1.epilog == PERCENT_OF_TOTAL_AMOUNT) total[levelID].numMissVals = gridsize;
          if (IS_SET(request.var1.f5)) var15[levelID].numMissVals = gridsize;
          if (IS_SET(request.var2.h2)) var22[levelID].numMissVals = gridsize;
        }

        cdo_read_field(istreamID1, field1);
        field1.grid = gridID;
        field1.missval = missval1;
        field2.grid = varsData2[dayOfYear][0][levelID].grid;
        field2.numMissVals = varsData2[dayOfYear][0][levelID].numMissVals;
        field2.missval = missval2;
        field_copy(varsData2[dayOfYear][0][levelID], field2);

        vfarnum(samp1[levelID], field1);
        vfarnum(samp2[levelID], field2);

        if (request.var1.epilog == PERCENT_OF_TOTAL_AMOUNT) field2_sum(total[levelID], field1);

        if (IS_SET(request.var1.f1)) request.var1.f1(field1, request.var1.f1arg);
        if (IS_SET(request.var1.f2)) request.var1.f2(field2, request.var1.f2arg);

        if (field1.numMissVals || !samp3[levelID].empty())
        {
          if (samp3[levelID].empty())
          {
            samp3[levelID].resize(gridsize);
            field_fill(samp3[levelID], numSets);
          }
          for (size_t i = 0; i < gridsize; ++i)
          {
            if (fp_is_equal(field1.vec_d[i], field1.missval)) continue;
            samp3[levelID].vec_d[i]++;
          }
        }

        if (IS_SET(request.var1.f5) && request.var1.f4 == &vfarnum2)
        {
          varray_copy(gridsize, var14[levelID].vec_d, field3.vec_d);
          field3.numMissVals = var14[levelID].numMissVals;
          field3.grid = var14[levelID].grid;
          field3.missval = var14[levelID].missval;
        }

        request.var1.f3(field1, field2);
        request.var1.f4(var14[levelID], field1);

        if (IS_SET(request.var2.h2))
        {
          varray_copy(gridsize, var14[levelID].vec_d, field2.vec_d);
          field2.numMissVals = var14[levelID].numMissVals;
          field2.grid = var14[levelID].grid;
          field2.missval = var14[levelID].missval;

          if (IS_SET(request.var2.h1)) request.var2.h1(field2, request.var2.h1arg);

          request.var2.h2(var22[levelID], field2);
        }

        if (IS_SET(request.var1.f5))
        {
          if (request.var1.f4 == &vfarnum2)
          {
            auto &array1 = field3.vec_d;
            auto const &array2 = var14[levelID].vec_d;
            auto missvaltemp = field1.missval;

            auto len = field1.size;
            if (len != field3.size) cdo_abort("Fields have different size (%s)", __func__);

            if (field1.numMissVals)
            {
              for (size_t i = 0; i < len; ++i)
              {
                if (fp_is_not_equal(array2[i], missvaltemp) && fp_is_not_equal(array2[i], 0.0)) array1[i] = 0.0;
              }
            }
            else
            {
              for (size_t i = 0; i < len; ++i)
              {
                if (fp_is_not_equal(array2[i], 0.0)) array1[i] = 0.0;
              }
            }
            request.var1.f5(var15[levelID], field3, request.var1.f5arg);
          }
          else
            request.var1.f5(var15[levelID], var14[levelID], request.var1.f5arg);
        }
      }

      ovDateTime = ivDateTime;
      numSets++;
      itsID++;
    }

    if (request.var1.f4 == &vfarnum2)
    {
      for (int levelID = 0; levelID < nlevels; ++levelID)
      {
        if (numFields == 0) { request.var1.f5(var15[levelID], var14[levelID], request.var1.f5arg); }
        else
        {
          auto const &array2 = var14[levelID].vec_d;
          auto &array1 = var15[levelID].vec_d;
          auto len = field1.size;
          auto missvaltemp = field1.missval;
          for (size_t i = 0; i < len; ++i)
          {
            if (fp_is_equal(array2[i], numSets) || array2[i] > numSets) array1[i] = missvaltemp;
          }
        }
      }
    }

    if (numFields == 0 && numSets == 0) break;

    if (request.var1.epilog == MEAN || request.var1.epilog == PERCENT_OF_TIME)
      for (int levelID = 0; levelID < nlevels; ++levelID)
      {
        auto &rvar = IS_SET(request.var1.f5) ? var15[levelID] : var14[levelID];

        if (samp3[levelID].empty())
          fieldc_div(rvar, numSets);
        else
          field2_div(rvar, samp3[levelID]);

        if (request.var1.epilog == PERCENT_OF_TIME) fieldc_mul(rvar, 100.0);
      }
    else if (request.var1.epilog == PERCENT_OF_TOTAL_AMOUNT)
      for (int levelID = 0; levelID < nlevels; ++levelID)
      {
        Field &rvar = IS_SET(request.var1.f5) ? var15[levelID] : var14[levelID];

        field2_div(rvar, total[levelID]);
        fieldc_mul(rvar, 100.0);
      }

    if (request.var1.refdate == 19550101)
    {
      taxisDefVdatetime(otaxisID, ovDateTime);
      cdo_def_timestep(ostreamID, otsID);
    }
    else
    {
      int year, month, day;
      cdiDate_decode(inDateTime21.date, &year, &month, &day);
      define_mid_of_time(request.compare_type, otaxisID, year, month, 12);
      cdo_def_timestep(ostreamID, otsID);
    }

    if (otsID && var1First.timeType == TIME_CONSTANT) continue;

    for (int levelID = 0, varID = 0; levelID < nlevels; ++levelID)
    {
      Field &rvar = IS_SET(request.var1.f5) ? var15[levelID] : var14[levelID];

      vfarsel(rvar, samp1[levelID]);
      vfarsel(rvar, samp2[levelID]);

      cdo_def_field(ostreamID, varID, levelID);
      cdo_write_field(ostreamID, rvar);
    }

    if (IS_SET(request.var2.h2))
    {
      for (int levelID = 0, varID = 1; levelID < nlevels; ++levelID)
      {
        auto &rvar = var22[levelID];

        vfarsel(rvar, samp1[levelID]);
        vfarsel(rvar, samp2[levelID]);

        cdo_def_field(ostreamID, varID, levelID);
        cdo_write_field(ostreamID, rvar);
      }
    }

    if (numFields == 0) break;

    otsID++;
  }

  cdo_stream_close(ostreamID);
  cdo_stream_close(istreamID2);
  cdo_stream_close(istreamID1);
}

void
eca3(const ECA_REQUEST_3 &request)
{
  CdiDateTime ovDateTime{};
  CdiDateTime inDateTime21{};
  int itsID;
  int otsID;

  auto compareDate = request.compare_type;

  auto istreamID1 = cdo_open_read(0);
  auto istreamID2 = cdo_open_read(1);

  auto ivlistID1 = cdo_stream_inq_vlist(istreamID1);
  auto ivlistID2 = cdo_stream_inq_vlist(istreamID2);
  auto ovlistID = vlistCreate();

  VarList varList1(ivlistID1);
  VarList varList2(ivlistID2);
  auto const &var1First = varList1.vars[FIRST_VAR_ID];

  varList_compare(varList1, varList2);

  auto gridID = var1First.gridID;
  auto zaxisID = var1First.zaxisID;
  auto missval = var1First.missval;

  auto ovarID1 = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);
  vlistDefVarMissval(ovlistID, ovarID1, missval);
  if (IS_SET(request.name)) cdiDefKeyString(ovlistID, ovarID1, CDI_KEY_NAME, request.name);
  if (IS_SET(request.longname)) cdiDefKeyString(ovlistID, ovarID1, CDI_KEY_LONGNAME, request.longname);
  if (IS_SET(request.units)) cdiDefKeyString(ovlistID, ovarID1, CDI_KEY_UNITS, request.units);

  if (request.compare_type == 16) vlistDefNtsteps(ovlistID, 1);

  auto itaxisID1 = vlistInqTaxis(ivlistID1);
  auto itaxisID2 = vlistInqTaxis(ivlistID2);
  auto otaxisID = cdo_taxis_create(TAXIS_RELATIVE);
  taxisDefTunit(otaxisID, TUNIT_DAY);
  //  taxisDefTunit(otaxisID, TUNIT_MINUTE);
  //  taxisDefCalendar(otaxisID, CALENDAR_PROLEPTIC);
  taxisDefCalendar(otaxisID, taxisInqCalendar(itaxisID1));
  taxisDefRdatetime(otaxisID, cdiDateTime_set(request.refdate, 0));
  vlistDefTaxis(ovlistID, otaxisID);

  auto ostreamID = cdo_open_write(2);
  cdo_def_vlist(ostreamID, ovlistID);

  auto gridsize = gridInqSize(gridID);

  Field field1, field2;
  field1.resize(gridsize);
  field2.resize(gridsize);

  auto numLevels = zaxisInqSize(zaxisID);

  FieldVector var1(numLevels), var2(numLevels);

  for (int levelID = 0; levelID < numLevels; ++levelID)
  {
    init_field(var1[levelID], gridID, missval, gridsize);
    init_field(var2[levelID], gridID, missval, gridsize);
  }

  itsID = 0;
  otsID = 0;
  while (true)
  {
    int numFields = 0;
    long numSets = 0;
    while (true)
    {
      numFields = cdo_stream_inq_timestep(istreamID1, itsID);
      if (numFields == 0) break;

      if (!cdo_stream_inq_timestep(istreamID2, itsID)) cdo_abort("Input streams have different number of time steps!");

      auto ivDateTime1 = taxisInqVdatetime(itaxisID1);
      auto ivdate1 = cdiDate_get(ivDateTime1.date);
      auto ivtime1 = cdiTime_get(ivDateTime1.time);

      auto ivDateTime2 = taxisInqVdatetime(itaxisID2);
      auto ivdate2 = cdiDate_get(ivDateTime2.date);
      auto ivtime2 = cdiTime_get(ivDateTime2.time);

      if (ivdate1 != ivdate2) cdo_abort("Input streams have different verification dates at time step %d!", itsID + 1);
      if (ivtime1 != ivtime2) cdo_abort("Input streams have different verification times at time step %d!", itsID + 1);

      if (numSets == 0) inDateTime21 = ivDateTime1;

      if (date_is_neq(ivDateTime1, inDateTime21, compareDate))
      {
        cdo_add_steps(-1);
        break;
      }

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        (void) cdo_inq_field(istreamID1);
        auto [varID, levelID] = cdo_inq_field(istreamID2);
        if (varID != FIRST_VAR_ID) continue;

        if (numSets == 0)
        {
          for (size_t i = 0; i < gridsize; ++i)
          {
            var1[levelID].vec_d[i] = missval;
            var2[levelID].vec_d[i] = missval;
          }
          var1[levelID].numMissVals = gridsize;
          var2[levelID].numMissVals = gridsize;
        }

        cdo_read_field(istreamID1, field1);
        field1.grid = var1[levelID].grid;
        field1.missval = var1[levelID].missval;

        cdo_read_field(istreamID2, field2);
        field2.grid = var1[levelID].grid;
        field2.missval = var1[levelID].missval;

        request.f1(var1[levelID], field1);
        request.f2(var2[levelID], field2);
      }

      ovDateTime = ivDateTime1;
      numSets++;
      itsID++;
    }

    if (numFields == 0 && numSets == 0) break;

    for (int levelID = 0; levelID < numLevels; ++levelID) request.f3(var1[levelID], var2[levelID]);

    if (request.refdate == 19550101)
    {
      taxisDefVdatetime(otaxisID, ovDateTime);
      cdo_def_timestep(ostreamID, otsID);
    }
    else
    {
      int year, month, day;
      cdiDate_decode(inDateTime21.date, &year, &month, &day);
      define_mid_of_time(request.compare_type, otaxisID, year, month, 12);
      cdo_def_timestep(ostreamID, otsID);
    }

    if (otsID && var1First.timeType == TIME_CONSTANT) continue;

    for (int levelID = 0, varID = 0; levelID < numLevels; ++levelID)
    {
      cdo_def_field(ostreamID, varID, levelID);
      cdo_write_field(ostreamID, var1[levelID]);
    }

    if (numFields == 0) break;
    otsID++;
  }

  cdo_stream_close(ostreamID);
  cdo_stream_close(istreamID2);
  cdo_stream_close(istreamID1);
}

// check for non missval values
static bool
fldhvs(FieldVector const &fieldVector, size_t nlevels)
{
  for (size_t level = 0; level < nlevels; level++)
  {
    if (fieldVector[level].numMissVals != fieldVector[level].size) return true;
  }

  return false;
}

void
eca4(const ECA_REQUEST_4 &request)
{
  int yearcnt = 0;
  bool resetAtJan = false, resetAtJul = false;
  bool isFirstYear = true;
  CdiDateTime ovDateTime{};
  CdiDateTime inDateTime21{};
  int64_t ivdate = 0, ovdate = 0;

  auto compareDate = request.compare_type;

  auto istreamID1 = cdo_open_read(0);
  auto istreamID2 = cdo_open_read(1);

  auto ivlistID1 = cdo_stream_inq_vlist(istreamID1);
  auto ivlistID2 = cdo_stream_inq_vlist(istreamID2);
  auto ovlistID = vlistCreate();

  VarList varList1(ivlistID1);
  VarList varList2(ivlistID2);
  auto const &var1First = varList1.vars[FIRST_VAR_ID];
  auto const &var2First = varList2.vars[FIRST_VAR_ID];

  auto gridID = var1First.gridID;
  if (var1First.gridsize != var2First.gridsize) cdo_abort("Grid sizes of the input fields do not match!");

  auto zaxisID = var1First.zaxisID;
  auto missval = var1First.missval;

  auto ovarID1 = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);
  vlistDefVarMissval(ovlistID, ovarID1, missval);
  if (IS_SET(request.name)) cdiDefKeyString(ovlistID, ovarID1, CDI_KEY_NAME, request.name);
  if (IS_SET(request.longname)) cdiDefKeyString(ovlistID, ovarID1, CDI_KEY_LONGNAME, request.longname);
  if (IS_SET(request.units)) cdiDefKeyString(ovlistID, ovarID1, CDI_KEY_UNITS, request.units);

  auto ovarID2 = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);
  vlistDefVarMissval(ovlistID, ovarID2, missval);
  if (IS_SET(request.name2)) cdiDefKeyString(ovlistID, ovarID2, CDI_KEY_NAME, request.name2);
  if (IS_SET(request.longname2)) cdiDefKeyString(ovlistID, ovarID2, CDI_KEY_LONGNAME, request.longname2);
  if (IS_SET(request.units2)) cdiDefKeyString(ovlistID, ovarID2, CDI_KEY_UNITS, request.units2);

  if (request.compare_type == 16) vlistDefNtsteps(ovlistID, 1);

  auto itaxisID = vlistInqTaxis(ivlistID1);
  auto otaxisID = cdo_taxis_create(TAXIS_RELATIVE);
  taxisDefTunit(otaxisID, TUNIT_DAY);
  //  taxisDefTunit(otaxisID, TUNIT_MINUTE);
  //  taxisDefCalendar(otaxisID, CALENDAR_PROLEPTIC);
  taxisDefCalendar(otaxisID, taxisInqCalendar(itaxisID));
  taxisDefRdatetime(otaxisID, cdiDateTime_set(19550101, 0));
  vlistDefTaxis(ovlistID, otaxisID);

  auto ostreamID = cdo_open_write(2);
  cdo_def_vlist(ostreamID, ovlistID);

  bool lyvals = true;
  auto gridtype = gridInqType(gridID);
  if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_PROJECTION)
  {
    gridID = gridToCurvilinear(gridID, NeedCorners::Yes);
  }
  else if (gridtype == GRID_GME) { gridID = gridToUnstructured(gridID, NeedCorners::Yes); }
  else { lyvals = false; }

  auto gridsize = gridInqSize(gridID);
  // for later check on northern\southern hemisphere
  std::vector<double> yvals(gridsize);
  if (lyvals) { gridInqYvals(gridID, yvals.data()); }
  else
  {
    for (size_t i = 0; i < gridsize; ++i) yvals[i] = 20;  // Northern hemisphere
  }

  // Two fields are needed because of the definition of gsl for northern and southern hemisphere
  Field fieldGt, fieldLt;
  fieldGt.resize(gridsize);
  fieldLt.resize(gridsize);

  // field for the land-water-distribution
  Field mask;
  mask.resize(gridsize);

  auto nlevels = zaxisInqSize(zaxisID);

  FieldVector startCount(nlevels), endCount(nlevels);
  FieldVector gslDuration(nlevels), gslFirstDay(nlevels);

  FieldVector2D startDateWithHist(2), endDateWithHist(2);
  /* because of the different definitions for northern and southern hemisphere,
   * the values of the last year have to be present THE LAST YEAR HAS THE INDEX 1 */
  for (int h = 0; h < 2; h++)
  {
    startDateWithHist[h].resize(nlevels);
    endDateWithHist[h].resize(nlevels);
  }

  for (int levelID = 0; levelID < nlevels; ++levelID)
  {
    init_field(startCount[levelID], gridID, missval, gridsize, 0.0);
    init_field(endCount[levelID], gridID, missval, gridsize, 0.0);
    init_field(gslDuration[levelID], gridID, missval, gridsize);
    init_field(gslFirstDay[levelID], gridID, missval, gridsize);

    for (int h = 0; h < 2; h++) init_field(startDateWithHist[h][levelID], gridID, missval, gridsize);
    for (int h = 0; h < 2; h++) init_field(endDateWithHist[h][levelID], gridID, missval, gridsize);
  }

  int itsID = 0;
  int otsID = 0;

  if (cdo_stream_inq_timestep(istreamID2, itsID))
  {
    (void) cdo_inq_field(istreamID2);
    cdo_read_field(istreamID2, mask);
    mask.grid = gridID;
    mask.missval = var2First.missval;

    request.s3(mask, request.s3arg);
  }
  else
    cdo_abort("Could not read land-water mask!");

  while (true)
  {
    int numFields = 0;
    long numSets = 0;
    while (true)
    {
      numFields = cdo_stream_inq_timestep(istreamID1, itsID);
      if (numFields == 0) break;

      auto ivDateTime = taxisInqVdatetime(itaxisID);
      ivdate = cdiDate_get(ivDateTime.date);

      int month = (ivdate % 10000) / 100;
      if (month < 1 || month > 12) cdo_abort("month %d out of range!", month);

      if (numSets == 0) inDateTime21 = ivDateTime;

      if (date_is_neq(ivDateTime, inDateTime21, compareDate))
      {
        resetAtJan = false;
        resetAtJul = false;
        cdo_add_steps(-1);
        break;
      }

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, ilevelID] = cdo_inq_field(istreamID1);
        auto levelID = ilevelID;
        if (varID != FIRST_VAR_ID) continue;

        if (numSets == 0)
        {
          field_fill(gslDuration[levelID], missval);
          field_fill(gslFirstDay[levelID], missval);
          // reinitialize the current year
          field_fill(startDateWithHist[0][levelID], missval);
          field_fill(endDateWithHist[0][levelID], missval);

          gslDuration[levelID].numMissVals = 0;
          gslFirstDay[levelID].numMissVals = 0;
          // reinitialize the current year
          startDateWithHist[0][levelID].numMissVals = gridsize;
          endDateWithHist[0][levelID].numMissVals = gridsize;
        }
        // init the history ONCE
        if (0 == itsID)
        {
          field_fill(startDateWithHist[1][levelID], missval);
          field_fill(endDateWithHist[1][levelID], missval);

          startDateWithHist[1][levelID].numMissVals = gridsize;
          endDateWithHist[1][levelID].numMissVals = gridsize;
        }

        cdo_read_field(istreamID1, fieldGt);
        fieldLt.vec_d = fieldGt.vec_d;
        fieldLt.numMissVals = fieldGt.numMissVals;
        fieldGt.grid = startCount[levelID].grid;
        fieldGt.missval = startCount[levelID].missval;
        fieldLt.grid = startCount[levelID].grid;
        fieldLt.missval = startCount[levelID].missval;

        // Reinitialization of (start|end)Count variables has to be done different for norther and southern hemisphere
        if (1 == month && !resetAtJan)
        {
          // reset northern startCount
          for (size_t i = 0; i < gridsize; ++i)
          {
            if (yvals[i] >= 0.0)
              if (fp_is_not_equal(startCount[levelID].vec_d[i], missval))
              {
                startCount[levelID].vec_d[i] = missval;
                startCount[levelID].numMissVals++;
              }
          }
          // reset southern endCount
          for (size_t i = 0; i < gridsize; ++i)
          {
            if (yvals[i] < 0.0)
              if (fp_is_not_equal(endCount[levelID].vec_d[i], missval))
              {
                endCount[levelID].vec_d[i] = missval;
                endCount[levelID].numMissVals++;
              }
          }

          resetAtJan = true;
        }
        if (7 == month && !resetAtJul)
        {
#ifdef _OPENMP
#pragma omp sections
#endif
          {
#ifdef _OPENMP
#pragma omp section
#endif
            {
              // reset northern endCount
              for (size_t i = 0; i < gridsize; ++i)
              {
                if (yvals[i] >= 0.0)
                {
                  if (fp_is_not_equal(endCount[levelID].vec_d[i], missval))
                  {
                    endCount[levelID].vec_d[i] = missval;
                    endCount[levelID].numMissVals++;
                  }
                }
              }
            }
#ifdef _OPENMP
#pragma omp section
#endif
            {
              // reset southern startCount
              for (size_t i = 0; i < gridsize; ++i)
              {
                if (yvals[i] < 0.0)
                {
                  if (fp_is_not_equal(startCount[levelID].vec_d[i], missval))
                  {
                    startCount[levelID].vec_d[i] = missval;
                    startCount[levelID].numMissVals++;
                  }
                }
              }
            }
          }
          resetAtJul = true;
        }

        // count the day with temperature larger/smaller than the given limit
#ifdef _OPENMP
#pragma omp sections
#endif
        {
#ifdef _OPENMP
#pragma omp section
#endif
          {
            vfarsel(fieldGt, mask);
            request.s1(fieldGt, request.s1arg);
            vfarnum2(startCount[levelID], fieldGt);
          }
#ifdef _OPENMP
#pragma omp section
#endif
          {
            vfarsel(fieldLt, mask);
            request.s2(fieldLt, request.s1arg);
            vfarnum2(endCount[levelID], fieldLt);
          }
        }

        if (month < 7)
        {
          for (size_t i = 0; i < gridsize; ++i)  // dictinct between northern and southern sphere
            // start with south
            if (yvals[i] < 0)
            {
              // south: periods can also start in the first half of the year, but this date has already gone into the
              // history
              if (fp_is_equal(startDateWithHist[1][levelID].vec_d[i], missval)
                  && is_equal(startCount[levelID].vec_d[i], request.consecutiveDays))
              {
                startDateWithHist[1][levelID].vec_d[i] = ivdate;
                // reset the endCount, because we are only interessted in the end of the eriod, if a start was found
                endCount[levelID].vec_d[i] = missval;
                endDateWithHist[0][levelID].vec_d[i] = missval;
              }
              if (fp_is_equal(endDateWithHist[0][levelID].vec_d[i], missval)
                  && is_equal(endCount[levelID].vec_d[i], request.consecutiveDays))
              {
                endDateWithHist[0][levelID].vec_d[i] = ivdate;
              }
            }
            else
            {
              if (fp_is_equal(startDateWithHist[0][levelID].vec_d[i], missval)
                  && is_equal(startCount[levelID].vec_d[i], request.consecutiveDays))
              {
                startDateWithHist[0][levelID].vec_d[i] = ivdate;
              }
            }
        }
        else
        {
          for (size_t i = 0; i < gridsize; ++i)
          {
            if (yvals[i] < 0)
            {
              if (fp_is_equal(startDateWithHist[0][levelID].vec_d[i], missval)
                  && is_equal(startCount[levelID].vec_d[i], request.consecutiveDays))
              {
                startDateWithHist[0][levelID].vec_d[i] = ivdate;
              }
            }
            else
            {
              // north: periods can also start in the second half of the year
              if (fp_is_equal(startDateWithHist[0][levelID].vec_d[i], missval)
                  && is_equal(startCount[levelID].vec_d[i], request.consecutiveDays))
              {
                startDateWithHist[0][levelID].vec_d[i] = ivdate;
                // reset the endCount, because we are only interessted in the end of the eriod, if a start was found
                endCount[levelID].vec_d[i] = missval;
                endDateWithHist[0][levelID].vec_d[i] = missval;
              }
              if (fp_is_equal(endDateWithHist[0][levelID].vec_d[i], missval)
                  && is_equal(endCount[levelID].vec_d[i], request.consecutiveDays))
              {
                endDateWithHist[0][levelID].vec_d[i] = ivdate;
              }
            }
          }
        }
        // update numMissVals for saving data in GRIB
        field_num_mv(startCount[levelID]);
        field_num_mv(endCount[levelID]);
        field_num_mv(startDateWithHist[1][levelID]);
        field_num_mv(startDateWithHist[0][levelID]);
        field_num_mv(endDateWithHist[1][levelID]);
        field_num_mv(endDateWithHist[0][levelID]);
      }

      ovDateTime = ivDateTime;
      ovdate = cdiDate_get(ovDateTime.date);
      numSets++;
      itsID++;
    }

    if (numFields == 0 && numSets == 0) break;

    adjust_end_date(nlevels, gridsize, yvals, missval, ovdate, startDateWithHist, endDateWithHist);

    /*  compute and write GSL for the previous year
     *  AND
     *  write the current start/end dates into the history
     *
     *  this is the default action if more than a year is available */
    if (yearcnt != 0)
    {
      compute_gsl(nlevels, gridsize, yvals, missval, startDateWithHist, endDateWithHist, gslDuration, gslFirstDay, false);

      // values of the privous year
      ovDateTime.date = cdiDate_encode(ovdate / 10000 - 1, 12, 31);
      write_gsl_stream(ostreamID, otaxisID, otsID, ovarID1, ovarID2, ivlistID1, FIRST_VAR_ID, gslDuration, gslFirstDay, ovDateTime,
                       nlevels);
      otsID++;
    }

    // if there is a previous year
    if (ovdate != ivdate)
    {
      /*  if the first year of data was processed, the history has to
       *  be checked befor it get's updated. This is necessary, if a
       *  growing period on the southern hemisphere was found. Otherwise,
       *  it would get overwritten. */
      if (isFirstYear)
      {
        // Check for non missing values, i.e. is there any data for the previous year?
        if (fldhvs(startDateWithHist[1], nlevels))
        {
          compute_gsl(nlevels, gridsize, yvals, missval, startDateWithHist, endDateWithHist, gslDuration, gslFirstDay, false);
          ovDateTime.date = cdiDate_encode(ovdate / 10000 - 1, 12, 31);
          write_gsl_stream(ostreamID, otaxisID, otsID, ovarID1, ovarID2, ivlistID1, FIRST_VAR_ID, gslDuration, gslFirstDay,
                           ovDateTime, nlevels);
          otsID++;
        }
        isFirstYear = false;
      }
#ifdef _OPENMP
#pragma omp sections
#endif
      {
        update_hist(startDateWithHist, nlevels, gridsize, yvals, false);
#ifdef _OPENMP
#pragma omp section
#endif
        update_hist(endDateWithHist, nlevels, gridsize, yvals, true);
      }
    }
    else  // process the current year, this only happens, if the last timestep is reached OR if data for only one year is present
    {
      compute_gsl(nlevels, gridsize, yvals, missval, startDateWithHist, endDateWithHist, gslDuration, gslFirstDay, true);

      write_gsl_stream(ostreamID, otaxisID, otsID, ovarID1, ovarID2, ivlistID1, FIRST_VAR_ID, gslDuration, gslFirstDay, ovDateTime,
                       nlevels);
      otsID++;
    }
    yearcnt++;

    if (numFields == 0) break;
  }

  cdo_stream_close(ostreamID);
  cdo_stream_close(istreamID2);
  cdo_stream_close(istreamID1);
}
