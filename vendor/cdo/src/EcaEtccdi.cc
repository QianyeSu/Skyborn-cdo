/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Fabian Wachsmann

*/

/*
   This module contains the following operators:

      EcaEtccdi    eca_etccdi         Etccdi conform indices
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_vlist.h"
#include "datetime.h"
#include "process_int.h"
#include "param_conversion.h"
#include "percentiles_hist.h"
#include "percentiles.h"
#include "util_date.h"
#include "ecautil.h"
#include "ecacore.h"
#include "field_functions.h"

enum functions_select
{
  func_selle = 1,
  func_selge = 2
};

static const char TX90P_UNITS[] = "%";
static const char TX90P_NAME[] = "tx90pETCCDI";
static const char TX90P_LONGNAME[] = "Percentage of Days when Daily Maximum Temperature is Above the 90th Percentile";

static const char TX10P_UNITS[] = "%";
static const char TX10P_NAME[] = "tx10pETCCDI";
static const char TX10P_LONGNAME[] = "Percentage of Days when Daily Maximum Temperature is Below the 10th Percentile";

static const char TN90P_UNITS[] = "%";
static const char TN90P_NAME[] = "tn90pETCCDI";
static const char TN90P_LONGNAME[] = "Percentage of Days when Daily Minimum Temperature is Above the 90th Percentile";

static const char TN10P_UNITS[] = "%";
static const char TN10P_NAME[] = "tn10pETCCDI";
static const char TN10P_LONGNAME[] = "Percentage of Days when Daily Minimum Temperature is Below the 10th Percentile";

static const char R99P_UNITS[] = "mm";
static const char R99P_NAME[] = "r99pETCCDI";
static const char R99P_LONGNAME[]
    = "Annual Total Precipitation when Daily Precipitation Exceeds the 99th Percentile of Wet Day Precipitation";

static const char R95P_UNITS[] = "mm";
static const char R95P_NAME[] = "r95pETCCDI";
static const char R95P_LONGNAME[]
    = "Annual Total Precipitation when Daily Precipitation Exceeds the 95th Percentile of Wet Day Precipitation";

/* windowDays()
   Saves for each bootstrap year
   for each day of year
   the day of each window day */
static void
windowDays(int dayOfYear, std::vector<int> &wdays, std::vector<bool> const &wdaysRead, int MaxDays, int ndates, int sumboot)
{
  int wdayOfYear = dayOfYear;
  for (int gobackDays = ceil(ndates / 2. - 1); gobackDays != 0; gobackDays--)
    {
      wdayOfYear--;
      if (wdayOfYear < 1) wdayOfYear += (MaxDays - 1) * sumboot;
      while (wdaysRead[wdayOfYear] == false)
        {
          wdayOfYear--;
          if (wdayOfYear == dayOfYear) cdo_abort("Too less timesteps!");
          if (wdayOfYear < 1) wdayOfYear += (MaxDays - 1) * sumboot;
        }
    }
  int base = (dayOfYear - 1) * ndates + 1;
  wdays[base] = wdayOfYear;
  int nndates = 1;
  while (nndates != ndates)
    {
      wdayOfYear++;
      if (wdayOfYear > sumboot * (MaxDays - 1)) wdayOfYear -= sumboot * (MaxDays - 1);
      if (wdaysRead[wdayOfYear] != false)
        {
          wdays[base + nndates] = wdayOfYear;
          nndates++;
        }
    }
}

static void
writeTimesteps(int MaxMonths, int recentYear, FieldVector3D &cei, int frequency, int taxisID4, const CdoStreamID &streamID4,
               int *otsID, VarList const &varList1, std::vector<FieldInfo> const &fieldInfoList, std::vector<int> const &tempdpm,
               int tempdpy, int func2)
{
  const int maxRecs = fieldInfoList.size();

  if (frequency == 8)
    {
      for (int loopmonth = 1; loopmonth < MaxMonths + 1; loopmonth++)
        {
          define_mid_of_time(frequency, taxisID4, recentYear, loopmonth, MaxMonths);
          cdo_def_timestep(streamID4, *otsID);
          for (int fieldID = 0; fieldID < maxRecs; ++fieldID)
            {
              auto [varIDo, levelIDo] = fieldInfoList[fieldID].get();
              if (*otsID && varList1.vars[varIDo].isConstant) continue;

              fieldc_div(cei[loopmonth - 1][0][levelIDo], (double) (tempdpm[loopmonth - 1] / 100.0));
              cdo_def_field(streamID4, varIDo, levelIDo);
              cdo_write_field(streamID4, cei[loopmonth - 1][0][levelIDo]);
            }
          (*otsID)++;
        }
    }
  else
    {
      define_mid_of_time(frequency, taxisID4, recentYear, 0, MaxMonths);
      cdo_def_timestep(streamID4, *otsID);
      for (int fieldID = 0; fieldID < maxRecs; ++fieldID)
        {
          auto [varIDo, levelIDo] = fieldInfoList[fieldID].get();
          if (*otsID && varList1.vars[varIDo].isConstant) continue;

          if (func2 == FieldFunc_Avg) fieldc_div(cei[0][0][levelIDo], (double) (tempdpy / 100.0));
          cdo_def_field(streamID4, varIDo, levelIDo);
          cdo_write_field(streamID4, cei[0][0][levelIDo]);
        }
      (*otsID)++;
    }
}

static void
calculateOuterPeriod(Field &field, int MaxMonths, int recentYear, int endOfCalc, FieldVector3D &cei, FieldVector3D &varsPtemp,
                     int frequency, int taxisID4, const CdoStreamID &streamID4, int *otsID, VarList const &varList1,
                     std::vector<FieldInfo> const &fieldInfoList, int selection, int func2)
{
  if (cdo_assert_files_only() == false) cdo_abort("This operator can't be combined with other operators!");

  auto cdiStream = streamOpenRead(cdo_get_stream_name(0));

  auto cdiVlistID = streamInqVlist(cdiStream);
  auto cdiTaxisID = vlistInqTaxis(cdiVlistID);
  std::vector<int> tempdpm(MaxMonths);
  int tempdpy = 0;
  for (int i = 0; i < MaxMonths; ++i) tempdpm[i] = 0;
  int year, month, day, tsID = 0, varID, levelID;
  bool lHasStarted = false;

  if (Options::cdoVerbose) cdo_print("Start to process variables");

  while (true)
    {
      auto numFields = streamInqTimestep(cdiStream, tsID++);
      if (numFields == 0) break;

      auto vDateTime = taxisInqVdatetime(cdiTaxisID);
      cdiDate_decode(vDateTime.date, &year, &month, &day);
      if (!lHasStarted && year != recentYear)
        continue;
      else if (!lHasStarted)
        lHasStarted = true;

      if (year != recentYear)
        {
          writeTimesteps(MaxMonths, recentYear, cei, frequency, taxisID4, streamID4, otsID, varList1, fieldInfoList, tempdpm,
                         tempdpy, func2);
          recentYear = year;
          tempdpy = 0;
          tempdpm[0] = 0;
          field_fill(cei[0][0][0], 0.);
          if (frequency == 8)
            for (int loopmonth = 1; loopmonth < MaxMonths; loopmonth++)
              {
                tempdpm[loopmonth] = 0;
                field_fill(cei[loopmonth][0][0], 0.);
              }
        }
      if (year == endOfCalc && func2 == FieldFunc_Avg) break;
      tempdpy++;
      auto dayOfYear = decode_day_of_year(vDateTime.date);
      tempdpm[month - 1]++;

      if (func2 == FieldFunc_Sum) dayOfYear = 1;

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          streamInqField(cdiStream, &varID, &levelID);
          streamReadField(cdiStream, field.vec_d.data(), &field.numMissVals);

          Field &pctls = varsPtemp[dayOfYear][0][levelID];
          if (selection == func_selle)
            vfarselle(field, pctls);
          else if (selection == func_selge)
            vfarselge(field, pctls);

          auto &array = field.vec_d;
          if (func2 == FieldFunc_Avg)
            for (size_t i = 0; i < field.size; ++i) array[i] = (fp_is_equal(array[i], field.missval)) ? 0.0 : 1.0;
          else
            for (size_t i = 0; i < field.size; ++i) array[i] = (fp_is_equal(array[i], field.missval)) ? 0.0 : array[i];
          if (frequency == 8)
            field2_add(cei[(int) ((dayOfYear - 1) / 31.)][0][levelID], field);
          else
            field2_add(cei[0][0][levelID], field);
        }
    }
  if (Options::cdoVerbose) cdo_print("Finished Processing variables");
  if (year != endOfCalc)
    writeTimesteps(MaxMonths, year, cei, frequency, taxisID4, streamID4, otsID, varList1, fieldInfoList, tempdpm, tempdpy, func2);

  field_fill(cei[0][0][0], 0.);

  if (frequency == 8)
    for (int loopmonth = 1; loopmonth < MaxMonths; loopmonth++)
      {
        tempdpm[loopmonth] = 0;
        field_fill(cei[loopmonth][0][0], 0.);
      }

  streamClose(cdiStream);
}

void
etccdi_op(ETCCDI_REQUEST &request)
{
  constexpr int MaxDays = 373;
  constexpr int MaxMonths = 12;
  FieldVector2D varsData2[MaxDays];
  HistogramSet hsets[MaxDays];

  const int operatorID = cdo_operator_id();
  auto selection = cdo_operator_f1(operatorID);
  auto frequency = request.compare_type;

  percentile_set_method("rtype8");

  int FIELD_MEMTYPE = (Options::CDO_Memtype == MemType::Float) ? FIELD_FLT : 0;

  bool wdaysSrc[MaxDays];

  if (request.endboot < request.startboot)
    {
      cdo_warning("Your interval end '%d' is before the interval start '%d'. Switched interval years.", request.endboot,
                  request.startboot);
      request.startboot = request.endboot;
      request.endboot = request.startboot;
    }
  int sumboot = request.endboot - request.startboot + 1;
  std::vector<bool> wdaysRead((MaxDays - 1) * sumboot + 1);
  std::vector<int> wdays(request.ndates * (MaxDays - 1) * sumboot + 1);
  std::vector<int> dpy(sumboot), dpm(MaxMonths * sumboot);

  for (int year = 0; year < sumboot; year++)
    {
      dpy[year] = 0;
      for (int month = 0; month < MaxMonths; ++month) dpm[month + year * MaxMonths] = 0;
    }

  auto streamID1 = cdo_open_read(0);
  auto streamID2 = cdo_open_read(1);
  auto streamID3 = cdo_open_read(2);

  auto vlistID1 = cdo_stream_inq_vlist(streamID1);
  auto vlistID2 = cdo_stream_inq_vlist(streamID2);
  auto vlistID3 = cdo_stream_inq_vlist(streamID3);
  auto vlistID4 = vlistDuplicate(vlistID1);

  vlist_unpack(vlistID4);

  VarList varList1(vlistID1);
  VarList varList2(vlistID2);

  varList_compare(varList1, varList2);
  varList_compare(varList1, VarList(vlistID3));

  auto taxisID1 = vlistInqTaxis(vlistID1);
  auto taxisID2 = vlistInqTaxis(vlistID2);
  auto taxisID3 = vlistInqTaxis(vlistID3);
  // TODO - check that time axes 2 and 3 are equal

  auto taxisID4 = taxisDuplicate(taxisID1);
  if (taxisHasBounds(taxisID4)) taxisDeleteBounds(taxisID4);
  vlistDefTaxis(vlistID4, taxisID4);

  auto streamID4 = cdo_open_write(3);

  auto numVars = varList1.numVars();
  auto lOnlyOneVar = true;

  if (numVars == 1)
    {
      cdiDefKeyString(vlistID4, 0, CDI_KEY_NAME, request.name);
      cdiDefKeyString(vlistID4, 0, CDI_KEY_LONGNAME, request.longname);
      cdiDefKeyString(vlistID4, 0, CDI_KEY_UNITS, request.units);
      cdiDefAttTxt(vlistID4, 0, "cell_methods", (int) std::strlen("time: maximum"), "time: maximum");
    }
  else
    {
      lOnlyOneVar = false;
      cdo_warning("Your input file has more than one variable. No attributes can be set.");
    }

  cdo_def_vlist(streamID4, vlistID4);

  auto maxFields = varList1.maxFields();
  std::vector<FieldInfo> fieldInfoList(maxFields);

  auto gridsizeMax = varList1.gridsizeMax();

  Field field;

  FieldVector3D varsData1(sumboot * (MaxDays - 1) + 1);
  FieldVector3D cei(sumboot * MaxMonths);
  FieldVector3D varsPtemp(MaxDays);

  for (int dayOfYear = 0; dayOfYear < MaxDays; dayOfYear++)
    {
      field2D_init(varsData1[dayOfYear], varList1, FIELD_VEC | FIELD_MEMTYPE);
      wdaysSrc[dayOfYear] = false;
      for (int year = 0; year < sumboot; year++) wdaysRead[dayOfYear + year * (MaxDays - 1)] = false;
    }

  for (int dayOfYear = MaxDays; dayOfYear < sumboot * (MaxDays - 1) + 1; dayOfYear++)
    field2D_init(varsData1[dayOfYear], varList1, FIELD_VEC | FIELD_MEMTYPE);

  int tsID = 0;

  while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID2, tsID);
      if (numFields == 0) break;

      if (numFields != cdo_stream_inq_timestep(streamID3, tsID))
        cdo_abort("Number of fields at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(1),
                  cdo_get_stream_name(2));

      auto vDateTime2 = taxisInqVdatetime(taxisID2);
      auto vDateTime3 = taxisInqVdatetime(taxisID3);

      if (cdiDate_get(vDateTime2.date) != cdiDate_get(vDateTime3.date))
        cdo_abort("Verification dates at time step %d of %s and %s differ!", tsID + 1, cdo_get_stream_name(1),
                  cdo_get_stream_name(2));

      auto dayOfYear = decode_day_of_year(vDateTime2.date);

      if (request.func2 == FieldFunc_Sum) dayOfYear = 1;

      if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);

      if (!varsData2[dayOfYear].size())
        {
          wdaysSrc[dayOfYear] = true;
          field2D_init(varsData2[dayOfYear], varList2, FIELD_VEC | FIELD_MEMTYPE);
          field2D_init(varsPtemp[dayOfYear], varList2, FIELD_VEC);
          hsets[dayOfYear].create(numVars);

          for (auto const &var : varList1.vars) hsets[dayOfYear].createVarLevels(var.ID, var.nlevels, var.gridsize);
        }

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID2);
          cdo_read_field(streamID2, varsData2[dayOfYear][varID][levelID]);
          varsPtemp[dayOfYear][varID][levelID].numMissVals = varsData2[dayOfYear][varID][levelID].numMissVals;
        }

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
        {
          auto [varID, levelID] = cdo_inq_field(streamID3);
          field.init(varList1.vars[varID]);
          cdo_read_field(streamID3, field);

          hsets[dayOfYear].defVarLevelBounds(varID, levelID, varsData2[dayOfYear][varID][levelID], field);
        }
      // fieldsFree(vlistID2, vars2[dayOfYear]);
      // field2D_init(vars2[dayOfYear], varList2, FIELD_VEC);

      tsID++;
    }

  if (Options::cdoVerbose) cdo_print("Defined the boundaries for the histograms");
  tsID = 0;
  bool lOnlyRefPeriod = true;
  int firstYear = 0, lastYear = 0;
  int year;
  while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID++);
      if (numFields == 0) break;

      auto vDateTime = taxisInqVdatetime(taxisID1);

      int month, day;
      cdiDate_decode(vDateTime.date, &year, &month, &day);
      if (tsID == 1)
        {
          if (year > request.startboot)
            cdo_abort("The interval start year '%d' is before infile start year '%d'.", request.startboot, year);
          firstYear = year;
        }
      lastYear = year;

      if (year >= request.startboot && year <= request.endboot)
        {
          auto dayOfYear = decode_day_of_year(vDateTime.date);
          if (dayOfYear < 0 || dayOfYear >= MaxDays) cdo_abort("Day %d out of range!", dayOfYear);
          if (wdaysSrc[dayOfYear] || request.func2 == FieldFunc_Sum)
            {
              // Variable independent ?
              wdaysRead[dayOfYear + (year - request.startboot) * (MaxDays - 1)]
                  = dayOfYear + (year - request.startboot) * (MaxDays - 1);
              dpy[year - request.startboot]++;
              dpm[(year - request.startboot) * MaxMonths + (int) ((dayOfYear - 1) / 31.)]++;
              for (int fieldID = 0; fieldID < numFields; ++fieldID)
                {
                  auto [varID, levelID] = cdo_inq_field(streamID1);
                  cdo_read_field(streamID1, varsData1[dayOfYear + (year - request.startboot) * (MaxDays - 1)][varID][levelID]);

                  if (tsID == 0) fieldInfoList[fieldID].set(varID, levelID);
                }
            }
          else
            cdo_warning("Could not find histogram minimum or maximum for day of year: '%d'", dayOfYear);
        }
      else
        lOnlyRefPeriod = false;
    }
  if (Options::cdoVerbose) cdo_print("Read in variables");

  if (year < request.endboot) cdo_abort("The interval end year '%d' is after infile end year '%d'.", request.endboot, year);

  for (year = 0; year < sumboot; year++)
    {
      field2D_init(cei[year * MaxMonths], varList1, FIELD_VEC);
      if (frequency == 8)
        for (int month = 1; month < MaxMonths; ++month) field2D_init(cei[year * MaxMonths + month], varList1, FIELD_VEC);
    }

  // printf("Wir beginnen nun mit der Schleife.\n");
  int bootsyear = 0;
  int subyear = 0;
  int otsID = 0;

  for (int loopdoy = 1; loopdoy < MaxDays; loopdoy++)
    {
      for (int ytoadd = 0; ytoadd < sumboot; ytoadd++)
        {
          if (wdaysRead[loopdoy + ytoadd * (MaxDays - 1)])
            {
              windowDays(loopdoy + ytoadd * (MaxDays - 1), wdays, wdaysRead, MaxDays, request.ndates, sumboot);
            }
        }
    }
  if (Options::cdoVerbose) cdo_print("Calculated window days");

  for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var1 = varList1.vars[varID];
      if (var1.isConstant) continue;
      for (int levelID = 0; levelID < var1.nlevels; ++levelID)
        {
          if (request.func2 == FieldFunc_Sum)
            {
              for (int loopdoy = 1; loopdoy < MaxDays; loopdoy++)
                {
                  for (int ytoadd = 0; ytoadd < sumboot; ytoadd++)
                    {
                      if (wdaysRead[loopdoy + ytoadd * (MaxDays - 1)])
                        {
                          auto &source = varsData1[loopdoy + ytoadd * (MaxDays - 1)][varID][levelID];
                          auto &hset = hsets[1];
                          hset.addVarLevelValues(varID, levelID, source);
                        }
                    }
                }
            }
          else
            {
#ifdef _OPENMP
#pragma omp parallel for shared(sumboot, varsData1, request, varID, levelID, hsets, wdays, wdaysRead) schedule(dynamic)
#endif
              for (int loopdoy = 1; loopdoy < MaxDays; loopdoy++)
                {
                  for (int ytoadd = 0; ytoadd < sumboot; ytoadd++)
                    {
                      if (wdaysRead[loopdoy + ytoadd * (MaxDays - 1)])
                        {
                          for (int ano = 0; ano < request.ndates; ano++)
                            {
                              auto &source = varsData1[wdays[ytoadd * request.ndates * (MaxDays - 1)
                                                             + (loopdoy - 1) * request.ndates + ano + 1]][varID][levelID];
                              auto &hset = hsets[loopdoy];
                              hset.addVarLevelValues(varID, levelID, source);
                            }
                        }
                    }
                }
            }
        }
    }

  tsID = 0;
  if (Options::cdoVerbose) cdo_print("Added 30 years to histograms");
  if (lOnlyOneVar && ((!lOnlyRefPeriod && firstYear != request.startboot) || request.func2 == FieldFunc_Sum))
    {
      for (int varID = 0; varID < numVars; ++varID)
        {
          auto const &var1 = varList1.vars[varID];
          if (var1.isConstant) continue;
          for (int levelID = 0; levelID < var1.nlevels; ++levelID)
            {
              if (request.func2 == FieldFunc_Sum)
                {
                  auto &pctls = varsPtemp[1][varID][levelID];
                  hsets[1].getVarLevelPercentiles(pctls, varID, levelID, request.pn);
                }
              else
                {
#ifdef _OPENMP
#pragma omp parallel for shared(request, wdaysSrc, varID, levelID, hsets, varsPtemp) schedule(dynamic)
#endif
                  for (int loopdoy = 1; loopdoy < MaxDays; loopdoy++)
                    {
                      if (wdaysSrc[loopdoy])
                        {
                          auto &pctls = varsPtemp[loopdoy][varID][levelID];
                          hsets[loopdoy].getVarLevelPercentiles(pctls, varID, levelID, request.pn);
                        }
                    }
                }
            }
        }
      field.resize(gridsizeMax);
      calculateOuterPeriod(field, MaxMonths, firstYear, request.startboot, cei, varsPtemp, frequency, taxisID4, streamID4, &otsID,
                           varList1, fieldInfoList, selection, request.func2);
    }
  else if (!lOnlyRefPeriod && firstYear != request.startboot)
    cdo_warning("Since you have more than one variable in the input file, only the bootstrapping period can be calculated");

  if (request.func2 == FieldFunc_Avg)
    {
      for (bootsyear = request.startboot; bootsyear < request.endboot + 1; bootsyear++)
        {
          if (Options::cdoVerbose) cdo_print("Bootsyear: %d", bootsyear);
          for (int varID = 0; varID < numVars; ++varID)
            {
              if (varList1.vars[varID].isConstant) continue;
              for (int levelID = 0; levelID < varList1.vars[varID].nlevels; ++levelID)
                {
#ifdef _OPENMP
#pragma omp parallel for shared(sumboot, wdaysRead, request, varsData1, varID, levelID, hsets, wdays, cei) schedule(dynamic)
#endif
                  for (int loopdoy = 1; loopdoy < MaxDays; loopdoy++)
                    {
                      for (int ytoadd = 0; ytoadd < sumboot; ytoadd++)
                        {
                          if (wdaysRead[loopdoy + ytoadd * (MaxDays - 1)])
                            {
                              for (int ano = 0; ano < request.ndates; ano++)
                                {
                                  int recentWday
                                      = ytoadd * request.ndates * (MaxDays - 1) + (loopdoy - 1) * request.ndates + ano + 1;
                                  if ((int((wdays[recentWday] - 1) / (MaxDays - 1)) + request.startboot) == bootsyear)
                                    {
                                      auto &source = varsData1[wdays[recentWday]][varID][levelID];
                                      auto &hset = hsets[loopdoy];
                                      hset.subVarLevelValues(varID, levelID, source);
                                    }
                                  // percyear cannot be smaller than request.startboot
                                  if ((int((wdays[recentWday] - 1) / (MaxDays - 1)) + request.startboot) == bootsyear - 1)
                                    {
                                      auto &source = varsData1[wdays[recentWday]][varID][levelID];
                                      auto &hset = hsets[loopdoy];
                                      hset.addVarLevelValues(varID, levelID, source);
                                    }
                                }
                            }
                        }
                    }
                  for (subyear = request.startboot; subyear < request.endboot + 1; subyear++)
                    {
                      if (Options::cdoVerbose) cdo_print("Subyear: %d", subyear);
                      if (subyear != bootsyear)
                        {
#ifdef _OPENMP
#pragma omp parallel for shared(sumboot, request, varsData1, varID, levelID, hsets, wdaysRead, varsPtemp, varsData2, cei, subyear, \
                                    bootsyear, wdays, frequency) schedule(dynamic)
#endif
                          for (int loopdoy = 1; loopdoy < MaxDays; loopdoy++)
                            {
                              for (int ytoadd = 0; ytoadd < sumboot; ytoadd++)
                                {
                                  if (wdaysRead[loopdoy + ytoadd * (MaxDays - 1)])
                                    {
                                      for (int ano = 0; ano < request.ndates; ano++)
                                        {
                                          int recentWday
                                              = ytoadd * request.ndates * (MaxDays - 1) + (loopdoy - 1) * request.ndates + ano + 1;
                                          if ((int((wdays[recentWday] - 1) / (MaxDays - 1)) + request.startboot) == subyear)
                                            {
                                              auto &source = varsData1[wdays[recentWday]][varID][levelID];
                                              auto &hset = hsets[loopdoy];
                                              if (hset.addVarLevelValues(varID, levelID, source) == 1)
                                                cdo_print("'%d', '%d", loopdoy, wdays[recentWday]);
                                            }
                                        }
                                    }
                                }
                              // printf("Haben es zum temp array addiert.\n");

                              /*** Calculate percentile  ***/
                              if (wdaysRead[loopdoy + (bootsyear - request.startboot) * (MaxDays - 1)])
                                {
                                  auto &pctls = varsPtemp[loopdoy][varID][levelID];
                                  hsets[loopdoy].getVarLevelPercentiles(pctls, varID, levelID, request.pn);
                                  /*** Compare data with percentile ***/
                                  auto &source
                                      = varsData1[loopdoy + (bootsyear - request.startboot) * (MaxDays - 1)][varID][levelID];
                                  auto &toCompare = varsData2[loopdoy][varID][levelID];
                                  field_copy(source, toCompare);
                                  if (selection == func_selle)
                                    vfarselle(toCompare, pctls);
                                  else if (selection == func_selge)
                                    vfarselge(toCompare, pctls);
                                  if (request.func2 == FieldFunc_Avg)
                                    {
                                      auto &array = toCompare.vec_d;
                                      for (size_t i = 0; i < toCompare.size; ++i)
                                        array[i] = (fp_is_equal(array[i], toCompare.missval)) ? 0.0 : 1.0;
                                    }
                                  else
                                    {
                                      auto &array = toCompare.vec_d;
                                      for (size_t i = 0; i < toCompare.size; ++i)
                                        array[i] = (fp_is_equal(array[i], toCompare.missval)) ? 0.0 : array[i];
                                    }
                                  // printf("Haben ein Percentil berechnet.\n");
                                  // Year sum
                                  if (frequency == 8)
#ifdef _OPENMP
#pragma omp critical
#endif
                                    field2_add(cei[(bootsyear - request.startboot) * MaxMonths + (int) ((loopdoy - 1) / 31.)][varID]
                                                  [levelID],
                                               toCompare);
                                  else
#ifdef _OPENMP
#pragma omp critical
#endif
                                    field2_add(cei[(bootsyear - request.startboot) * MaxMonths][varID][levelID], toCompare);
                                }
                              for (int ytoadd = 0; ytoadd < sumboot; ytoadd++)
                                {
                                  if (wdaysRead[loopdoy + ytoadd * (MaxDays - 1)])
                                    {
                                      for (int ano = 0; ano < request.ndates; ano++)
                                        {
                                          int recentWday
                                              = ytoadd * request.ndates * (MaxDays - 1) + (loopdoy - 1) * request.ndates + ano + 1;
                                          if ((int((wdays[recentWday] - 1) / (MaxDays - 1)) + request.startboot) == subyear)
                                            {
                                              auto &source = varsData1[wdays[recentWday]][varID][levelID];
                                              auto &hset = hsets[loopdoy];
                                              if (hset.subVarLevelValues(varID, levelID, source) == 1)
                                                cdo_print("'%d', '%d", loopdoy, wdays[recentWday]);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                  if (frequency == 8)
                    {
                      for (int month = 0; month < MaxMonths; ++month)
                        if (cei[(bootsyear - request.startboot) * MaxMonths + month][varID][levelID].vec_d.data())
                          {
                            // Divide vars2 to receive average
                            fieldc_div(cei[(bootsyear - request.startboot) * MaxMonths + month][varID][levelID],
                                       (double) ((sumboot - 1) * dpm[(bootsyear - request.startboot) * MaxMonths + month] / 100.0));
                          }
                    }
                  else if (cei[(bootsyear - request.startboot) * MaxMonths][varID][levelID].vec_d.data())
                    {
                      fieldc_div(cei[(bootsyear - request.startboot) * MaxMonths][varID][levelID],
                                 (double) ((sumboot - 1) * dpy[bootsyear - request.startboot] / 100.0));
                    }
                }
            }
          if (frequency == 8)
            {
              for (int month = 1; month < MaxMonths + 1; ++month)
                {
                  define_mid_of_time(frequency, taxisID4, bootsyear, month, MaxMonths);
                  cdo_def_timestep(streamID4, otsID);

                  for (int fieldID = 0; fieldID < maxFields; ++fieldID)
                    {
                      auto [varIDo, levelIDo] = fieldInfoList[fieldID].get();
                      if (otsID && varList1.vars[varIDo].isConstant) continue;

                      cdo_def_field(streamID4, varIDo, levelIDo);
                      cdo_write_field(streamID4, cei[(bootsyear - request.startboot) * MaxMonths + (month - 1)][varIDo][levelIDo]);
                    }
                  otsID++;
                }
            }
          else
            {
              define_mid_of_time(frequency, taxisID4, bootsyear, 0, MaxMonths);
              cdo_def_timestep(streamID4, otsID);

              for (int fieldID = 0; fieldID < maxFields; ++fieldID)
                {
                  auto [varIDo, levelIDo] = fieldInfoList[fieldID].get();
                  if (otsID && varList1.vars[varIDo].isConstant) continue;

                  cdo_def_field(streamID4, varIDo, levelIDo);
                  cdo_write_field(streamID4, cei[(bootsyear - request.startboot) * MaxMonths][varIDo][levelIDo]);
                }
              otsID++;
            }
          // printf("Haben ein Mittel für Jahr '%d' berechnet.\n", bootsyear);
        }
    }
  if (!lOnlyRefPeriod && lOnlyOneVar && lastYear != request.endboot && request.func2 == FieldFunc_Avg)
    {
      field_fill(cei[0][0][0], 0.);
      if (frequency == 8)
        for (int loopmonth = 1; loopmonth < MaxMonths; loopmonth++) field_fill(cei[loopmonth][0][0], 0.);

      for (int varID = 0; varID < numVars; ++varID)
        {
          if (varList1.vars[varID].isConstant) continue;
          for (int levelID = 0; levelID < varList1.vars[varID].nlevels; ++levelID)
            {
#ifdef _OPENMP
#pragma omp parallel for shared(request, wdaysRead, varID, levelID, hsets, varsPtemp) schedule(dynamic)
#endif
              for (int loopdoy = 1; loopdoy < MaxDays; loopdoy++)
                {
                  for (int ytoadd = 0; ytoadd < sumboot; ytoadd++)
                    {
                      if (wdaysRead[loopdoy + ytoadd * (MaxDays - 1)])
                        {
                          for (int ano = 0; ano < request.ndates; ano++)
                            {
                              int recentWday = ytoadd * request.ndates * (MaxDays - 1) + (loopdoy - 1) * request.ndates + ano + 1;
                              // percyear cannot be smaller than request.startboot
                              if ((int((wdays[recentWday] - 1) / (MaxDays - 1)) + request.startboot) == bootsyear - 1)
                                {
                                  auto &source = varsData1[wdays[recentWday]][varID][levelID];
                                  auto &hset = hsets[loopdoy];
                                  hset.addVarLevelValues(varID, levelID, source);
                                }
                            }
                        }
                    }
                  if (wdaysSrc[loopdoy])
                    {
                      auto &pctls = varsPtemp[loopdoy][varID][levelID];
                      hsets[loopdoy].getVarLevelPercentiles(pctls, varID, levelID, request.pn);
                    }
                }
            }
        }
      field.resize(gridsizeMax);
      field.missval = varsData1[1][0][0].missval;
      field.size = varsData1[1][0][0].size;
      field.grid = varsData1[1][0][0].grid;
      calculateOuterPeriod(field, MaxMonths, request.endboot + 1, lastYear + 1, cei, varsPtemp, frequency, taxisID4, streamID4,
                           &otsID, varList1, fieldInfoList, selection, request.func2);
    }

  cdo_stream_close(streamID4);
  cdo_stream_close(streamID3);
  cdo_stream_close(streamID2);
  cdo_stream_close(streamID1);
}

class EcaEtccdi : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "EcaEtccdi",
    .operators = { { "etccdi_tx90p", func_selge, CMP_DATE, EcaEtccdiHelp },
                   { "etccdi_tx10p", func_selle, CMP_DATE, EcaEtccdiHelp },
                   { "etccdi_tn90p", func_selge, CMP_DATE, EcaEtccdiHelp },
                   { "etccdi_tn10p", func_selle, CMP_DATE, EcaEtccdiHelp },
                   { "etccdi_r95p", func_selge, CMP_DATE, EcaEtccdiHelp },
                   { "etccdi_r99p", func_selge, CMP_DATE, EcaEtccdiHelp },
                   { "etccdi", 0, CMP_DATE, EcaEtccdiHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 3, 1, FilesOnly },
  };
  inline static RegisterEntry<EcaEtccdi> registration = RegisterEntry<EcaEtccdi>(module);

  int ETCCDI_TX90P, ETCCDI_R99P, ETCCDI_R95P, ETCCDI_TX10P, ETCCDI_TN90P, ETCCDI_TN10P;
  ETCCDI_REQUEST request;

public:
  void
  init() override
  {
    if (cdo_operator_argc() < 3) cdo_abort("Too few arguments!");
    if (cdo_operator_argc() > 4) cdo_abort("Too many arguments!");

    ETCCDI_TX90P = module.get_id("etccdi_tx90p");
    ETCCDI_R99P = module.get_id("etccdi_r99p");
    ETCCDI_R95P = module.get_id("etccdi_r95p");
    ETCCDI_TX10P = module.get_id("etccdi_tx10p");
    ETCCDI_TN90P = module.get_id("etccdi_tn90p");
    ETCCDI_TN10P = module.get_id("etccdi_tn10p");

    request.ndates = parameter_to_int(cdo_operator_argv(0));
    request.startboot = parameter_to_int(cdo_operator_argv(1));

    auto operatorID = cdo_operator_id();
    request.compare_type = cdo_operator_f2(cdo_operator_id());
    if (cdo_operator_argc() == 4 && 'm' == cdo_operator_argv(3)[0]) { request.compare_type = CMP_MONTH; }

    if (operatorID == ETCCDI_TX90P || operatorID == ETCCDI_TN90P || operatorID == ETCCDI_R95P || operatorID == ETCCDI_R99P)
      {
        if (operatorID == ETCCDI_TX90P || operatorID == ETCCDI_TN90P)
          {
            if (cdo_operator_argc() < 3)
              cdo_abort("Operator requires at least 3 parameter values, you provided '%d'!", cdo_operator_argc());
            request.endboot = parameter_to_int(cdo_operator_argv(2));
            if (operatorID == ETCCDI_TX90P)
              {
                request.name = TX90P_NAME;
                request.longname = TX90P_LONGNAME;
                request.units = TX90P_UNITS;
                request.func2 = FieldFunc_Avg;
              }
            else if (operatorID == ETCCDI_TN90P)
              {
                request.name = TN90P_NAME;
                request.longname = TN90P_LONGNAME;
                request.units = TN90P_UNITS;
                request.func2 = FieldFunc_Avg;
              }
            request.pn = 90;
          }
        else
          {
            if (cdo_operator_argc() < 2)
              cdo_abort("Operator requires at least 2 parameter values, you provided '%d'!", cdo_operator_argc());
            request.ndates = 1;
            request.startboot = parameter_to_int(cdo_operator_argv(0));
            request.endboot = parameter_to_int(cdo_operator_argv(1));
            if (operatorID == ETCCDI_R95P)
              {
                request.name = R95P_NAME;
                request.longname = R95P_LONGNAME;
                request.units = R95P_UNITS;
                request.pn = 95;
                request.func2 = FieldFunc_Sum;
              }
            else if (operatorID == ETCCDI_R99P)
              {
                request.name = R99P_NAME;
                request.longname = R99P_LONGNAME;
                request.units = R99P_UNITS;
                request.pn = 99;
                request.func2 = FieldFunc_Sum;
              }
          }
      }
    else if (operatorID == ETCCDI_TX10P || operatorID == ETCCDI_TN10P)
      {
        if (cdo_operator_argc() < 3)
          cdo_abort("Operator requires at least 3 parameter values, you provided '%d'!", cdo_operator_argc());
        request.endboot = parameter_to_int(cdo_operator_argv(2));
        if (operatorID == ETCCDI_TX10P)
          {
            request.name = TX10P_NAME;
            request.longname = TX10P_LONGNAME;
            request.units = TX10P_UNITS;
            request.func2 = FieldFunc_Avg;
          }
        else
          {
            request.name = TN10P_NAME;
            request.longname = TN10P_LONGNAME;
            request.units = TN10P_UNITS;
            request.func2 = FieldFunc_Avg;
          }
        request.pn = 10;
      }
  }

  void
  run() override
  {
    etccdi_op(request);
    /*  else
        EcaEtccdi(-1, ndates, startboot, endboot); */
  }

  void
  close() override
  {
  }
};
