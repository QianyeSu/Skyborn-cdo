/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Vertstat   vertrange       Vertical range
      Vertstat   vertmin         Vertical minimum
      Vertstat   vertmax         Vertical maximum
      Vertstat   vertsum         Vertical sum
      Vertstat   vertint         Vertical integral
      Vertstat   vertmean        Vertical mean
      Vertstat   vertavg         Vertical average
      Vertstat   vertvar         Vertical variance
      Vertstat   vertvar1        Vertical variance [Normalize by (n-1)]
      Vertstat   vertstd         Vertical standard deviation
      Vertstat   vertstd1        Vertical standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_stepstat.h"
#include "process_int.h"
#include "cdo_zaxis.h"
#include "param_conversion.h"
#include "pmlist.h"
#include "cdi_lockedIO.h"
#include "field_functions.h"

#define IS_SURFACE_LEVEL(zaxisID) (zaxisInqType(zaxisID) == ZAXIS_SURFACE && zaxisInqSize(zaxisID) == 1)

int
get_surface_ID(int vlistID)
{
  int surfaceID = -1;

  auto numZaxes = vlistNumZaxis(vlistID);
  for (int index = 0; index < numZaxes; ++index)
  {
    auto zaxisID = vlistZaxis(vlistID, index);
    if (IS_SURFACE_LEVEL(zaxisID))
    {
      surfaceID = vlistZaxis(vlistID, index);
      break;
    }
  }

  if (surfaceID == -1) surfaceID = zaxis_from_name("surface");

  return surfaceID;
}

static void
set_surface_ID(int vlistID, int surfaceID)
{
  auto numZaxes = vlistNumZaxis(vlistID);
  for (int index = 0; index < numZaxes; ++index)
  {
    auto zaxisID = vlistZaxis(vlistID, index);
    if (zaxisID != surfaceID || !IS_SURFACE_LEVEL(zaxisID)) vlistChangeZaxisIndex(vlistID, index, surfaceID);
  }
}

static void
vertstat_get_parameter(bool &weights, bool &genbounds)
{
  auto numArgs = cdo_operator_argc();
  if (numArgs)
  {
    auto const &argList = cdo_get_oper_argv();

    KVList kvlist;
    kvlist.name = cdo_module_name();
    if (kvlist.parse_arguments(argList) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);
      auto const &value = kv.values[0];

      // clang-format off
      if      (key == "weights")   weights = parameter_to_bool(value);
      else if (key == "genbounds") genbounds = parameter_to_bool(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }
}

class Vertstat : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Vertstat",
    .operators = { { "vertrange", FieldFunc_Range, 0, VertstatHelp },
                   { "vertmin", FieldFunc_Min, 0, VertstatHelp },
                   { "vertmax", FieldFunc_Max, 0, VertstatHelp },
                   { "vertsum", FieldFunc_Sum, 0, VertstatHelp },
                   { "vertint", FieldFunc_Sum, 1, VertstatHelp },
                   { "vertmean", FieldFunc_Mean, 1, VertstatHelp },
                   { "vertavg", FieldFunc_Avg, 1, VertstatHelp },
                   { "vertstd", FieldFunc_Std, 1, VertstatHelp },
                   { "vertstd1", FieldFunc_Std1, 1, VertstatHelp },
                   { "vertvar", FieldFunc_Var, 1, VertstatHelp },
                   { "vertvar1", FieldFunc_Var1, 1, VertstatHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Vertstat> registration = RegisterEntry<Vertstat>();

private:
  struct VertInfo
  {
    int zaxisID = -1;
    int status = -1;
    int numLevels = 0;
    Varray<double> thickness;
    Varray<double> weights;
  };

  int VERTINT{};
  int operatorID{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID2{ CDI_UNDEFID };

  int numVars{};

  VarList varList1{};

  std::vector<VertInfo> vert{};

  bool needWeights{};
  cdo::StepStat1Dvars stepStat;

public:
  void
  init() override
  {
    VERTINT = module.get_id("vertint");

    operatorID = cdo_operator_id();
    auto operfunc = cdo_operator_f1(operatorID);
    needWeights = cdo_operator_f2(operatorID);

    stepStat.init(operfunc);

    // int applyWeights = lmean;

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistClearFlag(vlistID1);

    varList1 = VarList(vlistID1);

    numVars = varList1.numVars();
    for (int varID = 0; varID < numVars; ++varID) vlistDefFlag(vlistID1, varID, 0, true);

    vlistID2 = vlistCreate();
    cdo_vlist_copy_flag(vlistID2, vlistID1);
    vlistDefNtsteps(vlistID2, varList1.numSteps());

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    auto surfID = get_surface_ID(vlistID1);
    set_surface_ID(vlistID2, surfID);

    auto numZaxes = varList1.numZaxes();
    vert.resize(numZaxes);
    if (needWeights)
    {
      auto useweights = true;
      auto genbounds = false;
      vertstat_get_parameter(useweights, genbounds);

      if (!useweights)
      {
        genbounds = false;
        cdo_print("Using constant vertical weights!");
      }

      for (int index = 0; index < numZaxes; ++index)
      {
        auto zaxisID = vlistZaxis(vlistID1, index);
        auto nlev = zaxisInqSize(zaxisID);
        vert[index].numLevels = 0;
        vert[index].status = 0;
        vert[index].zaxisID = zaxisID;
        // if (nlev > 1)
        {
          vert[index].numLevels = nlev;
          vert[index].thickness.resize(nlev);
          vert[index].weights.resize(nlev);
          vert[index].status
              = get_layer_thickness(useweights, genbounds, index, zaxisID, nlev, vert[index].thickness, vert[index].weights);
        }
        if (!useweights) vert[index].status = 3;
      }
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    int VARS_MEMTYPE = stepStat.lminmax ? FIELD_NAT : 0;
    stepStat.alloc(varList1, VARS_MEMTYPE);
  }

  void
  run() override
  {
    Field field;

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      std::vector<bool> varsLevelInit(numVars, false);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);

        auto const &var = varList1.vars[varID];

        auto &rsamp1 = stepStat.samp(varID);
        auto &rvar1 = stepStat.var1(varID);
        auto &rvar2 = stepStat.var2(varID);

        rvar1.nsamp++;
        if (stepStat.lrange) rvar2.nsamp++;

        auto gridsize = var.gridsize;

        auto layerWeight = 1.0;
        auto layerThickness = 1.0;
        if (needWeights)
        {
          for (int index = 0, n = vert.size(); index < n; ++index)
            if (vert[index].zaxisID == var.zaxisID)
            {
              if (vert[index].status == 0 && tsID == 0 && levelID == 0 && var.nlevels > 1)
              {
                cdo_warning("Layer bounds not available, using constant vertical weights for variable %s!", var.name);
              }
              else
              {
                layerWeight = vert[index].weights[levelID];
                layerThickness = vert[index].thickness[levelID];
              }

              break;
            }
        }

        field.init(var);
        cdo_read_field(streamID1, field);

        if (varsLevelInit[varID] == false)
        {
          varsLevelInit[varID] = true;
          field_copy(field, rvar1);

          if (stepStat.lrange) field_copy(field, rvar2);

          if (operatorID == VERTINT && is_not_equal(layerThickness, 1.0)) fieldc_mul(rvar1, layerThickness);
          if (stepStat.lmean && is_not_equal(layerWeight, 1.0)) fieldc_mul(rvar1, layerWeight);

          if (stepStat.lvarstd)
          {
            if (is_not_equal(layerWeight, 1.0))
            {
              field2_moqw(rvar2, rvar1, layerWeight);
              fieldc_mul(rvar1, layerWeight);
            }
            else { field2_moq(rvar2, rvar1); }
          }

          if (rvar1.numMissVals || !rsamp1.empty() || needWeights)
          {
            if (rsamp1.empty()) rsamp1.resize(gridsize);

            for (size_t i = 0; i < gridsize; ++i)
              rsamp1.vec_d[i] = (fp_is_equal(rvar1.vec_d[i], rvar1.missval)) ? 0.0 : layerWeight;
          }
        }
        else
        {
          if (operatorID == VERTINT && is_not_equal(layerThickness, 1.0)) fieldc_mul(field, layerThickness);
          if (stepStat.lmean && is_not_equal(layerWeight, 1.0)) fieldc_mul(field, layerWeight);

          if (field.numMissVals || !rsamp1.empty())
          {
            if (rsamp1.empty()) rsamp1.resize(gridsize, rvar1.nsamp);

            auto func = [&](auto const &v1, auto &v2, std::decay_t<decltype(v1[0])> missval)
            {
              for (size_t i = 0; i < gridsize; ++i)
                if (fp_is_not_equal(v1[i], missval)) { v2[i] += layerWeight; }
            };
            field_operation2(func, field, rsamp1, rvar1.missval);
          }

          if (stepStat.lvarstd)
          {
            if (is_not_equal(layerWeight, 1.0))
            {
              field2_sumqw(rvar2, field, layerWeight);
              field2_sumw(rvar1, field, layerWeight);
            }
            else { field2_sumsumq(rvar1, rvar2, field); }
          }
          else if (stepStat.lrange) { field2_maxmin(rvar1, rvar2, field); }
          else { field2_function(rvar1, field, stepStat.operfunc); }
        }
      }

      for (int varID = 0; varID < numVars; ++varID)
      {
        auto numSets = stepStat.var1(varID).nsamp;
        if (numSets)
        {
          stepStat.process(varID, numSets);

          cdo_def_field(streamID2, varID, 0);
          cdo_write_field(streamID2, stepStat.var1(varID));
          stepStat.var1(varID).nsamp = 0;
        }
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    vlistDestroy(vlistID2);
  }
};
