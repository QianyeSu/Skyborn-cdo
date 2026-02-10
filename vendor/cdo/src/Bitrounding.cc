/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_options.h"
#include "pmlist.h"
#include "param_conversion.h"
#include "progress.h"
#include "bitinformation.h"

namespace
{
struct Parameter
{
  double infLevel = 0.9999;
  int minBits = 1;
  int maxBits = 23;
  int addBits = 0;
  int numBits = -1;
  int numSteps = -1;
  std::string filename{};
  bool printBits = false;
};

struct VarStat
{
  int nsbMin = 1000;
  int nsbMax = -1000;
};
}  // namespace

static int
get_keepbits(const MutualInformation &bitInfo, double inflevel)
{
  // xbitinfo::get_keepbits v0.0.1 (https://github.com/observingClouds/xbitinfo)
  // Converted from Python to C++ by Uwe Schulzweida

  constexpr int floatNMBITS = 9;  // number of non mantissa bits for float
  int keepMantissaBits = 23;

  double bitInfoMax = -9.e33;
  for (int i = 0; i < NBITS; ++i) bitInfoMax = std::max(bitInfoMax, bitInfo.M[i]);
  // printf("bitInfoMax %g\n", bitInfoMax);

  double bitInfoMaxLast4 = -9.e33;
  for (int i = NBITS - 4; i < NBITS; ++i) bitInfoMaxLast4 = std::max(bitInfoMaxLast4, bitInfo.M[i]);
  // printf("bitInfoMax/bitInfoMaxLast4 %g\n", bitInfoMax/bitInfoMaxLast4);
  bitInfoMaxLast4 *= 1.5;
  // printf("bitInfoMaxLast4 %g\n", bitInfoMaxLast4);

  MutualInformation infoPerBitCleaned;
  for (int i = 0; i < NBITS; ++i) infoPerBitCleaned.M[i] = (bitInfo.M[i] > bitInfoMaxLast4) ? bitInfo.M[i] : 0.0;
  // for (int i = 0; i < NBITS; ++i) printf("cleaned[%d] %g\n", i + 1, infoPerBitCleaned.M[i]);

  for (int i = 1; i < NBITS; ++i) infoPerBitCleaned.M[i] += infoPerBitCleaned.M[i - 1];
  // for (int i = 0; i < NBITS; ++i) printf("cumsum[%d] %g\n", i + 1, infoPerBitCleaned.M[i]);

  auto lastBit = infoPerBitCleaned.M[NBITS - 1];
  if (lastBit > 0.0)
  {
    MutualInformation cdf;
    for (int i = 0; i < NBITS; ++i) cdf.M[i] = infoPerBitCleaned.M[i] / lastBit;
    // for (int i = 0; i < NBITS; ++i) printf("cdf[%d] %g\n", i + 1, infoPerBitCleaned.M[i]);

    constexpr int nonMantissaBits = floatNMBITS;

    for (int i = 0; i < NBITS; ++i)
      if (cdf.M[i] > inflevel)
      {
        keepMantissaBits = i + 1 - nonMantissaBits;
        break;
      }
  }

  // printf("keepMantissaBits: %d\n", keepMantissaBits);

  int nsb = std::clamp(keepMantissaBits, 1, 23);

  return nsb;
}

static int
bit_rounding(size_t len, Varray<float> v, double infLevel)  // copy v!
{
  signed_exponent(v.data(), len);

  auto bitInfo = bitinformation(v.data(), len);
  // if (Options::cdoVerbose) for (int i = 0; i < NBITS; ++i) fprintf(stderr, "bitInfo[%d] %.8e %g\n", i+1, bitInfo.M[i],
  // bitInfo.M[i]);

  return get_keepbits(bitInfo, infLevel);
}

static void
bitround(int nsb, size_t len, Varray<float> &v, float missval)
{
  // BitRound from NetCDF 4.9.0; routine nv4var.c

  constexpr uint32_t BIT_XPL_NBR_SGN_FLT = 23;

  // BitRound interprets nsb as number of significant binary digits (bits)
  uint32_t prc_bnr_xpl_rqr = nsb;

  uint32_t bit_xpl_nbr_zro = BIT_XPL_NBR_SGN_FLT - prc_bnr_xpl_rqr;

  // Create mask
  uint32_t msk_f32_u32_zro = 0u;       // Zero all bits
  msk_f32_u32_zro = ~msk_f32_u32_zro;  // Turn all bits to ones

  // BitShave mask for AND: Left shift zeros into bits to be rounded, leave ones in untouched bits.
  msk_f32_u32_zro <<= bit_xpl_nbr_zro;

  // BitSet mask for OR: Put ones into bits to be set, zeros in untouched bits.
  uint32_t msk_f32_u32_one = ~msk_f32_u32_zro;

  // BitRound mask for ADD: Set one bit: the MSB of LSBs
  uint32_t msk_f32_u32_hshv = msk_f32_u32_one & (msk_f32_u32_zro >> 1);

  // BitRound: Quantize to user-specified NSB with IEEE-rounding
  uint32_t *u32_ptr = (uint32_t *) v.data();

  for (size_t idx = 0; idx < len; idx++)
  {
    if (is_not_equal(v[idx], missval))
    {
      u32_ptr[idx] += msk_f32_u32_hshv;  // Add 1 to the MSB of LSBs, carry 1 to mantissa or even exponent
      u32_ptr[idx] &= msk_f32_u32_zro;   // Shave it
    }
  }
}

static void
check_range(double value, double minVal, double maxVal, std::string const &key)
{
  if (value < minVal || value > maxVal) cdo_abort("Parameter %s=%g out of range (min=%g/max=%g)!", key, value, minVal, maxVal);
}

static void
check_range(int value, int minVal, int maxVal, std::string const &key)
{
  if (value < minVal || value > maxVal) cdo_abort("Parameter %s=%d out of range (min=%d/max=%d)!", key, value, minVal, maxVal);
}

static Parameter
get_parameter()
{
  Parameter params;

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
      if      (key == "inflevel")  check_range(params.infLevel = parameter_to_double(value), 0.0, 1.0, key);
      else if (key == "minbits")   check_range(params.minBits = parameter_to_int(value), 1, 23, key);
      else if (key == "maxbits")   check_range(params.maxBits = parameter_to_int(value), 1, 23, key);
      else if (key == "addbits")   check_range(params.addBits = parameter_to_int(value), 0, 22, key);
      else if (key == "numbits")   check_range(params.numBits = parameter_to_int(value), 1, 23, key);
      else if (key == "numsteps")  check_range(params.numSteps = parameter_to_int(value), 1, 1, key);
      else if (key == "printbits") params.printBits = parameter_to_bool(value);
      else if (key == "filename")  params.filename = parameter_to_word(value);
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }
  }

  return params;
}

static void
print_parameter(Parameter const &params)
{
  std::stringstream outbuffer;

  outbuffer << "inflevel=" << params.infLevel;
  outbuffer << ", minbits=" << params.minBits;
  outbuffer << ", maxbits=" << params.maxBits;
  outbuffer << ", addbits=" << params.addBits;
  outbuffer << ", numbits=" << params.numBits;
  outbuffer << ", numsteps=" << params.numSteps;
  outbuffer << ", printbits=" << params.printBits;
  outbuffer << ", filename=" << params.filename;

  cdo_verbose("%s", outbuffer.str());
}

static void
check_attributes(int vlistID)
{
  int numBits = -1;
  auto status1 = cdiInqAttInt(vlistID, CDI_GLOBAL, "cdo_bitrounding_numbits", 1, &numBits);
  double infLevel = -1.0;
  auto status2 = cdiInqAttFlt(vlistID, CDI_GLOBAL, "cdo_bitrounding_inflevel", 1, &infLevel);
  char filename[2];
  auto status3 = cdiInqAttTxt(vlistID, CDI_GLOBAL, "cdo_bitrounding_filename", 1, filename);

  if ((status1 == 0 && numBits != -1) || (status2 == 0 && infLevel > 0.0) || status3 == 0)
    cdo_warning("It looks like CDO bitrounding has been applied to the input data before!");
}

static void
set_local_attributes(int vlistID, int varID, int numBits)
{
  cdiDefAttInt(vlistID, varID, "_QuantizeBitRoundNumberOfSignificantBits", CDI_DATATYPE_INT32, 1, &numBits);
}

static void
set_global_attributes(int vlistID, const Parameter &params, int numVarsHaveNumbits)
{
  if (params.filename.size() && numVarsHaveNumbits > 0)
    cdiDefAttTxt(vlistID, CDI_GLOBAL, "cdo_bitrounding_filename", (int) params.filename.size(), params.filename.c_str());

  if (numVarsHaveNumbits == vlistNvars(vlistID)) return;

  if (params.numBits != -1)
  {
    cdiDefAttInt(vlistID, CDI_GLOBAL, "cdo_bitrounding_numbits", CDI_DATATYPE_INT32, 1, &params.numBits);
  }
  else
  {
    cdiDefAttFlt(vlistID, CDI_GLOBAL, "cdo_bitrounding_inflevel", CDI_DATATYPE_FLT64, 1, &params.infLevel);
    if (params.addBits) cdiDefAttInt(vlistID, CDI_GLOBAL, "cdo_bitrounding_addbits", CDI_DATATYPE_INT32, 1, &params.addBits);
    if (params.minBits > 1) cdiDefAttInt(vlistID, CDI_GLOBAL, "cdo_bitrounding_minbits", CDI_DATATYPE_INT32, 1, &params.minBits);
    if (params.maxBits < 23) cdiDefAttInt(vlistID, CDI_GLOBAL, "cdo_bitrounding_maxbits", CDI_DATATYPE_INT32, 1, &params.maxBits);
    if (params.numSteps != -1)
      cdiDefAttInt(vlistID, CDI_GLOBAL, "cdo_bitrounding_numsteps", CDI_DATATYPE_INT32, 1, &params.numSteps);
  }
}

static std::vector<int>
get_vars_numbits(VarList const &varList, std::string const &filename)
{
  auto numVars = varList.numVars();
  std::vector<int> varsNumbits(numVars, -1);

  if (filename.size())
  {
    auto fp = std::fopen(filename.c_str(), "r");
    if (fp == nullptr) cdo_abort("Open failed on: %s\n", filename);
    PMList pmlist;
    pmlist.read_namelist(fp, filename.c_str());
    auto &kvlist = pmlist.front();
    std::fclose(fp);
    if (Options::cdoVerbose) kvlist.print();

    for (auto const &kv : kvlist)
    {
      auto const &key = kv.key;
      if (kv.nvalues > 1) cdo_abort("Too many values for parameter key >%s<!", key);
      if (kv.nvalues < 1) cdo_abort("Missing value for parameter key >%s<!", key);

      for (auto const &var : varList.vars)
      {
        if (key == var.name)
        {
          auto const &value = kv.values[0];
          auto numBits = parameter_to_int(value);
          check_range(numBits, 1, 23, key);
          varsNumbits[var.ID] = numBits;
        }
      }
    }
  }

  return varsNumbits;
}

static int
num_vars_have_numbits(std::vector<int> const &varsNumbits)
{
  int numVarsHaveNumbits = 0;
  for (size_t i = 0, n = varsNumbits.size(); i < n; ++i)
    if (varsNumbits[i] != -1) numVarsHaveNumbits++;

  return numVarsHaveNumbits;
}

class Bitrounding : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Bitrounding",
    .operators = { { "bitrounding", BitroundingHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Bitrounding> registration = RegisterEntry<Bitrounding>(module);

  CdoStreamID streamID1;
  int taxisID1{ CDI_UNDEFID };

  CdoStreamID streamID2;
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  Parameter params;

  VarList varList1;

  std::vector<VarStat> varsStatGlob;
  std::vector<bool> varsCheckMiss;
  std::vector<bool> varsCheckFloat;
  std::vector<int> varsNumbits;
  std::vector<std::vector<int>> nsbVarLevels;

public:
  void
  init() override
  {
    params = get_parameter();
    if (Options::cdoVerbose) print_parameter(params);

    streamID1 = cdo_open_read(0);
    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    taxisID1 = vlistInqTaxis(vlistID1);

    varList1 = VarList(vlistID1);

    varsNumbits = get_vars_numbits(varList1, params.filename);
    auto numVarsHaveNumbits = num_vars_have_numbits(varsNumbits);

    vlistID2 = vlistDuplicate(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    check_attributes(vlistID1);
    set_global_attributes(vlistID2, params, numVarsHaveNumbits);

    auto numVars = varList1.numVars();
    for (auto const &var : varList1.vars)
    {
      if (var.memType == MemType::Float)
      {
        int nsb = (varsNumbits[var.ID] != -1) ? varsNumbits[var.ID] : params.numBits;
        if (nsb >= 1 && nsb <= 23) set_local_attributes(vlistID2, var.ID, nsb);
      }
    }

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varsStatGlob = std::vector<VarStat>(numVars);
    varsCheckMiss = std::vector<bool>(numVars, true);
    varsCheckFloat = std::vector<bool>(numVars, true);

    nsbVarLevels = std::vector<std::vector<int>>(numVars);
    for (int varID = 0; varID < numVars; ++varID) nsbVarLevels[varID].resize(varList1.vars[varID].nlevels, 0);
  }

  void
  run() override
  {
    Field field;
    auto numVars = varList1.numVars();
    auto numSteps = varList1.numSteps();
    cdo::Progress progress(get_id());

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      std::vector<VarStat> varsStat(numVars);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto fstatus = (tsID + (fieldID + 1.0) / numFields) / numSteps;
        if (numSteps > 0) progress.update(fstatus);

        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_def_field(streamID2, varID, levelID);

        auto const &var = varList1.vars[varID];
        field.init(var);
        cdo_read_field(streamID1, field);

        if (field.memType == MemType::Double)
        {
          if (varsCheckFloat[varID])
          {
            varsCheckFloat[varID] = false;
            cdo_warning("64-bit floats unsupported, bitrounding disabled for %s!", var.name);
          }
        }
        else if (field.memType == MemType::Float)
        {
          int nsb = (varsNumbits[varID] != -1) ? varsNumbits[varID] : params.numBits;

          if (field.numMissVals == 0)
          {
            if (nsb == -1 && (tsID == 0 || params.numSteps == -1))
            {
              nsb = bit_rounding(field.size, field.vec_f, params.infLevel);
              // printf("nsb=%d\n", nsb);
              if (params.addBits) nsb += params.addBits;
              nsb = std::clamp(nsb, params.minBits, params.maxBits);
            }

            if (tsID == 0) { nsbVarLevels[varID][levelID] = nsb; }
            else if (params.numSteps == 1) { nsb = nsbVarLevels[varID][levelID]; }
          }

          auto &varStat = varsStat[varID];
          varStat.nsbMin = std::min(varStat.nsbMin, nsb);
          varStat.nsbMax = std::max(varStat.nsbMax, nsb);

          if (nsb >= 1 && nsb <= 23) bitround(nsb, field.size, field.vec_f, var.missval);

          if (nsb == -1 && field.numMissVals > 0 && varsCheckMiss[varID])
          {
            varsCheckMiss[varID] = false;
            cdo_warning("Missing values unsupported, bitrounding disabled for %s!", var.name);
          }
        }

        cdo_write_field(streamID2, field);
      }

      if (Options::cdoVerbose && params.numBits == -1)
      {
        fprintf(stderr, "NSB: step=%d:", tsID + 1);
        for (auto const &var1 : varList1.vars)
        {
          auto const &varStat = varsStat[var1.ID];
          if (varStat.nsbMin >= 1 && varStat.nsbMin <= 23)
          {
            fprintf(stderr, " %s=%d", var1.name.c_str(), varStat.nsbMin);
            if (var1.nlevels > 1) fprintf(stderr, "-%d ", varStat.nsbMax);
          }
        }
        fprintf(stderr, "\n");
      }

      for (int varID = 0; varID < numVars; ++varID)
      {
        auto &varStatGlob = varsStatGlob[varID];
        varStatGlob.nsbMin = std::min(varStatGlob.nsbMin, varsStat[varID].nsbMin);
        varStatGlob.nsbMax = std::max(varStatGlob.nsbMax, varsStat[varID].nsbMax);
      }

      if (params.printBits) break;

      tsID++;
    }

    if (params.printBits)
    {
      for (auto const &var1 : varList1.vars)
      {
        auto const &varStatGlob = varsStatGlob[var1.ID];
        if (varStatGlob.nsbMin >= 1 && varStatGlob.nsbMin <= 23) fprintf(stdout, "%s=%d\n", var1.name.c_str(), varStatGlob.nsbMax);
      }
    }
    else if (Options::cdoVerbose && params.numBits == -1)
    {
      fprintf(stderr, "NSB: step=all:");
      for (auto const &var1 : varList1.vars)
      {
        auto const &varStatGlob = varsStatGlob[var1.ID];
        if (varStatGlob.nsbMin >= 1 && varStatGlob.nsbMin <= 23)
        {
          fprintf(stderr, " %s=%d", var1.name.c_str(), varStatGlob.nsbMin);
          if (varStatGlob.nsbMin != varStatGlob.nsbMax) fprintf(stderr, "-%d", varStatGlob.nsbMax);
        }
        fprintf(stderr, "\n");
      }
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);

    vlistDestroy(vlistID2);
  }
};
