#ifndef CDO_STEPSTAT_H
#define CDO_STEPSTAT_H

#include "process_int.h"
#include "field.h"
#include "field_functions.h"

namespace cdo
{
class StepStatBase
{
public:
  int operfunc{};
  bool lminmax{ false };
  bool lminidx{ false };
  bool lmaxidx{ false };
  bool lrange{ false };
  bool lmean{ false };
  bool lmeanavg{ false };
  bool lstd{ false };
  bool lvarstd{ false };
  double divisor{};

  void
  init(int _operfunc)
  {
    operfunc = _operfunc;
    lminmax = (operfunc == FieldFunc_Min || operfunc == FieldFunc_Max);
    lminidx = (operfunc == FieldFunc_Minidx);
    lmaxidx = (operfunc == FieldFunc_Maxidx);
    lrange = (operfunc == FieldFunc_Range);
    lmean = (operfunc == FieldFunc_Mean);
    lmeanavg = (operfunc == FieldFunc_Mean || operfunc == FieldFunc_Avg);
    lstd = (operfunc == FieldFunc_Std || operfunc == FieldFunc_Std1);
    lvarstd = (lstd || operfunc == FieldFunc_Var || operfunc == FieldFunc_Var1);
    divisor = (operfunc == FieldFunc_Std1 || operfunc == FieldFunc_Var1);
  }

  void
  add_field_kernel(Field const &field, Field &sampData, Field &varData1, Field &varData2, int numSets)
  {
    if (numSets == 0)
      {
        if (lminidx || lmaxidx)
          field_fill(varData1, 0.0);
        else
          field_copy(field, varData1);

        if (lrange || lminidx || lmaxidx) field_copy(field, varData2);

        if (lvarstd) field2_moq(varData2, varData1);

        if (field.numMissVals || !sampData.empty())
          {
            if (sampData.empty()) sampData.resize(varData1.size);
            field2_vinit(sampData, field);
          }
      }
    else
      {
        if (field.numMissVals || !sampData.empty())
          {
            if (sampData.empty()) sampData.resize(varData1.size, numSets);
            field2_vincr(sampData, field);
          }

        // clang-format off
        if      (lvarstd) field2_sumsumq(varData1, varData2, field);
        else if (lrange)  field2_maxmin(varData1, varData2, field);
        else if (lminidx) field2_minidx(varData1, varData2, field, numSets);
        else if (lmaxidx) field2_maxidx(varData1, varData2, field, numSets);
        else              field2_function(varData1, field, operfunc);
        // clang-format on
      }
  }

  void
  process_kernel(Field const &sampData, Field &varData1, Field const &varData2, int numSets)
  {
    auto field2_stdvar_func = lstd ? field2_std : field2_var;
    auto fieldc_stdvar_func = lstd ? fieldc_std : fieldc_var;

    if (lmeanavg)
      {
        if (!sampData.empty())
          field2_div(varData1, sampData);
        else
          fieldc_div(varData1, (double) numSets);
      }
    else if (lvarstd)
      {
        if (!sampData.empty())
          field2_stdvar_func(varData1, varData2, sampData, divisor);
        else
          fieldc_stdvar_func(varData1, varData2, numSets, divisor);
      }
    else if (lrange) { field2_sub(varData1, varData2); }
  }
};

class StepStat1Dvars : public StepStatBase
{
private:
  FieldVector sampsData;
  FieldVector varsData1;
  FieldVector varsData2;

public:
  void
  alloc(VarList const &varList, int VARS_MEMTYPE)
  {
    auto var2needed = (lvarstd || lrange || lminidx || lmaxidx);
    field1Dvars_init(sampsData, varList);
    field1Dvars_init(varsData1, varList, FIELD_VEC | VARS_MEMTYPE);
    field1Dvars_init(varsData2, varList, var2needed ? FIELD_VEC : 0);
  }

  Field &
  var1(int varID)
  {
    return varsData1[varID];
  }

  Field &
  var2(int varID)
  {
    return varsData2[varID];
  }

  Field &
  samp(int varID)
  {
    return sampsData[varID];
  }

  void
  process(int varID, int numSets)
  {
    process_kernel(sampsData[varID], varsData1[varID], varsData2[varID], numSets);
  }
};

class StepStat1Dlevels : public StepStatBase
{
private:
  FieldVector sampsData;
  FieldVector varsData1;
  FieldVector varsData2;

public:
  void
  alloc(VarList const &varList, int VARS_MEMTYPE)
  {
    auto var2needed = (lvarstd || lrange || lminidx || lmaxidx);
    field1Dlevels_init(sampsData, varList);
    field1Dlevels_init(varsData1, varList, FIELD_VEC | VARS_MEMTYPE);
    field1Dlevels_init(varsData2, varList, var2needed ? FIELD_VEC : 0);
  }

  Field &
  var1(int levelID)
  {
    return varsData1[levelID];
  }

  void
  add_field(Field const &field, int levelID, int numSets)
  {
    auto &sampData = sampsData[levelID];
    auto &varData1 = varsData1[levelID];
    auto &varData2 = varsData2[levelID];

    varData1.nsamp++;
    if (lrange) varData2.nsamp++;
    add_field_kernel(field, sampData, varData1, varData2, numSets);
  }

  void
  moq(int levelID)
  {
    field2_moq(varsData2[levelID], varsData1[levelID]);
  }

  void
  process(int levelID, int numSets)
  {
    process_kernel(sampsData[levelID], varsData1[levelID], varsData2[levelID], numSets);
  }
};

class StepStat2D : public StepStatBase
{
private:
  Varray<double> vsamp;
  FieldVector2D sampsData;
  FieldVector2D varsData1;
  FieldVector2D varsData2;

  static void
  set_missval(Field &field, Field const &sampData, int numSets, double vfraction)
  {
    auto fieldsize = field.size;
    auto missval = field.missval;

    size_t irun = 0;
    for (size_t i = 0; i < fieldsize; ++i)
      {
        if ((sampData.vec_d[i] / numSets) < vfraction)
          {
            field.vec_d[i] = missval;
            irun++;
          }
      }

    if (irun) field_num_mv(field);
  }

public:
  void
  alloc(VarList const &varList, int VARS_MEMTYPE)
  {
    auto var2needed = (lvarstd || lrange || lminidx || lmaxidx);
    field2D_init(sampsData, varList);
    field2D_init(varsData1, varList, FIELD_VEC | VARS_MEMTYPE);
    field2D_init(varsData2, varList, var2needed ? FIELD_VEC : 0);
  }

  Field &
  var1(int varID, int levelID)
  {
    return varsData1[varID][levelID];
  }

  Varray<double> &
  samp(int varID, int levelID, int numSets)
  {
    auto const &sampData = sampsData[varID][levelID];
    auto const &varData1 = varsData1[varID][levelID];

    vsamp.resize(varData1.size);
    if (!sampData.empty())
      vsamp = sampData.vec_d;
    else
      std::ranges::fill(vsamp, (double) numSets);

    return vsamp;
  }

  void
  add_field(Field const &field, int varID, int levelID, int numSets)
  {
    auto &sampData = sampsData[varID][levelID];
    auto &varData1 = varsData1[varID][levelID];
    auto &varData2 = varsData2[varID][levelID];

    add_field_kernel(field, sampData, varData1, varData2, numSets);
  }

  void
  set_missval(int varID, int levelID, int numSets, double vfraction)
  {
    auto const &sampData = sampsData[varID][levelID];
    if (!sampData.empty()) set_missval(varsData2[varID][levelID], sampData, numSets, vfraction);
  }

  void
  process(int varID, int levelID, int numSets)
  {
    process_kernel(sampsData[varID][levelID], varsData1[varID][levelID], varsData2[varID][levelID], numSets);
  }
};

class StepStat3D : public StepStatBase
{
private:
  FieldVector3D sampsData;
  FieldVector3D varsData1;
  FieldVector3D varsData2;
  int m_dimlen0{ 0 };

public:
  void
  set_dimlen0(int dimlen0)
  {
    m_dimlen0 = dimlen0;
    sampsData.resize(dimlen0);
    varsData1.resize(dimlen0);
    varsData2.resize(dimlen0);
  }

  void
  alloc(int dim0, VarList const &varList, int VARS_MEMTYPE)
  {
    auto var2needed = (lvarstd || lrange || lminidx || lmaxidx);
    field2D_init(sampsData[dim0], varList);
    field2D_init(varsData1[dim0], varList, FIELD_VEC | VARS_MEMTYPE);
    field2D_init(varsData2[dim0], varList, var2needed ? FIELD_VEC : 0);
  }

  FieldVector2D &
  samp(int dim0)
  {
    return sampsData[dim0];
  }

  Field &
  samp(int dim0, int varID, int levelID)
  {
    return sampsData[dim0][varID][levelID];
  }

  FieldVector2D &
  var1(int dim0)
  {
    return varsData1[dim0];
  }

  Field &
  var1(int dim0, int varID, int levelID)
  {
    return varsData1[dim0][varID][levelID];
  }

  FieldVector2D &
  var2(int dim0)
  {
    return varsData2[dim0];
  }

  Field &
  var2(int dim0, int varID, int levelID)
  {
    return varsData2[dim0][varID][levelID];
  }

  void
  add_field(Field const &field, int dim0, int varID, int levelID, int numSets)
  {
    auto &sampData = sampsData[dim0][varID][levelID];
    auto &varData1 = varsData1[dim0][varID][levelID];
    auto &varData2 = varsData2[dim0][varID][levelID];

    add_field_kernel(field, sampData, varData1, varData2, numSets);
  }

  void
  process(int dim0, int varID, int levelID, int numSets)
  {
    process_kernel(sampsData[dim0][varID][levelID], varsData1[dim0][varID][levelID], varsData2[dim0][varID][levelID], numSets);
  }
};

const auto write_out_stream = [](CdoStreamID streamID2, std::vector<FieldInfo> const &fieldInfoList, VarList const &varList1,
                                 cdo::StepStat2D &stepStat, int otsID) noexcept {
  cdo_def_timestep(streamID2, otsID);

  for (auto const &fieldInfo : fieldInfoList)
    {
      auto [varID, levelID] = fieldInfo.get();
      if (otsID && varList1.vars[varID].isConstant) continue;

      cdo_def_field(streamID2, varID, levelID);
      cdo_write_field(streamID2, stepStat.var1(varID, levelID));
    }
};

const auto write_diag_stream = [](CdoStreamID streamID3, std::vector<FieldInfo> const &fieldInfoList, VarList const &varList1,
                                  cdo::StepStat2D &stepStat, int otsID, int numSets) noexcept {
  cdo_def_timestep(streamID3, otsID);

  for (auto const &fieldInfo : fieldInfoList)
    {
      auto [varID, levelID] = fieldInfo.get();
      if (otsID && varList1.vars[varID].isConstant) continue;

      auto &vsamp = stepStat.samp(varID, levelID, numSets);

      cdo_def_field(streamID3, varID, levelID);
      cdo_write_field(streamID3, vsamp.data(), 0);
    }
};

const auto fields_process
    = [](std::vector<FieldInfo> const &fieldInfoList, VarList const &varList1, cdo::StepStat2D &stepStat, int numSets) noexcept {
        for (auto const &fieldInfo : fieldInfoList)
          {
            auto [varID, levelID] = fieldInfo.get();
            if (varList1.vars[varID].isConstant) continue;

            stepStat.process(varID, levelID, numSets);
          }
      };

const auto fields_set_missval = [](std::vector<FieldInfo> const &fieldInfoList, VarList const &varList1, cdo::StepStat2D &stepStat,
                                   int numSets, double vfraction) noexcept {
  for (auto const &fieldInfo : fieldInfoList)
    {
      auto [varID, levelID] = fieldInfo.get();
      if (varList1.vars[varID].isConstant) continue;

      stepStat.set_missval(varID, levelID, numSets, vfraction);
    }
};

const auto fields_process_3D = [](int dim0, std::vector<FieldInfo> const &fieldInfoList, VarList const &varList1,
                                  cdo::StepStat3D &stepStat, int numSets) noexcept {
  for (auto const &fieldInfo : fieldInfoList)
    {
      auto [varID, levelID] = fieldInfo.get();
      if (varList1.vars[varID].isConstant) continue;

      stepStat.process(dim0, varID, levelID, numSets);
    }
};

};  // namespace cdo

#endif
