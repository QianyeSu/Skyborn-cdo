/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "arithmetic.h"
#include "process_int.h"

class Complextorect : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Complextorect",
    .operators = { { "complextorect" }, { "complextopol" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_COMP,  // Allowed number type
    .constraints = { 1, 2, OnlyFirst },
  };
  inline static RegisterEntry<Complextorect> registration = RegisterEntry<Complextorect>();

private:
  int COMPLEXTORECT{}, COMPLEXTOPOL{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};
  CdoStreamID streamID3{};
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int taxisID3{};
  int vlistID2{ CDI_UNDEFID };
  int vlistID3{};

  int operatorID{};
  VarList varList1;

public:
  void
  init() override
  {
    COMPLEXTORECT = module.get_id("complextorect");
    COMPLEXTOPOL = module.get_id("complextopol");

    operatorID = cdo_operator_id();

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    int vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);
    vlistID3 = vlistDuplicate(vlistID1);

    auto nvars = vlistNvars(vlistID2);
    for (int varID = 0; varID < nvars; ++varID)
    {
      auto datatype = vlistInqVarDatatype(vlistID2, varID);
      datatype = (datatype == CDI_DATATYPE_CPX64) ? CDI_DATATYPE_FLT64 : CDI_DATATYPE_FLT32;
      vlistDefVarDatatype(vlistID2, varID, datatype);
      vlistDefVarDatatype(vlistID3, varID, datatype);
    }

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    taxisID3 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);
    vlistDefTaxis(vlistID3, taxisID3);

    streamID2 = cdo_open_write(1);
    streamID3 = cdo_open_write(2);

    cdo_def_vlist(streamID2, vlistID2);
    cdo_def_vlist(streamID3, vlistID3);
    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    auto gridsizeMax = varList1.gridsizeMax();
    Varray<double> array1(2 * gridsizeMax);
    Varray<double> array2(gridsizeMax);
    Varray<double> array3(gridsizeMax);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_taxis_copy_timestep(taxisID3, taxisID1);

      cdo_def_timestep(streamID2, tsID);
      cdo_def_timestep(streamID3, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_def_field(streamID2, varID, levelID);
        cdo_def_field(streamID3, varID, levelID);

        auto const &var1 = varList1.vars[varID];
        auto gridsize = var1.gridsize;

        size_t numMissVals;
        cdo_read_field(streamID1, array1.data(), &numMissVals);

        auto is_EQ = fp_is_equal;
        auto missval1 = var1.missval;
        auto missval2 = missval1;

        if (operatorID == COMPLEXTORECT)
        {
          for (size_t i = 0; i < gridsize; ++i)
          {
            array2[i] = array1[2 * i];
            array3[i] = array1[2 * i + 1];
          }
        }
        else if (operatorID == COMPLEXTOPOL)
        {
          for (size_t i = 0; i < gridsize; ++i)
          {
            array2[i] = SQRTM(ADDM(MULM(array1[2 * i], array1[2 * i]), MULM(array1[2 * i + 1], array1[2 * i + 1])));
            array3[i] = (is_EQ(array1[2 * i], missval1) || is_EQ(array1[2 * i + 1], missval1))
                            ? missval1
                            : std::atan2(array1[2 * i + 1], array1[2 * i]);
          }
        }

        cdo_write_field(streamID2, array2.data(), numMissVals);
        cdo_write_field(streamID3, array3.data(), numMissVals);
      }

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID3);
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);

    vlistDestroy(vlistID2);
    vlistDestroy(vlistID3);
  }
};
