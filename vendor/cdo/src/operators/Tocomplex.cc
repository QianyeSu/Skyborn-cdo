/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"
#include "cdo_default_values.h"  // Namespace CdoDefault

class Tocomplex : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Tocomplex",
    .operators = { { "retocomplex" }, { "imtocomplex" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Tocomplex> registration = RegisterEntry<Tocomplex>();

  int RETOCOMPLEX, IMTOCOMPLEX;
  CdoStreamID streamID1;
  CdoStreamID streamID2;
  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int operatorID;

  VarList varList1;

public:
  void
  init() override
  {
    RETOCOMPLEX = module.get_id("retocomplex");
    IMTOCOMPLEX = module.get_id("imtocomplex");

    operatorID = cdo_operator_id();

    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    auto nvars = vlistNvars(vlistID2);
    for (int varID = 0; varID < nvars; ++varID)
    {
      auto datatype = vlistInqVarDatatype(vlistID2, varID);
      datatype = (datatype == CDI_DATATYPE_FLT64) ? CDI_DATATYPE_CPX64 : CDI_DATATYPE_CPX32;
      vlistDefVarDatatype(vlistID2, varID, datatype);
    }

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    // if (CdoDefault::FileType != CDI_FILETYPE_EXT) cdo_abort("Complex numbers need EXTRA format; used CDO option -f ext!");
    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    Varray<double> array1(varList1.gridsizeMax());
    Varray<double> array2(2 * varList1.gridsizeMax());

    int tsID = 0;
    int tsID2 = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID2++);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_def_field(streamID2, varID, levelID);

        size_t numMissVals;
        cdo_read_field(streamID1, array1.data(), &numMissVals);

        auto gridsize = varList1.vars[varID].gridsize;
        if (operatorID == RETOCOMPLEX)
        {
          for (size_t i = 0; i < gridsize; ++i)
          {
            array2[2 * i] = array1[i];
            array2[2 * i + 1] = 0;
          }
        }
        else if (operatorID == IMTOCOMPLEX)
        {
          for (size_t i = 0; i < gridsize; ++i)
          {
            array2[2 * i] = 0;
            array2[2 * i + 1] = array1[i];
          }
        }

        cdo_write_field(streamID2, array2.data(), numMissVals);
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
