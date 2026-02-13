/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "c_wrapper.h"
#include "process_int.h"
#include "varray.h"

class Testdata : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Testdata",
    .operators = { { "testdata" } },
    .aliases = {},
    .mode = INTERNAL,    // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Testdata> registration = RegisterEntry<Testdata>();

  CdoStreamID streamID1;
  CdoStreamID streamID2;
  int taxisID1{ CDI_UNDEFID };
  int vlistID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  void
  init() override
  {
    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);
    taxisID1 = vlistInqTaxis(vlistID1);

    streamID2 = cdo_open_write(1);

    int vlistID2 = vlistDuplicate(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    VarList varList1(vlistID1);

    auto gridsizeMax = varList1.gridsizeMax();
    Varray<double> array(gridsizeMax);
    Varray<float> fval(gridsizeMax);
    Varray<int> ival(gridsizeMax);
    Varray<unsigned char> cval(gridsizeMax * 4);
    Varray<unsigned char> cval2(gridsizeMax * 4);

    auto fobj = c_fopen("testdata", "w");

    int tsID2 = 0;
    int tsID1 = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID1);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID2);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        cdo_def_field(streamID2, varID, levelID);

        size_t numMissVals;
        cdo_read_field(streamID1, array.data(), &numMissVals);

        gridsizeMax = varList1.vars[varID].gridsize;
        for (size_t i = 0; i < gridsizeMax; ++i)
        {
          fval[i] = (float) array[i];

          std::memcpy(&ival[i], &fval[i], 4);
          std::memcpy(&cval[i * 4], &fval[i], 4);

          cval2[i + gridsizeMax * 0] = cval[i * 4 + 0];
          cval2[i + gridsizeMax * 1] = cval[i * 4 + 1];
          cval2[i + gridsizeMax * 2] = cval[i * 4 + 2];
          cval2[i + gridsizeMax * 3] = cval[i * 4 + 3];

          if (tsID1 == 0 && fieldID == 0)
            printf("%4zu %3u %3u %3u %3u %d %g\n", i, (unsigned int) cval[4 * i + 0], (unsigned int) cval[4 * i + 1],
                   (unsigned int) cval[4 * i + 2], (unsigned int) cval[4 * i + 3], ival[i], fval[i]);
        }

        cdo_write_field(streamID2, array.data(), numMissVals);

        std::fwrite(cval.data(), 4, gridsizeMax, fobj.get());
      }

      tsID1++;
      tsID2++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);
  }
};
