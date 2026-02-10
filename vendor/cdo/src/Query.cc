/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "cdo_options.h"
#include "cdo_query.h"
#include "process_int.h"
#include "param_conversion.h"
#include "pmlist.h"
#include "cdo_default_values.h"

class Query : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Query",
    .operators = { { "query", 0, 0, "queryentries" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, FilesOnly },
  };
  inline static RegisterEntry<Query> registration = RegisterEntry<Query>(module);

  int streamID1{};  // QueryStream
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  VarList varList1{};

public:
  void
  init() override
  {
    // auto dataIsUnchanged = data_is_unchanged();

    auto operatorID = cdo_operator_id();

    operator_input_arg(cdo_operator_enter(operatorID));

    auto natts = cdo_operator_argc();
    if (natts == 0) cdo_abort("Parameter missing!");

    if (cdo_assert_files_only() == false) cdo_abort("This operator can't be combined with other operators!");

    PMList pmlist;
    KVList kvlist;
    kvlist.name = cdo_module_name();
    if (kvlist.parse_arguments(cdo_get_oper_argv()) != 0) cdo_abort("Parse error!");
    if (Options::cdoVerbose) kvlist.print();

    auto pkvlist = &kvlist;
    if (natts == 1)
    {
      KeyValues &kv = kvlist.front();
      if (kv.key == "FILE")
      {
        if (Options::cdoVerbose) cdo_print("Reading query from: %s", kv.values[0]);
        auto filename = parameter_to_word(kv.values[0].c_str());
        auto fp = std::fopen(filename, "r");
        if (fp == nullptr) cdo_abort("Open failed on: %s\n", filename);
        pmlist.read_namelist(fp, filename);
        pkvlist = &pmlist.front();
        std::fclose(fp);
        if (Options::cdoVerbose) pkvlist->print();
      }
    }

    auto query = cdiQueryCreate();
    set_query_parameter(*pkvlist, query);
    if (Options::cdoVerbose) cdiQueryPrint(query);

    streamID1 = streamOpenReadQuery(cdo_get_stream_name(0), query);
    if (streamID1 < 0) cdi_open_error(streamID1, "Open failed on >%s<", cdo_get_stream_name(0));

    cdiQueryPrintEntriesNotFound(query);
    cdiQueryDelete(query);

    auto filetype = streamInqFiletype(streamID1);
    if (CdoDefault::FileType == CDI_UNDEFID) CdoDefault::FileType = filetype;

    auto vlistID1 = streamInqVlist(streamID1);
    auto vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);

    varList1 = VarList(vlistID1);
  }

  void
  run() override
  {
    Field field;

    int tsID = 0;
    while (true)
    {
      auto numFields = streamInqTimestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        int varID, levelID;
        streamInqField(streamID1, &varID, &levelID);
        cdo_def_field(streamID2, varID, levelID);
        /*
        if (dataIsUnchanged)
          {
            streamCopyField(streamID2, streamID1);
          }
        else
        */
        {
          field.init(varList1.vars[varID]);
          if (field.memType == MemType::Float)
            streamReadFieldF(streamID1, field.vec_f.data(), &field.numMissVals);
          else
            streamReadField(streamID1, field.vec_d.data(), &field.numMissVals);
          cdo_write_field(streamID2, field);
        }
      }

      tsID++;
    }
  }

  void
  close() override
  {
    streamClose(streamID1);
    cdo_stream_close(streamID2);
  }
};
