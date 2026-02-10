/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:

      Filedes    codetab         Parameter code table
      Filedes    griddes         Grid description
      Filedes    vct             Vertical coordinate table
*/

#include <cdi.h>

#include "cdo_options.h"
#include "process_int.h"
#include <mpim_grid.h>

void cdoPrintZaxis(int zaxisID);

static void
printSource(FILE *fp, int vlistID, int varID)
{
  // institute info
  auto instptr = institutInqLongnamePtr(vlistInqVarInstitut(vlistID, varID));
  if (instptr) fprintf(fp, "  institution=\"%s\"\n", instptr);

  // source info
  auto modelptr = modelInqNamePtr(vlistInqVarModel(vlistID, varID));
  if (modelptr) fprintf(fp, "  source=\"%s\"\n", modelptr);
}

static void
printVCT(int vlistID, bool lvct)
{
  auto numZaxes = vlistNumZaxis(vlistID);
  for (int index = 0; index < numZaxes; ++index)
  {
    auto zaxisID = vlistZaxis(vlistID, index);
    auto type = zaxisInqType(zaxisID);
    if (type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF)
    {
      auto vctsize = zaxisInqVctSize(zaxisID);
      auto vct = zaxisInqVctPtr(zaxisID);

      if (vctsize % 2 == 0)
      {
        if (lvct)
        {
          fprintf(stdout, "#   k         vct_a(k) [Pa]             vct_b(k) []\n");
          for (int i = 0; i < vctsize / 2; ++i) fprintf(stdout, "%5d %25.17f %25.17f\n", i, vct[i], vct[vctsize / 2 + i]);
        }
        else
        {
          fprintf(stdout, "vctsize   = %d\n", vctsize);
          int nbyte0 = fprintf(stdout, "vct       = ");
          int nbyte = nbyte0;
          for (int i = 0; i < vctsize; ++i)
          {
            if (nbyte > 70 || i == vctsize / 2)
            {
              fprintf(stdout, "\n%*s", nbyte0, "");
              nbyte = nbyte0;
            }
            nbyte += fprintf(stdout, "%.9g ", vct[i]);
          }
          fprintf(stdout, "\n");
        }
      }
      else
        for (int i = 0; i < vctsize; ++i) fprintf(stdout, "%5d %25.17f\n", i, vct[i]);

      break;
    }
  }
}

static void
printCodeTable(VarList const &varList)
{
  int nvars = varList.numVars();
  for (int varID = 0; varID < nvars; ++varID)
  {
    auto const &var = varList.vars[varID];
    fprintf(stdout, "%4d  %-12s", var.code, var.name.c_str());
    if (var.longname.size())
    {
      fprintf(stdout, "  %s", var.longname.c_str());
      if (var.units.size()) fprintf(stdout, " [%s]", var.units.c_str());
    }
    fprintf(stdout, "\n");
  }
}

static void
partab(FILE *fp, int vlistID, VarList const &varList, int option)
{
  int varID, datatype = -1;
  char paramstr[32];

  int numVars = varList.numVars();
  auto linebreak = (option != 4);

  if (option == 2)
  {
    int natts;
    cdiInqNatts(vlistID, CDI_GLOBAL, &natts);
    if (natts > 0)
    {
      fprintf(fp, "&parameter\n");
      fprintf(fp, "  name=_GLOBAL_\n");
      printSource(fp, vlistID, 0);
      cdo_print_attributes(fp, vlistID, CDI_GLOBAL, 2);
      fprintf(fp, "/\n");
    }
  }

  if (numVars > 1)
  {
    datatype = varList.vars[0].dataType;
    for (varID = 1; varID < numVars; ++varID)
    {
      if (datatype != varList.vars[varID].dataType)
      {
        datatype = -1;
        break;
      }
    }

    if (datatype != -1)
    {
      fprintf(fp, "&parameter");
      if (linebreak) fprintf(fp, "\n");
      fprintf(fp, "  name=_default_");
      if (linebreak) fprintf(fp, "\n");
      auto datatypestr = cdo::datatype_to_cstr(datatype);
      if (*datatypestr)
      {
        fprintf(fp, "  datatype=%s", datatypestr);
        if (linebreak) fprintf(fp, "\n");
      }
      fprintf(fp, "/\n");
    }
  }

  for (varID = 0; varID < numVars; ++varID)
  {
    auto const &var = varList.vars[varID];

    fprintf(fp, "&parameter");
    if (linebreak) fprintf(fp, "\n");

    fprintf(fp, "  name=%s", var.name.c_str());
    if (linebreak) fprintf(fp, "\n");

    if (var.param >= 0)
    {
      cdiParamToString(var.param, paramstr, sizeof(paramstr));
      fprintf(fp, "  param=%s", paramstr);
      if (linebreak) fprintf(fp, "\n");
    }
    if (var.stdname.size())
    {
      fprintf(fp, "  standard_name=%s", var.stdname.c_str());
      if (linebreak) fprintf(fp, "\n");
    }
    if (var.longname.size())
    {
      fprintf(fp, "  long_name=\"%s\"", var.longname.c_str());
      if (linebreak) fprintf(fp, "\n");
    }
    if (var.units.size())
    {
      fprintf(fp, "  units=\"%s\"", var.units.c_str());
      if (linebreak) fprintf(fp, "\n");
    }

    if (datatype == -1)
    {
      auto datatypestr = cdo::datatype_to_cstr(var.dataType);
      if (*datatypestr)
      {
        fprintf(fp, "  datatype=%s", datatypestr);
        if (linebreak) fprintf(fp, "\n");
      }
    }

    int uvRelativeToGrid = 0;
    if (cdiInqKeyInt(vlistID, varID, CDI_KEY_UVRELATIVETOGRID, &uvRelativeToGrid) == CDI_NOERR)
    {
      fprintf(fp, "  uvRelativeToGrid=%d", uvRelativeToGrid);
      if (linebreak) fprintf(fp, "\n");
    }

    int chunkType = -1;
    cdiInqKeyInt(vlistID, varID, CDI_KEY_CHUNKTYPE, &chunkType);
    const char *chunkName = (chunkType == CDI_CHUNK_AUTO)
                                ? "auto"
                                : ((chunkType == CDI_CHUNK_GRID) ? "grid" : ((chunkType == CDI_CHUNK_LINES) ? "lines" : nullptr));

    if (chunkName)
    {
      fprintf(fp, "  chunkType=%s", chunkName);
      if (linebreak) fprintf(fp, "\n");
    }

    if (option == 2)
    {
      fprintf(fp, "  missing_value=%g\n", var.missval);
      cdo_print_attributes(fp, vlistID, varID, 2);
    }

    if (!linebreak) fprintf(fp, "  ");
    fprintf(fp, "/\n");
  }
}

static void
filedes(CdoStreamID streamID)
{
  printf("\n");
  auto filetype = cdo_inq_filetype(streamID);

  auto filetypestr = cdo::filetype_to_cstr(filetype);
  if (filetypestr == nullptr || *filetypestr == 0)
    printf("  unsupported filetype %d\n", filetype);
  else
    printf("  %s data\n", filetypestr);

  if (filetype == CDI_FILETYPE_SRV || filetype == CDI_FILETYPE_EXT || filetype == CDI_FILETYPE_IEG)
  {
    auto byteorder = cdo_inq_byteorder(streamID);
    switch (byteorder)
    {
      case CDI_BIGENDIAN: printf("  byteorder is BIGENDIAN\n"); break;
      case CDI_LITTLEENDIAN: printf("  byteorder is LITTLEENDIAN\n"); break;
      default: printf("  byteorder %d undefined\n", byteorder); break;
    }
  }

  printf("\n");
}

class Filedes : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Filedes",
    .operators = { { "filedes", FiledesHelp },
                   { "griddes", FiledesHelp },
                   { "griddes2", FiledesHelp },
                   { "zaxisdes", FiledesHelp },
                   { "vct", FiledesHelp },
                   { "vct2", FiledesHelp },
                   { "codetab", FiledesHelp },
                   { "vlist", FiledesHelp },
                   { "partab", FiledesHelp },
                   { "partab2", FiledesHelp },
                   { "spartab", FiledesHelp } },
    .aliases = { { "vardes", "codetab" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_BOTH,  // Allowed number type
    .constraints = { 1, 0, NoRestriction },
  };
  inline static RegisterEntry<Filedes> registration = RegisterEntry<Filedes>(module);

public:
  void
  init() override
  {
    auto GRIDDES = module.get_id("griddes");
    auto GRIDDES2 = module.get_id("griddes2");
    auto ZAXISDES = module.get_id("zaxisdes");
    auto VCT = module.get_id("vct");
    auto VCT2 = module.get_id("vct2");
    auto CODETAB = module.get_id("codetab");
    auto FILEDES = module.get_id("filedes");
    auto VLIST = module.get_id("vlist");
    auto SPARTAB = module.get_id("spartab");
    auto PARTAB = module.get_id("partab");
    auto PARTAB2 = module.get_id("partab2");

    auto operatorID = cdo_operator_id();

    auto loadGrid = (operatorID == GRIDDES || operatorID == GRIDDES2);
    if (Options::lazyGridLoad && this_is_the_only_process()) { cdiDefGlobal("NETCDF_LAZY_GRID_LOAD", true); }
    if (not loadGrid && this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CORNERS", false); }
    if (not loadGrid && this_is_the_only_process()) { cdiDefGlobal("READ_CELL_CENTER", false); }

    operator_check_argc(0);

    auto streamID = cdo_open_read(0);
    auto vlistID = cdo_stream_inq_vlist(streamID);

    VarList varList(vlistID);

    if (operatorID == GRIDDES || operatorID == GRIDDES2)
    {
      auto opt = (operatorID == GRIDDES) ? 1 : 0;
      auto numGrids = vlistNumGrids(vlistID);
      for (int index = 0; index < numGrids; ++index)
      {
        printf("#\n# gridID %d\n#\n", index + 1);
        cdo_print_griddes(vlistGrid(vlistID, index), opt);
        auto nsubtypes = vlistNsubtypes(vlistID);
        for (int i = 0; i < nsubtypes; ++i) subtypePrint(vlistSubtype(vlistID, i));
      }
    }
    else if (operatorID == ZAXISDES)
    {
      auto numZaxes = vlistNumZaxis(vlistID);
      for (int index = 0; index < numZaxes; ++index)
      {
        printf("#\n# zaxisID %d\n#\n", index + 1);
        cdoPrintZaxis(vlistZaxis(vlistID, index));
      }
    }
    else if (operatorID == VCT || operatorID == VCT2) { printVCT(vlistID, operatorID == VCT); }
    else if (operatorID == VLIST) { vlistPrint(vlistID); }
    else if (operatorID == CODETAB) { printCodeTable(varList); }
    else if (operatorID == PARTAB || operatorID == SPARTAB || operatorID == PARTAB2)
    {
      auto option = (operatorID == SPARTAB) ? 4 : ((operatorID == PARTAB2) ? 2 : 1);
      partab(stdout, vlistID, varList, option);
    }
    else if (operatorID == FILEDES) { filedes(streamID); }

    cdo_stream_close(streamID);

    if (!cdo::stdoutIsTerminal) Options::silentMode = true;
  }

  void
  run() override
  {
  }

  void
  close() override
  {
  }
};
