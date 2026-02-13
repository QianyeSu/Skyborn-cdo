/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "process_int.h"

class Importfv3grid : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Importfv3grid",
    .operators = { { "import_fv3grid" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Importfv3grid> registration = RegisterEntry<Importfv3grid>();

  size_t gridsize1{};
  size_t gridsize2{};

  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  size_t nx{};
  size_t ny{};

  int gridIDo{};

  int vlistID2{ CDI_UNDEFID };
  int surfaceID{};

public:
  void
  init() override
  {
    operator_check_argc(0);

    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    auto nvars = vlistNvars(vlistID1);
    if (nvars != 5) cdo_abort("Found %d variables, expected 5 variables!", nvars);

    std::vector<std::string> vars(nvars);
    vars[0] = "grid_lon";
    vars[1] = "grid_lat";
    vars[2] = "grid_lont";
    vars[3] = "grid_latt";
    vars[4] = "area";

    for (int varID = 0; varID < nvars; ++varID)
    {
      auto varname = cdo::inq_var_name(vlistID1, varID);
      if (varname == vars[varID]) cdo_abort("Found variable %s, expected variable %s!", varname, vars[varID]);
    }

    auto numGrids = vlistNumGrids(vlistID1);
    if (numGrids != 2) cdo_abort("Found %d grids, expected 2 grids!", nvars);

    auto gridIDi1 = vlistGrid(vlistID1, 0);
    auto gridIDi2 = vlistGrid(vlistID1, 1);

    nx = gridInqXsize(gridIDi1);
    ny = gridInqYsize(gridIDi1);
    gridsize1 = gridInqSize(gridIDi1);
    gridsize2 = gridInqSize(gridIDi2);

    cdo_stream_inq_timestep(streamID1, 0);

    gridIDo = gridCreate(GRID_UNSTRUCTURED, gridsize2);
    gridDefNvertex(gridIDo, 4);
  }

  void
  run() override
  {
    Varray<double> buffer(gridsize1);
    size_t numMissVals;

    {
      Varray<double> grid_corner(4 * gridsize2);

      (void) cdo_inq_field(streamID1);  // grid_lon
      cdo_read_field(streamID1, buffer.data(), &numMissVals);

      for (size_t j = 1; j < ny; ++j)
        for (size_t i = 1; i < nx; ++i)
        {
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 0] = buffer[j * nx + i];
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 1] = buffer[j * nx + (i - 1)];
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 2] = buffer[(j - 1) * nx + (i - 1)];
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 3] = buffer[(j - 1) * nx + i];
        }

      gridDefXbounds(gridIDo, grid_corner.data());

      (void) cdo_inq_field(streamID1);  // grid_lat
      cdo_read_field(streamID1, buffer.data(), &numMissVals);

      for (size_t j = 1; j < ny; ++j)
        for (size_t i = 1; i < nx; ++i)
        {
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 0] = buffer[j * nx + i];
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 1] = buffer[j * nx + (i - 1)];
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 2] = buffer[(j - 1) * nx + (i - 1)];
          grid_corner[((j - 1) * (nx - 1) + (i - 1)) * 4 + 3] = buffer[(j - 1) * nx + i];
        }

      gridDefYbounds(gridIDo, grid_corner.data());
    }

    (void) cdo_inq_field(streamID1);  // grid_lont
    cdo_read_field(streamID1, buffer.data(), &numMissVals);
    gridDefXvals(gridIDo, buffer.data());

    (void) cdo_inq_field(streamID1);  // grid_latt
    cdo_read_field(streamID1, buffer.data(), &numMissVals);
    gridDefYvals(gridIDo, buffer.data());

    (void) cdo_inq_field(streamID1);  // area
    cdo_read_field(streamID1, buffer.data(), &numMissVals);

    double sfclevel = 0;
    surfaceID = zaxisCreate(ZAXIS_SURFACE, 1);
    zaxisDefLevels(surfaceID, &sfclevel);

    vlistID2 = vlistCreate();
    int varID = vlistDefVar(vlistID2, gridIDo, surfaceID, TIME_CONSTANT);
    cdiDefKeyString(vlistID2, varID, CDI_KEY_NAME, "area");
    cdiDefKeyString(vlistID2, varID, CDI_KEY_STDNAME, "area");
    cdiDefKeyString(vlistID2, varID, CDI_KEY_LONGNAME, "cell area");
    cdiDefKeyString(vlistID2, varID, CDI_KEY_UNITS, "m2");
    vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT32);

    auto taxisID = cdo_taxis_create(TAXIS_ABSOLUTE);
    vlistDefTaxis(vlistID2, taxisID);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
    cdo_def_timestep(streamID2, 0);
    cdo_def_field(streamID2, 0, 0);
    cdo_write_field(streamID2, buffer.data(), 0);
  }

  void
  close() override
  {
    cdo_stream_close(streamID2);
    cdo_stream_close(streamID1);
    vlistDestroy(vlistID2);
  }
};
