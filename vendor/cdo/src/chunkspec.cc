#include "cdi.h"

#include "chunkspec.h"
#include "cdo_output.h"
#include "util_string.h"
#include "param_conversion.h"

namespace cdo
{

ChunkSpec
get_chunkspec(int vlistID, int varID)
{
  ChunkSpec chunkSpec;
  std::vector<std::tuple<int, int &>> chunkList = { { CDI_KEY_CHUNKSIZE_DIMX, chunkSpec.x },
                                                    { CDI_KEY_CHUNKSIZE_DIMY, chunkSpec.y },
                                                    { CDI_KEY_CHUNKSIZE_DIMZ, chunkSpec.z },
                                                    { CDI_KEY_CHUNKSIZE_DIMT, chunkSpec.t } };
  int chunkSize;
  for (auto [key, dimSpec] : chunkList)
    {
      if (cdiInqKeyInt(vlistID, varID, key, &chunkSize) == 0)
        {
          if (chunkSize) dimSpec = chunkSize;
        }
    }

  return chunkSpec;
}

ChunkSpec
parse_chunkspec_parameter(std::string const &argument)
{
  ChunkSpec chunkSpec;
  auto pargv = split_string(argument, ",");
  for (auto const &parg : pargv)
    {
      auto const keyValue = split_string(parg, "=");
      if (keyValue.size() != 2) cdo_abort("Invalid chunkspec parameter %s", argument);
      auto const &key = keyValue[0];
      auto const &value = keyValue[1];
      if (key == "t")
        {
          int chunkSize = parameter_to_bytes(value);
          if (chunkSize > 0) chunkSpec.t = chunkSize;
        }
      else if (key == "z")
        {
          int chunkSize = parameter_to_bytes(value);
          if (chunkSize > 0 || chunkSize == -1) chunkSpec.z = chunkSize;
        }
      else if (key == "y")
        {
          int chunkSize = parameter_to_bytes(value);
          if (chunkSize > 0 || chunkSize == -1) chunkSpec.y = chunkSize;
        }
      else if (key == "x")
        {
          int chunkSize = parameter_to_bytes(value);
          if (chunkSize > 0 || chunkSize == -1) chunkSpec.x = chunkSize;
        }
      else { cdo_abort("Invalid chunkspec parameter %s (dim=%s not available)", argument, key); }
    }

  return chunkSpec;
}

}  // namespace cdo
