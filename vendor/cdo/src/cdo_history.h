#ifndef CDO_HISTORY_H
#define CDO_HISTORY_H

#include <string>

std::string cdo_inq_history(int vlistID);
void cdo_append_history(int vlistID, const char *histstring);
void cdo_def_creation_date(int vlistID);
void cdo_def_tracking_id(int vlistID, const char *uuid_attribute);

#endif
