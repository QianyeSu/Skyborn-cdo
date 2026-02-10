#ifndef PRINTINFO_H
#define PRINTINFO_H

#include <string>
#include <cdoStream.h>

std::string date_to_string(CdiDate date);
std::string time_to_string(CdiTime time);
std::string datetime_to_string(CdiDateTime dt);

const char *comptype_to_name(int comptype);

void print_filetype(CdoStreamID p_streamID, int vlistID);
void print_grid_info(int vlistID);
void print_zaxis_info(int vlistID);
void print_subtype_info(int vlistID);
void print_timesteps(CdoStreamID streamID, int taxisID, int verbose);

#endif
