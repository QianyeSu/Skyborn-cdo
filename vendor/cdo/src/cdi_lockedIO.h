#ifndef CDI_LOCKEDIO_H
#define CDI_LOCKEDIO_H

#include <cstddef>
#include <cstdbool>

int stream_open_read_locked(const char *filename);
void stream_close_locked(int p_fileID);
void stream_inq_field_locked(int p_fileID, int *p_varID, int *p_levelID);
void stream_def_field_locked(int p_fileID, int p_varID, int levelID);
void stream_read_field_float_locked(int p_fileID, float *p_data, size_t *p_numMissVals);
void stream_read_field_double_locked(int p_fileID, double *p_data, size_t *p_numMissVals);
void stream_def_vlist_locked(int p_fileID, int p_vlistID);
int stream_inq_vlist_locked(int p_fileID);
void stream_write_field_double_locked(int p_fileID, const double *const p_data, size_t p_numMissVals);
void stream_write_field_float_locked(int p_fileID, const float *const p_data, size_t p_numMissVals);
int stream_inq_time_step_locked(int p_fileID, int p_tsID);
int stream_def_time_step_locked(int p_fileID, int p_tsID);
int stream_copy_field_locked(int p_fileID, int p_targetFileID);
void vlist_copy_flag_locked(int p_vlistID2, int p_vlistID1);
void open_lock();
void open_unlock();
void cdo_vlist_copy_flag(int vlistID2, int vlistID1);

#endif
