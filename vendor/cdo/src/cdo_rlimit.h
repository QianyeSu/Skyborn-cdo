#ifndef CDO_RLIMIT_H
#define CDO_RLIMIT_H

namespace cdo
{
void print_rlimits(void);
void set_numfiles(long numfiles);
void set_stacksize(long stacksize);
void set_coresize(long coresize);
long get_rss_cur();
}  // namespace cdo

#endif
