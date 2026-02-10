/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <mutex>
#include <condition_variable>

#include "cdo_output.h"

void
Cthread_mutex_lock(const char *caller, std::mutex &p_mutex)
{
  try
    {
      p_mutex.lock();
    }
  catch (const std::system_error &e)
    {
      std::cout << "locking failed in " << caller << ". ErrorCode:" << e.code() << " " << e.what() << std::endl;
    }
}

void
Cthread_mutex_unlock(std::mutex &p_mutex)
{
  p_mutex.unlock();
}

void
Cthread_cond_signal(const char *caller, std::condition_variable &p_cond_var)
{
  Debug(PTHREAD, "+%s (cond = %p)", caller, (void *) &p_cond_var);
  p_cond_var.notify_all();
  Debug(PTHREAD, "-%s (cond = %p)", caller, (void *) &p_cond_var);
}
