#ifndef CTHREAD_DEBUG_H
#define CTHREAD_DEBUG_H

#include <mutex>
#include <condition_variable>

void Cthread_mutex_lock(const char *caller, std::mutex &p_mutex);
void Cthread_mutex_unlock(std::mutex &p_mutex);

void Cthread_cond_signal(const char *caller, std::condition_variable &p_cond_var);

#define cthread_mutex_lock(a) Cthread_mutex_lock(__func__, a)
#define cthread_mutex_unlock(a) Cthread_mutex_unlock(a)

#define cthread_cond_signal(a) Cthread_cond_signal(__func__, a)

#endif /* CTHREAD_DEBUG_H */
