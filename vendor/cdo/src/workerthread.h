/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef WORKERTHREAD_H
#define WORKERTHREAD_H

#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>

class WorkerThread
{
private:
  enum struct State
  {
    SETUP,
    IDLE,
    JOB,
    DIE
  };

  std::function<void()> function;

  State state{ State::SETUP };
  std::thread thread;
  std::mutex workMutex;
  std::mutex bossMutex;
  std::condition_variable workCond;
  std::condition_variable_any bossCond;
  static void task(WorkerThread *taskInfo);

public:
  WorkerThread();
  ~WorkerThread();

  void doAsync(const std::function<void()> &_function);
  void wait();
};

#endif /* WORKERTHREAD_H */
