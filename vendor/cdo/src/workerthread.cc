/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#if defined(_OPENMP)
#include <omp.h>
#include "cdo_options.h"
#endif

#include "workerthread.h"

void
WorkerThread::task(WorkerThread *taskInfo)
{
#if defined(_OPENMP)
  omp_set_num_threads(Threading::ompNumMaxThreads);  // Has to be called for every thread!
#endif

  // cond.wait mutex must be locked before we can wait
  std::unique_lock<std::mutex> workLock(taskInfo->workMutex);
  // ensure boss is waiting
  taskInfo->bossMutex.lock();
  // signal to boss that setup is complete
  taskInfo->state = State::IDLE;
  // wake-up signal
  taskInfo->bossCond.notify_one();
  taskInfo->bossMutex.unlock();

  while (1)
  {
    taskInfo->workCond.wait(workLock);

    if (State::DIE == taskInfo->state) break;      // kill thread
    if (State::IDLE == taskInfo->state) continue;  // accidental wake-up

    // do blocking task
    // printf("<worker> JOB start\n");
    taskInfo->function();
    // printf("<worker> JOB end\n");
    // ensure boss is waiting
    taskInfo->bossMutex.lock();
    // indicate that job is done
    taskInfo->state = State::IDLE;
    // wake-up signal
    taskInfo->bossCond.notify_one();
    taskInfo->bossMutex.unlock();
  }
}

void
WorkerThread::doAsync(const std::function<void()> &_function)
{
  // ensure worker is waiting
  std::lock_guard<std::mutex> _(workMutex);
  // set job information & state
  this->function = _function;
  this->state = State::JOB;
  // wake-up signal
  workCond.notify_one();
}

void
WorkerThread::wait()
{
  while (1)
  {
    if (State::IDLE == this->state) break;
    bossCond.wait(bossMutex);
  }
}

WorkerThread::WorkerThread()
{
  bossMutex.lock();
  this->thread = std::thread(this->task, this);
  this->wait();
}

WorkerThread::~WorkerThread()
{
  // ensure the worker is waiting
  workMutex.lock();
  // printf("Task::delete: send DIE to <worker>\n");
  this->state = State::DIE;
  // wake-up signal
  workCond.notify_one();
  workMutex.unlock();
  // wait for thread to exit
  this->thread.join();
  bossMutex.unlock();
}

#ifdef TEST_WORKERTHREAD
// g++ -g -Wall -O2 -DTEST_WORKERTHREAD workerthread.cc

void
func(int &intArg)
{
  intArg = -1;
  printf("run myfunc\n");
}

void
worker1(void)
{
  WorkerThread workerThread;

  int ivalue = 0;
  std::function<void()> my_task = std::bind(func, std::ref(ivalue));

  workerThread.doAsync(my_task);
  workerThread.wait();

  printf("worker1: %d\n", ivalue);
}

void
worker2(void)
{
  bool useWorkerThread = true;
  auto workerThread = useWorkerThread ? std::make_unique<WorkerThread>() : nullptr;

  int ivalue = 0;
  std::function<void()> my_task = std::bind(func, std::ref(ivalue));

  useWorkerThread ? workerThread->doAsync(my_task) : my_task();

  if (useWorkerThread) workerThread->wait();

  printf("worker2: %d\n", ivalue);
}

int
main(int argc, char **argv)
{
  worker1();
  worker2();

  return 0;
}
#endif
