#ifndef PERCENTILES_HIST_H_
#define PERCENTILES_HIST_H_

#include <cstdlib>
#include <cassert>

#include "field.h"

struct HistogramEntry
{
  void *ptr = nullptr;
  float min = 0.0f;
  float max = 0.0f;
  float step = 0.0f;
  int nsamp = 0;
  int capacity = 0;
  int numBins = 0;
  bool isUint32 = false;
};

// clang-format off
class  // HistogramSet
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
HistogramSet
// clang-format on
{
private:
  int numVars = 0;
  int numSteps = 0;
  std::vector<int> var_numLevels;
  std::vector<size_t> var_numHists;
  std::vector<std::vector<std::vector<HistogramEntry>>> histograms;

  void
  init()
  {
    var_numLevels.resize(numVars, 0);
    var_numHists.resize(numVars, 0);
    histograms.resize(numVars);
  }

public:
  HistogramSet() {}

  explicit HistogramSet(int _numVars) : numVars(_numVars)
  {
    assert(numVars > 0);
    init();
  }

  HistogramSet(int _numVars, int _numSteps) : numVars(_numVars), numSteps(_numSteps)
  {
    assert(numVars > 0);
    init();
  }

  void
  create(int _numVars, int _numSteps = 0)
  {
    numVars = _numVars;
    numSteps = _numSteps;
    assert(numVars > 0);
    init();
  }

  ~HistogramSet()
  {
    for (auto varID = numVars; varID-- > 0;)
      {
        auto numHists = this->var_numHists[varID];
        for (auto levelID = this->var_numLevels[varID]; levelID-- > 0;)
          {
            for (auto histID = numHists; histID-- > 0;) std::free(this->histograms[varID][levelID][histID].ptr);
          }
      }
  }

  void createVarLevels(int varID, int numLevels, size_t numHists);
  void defVarLevelBounds(int varID, int levelID, Field const &field1, Field const &field2);
  int addVarLevelValues(int varID, int levelID, Field const &field);
  int subVarLevelValues(int varID, int levelID, Field const &field);
  void getVarLevelPercentiles(Field &field, int varID, int levelID, double p);
  // void reset(int varID, int levelID); // unused
};

#endif /* PERCENTILES_HIST_H_ */
