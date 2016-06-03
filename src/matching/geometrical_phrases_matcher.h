#ifndef GFP_GEOMETRICAL_PHRASES_MATCHER_H
#define GFP_GEOMETRICAL_PHRASES_MATCHER_H

#include "base_matcher.h"
#include <boost/make_shared.hpp>

class GeometricalPhrasesMatcher : public BaseMatcher
{

public:

  GeometricalPhrasesMatcher(const size_t  vocabulary_size, unsigned int kernel_size = 2)
    : BaseMatcher(vocabulary_size)
    , gp_kernel_size_(kernel_size)
  {
    cache_binomial_coeff();
  }

  virtual void updateLaserScanData();

  virtual void getCandidates(const LaserScanData& scan, ScoreSet& score_set);

  virtual LaserScanData::Ptr getLaserScanData(const std::vector<int>& scan_words) const
  {
    return boost::make_shared<LaserScanData>(scan_words, LaserScanData::GEOMETRICAL_PHRASES);
  }

private:

  void cache_binomial_coeff();

  unsigned int gp_kernel_size_;

  std::vector<double> cached_binomial_coeff;

};



#endif // GFP_GEOMETRICAL_PHRASES_MATCHER_H
