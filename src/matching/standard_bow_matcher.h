#ifndef GFP_STANDARD_BOW_MATCHER_H
#define GFP_STANDARD_BOW_MATCHER_H

#include "base_matcher.h"
#include <boost/make_shared.hpp>

class StandardBowMatcher : public BaseMatcher
{

public:

  /**
   * @brief Represent the different ways to compute the term frequency in a document
   */
  enum BowType
  {
    STANDARD,
    SUBLINEAR,
    LENGTH_SMOOTHING
  };


  StandardBowMatcher(const size_t  vocabulary_size, BowType type = STANDARD, double alpha_smoothing = 0.4 )
    : BaseMatcher(vocabulary_size)
    , bow_type_(bow_type_)
    , alpha_smoothing_(alpha_smoothing)
  {}

  virtual void getCandidates(const LaserScanData& scan, ScoreSet& score_set);

  virtual void updateLaserScanData();

  virtual LaserScanData::Ptr getLaserScanData(const std::vector<int>& scan_words) const
  {
    return boost::make_shared<LaserScanData>(scan_words, LaserScanData::STANDARD_BOW);
  }

private:

  double alpha_smoothing_;
  BowType bow_type_;

};

#endif // GFP_STANDARD_BOW_MATCHER_H
