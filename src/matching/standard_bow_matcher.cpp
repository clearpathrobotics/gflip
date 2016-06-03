#include "standard_bow_matcher.h"


/**
 * @brief Simple Functor class to help normalize and transform
 * a map of id score, into a ScoreSet
 */
class ScoreNormalizer
{
public:
  ScoreNormalizer(const double norm) : norm_ (norm) {}

  std::pair<double, size_t> operator()(const std::pair<size_t, double>& score)
  {
    return std::make_pair(1 - score.second/norm_, score.first);
  }

private:
  double norm_;
};


void StandardBowMatcher::updateLaserScanData()
{
  BOOST_FOREACH (LaserScanData &bow, scans_ | boost::adaptors::map_values)
  {
    double norm = 0;

    const double inverse_max_histogram = 1.0 / bow.getMaxHistogram();
    BOOST_FOREACH (const int word, bow.uniqueWords())
    {
      const double weight_unnormalized = bow.getWeight(word).first;
      const double weight = bow.getWeight(word).second;
      double & tf_idf = bow.tf_idf(word);

      switch (bow_type_)
      {
        case STANDARD:
          tf_idf = weight * idf(word);
          break;
        case SUBLINEAR:
          tf_idf = (1 + log(weight_unnormalized) ) * idf(word);
          break;
        case LENGTH_SMOOTHING:
          tf_idf =
              (alpha_smoothing_ + ((1.0 - alpha_smoothing_) * weight_unnormalized) * inverse_max_histogram) * idf(word);
      }

      norm += tf_idf * tf_idf;
    }

    norm /= sqrt(norm);

    BOOST_FOREACH (const int word, bow.uniqueWords())
    {
      bow.tf_idf(word) *= norm;
    }
  }
}

void StandardBowMatcher::getCandidates(const LaserScanData &scan, ScoreSet &score_set)
{
  double norm = 0;
  std::map<size_t, double> candidates;

  BOOST_FOREACH (const int word, scan.uniqueWords())
  {
    norm += scan.getWeight(word).first * scan.getWeight(word).first;

    BOOST_FOREACH (size_t scan_id, occurence_list_[word])
    {
      candidates[scan_id] += scans_[scan_id].tf_idf(word);
    }
  }

  norm = std::sqrt(norm);

  ScoreNormalizer score_normalizer(norm);

  score_set.resize(candidates.size());

  std::transform(candidates.begin(), candidates.end(), score_set.begin(), score_normalizer);

  std::sort(score_set.begin(), score_set.end());

}





