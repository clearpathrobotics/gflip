#include "geometrical_phrases_matcher.h"

#include <boost/math/special_functions/binomial.hpp>

#include <boost/range/adaptor/map.hpp>


void GeometricalPhrasesMatcher::updateLaserScanData()
{
  BOOST_FOREACH (LaserScanData scan, scans_ | boost::adaptors::map_values)
  {
    double& norm_wgv = scan.getGpNorm();
    norm_wgv = 0;

    BOOST_FOREACH (const std::vector<int>& words, scan.getOffsetHistograms())
    {
      const size_t size = words.size();
      if (size < gp_kernel_size_)
      {
        continue;
      }

      double sum_idf = 0;

      BOOST_FOREACH (const int & word, words)
      {
        sum_idf += idf(word);
      }

      norm_wgv += cached_binomial_coeff[size - 1] * sum_idf;
    }

    if (norm_wgv > 0)
    {
      norm_wgv = sqrt(norm_wgv);
    }
  }
}

void GeometricalPhrasesMatcher::getCandidates(const LaserScanData &scan, ScoreSet &score_set)
{
  typedef std::map<int, std::pair<int, double> > OffsetHistogram;
  typedef std::map<size_t, OffsetHistogram > OffsetHistograms;

  OffsetHistograms offset_hitograms;

  for (unsigned int word_order = 0; word_order < scan.getNumberOfWords(); ++word_order)
  {
    const int word = scan.word(word_order);
    const double word_idf = idf(word);

    BOOST_FOREACH (size_t scan_id, occurence_list_[word])
    {
      BOOST_FOREACH (const int other_order, scans_[scan_id].orders(word))
      {
        offset_hitograms[scan_id][word_order - other_order].second += word_idf;
        offset_hitograms[scan_id][word_order - other_order].first++;
      }
    }
  }

  OffsetHistograms::const_iterator iter = offset_hitograms.begin();

  for (; iter != offset_hitograms.end(); ++iter)
  {
    double score = 0;
    BOOST_FOREACH (const OffsetHistogram::mapped_type& offset_score, iter->second | boost::adaptors::map_values)
    {
      if (offset_score.first < gp_kernel_size_)
      {
        continue;
      }

      score += offset_score.second * cached_binomial_coeff[offset_score.first - 1];
    }

    score /= (scans_[iter->first].getGpNorm() * scan.getGpNorm());
    score_set.push_back(std::make_pair(1 - score, iter->first));
  }

  std::sort(score_set.begin(), score_set.end());
}


void GeometricalPhrasesMatcher::cache_binomial_coeff()
{
  cached_binomial_coeff.resize (DEFAULT_CACHEBINOMIAL,0);

  const double kernel_size = static_cast<double>(gp_kernel_size_ - 1);

  for(unsigned int i = gp_kernel_size_ - 1 ; i < DEFAULT_CACHEBINOMIAL; i++)
  {
    cached_binomial_coeff[i] = boost::math::binomial_coefficient <double>(i, kernel_size);
  }
}


