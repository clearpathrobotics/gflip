#include "laser_scan_data.h"
#include <boost/bind.hpp>


void LaserScanData::init()
{
  if (type_ == GEOMETRICAL_PHRASES)
  {
    int i = 0;
    BOOST_FOREACH (int word, words_)
    {
      int j=0;
      BOOST_FOREACH(int other, words_)
      {
        if (other==word)
        {
          offsets_histogram_[i-j].push_back(word);

          orders_[word].insert(j);
        }
        ++j;
      }
    }
    // The unique words are the keys of orders_, cache them
    std::transform(orders_.begin(), orders_.end(), std::back_inserter(unique_words_),
                       boost::bind(&std::map<int, std::set<size_t> >::value_type::first,_1) );
  }

  else if (type_ == STANDARD_BOW)
  {
    BOOST_FOREACH (int word, words_)
    {
      BOOST_FOREACH(int other, words_)
      {
        if (other==word)
        {
          weights_[word].first++;
        }
      }
      if (weights_[word].first > max_histogram_)
      {
        max_histogram_ = weights_[word].first;
      }
    }

    BOOST_FOREACH (Weights& weights , weights_ | boost::adaptors::map_values)
    {
      weights.second = weights.first / words_.size();
    }

    // The unique words are the keys of weights_, cache them
    std::transform(weights_.begin(), weights_.end(), std::back_inserter(unique_words_),
                       boost::bind(&std::map<int, std::pair<double, double> >::value_type::first,_1) );
  }
}

