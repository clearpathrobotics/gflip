#ifndef GFP_LASER_SCAN_DATA_H
#define GFP_LASER_SCAN_DATA_H

#include <algorithm>
#include <map>
#include <vector>
#include <set>
#include <limits.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/any_range.hpp>
#include <boost/foreach.hpp>


/** Class to generate and store bag of words data
  * pertaining to a single scan
  */
class LaserScanData
{
public:

  typedef boost::shared_ptr<LaserScanData> Ptr;

  typedef std::pair<double, double> Weights;

  enum Type
  {
    GEOMETRICAL_PHRASES,
    STANDARD_BOW
  };

  /**
   * @brief Constructor
   * @param type GEOMETRICAL_PHRASES or STANDARD_BOW
   */
  LaserScanData(const Type type = GEOMETRICAL_PHRASES )
    : max_histogram_(0)
    , type_(type)
  {}

  /**
   * @brief LaserScanData
   * @param words list of features in the scan mapped into the vocabulary space
   * @param type
   */
  LaserScanData(const std::vector<int>& words, const Type type = GEOMETRICAL_PHRASES )
    : words_(words)
    , max_histogram_(0)
    , type_(type)
  {}

  /**
   * @brief init computes the initial statistics of the scan words
   */
  void init();

  /**
   * @brief Total number of words
   * @return the total number of words in the scan
   */
  size_t getNumberOfWords() const
  {
    return words_.size();
  }

  /**
   * @brief word
   * @param i
   * @return
   */
  int word(size_t i) const
  {
    return words_[i];
  }

  /**
   * @brief Geometrical Phrases norm computed from the offset histgram
   * @return
   */
  double & getGpNorm()
  {
    return gp_norm_;
  }

  /**
   * @brief Geometrical Phrases norm computed from the offset histgram
   * @return
   */
  double getGpNorm() const
  {
    return gp_norm_;
  }

  /**
   * @brief Returns a pair of weights
   * @param word
   * @return
   */
  Weights getWeight(const int word) const
  {
    // Note: operator [] can not be used as it returns a non-const ref
    return weights_.find(word)->second;
  }

  std::set<size_t>& orders(const int word)
  {
    return orders_[word];
  }

  double& tf_idf(int word)
  {
    return tf_idfs_[word];
  }

  int getMaxHistogram()
  {
    return max_histogram_;
  }


  boost::any_range<std::vector<int>, boost::forward_traversal_tag,  const std::vector<int>&, std::ptrdiff_t> getOffsetHistograms() const
  {
    return offsets_histogram_ | boost::adaptors::map_values;
  }


  /**
   * @brief Unique words in the scan
   * @return
   */
  std::vector<int> uniqueWords() const
  {
    return unique_words_;
  }

private:
  /**
   * Serialization stuff
   */
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive & ar, const unsigned int /* version */)
  {
    ar & orders_;
    ar & weights_;
    ar & offsets_histogram_;
    ar & tf_idfs_;
    ar & words_;
    ar & gp_norm_;
    ar & max_histogram_;
    ar & type_;
  }

  std::map<int, std::set<size_t> > orders_;
  std::map<int, std::pair<double, double> > weights_;

  std::map<int, std::vector<int> > offsets_histogram_;

  std::map<int, double> tf_idfs_;

  std::vector<int> words_;

  std::vector<int> unique_words_;

  double gp_norm_;
  size_t max_histogram_;
  Type type_;
};

#endif // GFP_LASER_SCAN_DATA_H
