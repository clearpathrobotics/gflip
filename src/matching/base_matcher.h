#ifndef GFP_BASE_MATCHER_H
#define GFP_BASE_MATCHER_H

#include <boost/serialization/set.hpp>
#include "laser_scan_data.h"

#define DEFAULT_CACHEBINOMIAL 10000

#define DEFAULT_ID std::numeric_limits<size_t>::max()

typedef std::vector<std::pair<double, size_t> > ScoreSet;

class BaseMatcher
{

public:
  BaseMatcher(const size_t  vocabulary_size)
    : occurence_list_(vocabulary_size)
    , addition_update_threshold_(1) /* update happens after every addition */
    , removal_update_threshold_(1)  /* update happens after every removal  */
    , removed_scans_(0)
  {}

  /**
   * @brief Add a scan to the database
   * @param[in] words features of the added scan mapped into vocabulary words
   * @param[in] type LaserScanData::GEOMETRICAL_PHRASES (default) or LaserScanData::STANDARD_BOW
   * @param[in] id id of the scan. Defaults to the addition order
   */
  void addScan(const std::vector <int>& words,
               const size_t id = DEFAULT_ID);

  /**
   * @brief Remove a scan from the database
   * @param[in] id id of the scan to be removed
   */
  void removeScan(const size_t id);

  /**
   * @brief Return the inverse document frequency of a word
   * @param[in] word  Word for which to get the idf
   * @return inverse document frequency of all the documents containing the word
   */
  double idf(const int word) const
  {
    return idfs_[word];
  }

  /**
   * @brief Generates the idf for all the stored scans
   */
  void updateIdfs()
  {
    const double corpus_size = static_cast<double>(scans_.size());

    for (int i=0; i<occurence_list_.size(); ++i)
    {
      idfs_[i] = log(corpus_size / occurence_list_[i].size());
    }
  }

  void getCandidates(const std::vector<int>& scan_words, ScoreSet& score_set,
                     const bool insert, const int id = DEFAULT_ID);

  /**
   * @brief The actual matching function
   * @param scan
   * @param score_set
   */
  virtual void getCandidates(const LaserScanData& scan, ScoreSet& score_set) = 0;

  virtual LaserScanData::Ptr getLaserScanData(const std::vector<int>& scan_words) const = 0;

  /**
   * @brief updateLaserScanData
   */
  virtual void updateLaserScanData() = 0;

  /**
   * @brief update updates the internal structures with the added and removed data
   */
  void update();

  /**
   * @brief Threshold on the number of scans to be removed before an update happens
   * @param removal_update_threshold  Defaults to 1. if 0 update won't happen unless
   * called explicitly
   */
  void setRemovalUpdateThreshold(const unsigned int removal_update_threshold)
  {
    removal_update_threshold_ = removal_update_threshold;
  }

  /**
   * @brief Threshold on the number of scans to be added before an update happens
   * @param addition_update_threshold  Defaults to 1. if 0 update won't happen unless
   * called explicitly
   */
  void setAdditionUpdateThreshold(const unsigned int addition_update_threshold)
  {
    addition_update_threshold_ = addition_update_threshold;
  }

  /** Serialize the matcher to an archive
  * @param o_archive output archive to serialize to
  */
  void save(boost::archive::binary_oarchive& o_archive) const;

  /**
   * Loads the matcher from a boost archive
   * @param i_archive input archive to load from
   */
  void load(boost::archive::binary_iarchive& i_archive);

protected:
  /**
   * @brief Information of the scans as bag of words
   */
  std::map<size_t, LaserScanData> scans_;

  /**
   * @brief Maintains a list of scans in which each word occurs
   */
  std::vector<std::vector<size_t> > occurence_list_;

  /**
   * @brief container to cache the idfs of words;
   */
  std::vector<double> idfs_;

  /**
   * @brief container to hold the added scans that are not yet being used
   */
  std::vector<size_t> added_scans_;


  /**
   * @brief System update will happen when addition_update_threshold_ scans are added
   */
  unsigned int addition_update_threshold_;

  /**
   * @brief number of removed scans
   */
  unsigned int removed_scans_;

  /**
   * @brief System update will happen when removal_update_threshold_ scans are removed
   */
  unsigned int removal_update_threshold_;
};

#endif //GFP_BASE_MATCHER_H
