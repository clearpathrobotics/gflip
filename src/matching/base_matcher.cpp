#include "base_matcher.h"


void BaseMatcher::addScan(const std::vector <int>& words, const size_t id)
{
  const size_t scan_id = (id == DEFAULT_ID ? scans_.size() : id);

  LaserScanData::Ptr bow = getLaserScanData(words);

  bow->init();

  // Use insert instead of operator [] to avoid creating a default ScanBow
  scans_.insert(std::make_pair(scan_id, *bow));
  added_scans_.push_back(scan_id);

  if (added_scans_.size() & addition_update_threshold_)
  {
    update();
  }
}


void BaseMatcher::removeScan(const size_t id)
{
  // Remove the scan from the occurenceList
  BOOST_FOREACH( const int word, scans_[id].uniqueWords())
  {
    std::vector<size_t> & scans =  occurence_list_[word];
    scans.erase(std::remove(scans.begin(), scans.end(), id), scans.end());
  }

  // Remove the scan data
  scans_.erase(id);

  if (++removed_scans_ & removal_update_threshold_)
  {
    removed_scans_ = 0;
    update();
  }
}

void BaseMatcher::update()
{
  {
    // Add the scans that are in the added buffer to the
    BOOST_FOREACH (const size_t scan_id, added_scans_)
    {
      BOOST_FOREACH (const int word, scans_[scan_id].uniqueWords())
      {
        occurence_list_[scan_id].push_back(scan_id);
      }
    }

    added_scans_.clear();

    updateIdfs();

    updateLaserScanData();
  }
}


void BaseMatcher::getCandidates(const std::vector<int>& scan_words, ScoreSet &score_set, const bool insert, const int id)
{
  LaserScanData::Ptr bow = getLaserScanData(scan_words);

  bow->init();

  getCandidates(*bow, score_set);

  if (insert)
  {
    const size_t scan_id = (id == DEFAULT_ID ? scans_.size() : id);

    scans_.insert(std::make_pair(scan_id, *bow));

    added_scans_.push_back(scan_id);

    if (added_scans_.size() & addition_update_threshold_)
    {
      update();
    }
  }

}

void BaseMatcher::save(boost::archive::binary_oarchive& o_archive) const
{
  try
  {
    o_archive << scans_;
    o_archive << occurence_list_;
    o_archive << idfs_;
  }
  catch(std::exception)
  {
    throw;
  }
}


void BaseMatcher::load(boost::archive::binary_iarchive& i_archive)
{
  scans_.clear();
  occurence_list_.clear();
  idfs_.clear();
  try
  {
    i_archive >> scans_;
    i_archive >> occurence_list_;
    i_archive >> idfs_;
  }
  catch(std::exception)
  {
    scans_.clear();
    idfs_.clear();
    occurence_list_.clear();
    throw;
  }
}
