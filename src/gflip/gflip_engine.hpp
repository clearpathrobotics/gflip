//
//
// GFLIP - Geometrical FLIRT Phrases for Large Scale Place Recognition
// Copyright (C) 2012-2013 Gian Diego Tipaldi and Luciano Spinello and Wolfram
// Burgard
//
// This file is part of GFLIP.
//
// GFLIP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GFLIP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GFLIP.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef GFP_NGN_H
#define GFP_NGN_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <set>
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/foreach.hpp>
#include <sys/time.h>  
#include <algorithm>
#include <map>


//~ Basic and bag of distances  defaults   
#define DEFAULT_BOWDST_START 0
#define	DEFAULT_BOWDST_INTERVAL 0.2
#define DEFAULT_BOWDST_END 15
#define DEFAULT_BOWSUBTYPE 0
#define DEFAULT_ALPHASMOOTH 0.4
#define DEFAULT_BAGDISTANCE 0
#define DEFAULT_CACHEBINOMIAL 10000

/**
 * Contains a 2D scan represented by FLIRT words identified by their index, their (TF-IDF) weights, their norm for GFP
 * 
 * @author Luciano Spinello
 */					 
 
class scan_bow
{
public:
  std::vector <int> w, word_weight_unnormalized;
  std::vector <double> w_x, w_y, word_weight, tfidf_w;
  std::map<int, std::vector<int> > histogram;
  std::map<int, std::vector<int> > normgfp_rc;
  //std::map<int, double> normgfp_rc;
  std::map<int, int> word_ref;
  double sum_weight,  norm_wgv;
  int max_histogram;
  scan_bow(int no)
  {
    w= std::vector <int> (no);
    w_x = std::vector <double> (no);
    w_y = std::vector <double> (no);
    max_histogram = 0;
  }
  void compute_word_weight()
  {
    word_weight_unnormalized = std::vector<int>(w.size(),0.0);
    word_weight = std::vector<double>(w.size(),0.0);
    int index = 0;
    sum_weight = 0.0;
    BOOST_FOREACH (int word, w)
    {
      int j=0;
      BOOST_FOREACH(int other, w)
      {

        if(other==word)
        {
          normgfp_rc[index-j].push_back(word);
          word_weight_unnormalized[index]++;
          if (std::find(histogram[word].begin(), histogram[word].end(), j) == histogram[word].end())
            histogram[word].push_back(j);
        }
        ++j;
      }
      if (word_weight_unnormalized[index] > max_histogram)
        max_histogram = word_weight_unnormalized[index];
      sum_weight += word_weight_unnormalized[index];
      word_weight[index] = word_weight_unnormalized[index];
      index++;
    }
    BOOST_FOREACH(double weight, word_weight)
        weight /= sum_weight;
  }

};

 

/**
 * Caches the word orders in a scan for GFP indexing
 * 
 * @author Luciano Spinello
 */	
class tf_idf_db_ordercache
{
	public:
		std::vector<int> pos;
    tf_idf_db_ordercache(std::vector<int> p = std::vector<int>()): pos(p){}
};


/**
 * Contains TF-IDF weight for a  FLIRT word
 * 
 * @author Luciano Spinello
 */	
class tf_idf_db
{
	public:
		//~ per doc
		std::vector <tf_idf_db_ordercache> word_order;
		std::vector <int> doc_id;
		std::vector <int> term_count_unnormalized;
		std::vector <double> tf_idf_doc_normed;
		std::vector <double> ntf_idf_doc_normed;
		std::vector <double> wf_idf_doc_normed;
		std::vector <int> num_words;
		std::vector <double> term_count;
		
		//~ per term
		int num_doc_containing_the_word, corpus_size;
		double idf;
		
		tf_idf_db()
		{
			num_doc_containing_the_word = 0;
			corpus_size = 0;
		}

                void compare(tf_idf_db other)
                {
                  std::vector<tf_idf_db_ordercache>::const_iterator db = word_order.begin(),
                      other_db = other.word_order.begin();
                  for(; db != word_order.end() && other_db != other.word_order.end(); ++db, ++other_db)
                    assert(db->pos == other_db->pos);
                  assert(doc_id == other.doc_id);
                  assert(term_count_unnormalized == other.term_count_unnormalized);
                  assert(term_count == other.term_count);
                  assert(num_words == other.num_words);
                  assert(corpus_size == other.corpus_size);
                  assert(idf == other.idf);
                  assert(ntf_idf_doc_normed == other.ntf_idf_doc_normed);
                  assert(tf_idf_doc_normed == other.tf_idf_doc_normed);
                  assert(wf_idf_doc_normed == other.wf_idf_doc_normed);
                  assert(num_doc_containing_the_word == other.num_doc_containing_the_word);
                }
};

/**
 * Geometrical FLIRT Phrases (GFP) for matching 2D laser scans represented FLIRT words
 * 
 * It includes methods for building the search index and matching
 * 
 * <a href="http://www.informatik.uni-freiburg.de/~spinello/tipaldiICRA13.pdf">"Geometrical FLIRT Phrases for Large Scale Place Recognition in 2D Range Data", G. D. Tipaldi, L. Spinello, W. Burgard -- Int. Conf. Robotics and Automation (ICRA) 2013</a> 
 * @author Luciano Spinello
 */		
class gflip_engine
{
	private:
		//~ vars
		std::vector <scan_bow> laserscan_bow;
		std::vector < std::pair <double, int> > scoreset;
		std::vector <tf_idf_db> tf_idf;
		std::string fileoutput_rootname;
		int dictionary_dimensions, start_l, stop_l, max_bow_len, wgv_kernel_size, bow_type, bow_subtype;
		double anglethres, bow_dst_start, bow_dst_interval, bow_dst_end, alpha_vss;
		uint number_of_scans, kbest;
		std::vector<double> cached_binomial_coeff, mtchgfp_rc_idf_sum, normgfp_rc_idf_sum;
		std::vector <int> mtchgfp_min_det_idx, mtchgfp_max_det_idx, mtchgfp_rc_weak_match, normgfp_rc_weak_match;
		std::vector<char> mtchgfp_used_doc_idx;

                // contains the word with the maximum frequency in each doc
                std::vector < double > mxtf_val;

		//~ functions
 		double norm_gfp(std::vector <int> & query_v);
                double norm_gfp(scan_bow & bow);
 		void matching_bow(std::vector <int> &query_v );
		void matching_gfp(std::vector <int> &query_v );
		void voting_tfidf_weak_verificationOLD(std::vector <int> &query_v );		
		void reformulate_to_bagofdistances(void);
		void cache_binomial_coeff(void);
                void add_doc_stats(int doc_id);
                void compute_idfs();
                void compute_tf_idfs();

	
	public:

		//~ functions		
		/**
		 * Reads file generated by FLIRTLIB in which each scan is described as a sequence of FLIRT words, represented each by a number
		 * 
		 * @author Luciano Spinello
		 */		 
 		int read_wordscan_file(std::string filename);

		/**
		 * Inserts a scan described as a sequence of FLIRT words, represented each by a number 
		 * @param wordscan scan identified as a sequence of ids
		 * @param xpos,ypos metric position of each word in \c wordscan
		 * @author Luciano Spinello
		 */		 
 		void insert_wordscan(std::vector <int> wordscan, std::vector <double> xpos, std::vector <double> ypos);
 		

		/**
		 * Builds TF-IDF index for standard and weak verification matching methods
		 * 
		 * It implements IDF, TF, TF-IDF, wordcount and improved TF-IDF models proposed in:
		 * <a href="http://comminfo.rutgers.edu/~muresan/IR/Docs/Articles/ipmSalton1988.pdf"> Gerard Salton and Christopher Buckley: "Term-weighting approaches in automatic text retrieval", Information processing & management, vol. 24, no. 5, 1988, Elsevier </a>
		 * @author Luciano Spinello
		 */		
 		void build_tfidf(void);

                /**
                 * @brief update_tfidf
                 * @param number of new scans added to laserscan_bow that are not yet
                 * considered in the tf_idf
                 */
                void update_tfidf(int num_scans);

                void build_tfidf_new();
                void build_tfidf_new1();

		/**
		 * Matches all the scans in the dataset vs all the scans in the dataset
		 * 
		 * Saves the k-best results on disk for each query along with computational time 
		 * 
		 * @param dtype  kind of matching method: 1 standard bag-of-words, 2 geometrical FLIRT phrases
		 * @author Luciano Spinello
		 */
		void run_evaluation(int dtype);
		

		/**
		 * Matches a query scan with the dataset
		 * 
		 * Saves the k-best results on disk for each query along with computational time 
		 * 
		 * @param dtype  kind of matching method: 1 standard bag-of-words, 2 geometrical FLIRT phrases
		 * @param query_v a query scan, composed by a vector of numbers, each indicating a FLIRT word
		 * @param scoreoutput a pointer to a sorted vector of pairs containing <scorematch, index of the scan in the dataset>. 
		 * @author Luciano Spinello
		 */
		void query(int dtype, std::vector <int>   &query_v, std::vector < std::pair <double, int> > **scoreoutput);


		/**
		 * Prepares indeces and cache for matching. Executed once at the beginning.
		 * Builds TF-IDF index for the dataset and norms all vectors on the dataset. Allocates also needed memory.
		 * Optionally it generates bag-of-distances
		 * @author Luciano Spinello
		 */
		void prepare(void);

		/**
		 * Constructor
		 * 
		 * @param krnl Kernel size
		 * @param kbt # of k-best results
		 * @param bt 1 for bag-of-distances, 0 otherwise
		 * @param bstype flavor of TF-IDF in case of standard bag-of-words: 0 standard TFIDF, 1 sublinear TFIDF scaling, 2 lenght smoothing TFIDF, see \link gflip_engine::build_tfidf\endlink
		 * @param a_vss alpha_smoothing in case of standard bag-of-words with lenght smoothing TFIDF (0.4 default)
		 * @author Luciano Spinello
		 */  

                gflip_engine (int krnl, int kbt, int bt=DEFAULT_BAGDISTANCE, int bstype=DEFAULT_BOWSUBTYPE, double a_vss=DEFAULT_ALPHASMOOTH)
		{
			bow_type = bt;
			kbest = kbt;
			wgv_kernel_size = krnl;
 			bow_subtype= bstype;
			alpha_vss = a_vss;

			
			//~ basic defaults for bag of distances
			bow_dst_start= DEFAULT_BOWDST_START;
		    bow_dst_interval=DEFAULT_BOWDST_INTERVAL;
		    bow_dst_end=DEFAULT_BOWDST_END;

		}

                void setDictionaryDimension(int dictionary_size)
                {
                  dictionary_dimensions = dictionary_size;
                  std::cout<<"Dictionary dimensions "<<dictionary_dimensions<<std::endl;
                  max_bow_len = (dictionary_dimensions+1)*2;
                  tf_idf = std::vector <tf_idf_db> (dictionary_dimensions);
                }
		
};

/**
 * Simple string tokenizer (avoids additional dependencies)
 * 
 * @author Luciano Spinello
 */	
void LSL_stringtoken(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters);
bool isBettermatched(std::pair <double, int> x, std::pair <double, int> y); 
#endif
