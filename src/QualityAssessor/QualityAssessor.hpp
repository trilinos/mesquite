/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      

    (2006) kraftche@cae.wisc.edu    
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file QualityAssessor.hpp

Header file for the Mesquite::QualityAssessor class

  \author Thomas Leurent
  \date   2002-05-01
  \author Jason Kraftcheck
  \date   2005-03-09
 */


#ifndef MSQ_QUALITYASSESSOR_HPP
#define MSQ_QUALITYASSESSOR_HPP

#include <math.h>

#include "Mesquite.hpp"
#include "Instruction.hpp"
#include "MeshInterface.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <list.h>
#  include <string.h>
#  include <vector.h>
#else
#  include <list>
#  include <string>
#  include <vector>
#endif

#ifdef MSQ_USE_OLD_IO_HEADERS
#  include <ostream.h>
#else
#  include <iosfwd>
#endif


namespace MESQUITE_NS 
{

   class QualityMetric;
   class MsqError;
   class ParallelHelper;

  /*! \class QualityAssessor

      \brief A QualityAssessor instance can be inserted into an 
      InstructionQueue to calculate and summarize registered
      QualityMetric or QualityMetrics for the mesh.

      QualityAssessor provides a summary of the mesh quality.  An 
      instance of QualityAssessor may be inserted in the InstructionQueue
      at any point to print a summary of the mesh quality at that
      time during the optimization.  The application is expected to
      register QualityMetric instances to be used in assessing the
      mesh quality.  If no QualityMetrics are registered, the only
      assessment that will be performed is a simple count of inverted
      elements.

      The "stopping assessor" and "stopping function", if set,
      determinte the value reported to Mesquite for the overall
      run of of the QualityAssessor.
      
      All summary data except the histogram and custom power-means is 
      accumulated for all registered metrics.  QualityAssessment data
      can be obtained through three different mechanisms:
        - QualityAssessor can print a summary to std::out or a specified
          output stream after they are calculated.
        - The get_results and get_all_results methods can be used to
          obtain the sumary data programatically.
        - Quality measures for element-based or vertex-based metrics can
          be stored on the corresponding entities using "tags".
      
  */
  class QualityAssessor : public Instruction
  {
  public:
    
    /**\brief Simple consturctor.  Metrics registered separately.
     *\param print_summary_to_std_out If true, summary of mesh quality
     *                will be written to std::out.  If false, quality
     *                assessment will be available via the get_results
     *                and get_all_results methods, but will not be printed.
     *\param free_elements_only If true, only quality values that depend
     *                on at least one free vertex will be uses.
     *\param  inverted_element_tag_name  If a non-null value is specified,
     *                an integer tag with the specified name will be used
     *                it store value of 0 for normal elements and 1 for
     *                inverted elements.
     *\param  name    Name to include in output.  Useful if several QualityAssessors
     *                are in use at the same time.
     */
    MESQUITE_EXPORT 
    QualityAssessor( bool print_summary_to_std_out = true,
                     bool free_elements_only = true,
                     const char* inverted_element_tag_name = 0,
                     msq_std::string name = msq_std::string() );

    /**\brief Simple consturctor.  Metrics registered separately.
     *\param output_stream IO stream to which to write a summary of the 
     *                mesh quality.
     *\param free_elements_only If true, only quality values that depend
     *                on at least one free vertex will be uses.
     *\param  inverted_element_tag_name  If a non-null value is specified,
     *                an integer tag with the specified name will be used
     *                it store value of 0 for normal elements and 1 for
     *                inverted elements.
     *\param  name    Name to include in output.  Useful if several QualityAssessors
     *                are in use at the same time.
     */
    MESQUITE_EXPORT 
    QualityAssessor( msq_stdio::ostream& output_stream,
                     bool free_elements_only = true,
                     const char* inverted_element_tag_name = 0,
                     msq_std::string name = msq_std::string() );
                     
    /**\brief Construct and register a QualityMetric
     *\param output_stream IO stream to which to write a summary of the 
     *                mesh quality.
     *\param metric   QualtiyMetric to register for use in assessing mesh
     *                quality.  Will also be used as the "stopping assessment".
     *\param historgram_intervals If non-zero, a histogram of quality metric
     *                values composed of the specified number of intervals
     *                will be generated.
     *\param power_mean If non-zero, in addition to the normal summary
     *                statistics for the quality metric, an additional
     *                general power mean with the specified power will
     *                be calculated.
     *\param metric_value_tag_name  If a non-null value is specified,
     *                a tag with the specified name will be populated
     *                with quality values for individual elements or
     *                vertices if metric is an element-based or vertex-
     *                based metric.  If metric is not element-based or
     *                vertex-based, this argument has no effect.
     *\param free_elements_only If true, only quality values that depend
     *                on at least one free vertex will be uses.
     *\param inverted_element_tag_name  If a non-null value is specified,
     *                an integer tag with the specified name will be used
     *                it store value of 0 for normal elements and 1 for
     *                inverted elements.
     *\param  name    Name to include in output.  Useful if several QualityAssessors
     *                are in use at the same time.
     */
    MESQUITE_EXPORT 
    QualityAssessor( msq_stdio::ostream& output_stream,
                     QualityMetric* metric, 
                     int histogram_intervals = 0,
                     double power_mean = 0.0,
                     bool free_elements_only = true,
                     const char* metric_value_tag_name = 0,
                     const char* inverted_element_tag_name = 0,
                     msq_std::string name = msq_std::string() );

                     
    /**\brief Construct and register a QualityMetric
     *\param print_summary_to_std_out If true, summary of mesh quality
     *                will be written to std::out.  If false, quality
     *                assessment will be available via the get_results
     *                and get_all_results methods, but will not be printed.
     *\param metric   QualtiyMetric to register for use in assessing mesh
     *                quality.  Will also be used as the "stopping assessment".
     *\param historgram_intervals If non-zero, a histogram of quality metric
     *                values composed of the specified number of intervals
     *                will be generated.
     *\param power_mean If non-zero, in addition to the normal summary
     *                statistics for the quality metric, an additional
     *                general power mean with the specified power will
     *                be calculated.
     *\param metric_value_tag_name  If a non-null value is specified,
     *                a tag with the specified name will be populated
     *                with quality values for individual elements or
     *                vertices if metric is an element-based or vertex-
     *                based metric.  If metric is not element-based or
     *                vertex-based, this argument has no effect.
     *\param free_elements_only If true, only quality values that depend
     *                on at least one free vertex will be uses.
     *\param  inverted_element_tag_name  If a non-null value is specified,
     *                an integer tag with the specified name will be used
     *                it store value of 0 for normal elements and 1 for
     *                inverted elements.
     *\param  name    Name to include in output.  Useful if several QualityAssessors
     *                are in use at the same time.
     */
    MESQUITE_EXPORT 
    QualityAssessor( QualityMetric* metric, 
                     int histogram_intervals = 0,
                     double power_mean = 0.0,
                     bool free_elements_only = true,
                     const char* metric_value_tag_name = 0,
                     bool print_summary_to_std_out = true,
                     const char* inverted_element_tag_name = 0,
                     msq_std::string name = msq_std::string() );

    MESQUITE_EXPORT virtual ~QualityAssessor();
    
      //! Provides a name to the QualityAssessor (use it for default name in constructor).
    MESQUITE_EXPORT void set_name(msq_std::string name) { qualityAssessorName = name; };
      //! Retrieves the QualityAssessor name. A default name should be set in the constructor.
    MESQUITE_EXPORT virtual msq_std::string get_name() const { return qualityAssessorName; }
    
      /**\brief All elements or only improvable ones.
       *
       * If set to true, the quality assessment results will include
       * quality values only for those elements (or more precisely metric
       * sample points) which are influenced by at least one free vertex.
       * That is, quality for elements (or sample points) that the sovler
       * cannot improve (e.g. an element with all vertices fixed) will not
       * be included in the quality assessment.
       *
       * If set to false, quality for all elements will be assessed.
       */
    MESQUITE_EXPORT void measure_free_samples_only( bool yesno )
      { skipFixedSamples = yesno; }
      
      /**\brief All elements or only improvable ones.
       *
       * If set to true, the quality assessment results will include
       * quality values only for those elements (or more precisely metric
       * sample points) which are influenced by at least one free vertex.
       * That is, quality for elements (or sample points) that the sovler
       * cannot improve (e.g. an element with all vertices fixed) will not
       * be included in the quality assessment.
       *
       * If set to false, quality for all elements will be assessed.
       */
    MESQUITE_EXPORT bool measuring_free_samples_only() const
      { return skipFixedSamples; }
    
    
    /**\brief Register a QualityMetric for use in quality assessment.
     *
     * Add a quality metric to the list of metrics used to assess
     * the quality of the mesh.
     *
     *\param metric   QualtiyMetric to register for use in assessing mesh
     *                quality.  Will also be used as the "stopping assessment".
     *\param historgram_intervals If non-zero, a histogram of quality metric
     *                values composed of the specified number of intervals
     *                will be generated.
     *\param power_mean If non-zero, in addition to the normal summary
     *                statistics for the quality metric, an additional
     *                general power mean with the specified power will
     *                be calculated.
     *\param metric_value_tag_name  If a non-null value is specified,
     *                a tag with the specified name will be populated
     *                with quality values for individual elements or
     *                vertices if metric is an element-based or vertex-
     *                based metric.  If metric is not element-based or
     *                vertex-based, this argument has no effect.
     *\param  inverted_element_tag_name  If a non-null value is specified,
     *                an integer tag with the specified name will be used
     *                it store value of 0 for normal elements and 1 for
     *                inverted elements.
     */
    MESQUITE_EXPORT 
    void add_quality_assessment( QualityMetric* metric,
                                 int histogram_intervals = 0,
                                 double power_mean = 0.0,
                                 const char* metric_value_tag_name = 0 );
    
    /**\brief Same as add_quality_assessment, except that the average
     *        metric value is also used as the return value from loop_over_mesh.
     *        Specify a power_mean value to control which average is used.
     */
    MESQUITE_EXPORT 
    void set_stopping_assessment( QualityMetric* metric,
                                  int histogram_intervals = 0,
                                  double power_mean = 0.0,
                                  const char* metric_value_tag_name = 0 );

    /**\brief Register a QualityMetric for use in quality assessment.
     *
     * Add a quality metric to the list of metrics used to assess
     * the quality of the mesh.  Specify more parameters controlling
     * histogram.
     *
     *\param metric   QualtiyMetric to register for use in assessing mesh
     *                quality.  Will also be used as the "stopping assessment".
     *\param min      Minimum of histogram rnage.
     *\param max      Maximum of histogram range.
     *\param intervals Histogram intervals.
     *\param power_mean If non-zero, in addition to the normal summary
     *                statistics for the quality metric, an additional
     *                general power mean with the specified power will
     *                be calculated.
     *\param metric_value_tag_name  If a non-null value is specified,
     *                a tag with the specified name will be populated
     *                with quality values for individual elements or
     *                vertices if metric is an element-based or vertex-
     *                based metric.  If metric is not element-based or
     *                vertex-based, this argument has no effect.
     *\param  inverted_element_tag_name  If a non-null value is specified,
     *                an integer tag with the specified name will be used
     *                it store value of 0 for normal elements and 1 for
     *                inverted elements.
     */
    MESQUITE_EXPORT 
    void add_histogram_assessment( QualityMetric* qm, 
                                   double min, 
                                   double max,
                                   int intervals,
                                   double power_mean = 0.0,
                                   const char* metric_value_tag_name = 0 );
    
      //! Does one sweep over the mesh and assess the quality with the metrics previously added.
    virtual double loop_over_mesh( Mesh* mesh,
                                   MeshDomain* domain,
                                   const Settings* settings,
                                   MsqError &err);

      //! Does one sweep over the mesh and assess the quality with the metrics previously added.
    virtual double loop_over_mesh( ParallelMesh* mesh,
                                   MeshDomain* domain,
                                   const Settings* settings,
                                   MsqError &err);

      //! Do not print results of assessment.
    MESQUITE_EXPORT void disable_printing_results()
       {
         printSummary = false;
       }
      
      //! Print accumulated summary data to specified stream. 
    MESQUITE_EXPORT void print_summary( msq_stdio::ostream& stream ) const;
    
      //! True if any metric evaluated to an invalid value
      //! for any element
    MESQUITE_EXPORT bool invalid_elements() const;

      //! Provides the number of inverted elements, inverted_elmes,
      //!  and the number of elements whose orientation can not be
      //!  determined, undefined_elems.
      //! Returns false if this information is not yet available.
      //! Returns true, otherwise.
    MESQUITE_EXPORT bool get_inverted_element_count(int &inverted_elems,
                                    int &undefined_elems,
                                    MsqError &err);
    
      //! Reset calculated data 
    MESQUITE_EXPORT void reset_data();
    
    MESQUITE_EXPORT void tag_inverted_elements( std::string tagname ) 
      { invertedTagName = tagname; }
    MESQUITE_EXPORT void dont_tag_inverted_elements()
      { invertedTagName.clear(); }
    MESQUITE_EXPORT bool tagging_inverted_elements() const
      { return !invertedTagName.empty(); }
    
    /** \brief Per-metric QualityAssessor data
     *
     * The Assessor class holds QualityAssessor data for
     * each metric added by the calling application, including
     * a pointer to the metric instance, QAFunction flags
     * dictating what is to be calculated and output, histogram
     * parameters, and the variables used to accumulate results as
     * the QualityAssessor is running.  It also provides 
     * methods to access the calculated data once the QualityAssessor
     * pass is completed.
     */
    class Assessor
    {
      public:
      
        Assessor( QualityMetric* metric );
        
        double get_average() const ;
        double get_maximum() const { return maximum; }
        double get_minimum() const { return minimum; }
        double get_rms()     const ;
        double get_stddev()  const ;
        double get_power_mean() const;
        double get_power() const   { return pMean; } 
        int get_count() const { return count; }
        
        bool have_power_mean() const { return 0.0 != pMean; }
        
        int get_invalid_element_count() const { return numInvalid; }
        
        /** Get historgram of data, if calculated.
         *\param lower_bound_out  The lower bound of the histogram
         *\param upper_bound_out  The upper bound of the histogram
         *\param counts_out       An array of counts of elements where
         *              the first entry is the number of elements for
         *              which the metric is below the lower bound, the
         *              last entry is the number of elements above the
         *              upper bound, and all other values are the counts
         *              for histogram intervals between the lower and
         *              upper bounds.
         */
        void get_histogram( double& lower_bound_out,
                            double& upper_bound_out,
                            msq_std::vector<int>& counts_out,
                            MsqError& err ) const;
                            
        /** Reset all calculated data */
        void reset_data();
       
        /** Print the histogram */
        void print_histogram( msq_stdio::ostream&, int width = 0 ) const;

        /** Get the QualityMetric */
        QualityMetric* get_metric() const { return qualMetric; }
        
        /** Add a value to the running counts */
        void add_value( double metric_value );
        
        /** Add a value to the hisogram data */
        void add_hist_value( double metric_value );
        
        /** Note invalid result */
        void add_invalid_value() ;
        
        bool have_histogram() const
          { return !histogram.empty(); }
        
        /** If range of histogram has not yet been determined,
          * calculate it from the min/max values.
          */
        void calculate_histogram_range();
        
        bool write_to_tag() const { return !tagName.empty(); }
        
        void set_stopping_function( bool value )
          { stoppingFunction = value; }
        
        bool stopping_function( ) const
          { return stoppingFunction; }
        
        double stopping_function_value() const;
        
      private:
      
        friend class QualityAssessor;
        
        QualityMetric *const qualMetric; //< The quality metric
        
        unsigned long count;  //< The total number of times the metric was evaluated
        
        double sum;       //< The sum of the metric over all elements
        double maximum;   //< The maximum of the metric
        double minimum;   //< The minimum value of the metric
        double sqrSum;    //< The sum of the square of the metric values
        double pSum;      //< The sum of the metric values raised to the pMean power
        unsigned long numInvalid;  //< Count of invalid metric values
        
        double pMean;     //< Power for general power-mean.
        
        
        /** The histogram counts, where the first and last values are
         * counts of values below the lower bound and above the upper
         * bound, respectively.  The remaining values are the histogram
         * counts.
         */
        bool haveHistRange;
        double histMin;   //< Lower bound of histogram
        double histMax;   //< Upper bound of histogram
        msq_std::vector<int> histogram;
        
        std::string tagName; //< Tag to which to write metric values.
        /** Cached tag handle */
        TagHandle tagHandle;
        
        /** Value is return value for all of QualityAssessor */
        bool stoppingFunction;
    };    
        
    typedef msq_std::list<Assessor> list_type;
        
    /** \brief Request summary data for a specific QualityMetric 
     * This method allows the application to request the summary
     * data for a metric it has registered with the QualityAssessor.
     * If the passed QualityMetric has not been registered with the
     * QualityAssessor instance, NULL is returned.
     */
    const Assessor* get_results( QualityMetric* metric ) const;
    
    /** \brief Get list of all summary data.
     *  Return a const reference to the internal list of 
     *  calculated data.
     */
   const list_type& get_all_results() const
      { return assessList; }
      
  private:

      //! Common code for serial and parallel loop_over_mesh
    double loop_over_mesh_internal( Mesh* mesh,
                                    MeshDomain* domain,
                                    const Settings* settings,
                                    ParallelHelper* helper,
                                    MsqError &err);
  
    /** Find an Assessor corresponding to the passed
     *  QualityMetric, or create it if is not found in
     *  the list.
     */
    list_type::iterator find_or_add( QualityMetric* qm );

    /** Find an Assessor corresponding to the passed
     *  QualityMetric, or create it if is not found in
     *  the list.
     */
    list_type::iterator find_stopping_assessment();
   
    /** Find or create tag */
    TagHandle get_tag( Mesh* mesh,
                       std::string name, 
                       Mesh::TagType type, 
                       unsigned size, 
                       MsqError& err );
   
    /** Try to determine if output stream is a terminal and if so, 
        the width of that terminal.  Returns zero if width cannot
        be determined. */
    int get_terminal_width() const;
   
    /** Name */
    msq_std::string qualityAssessorName;  
    
    /** List of quality metrics and corresponding data */
    list_type assessList;

      /** Count of inverted elements. */
    int invertedCount;

      /** Count of elements whose orientation can not be determined.*/
    int indeterminateCount;
    
      /** Number of elements */
    size_t elementCount;
    
      /** Number of elements with at least one free vertex */
    size_t freeElementCount;
   
    /** Stream to which to write summary of metric data */
    msq_stdio::ostream& outputStream;
    /** Disable printing */
    bool printSummary;
    
    std::string invertedTagName;
    
    bool skipFixedSamples;
  };

  
} //namespace


#endif // QualityAssessor_hpp
