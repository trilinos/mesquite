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
   
  ***************************************************************** */
/*!
  \file   QualityAssessor.cpp
  \brief  Member function of the Mesquite::QualityAssessor class

  \author Thomas Leurent
  \date   2002-05-23
*/

#include "QualityAssessor.hpp"
#include "QualityMetric.hpp"
#include "PatchData.hpp"
#include "MsqMeshEntity.hpp"
#include "MsqVertex.hpp"
#include "MsqDebug.hpp"
#include "MeshInterface.hpp"
#include "VertexPatches.hpp"
#include "ElementPatches.hpp"
#include "ParallelMeshInterface.hpp"
#include "ParallelHelperInterface.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <list.h>
#  include <vector.h>
#else
#  include <list>
#  include <vector>
#endif

#ifdef MSQ_USE_OLD_IO_HEADERS
#  include <iostream.h>
#  include <fstream.h>
#  include <iomanip.h>
#  include <set.h>
#else
#  include <iostream>
#  include <fstream>
#  include <iomanip>
#  include <set>
#endif

#ifdef HAVE_SYS_IOCTL_H
# include <sys/ioctl.h>
#endif
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif
#ifdef HAVE_TERMIOS_H
# include <termios.h>
#endif

namespace MESQUITE_NS {

const char* default_name( bool free_only )
{
  static const char all_name[] = "QualityAssessor";
  static const char free_name[] = "QualityAssessor(free only)";
  return free_only ? free_name : all_name;
}

QualityAssessor::QualityAssessor( bool print_summary,
                                  bool free_only,
                                  const char* inverted_tag_name,
                                  msq_std::string name) :
  qualityAssessorName(name),
  invertedElementCount(-1),
  invertedSampleCount(-1),
  outputStream( msq_stdio::cout ),
  printSummary( print_summary ),
  skipFixedSamples(free_only)
{
  if (inverted_tag_name)
    tag_inverted_elements( inverted_tag_name );
    
  if (qualityAssessorName.empty()) 
    qualityAssessorName = default_name( free_only );
}

QualityAssessor::QualityAssessor( msq_stdio::ostream& stream,
                                  bool free_only,
                                  const char* inverted_tag_name,
                                  msq_std::string name) :
  qualityAssessorName(name),
  invertedElementCount(-1),
  invertedSampleCount(-1),
  outputStream( stream ),
  printSummary( true ),
  skipFixedSamples(free_only)
{
  if (inverted_tag_name)
    tag_inverted_elements( inverted_tag_name );
    
  if (qualityAssessorName.empty()) 
    qualityAssessorName = default_name( free_only );
}

QualityAssessor::QualityAssessor( msq_stdio::ostream& output_stream,
                                  QualityMetric* metric, 
                                  int histogram_intervals,
                                  double power_mean,
                                  bool free_only,
                                  const char* metric_value_tag_name,
                                  const char* inverted_tag_name,
                                  msq_std::string name ) :
  qualityAssessorName(name),
  invertedElementCount(-1),
  invertedSampleCount(-1),
  outputStream( output_stream ),
  printSummary( true ),
  skipFixedSamples(free_only)
{
  if (inverted_tag_name)
    tag_inverted_elements( inverted_tag_name );
    
  if (qualityAssessorName.empty()) 
    qualityAssessorName = default_name( free_only );
    
  set_stopping_assessment( metric, histogram_intervals, power_mean, metric_value_tag_name );
}

QualityAssessor::QualityAssessor( QualityMetric* metric, 
                                  int histogram_intervals,
                                  double power_mean,
                                  bool free_only,
                                  const char* metric_value_tag_name,
                                  bool print_summary,
                                  const char* inverted_tag_name,
                                  msq_std::string name) :
  qualityAssessorName(name),
  invertedElementCount(-1),
  invertedSampleCount(-1),
  outputStream( msq_stdio::cout ),
  printSummary( print_summary ),
  skipFixedSamples(free_only)
{
  if (inverted_tag_name)
    tag_inverted_elements( inverted_tag_name );
    
  if (qualityAssessorName.empty()) 
    qualityAssessorName = default_name( free_only );
    
  set_stopping_assessment( metric, histogram_intervals, power_mean, metric_value_tag_name );
}

QualityAssessor::~QualityAssessor()
  { }

double QualityAssessor::Assessor::get_average() const
{
  return count ? sum/count : 0;
}

double QualityAssessor::Assessor::get_rms() const 
{
  return count ? sqrt(sqrSum/count) : 0;
}

double QualityAssessor::Assessor::get_stddev() const
{
  double sqr = count ? sqrSum/count - sum*sum/((double)count*count) : 0;
  return sqr < 0 ? 0 : sqrt(sqr);
}

double QualityAssessor::Assessor::get_power_mean() const
{
  return (count && pMean) ? pow(pSum/count, 1/pMean) : 0;
}

bool QualityAssessor::get_inverted_element_count(int &inverted_elems,
                                                 int &inverted_samples,
                                                 MsqError &err)
{
  if(invertedElementCount == -1){
    MSQ_SETERR(err)("Number of inverted elements has not yet been calculated.", MsqError::INVALID_STATE);
    return false;
  }
  inverted_elems = invertedElementCount;
  inverted_samples = invertedSampleCount;
  return true;
}



void QualityAssessor::add_quality_assessment( QualityMetric* metric,
                                              int histogram_intervals,
                                              double power_mean,
                                              const char* tag_name )
{ 
  list_type::iterator i = find_or_add( metric );
  i->pMean = power_mean;
  if (histogram_intervals > 0)
    i->histogram.resize( histogram_intervals+2 );
  else
    i->histogram.clear();
  i->haveHistRange = false;
  if (!tag_name)
    i->tagName.clear();
  else
    i->tagName = tag_name;
}

QualityAssessor::list_type::iterator QualityAssessor::find_or_add( QualityMetric* qm )
{
  list_type::iterator iter;
  
    // If metric is already in list, find it
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    if (iter->qualMetric == qm )
      break;
  
    // If metric not found in list, add it
  if (iter == assessList.end())
  {
    if (qm->get_metric_type() == QualityMetric::VERTEX_BASED)
    {
      assessList.push_back( Assessor(qm) );
      iter = --assessList.end();
    }
    else
    {
      assessList.push_front( Assessor(qm) );
      iter = assessList.begin();
    }
  }
  
  return iter;
}

QualityAssessor::list_type::iterator QualityAssessor::find_stopping_assessment()
{
  list_type::iterator iter;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    if (iter->stopping_function())
      break;
  return iter;
}


/*!Sets which QualityMetric and QAFunction
combination is used to determine the value return from assess_mesh_quality().
It first ensures that the inputed QAFunction was not HISTOGRAM.  It then
calls add_quality_assessment with the given QualityMetric and QAFunction,
to ensure that this combination will be computed.  Finally, it sets
the stoppingMetric pointer and the stoppingFunction data members.
\param qm Pointer to QualityMetric.     
\param func (QAFUNCTION) Wrapper function for qm (e.g. MINIMUM, MAXIMUM,...).
    */
void QualityAssessor::set_stopping_assessment( QualityMetric* metric,
                                               int histogram_intervals,
                                               double power_mean,
                                               const char* tag_name )
{
  list_type::iterator i = find_stopping_assessment();
  if (i != assessList.end())
    i->set_stopping_function(false);

  i = find_or_add( metric );
  i->pMean = power_mean;
  if (histogram_intervals > 0)
    i->histogram.resize( histogram_intervals+2 );
  else
    i->histogram.clear();
  i->haveHistRange = false;
  if (!tag_name)
    i->tagName.clear();
  else
    i->tagName = tag_name;
  i->set_stopping_function(true);
}


/*! 
Checks first to see if the QualityMetric, qm, has been added to this
QualityAssessor, and if it has not, adds it.  It then adds HISTOGRAM as a
QAFunciton for that metric.  It then sets the minimum and maximum values
for the histogram.
\param qm Pointer to the QualityMetric to be used in histogram.
\param min_val (double) Minimum range of histogram.
\param max_val (double) Maximum range of histogram.
\param intervals Number of histogram intervals
    */
void QualityAssessor::add_histogram_assessment( QualityMetric* metric,
                                                double min_val, 
                                                double max_val,
                                                int intervals,
                                                double power_mean,
                                                const char* tag_name )
{
  if (intervals < 1)
    intervals = 1;
  list_type::iterator i = find_or_add( metric );
  i->pMean = power_mean;
  i->histMin = min_val;
  i->histMax = max_val;
  i->haveHistRange = min_val < max_val;
  i->histogram.resize( intervals + 2 );
  if (!tag_name)
    i->tagName.clear();
  else
    i->tagName = tag_name;
} 


TagHandle QualityAssessor::get_tag( Mesh* mesh,
                                    std::string name, 
                                    Mesh::TagType type, 
                                    unsigned size, 
                                    MsqError& err )
{
  TagHandle tag  = mesh->tag_get( name, err );
  if (!err) {
    Mesh::TagType exist_type;
    msq_std::string junk;
    unsigned exist_size;
    mesh->tag_properties( tag, junk, exist_type, exist_size, err ); MSQ_ERRZERO(err);
    if (type != exist_type || size != exist_size) {
       MSQ_SETERR(err)( MsqError:: TAG_ALREADY_EXISTS,
       "Tag \"%s\" exists with incorrect type or length.",
       name.c_str());
    }
  }
  else if (err.error_code() == MsqError::TAG_NOT_FOUND) {
    err.clear();
    tag = mesh->tag_create( name, type, size, 0, err ); MSQ_ERRZERO(err);
  }
  else {
    MSQ_ERRZERO(err);
  }
  
  return tag;
}


/*! 
  Computes the quality data for a given
  MeshSet, ms. What quality information is calculated, depends
  on what has been requested through the use of the QualityAssessor
  constructor, add_quality_assessment(), and set_stopping_assessment().
  The resulting data is printed in a table unless disable_printing_results()
  has been called.  The double returned depends on the QualityMetric
  and QAFunction "return" combination, which can be set using
  set_stopping_assessemnt().
  \param ms (const MeshSet &) MeshSet used for quality assessment.
 */
double QualityAssessor::loop_over_mesh( Mesh* mesh,
                                        MeshDomain* domain,
                                        const Settings* settings,
                                        MsqError& err)
{
  return loop_over_mesh_internal( mesh, domain, settings, NULL, err );
}

/*! 
  Computes the quality data for a given
  MeshSet, ms. What quality information is calculated, depends
  on what has been requested through the use of the QualityAssessor
  constructor, add_quality_assessment(), and set_stopping_assessment().
  The resulting data is printed in a table unless disable_printing_results()
  has been called.  The double returned depends on the QualityMetric
  and QAFunction "return" combination, which can be set using
  set_stopping_assessemnt().
  \param ms (const MeshSet &) MeshSet used for quality assessment.
 */
double QualityAssessor::loop_over_mesh( ParallelMesh* mesh,
                                        MeshDomain* domain,
                                        const Settings* settings,
                                        MsqError& err)
{
  return loop_over_mesh_internal( mesh, domain, settings, 
                                  mesh->get_parallel_helper(), err );
}

double QualityAssessor::loop_over_mesh_internal( Mesh* mesh,
                                                 MeshDomain* domain,
                                                 const Settings* settings,
                                                 ParallelHelper* helper,
                                                 MsqError& err )
{
    // Clear out any previous data
  reset_data();

  PatchData patch;
  patch.set_mesh( mesh );
  patch.set_domain( domain );
  if (settings)
    patch.attach_settings( settings );
  
  ElementPatches elem_patches;
  elem_patches.set_mesh( mesh );
  VertexPatches vert_patches(1,false);
  vert_patches.set_mesh( mesh );
  
  msq_std::vector<PatchSet::PatchHandle> patches;
  msq_std::vector<PatchSet::PatchHandle>::iterator p;
  msq_std::vector<Mesh::VertexHandle> patch_verts;
  msq_std::vector<Mesh::ElementHandle> patch_elems;
  msq_std::vector<size_t> metric_handles;

    // Check if we really need the helper
  if (helper && helper->get_nprocs() == 1) helper = 0;

    // Get any necessary tag handles.
  TagHandle invertedTag = 0;
  if (tagging_inverted_elements()) {
    invertedTag = get_tag( mesh, invertedTagName, Mesh::INT, 1, err );
    MSQ_ERRZERO(err);
  }
  list_type::iterator iter;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter) {
    if (iter->write_to_tag()) {
      iter->tagHandle = get_tag( mesh, iter->tagName, Mesh::DOUBLE, 1, err );
      MSQ_ERRZERO(err);
    }
  }
  
  TagHandle fixedTag = 0;
  if (tagging_fixed_elements()) {
    fixedTag = get_tag( mesh, fixedTagName, Mesh::INT, 1, err );
    MSQ_ERRZERO(err);
  }
  
    // Check for any metrics for which a histogram is to be 
    // calculated and for which the user has not specified 
    // minimum and maximum values.  
    // Element-based metrics are first in list, followed
    // by vertex-based metrics.  Find first vertex-based
    // metric also such that element metrics go from
    // assessList.begin() to elem_end and vertex metrics
    // go from elem_end to assessList.end()
  list_type::iterator elem_end = assessList.end();
  bool need_second_pass_for_elements = false;
  bool need_second_pass_for_vertices = false;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
  {
    if (iter->get_metric()->get_metric_type() == QualityMetric::VERTEX_BASED)
      break;

    if (iter->have_histogram() && !iter->haveHistRange)
      need_second_pass_for_elements = true;
  }
  elem_end = iter;
  for ( ; iter != assessList.end(); ++iter)
  {
    if (iter->have_histogram() && !iter->haveHistRange)
      need_second_pass_for_vertices = true;
  }
  
    // Do element-based metrics
  elem_patches.get_patch_handles( patches, err ); MSQ_ERRZERO(err);

  invertedElementCount = invertedSampleCount = sampleCount = 0;
  elementCount = patches.size();
  freeElementCount = 0;
  
  int inverted, samples;
  bool first_pass = false;
  do { // might need to loop twice to calculate histograms
    first_pass = !first_pass;

      //until there are no more patches
      //there is another get_next_patch at
      //the end of this loop
    for (p = patches.begin(); p != patches.end(); ++p) {
      elem_patches.get_patch( *p, patch_elems, patch_verts, err ); MSQ_ERRZERO(err);
 
      if (helper) {
        bool ours = helper->is_our_element( patch_elems[0], err );
        MSQ_ERRZERO(err);
        if (!ours) {
          --elementCount;
          continue;
        }
      }
      patch.set_mesh_entities( patch_elems, patch_verts, err ); MSQ_ERRZERO(err);
      
      if (0 == patch.num_free_vertices()) {
        if (tagging_fixed_elements()) {
          Mesh::ElementHandle h = patch.get_element_handles_array()[0];
          int val = 1;
          mesh->tag_set_element_data( fixedTag, 1, &h, &val, err );
          MSQ_ERRZERO(err);
        }
        if (skipFixedSamples)
          continue;
      }
      
       //first check for inverted elements
      if (first_pass){
        ++freeElementCount;
        inverted = samples = 0;
        patch.element_by_index(0).check_element_orientation(patch, inverted, samples, err);
        MSQ_ERRZERO(err);
        invertedElementCount += (inverted != 0);
        invertedSampleCount += inverted;
        sampleCount += samples;

        if (tagging_inverted_elements()) {
          Mesh::ElementHandle h = patch.get_element_handles_array()[0];
          int val = (inverted != 0);
          mesh->tag_set_element_data( invertedTag, 1, &h, &val, err );
          MSQ_ERRZERO(err);
        }
      }

      // now process all element-based metrics
      for (iter = assessList.begin(); iter != elem_end; ++iter)
      {
          // If first pass, get values for all metrics
        if (first_pass)
        {
          double value = 0.0;
          metric_handles.clear();
          QualityMetric* qm = iter->get_metric();
          qm->get_single_pass( patch, metric_handles, skipFixedSamples, err ); MSQ_ERRZERO(err);
          for (msq_std::vector<size_t>::iterator j = metric_handles.begin(); 
               j != metric_handles.end(); ++j) 
          {
            bool valid = iter->get_metric()->evaluate( patch, *j, value, err ); MSQ_ERRZERO(err);
            iter->add_value(value);
            if (!valid) 
              iter->add_invalid_value();
          }
            // we don't do tag stuff unless metric is truely element-based
            // (only one value per element)
          if (iter->write_to_tag() && metric_handles.size() == 1) {
            Mesh::ElementHandle h = patch.get_element_handles_array()[0];
            mesh->tag_set_element_data( iter->tagHandle, 1, &h, &value, err ); MSQ_ERRZERO(err);
          }
        }
          // If second pass, only do metrics for which the
          // histogram hasn't been calculated yet.
        else if (iter->have_histogram() && !iter->haveHistRange)
        {
          metric_handles.clear();
          QualityMetric* qm = iter->get_metric();
          qm->get_evaluations( patch, metric_handles, skipFixedSamples, err ); MSQ_ERRZERO(err);
          for (msq_std::vector<size_t>::iterator j = metric_handles.begin(); 
               j != metric_handles.end(); ++j) 
          {
            double value;
            iter->get_metric()->evaluate( patch, *j, value, err ); MSQ_ERRZERO(err);
            iter->add_hist_value(value);
          }
        }
      }
    }
    if (MSQ_CHKERR(err)) return 0.0;

      // Fix up any histogram ranges which were calculated
    for (iter = assessList.begin(); iter != elem_end; ++iter)
      if (iter->have_histogram() && !iter->haveHistRange)
        if (first_pass) {
          if (helper) {
            helper->communicate_min_max_to_all(&(iter->minimum), &(iter->maximum), err); 
            MSQ_ERRZERO(err);
          }
        
          iter->calculate_histogram_range();
// Uncomment the following to have the QA keep the first
// calculated histogram range for all subsequent iterations.
//          else
//            iter->haveHistRange = true;
        }
  } while (first_pass && need_second_pass_for_elements);
      
    
      // Do vertex-based metrics
  if (assessList.end() != elem_end)
  {
    vert_patches.get_patch_handles( patches, err ); MSQ_ERRZERO(err);
    
    bool first_pass = false;
    do { // might need to loop twice to calculate histograms
      first_pass = !first_pass;
     
        //until there are no more patches
        //there is another get_next_patch at
        //the end of this loop
      for (p = patches.begin(); p != patches.end(); ++p) {
        vert_patches.get_patch( *p, patch_elems, patch_verts, err ); MSQ_ERRZERO(err);
        patch.set_mesh_entities( patch_elems, patch_verts, err ); MSQ_ERRZERO(err);
        if (skipFixedSamples && 0 == patch.num_free_vertices())
          continue;
        
        Mesh::VertexHandle vert_handle = reinterpret_cast<Mesh::VertexHandle>(*p);

        for (iter = elem_end; iter != assessList.end(); ++iter)
        {
            // If first pass, get values for all metrics
          if (first_pass)
          {
            double value = 0.0;
            metric_handles.clear();
            QualityMetric* qm = iter->get_metric();
            qm->get_single_pass( patch, metric_handles, skipFixedSamples, err ); MSQ_ERRZERO(err);
            for (msq_std::vector<size_t>::iterator j = metric_handles.begin(); 
                 j != metric_handles.end(); ++j) 
            {
              bool valid = iter->get_metric()->evaluate( patch, *j, value, err ); MSQ_ERRZERO(err);
              iter->add_value(value);
              if (!valid) 
                iter->add_invalid_value();
            }
              // we don't do tag stuff unless metric is truely vertex-based
              // (only one value per vertex)
            if (iter->write_to_tag() && metric_handles.size() == 1) {
              mesh->tag_set_vertex_data( iter->tagHandle, 1, &vert_handle, &value, err ); MSQ_ERRZERO(err);
            }
          }
            // If second pass, only do metrics for which the
            // histogram hasn't been calculated yet.
          else if (iter->have_histogram() && !iter->haveHistRange)
          {
            metric_handles.clear();
            QualityMetric* qm = iter->get_metric();
            qm->get_evaluations( patch, metric_handles, skipFixedSamples, err ); MSQ_ERRZERO(err);
            for (msq_std::vector<size_t>::iterator j = metric_handles.begin(); 
                 j != metric_handles.end(); ++j) 
            {
              double value;
              iter->get_metric()->evaluate( patch, *j, value, err ); MSQ_ERRZERO(err);
              iter->add_hist_value(value);
            }
          }
        }
      }
      if (MSQ_CHKERR(err)) return 0.0;
  
        // Fix up any histogram ranges which were calculated
      for (iter = elem_end; iter != assessList.end(); ++iter)
        if (iter->have_histogram() && !iter->haveHistRange)
          if (first_pass) {
            if (helper) {
              helper->communicate_min_max_to_all(&(iter->minimum), &(iter->maximum), err); 
              MSQ_ERRZERO(err);
            }
            iter->calculate_histogram_range();
// Uncomment the following to have the QA keep the first
// calculated histogram range for all subsequent iterations.
//          else
//            iter->haveHistRange = true;
          }
    } while (first_pass && need_second_pass_for_vertices);
  }  
  
  if (helper) {
    for (iter = assessList.begin(); iter != assessList.end(); ++iter) {

      helper->communicate_min_max_to_zero(&(iter->minimum), &(iter->maximum), err);
      MSQ_ERRZERO(err);

      helper->communicate_sums_to_zero(&freeElementCount, &invertedElementCount, &elementCount, &invertedSampleCount, &sampleCount, &(iter->count), &(iter->numInvalid), &(iter->sum), &(iter->sqrSum), err);
      MSQ_ERRZERO(err);

      if (iter->have_power_mean()) {
        helper->communicate_power_sum_to_zero( &(iter->pMean), err );
        MSQ_ERRZERO(err);
      }

      if (iter->have_histogram()) {
        helper->communicate_histogram_to_zero(iter->histogram, err);
        MSQ_ERRZERO(err);
      }
    }
  }
  
    // Print results, if requested
  if (printSummary && (!helper || helper->get_rank() == 0))
    print_summary( this->outputStream );
  
  list_type::iterator i = find_stopping_assessment();
  return i == assessList.end() ? 0 : i->stopping_function_value();
}

bool QualityAssessor::invalid_elements( ) const
{
  bool result = false;
  list_type::const_iterator iter;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    if (iter->get_invalid_element_count())
      result = true;
  return result;
}

void QualityAssessor::reset_data() 
{
  list_type::iterator iter;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    iter->reset_data();
  invertedElementCount = -1;
  invertedSampleCount = -1;
}

QualityAssessor::Assessor::Assessor( QualityMetric* metric )
  : qualMetric(metric),
    pMean(0.0),
    haveHistRange(false),
    histMin(1.0),
    histMax(0.0),
    tagHandle(0),
    stoppingFunction(false)
{
  reset_data();
}

const QualityAssessor::Assessor* QualityAssessor::get_results( QualityMetric* metric ) const
{
  list_type::const_iterator iter;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    if (iter->get_metric() == metric)
      return &*iter;
  return 0;
}


void QualityAssessor::Assessor:: get_histogram( double& lower_bound_out,
                                                double& upper_bound_out,
                                                msq_std::vector<int>& counts_out,
                                                MsqError& err ) const 
{
  if ( !have_histogram() )
  {
    MSQ_SETERR(err)("No histogram calculated.", MsqError::INVALID_STATE);
    return;
  }

  lower_bound_out = histMin;
  upper_bound_out = histMax;
  counts_out = histogram;
}

void QualityAssessor::Assessor::reset_data()
{
  count = 0;
  sum = 0;
  maximum = -HUGE_VAL;
  minimum = HUGE_VAL;
  sqrSum = 0;
  pSum = 0;
  numInvalid = 0;
    // zero histogram data
  size_t hist_size = histogram.size();
  histogram.clear();
  histogram.resize( hist_size, 0 );
}

void QualityAssessor::Assessor::add_value( double metric_value )
{
  sum += metric_value;
  sqrSum += metric_value*metric_value;
  if (metric_value > maximum)
    maximum = metric_value;
  if (metric_value < minimum)
    minimum = metric_value;
    // Only add value to histogram data from this function if
    // the user has specified the range.  If user has not 
    // specified the range, QualityAssessor will call add_hist_value()
    // directly once the range has been calculated.
  if (have_histogram() && haveHistRange)
    add_hist_value( metric_value );
  
  if (have_power_mean())
    pSum += pow( metric_value, pMean );
  
  ++count;
}

void QualityAssessor::Assessor::add_invalid_value()
{
  ++numInvalid;
}

void QualityAssessor::Assessor::add_hist_value( double metric_value )
{
    // First and last values in array are counts of values
    // outside the user-specified range of the histogram
    // (below and above, respectively.)
  if (metric_value < histMin)
    ++histogram[0];
  else if (metric_value > histMax)
    ++histogram[histogram.size()-1];
  else
  {
      // Calculate which interval the value is in.  Add one
      // because first entry is for values below user-specifed
      // minimum value for histogram.
    double fract = (metric_value - histMin) / (histMax - histMin);
    unsigned cell;
    if (fabs(fract - 1.0) < histMax*DBL_EPSILON)
      cell = histogram.size() - 1;
    else
      cell = 1 + (int)((fract * (histogram.size()-2)));

      // Add value to interval.
    ++histogram[cell];
  }
}

void QualityAssessor::Assessor::calculate_histogram_range()
{
  double step = (maximum - minimum) / (histogram.size() - 2);
  if (step == 0)
    step = 1.0;
  double size = pow( 10.0, ceil(log10(step)) );
  if (size < 1e-6) 
    size = 1.0;
  histMin = size * floor( minimum / size );
  histMax = size *  ceil( maximum / size );
}  

void QualityAssessor::print_summary( msq_stdio::ostream& stream ) const
{
  const int NAMEW = 19;  // Width of name column in table output
  const int NUMW = 12;   // Width of value columns in table output
  
    // Print title
  stream << msq_stdio::endl 
         << "************** " 
         << qualityAssessorName
         << " Summary **************"
         << msq_stdio::endl
         << msq_stdio::endl;
         
  if (freeElementCount != elementCount)
    stream << "  Evaluating quality for " << freeElementCount 
           << " free elements of " << elementCount 
           << " total elements." << msq_stdio::endl;
  else
    stream << "  Evaluating quality for " << elementCount 
           << " elements." << msq_stdio::endl;
  
  if (invertedElementCount) {
    stream << "  THERE ARE "
           << invertedElementCount
           << " INVERTED ELEMENTS. "
           << msq_stdio::endl
           << "  (Elements invalid at "
           << invertedSampleCount
           << " of " << sampleCount
           << " sample locations.)"
           << msq_stdio::endl
           << msq_stdio::endl;
  }
  else {
    stream << "  There were no inverted elements detected. "
           << msq_stdio::endl;
  }
  
    // Check if there are invalid values for any metric
  list_type::const_iterator iter;
  bool some_invalid = false;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter) {
    if (iter->get_invalid_element_count()) {
      some_invalid = true;
      stream << "  " << iter->get_invalid_element_count()
             << " OF " << iter->get_count()
             << " ENTITIES EVALUATED TO AN UNDEFINED VALUE FOR " 
             << iter->get_metric()->get_name()
             << msq_stdio::endl << msq_stdio::endl;
    }
  }
  if (!some_invalid) {
    stream << "  No entities had undefined values for any computed metric." 
           << msq_stdio::endl << msq_stdio::endl;
  }
         
    // Check if a user-define power-mean was calculated for any of the metrics
  msq_std::set<double> pmeans;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    if (iter->have_power_mean())
      pmeans.insert( iter->get_power() );
      
    // If power-mean of 1 or 2 was requested, change titles rather
    // than printing redundant values
  std::string average_str("average"), rms_str("rms");
  if (pmeans.find(1.0) != pmeans.end()) {
    pmeans.erase( pmeans.find(1.0) );
    average_str = "1-mean";
  }
  if (pmeans.find(2.0) != pmeans.end()) {
    pmeans.erase( pmeans.find(2.0) );
    rms_str = "2-mean";
  }
  
    // Number of values in table
  unsigned num_values = pmeans.size() + 5;
  
    // Decide how wide of a table field shoudl be used for the metric name
  unsigned twidth = get_terminal_width();
  unsigned maxnwidth = NAMEW;
  if (twidth) {
    unsigned valwidth = NUMW*num_values;
    maxnwidth = valwidth < twidth ? twidth - valwidth : 0;
  }
  unsigned namewidth = 0;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    if (iter->get_metric()->get_name().size() > namewidth)
      namewidth = iter->get_metric()->get_name().size();
  if (namewidth > maxnwidth)
    namewidth = maxnwidth;
  if (namewidth < 7)  // at least enough width for the column header
    namewidth = 7;
    
    // print comlumn label line
  std::set<double>::const_iterator piter;
  stream << msq_stdio::setw(namewidth) << "metric";
  stream << msq_stdio::setw(NUMW)      << "minimum";
  for (piter = pmeans.begin(); piter != pmeans.end() && *piter < 1.0; ++piter)
    stream << msq_stdio::setw(NUMW-6) << *piter << "-mean ";
  stream << msq_stdio::setw(NUMW)      << average_str;
  for (; piter != pmeans.end() && *piter < 2.0; ++piter)
    stream << msq_stdio::setw(NUMW-6) << *piter << "-mean ";
  stream << msq_stdio::setw(NUMW)      << rms_str;
  for (; piter != pmeans.end(); ++piter)
    stream << msq_stdio::setw(NUMW-6) << *piter << "-mean ";
  stream << msq_stdio::setw(NUMW)      << "maximum";
  stream << msq_stdio::setw(NUMW)      << "std.dev.";
  stream << msq_stdio::endl;
  
    // print metric values
  for (iter = assessList.begin(); iter != assessList.end(); ++iter) {
      // print name
    stream << msq_stdio::setw(namewidth) << iter->get_metric()->get_name();
    if (iter->get_metric()->get_name().size() > namewidth) 
      stream << std::endl << msq_stdio::setw(namewidth) << " ";
      // print minimum
    stream << msq_stdio::setw(NUMW) << iter->get_minimum();
      // print power-means with P less than 1.0
    for (piter = pmeans.begin(); piter != pmeans.end() && *piter < 1.0; ++piter) {
      if (iter->have_power_mean() && iter->get_power() == *piter) 
        stream << msq_stdio::setw(NUMW) << iter->get_power_mean();
      else
        stream << msq_stdio::setw(NUMW) << " ";
    }
      // print average
    stream << msq_stdio::setw(NUMW) << iter->get_average();
      // print power-means with P less than 2.0
    for ( ; piter != pmeans.end() && *piter < 2.0; ++piter) {
      if (iter->have_power_mean() && iter->get_power() == *piter) 
        stream << msq_stdio::setw(NUMW) << iter->get_power_mean();
      else
        stream << msq_stdio::setw(NUMW) << " ";
    }
      // print RMS
    stream << msq_stdio::setw(NUMW) << iter->get_rms();
      // print power-means with P greater than 2.0
    for ( ; piter != pmeans.end(); ++piter) {
      if (iter->have_power_mean() && iter->get_power() == *piter) 
        stream << msq_stdio::setw(NUMW) << iter->get_power_mean();
      else
        stream << msq_stdio::setw(NUMW) << " ";
    }
      // print maximum and standard deviation
    stream << msq_stdio::setw(NUMW) << iter->get_maximum();
    stream << msq_stdio::setw(NUMW) << iter->get_stddev();
    stream << msq_stdio::endl;
  }
  
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    if (iter->have_histogram())
      iter->print_histogram( stream, get_terminal_width() );
}


void QualityAssessor::Assessor::print_histogram( msq_stdio::ostream& stream,
                                                 int termwidth ) const
{
  // Portability notes:
  //  Use log10 rather than log10f because the float variations require
  //  including platform-dependent headers on some platforms.  
  //  Explicitly cast log10 argument to double because some platforms
  //  have overloaded float and double variations in C++ making an 
  //  implicit cast from an integer ambiguous.
  
  const char indent[] = "   ";
  const char GRAPH_CHAR = '=';  // Character used to create bar graphs
  const int FLOATW = 12;        // Width of floating-point output
  const int TOTAL_WIDTH = termwidth > 30 ? termwidth : 70;   // Width of histogram
  int GRAPHW = TOTAL_WIDTH - FLOATW - 4 - sizeof(indent);
  
    // range is either user-specified (histMin & histMax) or
    // calculated (minimum & maximum)
  double min, max;
  min = histMin;
  max = histMax;
    // Witdh of one interval of histogram
  double step = (max - min) / (histogram.size()-2);
  
    // Find maximum value for an interval of the histogram
  unsigned i;
  int max_interval = 1;
  for (i = 0; i < histogram.size(); ++i)
    if (histogram[i] > max_interval)
      max_interval = histogram[i];
  
  if (0 == max_interval)
    return; // no data 
  
    // Calculate width of field containing counts for 
    // histogram intervals (log10(max_interval)).
  int num_width = 1;
  for (int temp = max_interval; temp > 0; temp /= 10)
    ++num_width;
  GRAPHW -= num_width;

    // Create an array of bar graph characters for use in output
  msq_std::vector<char> graph_chars(GRAPHW+1, GRAPH_CHAR);
  
    // Check if bar-graph should be linear or log10 plot
    // Do log plot if standard deviation is less that 1.5
    // histogram intervals.
  bool log_plot = false;
  double stddev = get_stddev();
  if (stddev > 0 && stddev < 2.0*step)
  {
    int new_interval = (int)(log10((double)(1+max_interval)));
    if (new_interval > 0) {
      log_plot = true;
      max_interval = new_interval;
    }
  }

  
    // Write title
  stream << msq_stdio::endl << indent << get_metric()->get_name() << " histogram:";
  if (log_plot)
    stream << " (log10 plot)";
  stream << msq_stdio::endl;

  
    // For each interval of histogram
  for (i = 0; i < histogram.size(); ++i)
  {
      // First value is the count of the number of values that
      // were below the minimum value of the histogram.
    if (0 == i)
    {
      if (0 == histogram[i])
        continue;
      stream << indent << msq_stdio::setw(FLOATW) << "under min";
    }
      // Last value is the count of the number of values that
      // were above the maximum value of the histogram.
    else if (i+1 == histogram.size())
    {
      if (0 == histogram[i])
        continue;
      stream << indent << msq_stdio::setw(FLOATW) << "over max";
    }
      // Anything else is a valid interval of the histogram.
      // Print the lower bound for each interval.
    else
    {
      stream << indent << msq_stdio::setw(FLOATW) << min + (i-1)*step;
    }
    
      // Print interval count.
    stream << ": " << msq_stdio::setw(num_width) << histogram[i] << ": ";
    
      // Print bar graph
    
      // First calculate the number of characters to output
    int num_graph;
    if (log_plot)
      num_graph = GRAPHW * (int)log10((double)(1+histogram[i])) / max_interval;
    else
      num_graph = GRAPHW * histogram[i] / max_interval;
      
      // print num_graph characters using array of fill characters.
    graph_chars[num_graph] = '\0';
    stream << &graph_chars[0] << msq_stdio::endl;
    graph_chars[num_graph] = GRAPH_CHAR;
  }
  
  stream << msq_stdio::endl;
}

#ifdef _MSC_VER
# define fileno(A) _fileno( (A) )
#endif 
int QualityAssessor::get_terminal_width() const
{
#ifdef TIOCGWINSZ
  int fd = -1;
  if (&outputStream == &std::cout)
    fd = fileno(stdout);
  else if (&outputStream == &std::cerr)
    fd = fileno(stderr);
#ifdef FSTREAM_HAS_FD
  else if (std::ofstream* f = dynamic_cast<std::ofstream*>(&outputStream))
    fd = f->rdbuf()->fd();
#endif

  if (fd >= 0) {
    struct winsize ws;
    if (ioctl(fd, TIOCGWINSZ, &ws) >= 0)
      return ws.ws_col;
  }
#endif

  return 0;
}

double QualityAssessor::Assessor::stopping_function_value() const
{
  return have_power_mean() ? get_power_mean() : get_average();
}

  

} //namespace Mesquite
