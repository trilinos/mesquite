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
#include "PatchData.hpp"
#include "MsqMeshEntity.hpp"
#include "MsqVertex.hpp"
#include "MsqDebug.hpp"
#include "MeshSet.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <list.h>
#  include <vector.h>
#else
#  include <list>
#  include <vector>
   using namespace std;
#endif

#ifdef MSQ_USE_OLD_IO_HEADERS
#  include <iostream.h>
#  include <iomanip.h>
#else
#  include <iostream>
#  include <iomanip>
   using namespace std;
#endif

namespace Mesquite {

QualityAssessor::QualityAssessor(string name) :
  qualityAssessorName(name),
  outputStream( cout ),
  printSummary( true ),
  stoppingMetric( assessList.end() ),
  stoppingFunction( (QAFunction)0 )
{ 
  MsqError err;
  set_patch_type( PatchData::GLOBAL_PATCH, err, 0 );
}

QualityAssessor::QualityAssessor(ostream& stream, string name) :
  qualityAssessorName(name),
  outputStream( stream ),
  printSummary( true ),
  stoppingMetric( assessList.end() ),
  stoppingFunction( (QAFunction)0 )
{ 
  MsqError err;
  set_patch_type( PatchData::GLOBAL_PATCH, err, 0 );
}
 
QualityAssessor::QualityAssessor( QualityMetric* metric,
                                  QAFunction function,
                                  MsqError& err,
                                  string name ) :
  qualityAssessorName(name),
  outputStream( cout ),
  printSummary( true ),
  stoppingMetric( assessList.end() ),
  stoppingFunction( (QAFunction)0 )
{ 
  set_patch_type( PatchData::GLOBAL_PATCH, err, 0 );
  set_stopping_assessment( metric, function, err );
}

QualityAssessor::QualityAssessor( QualityMetric* metric,
                                  QAFunction function,
                                  ostream& stream, 
                                  MsqError& err,
                                  string name ) :
  qualityAssessorName(name),
  outputStream( stream ),
  printSummary( true ),
  stoppingMetric( assessList.end() ),
  stoppingFunction( (QAFunction)0 )
{ 
  set_patch_type( PatchData::GLOBAL_PATCH, err, 0 );
  set_stopping_assessment( metric, function, err );
}

QualityAssessor::~QualityAssessor()
  { }

string QualityAssessor::get_QAFunction_name(
                              enum QualityAssessor::QAFunction fun)
{
  switch(fun){
    case(AVERAGE):
      return "Average   ";
    case(HISTOGRAM):
      return "Histogram of metric values: ";
    case(MAXIMUM):
      return "Maximum   ";
    case(MINIMUM):
      return "Minimum   ";
    case(RMS):
      return "RMS       ";
    case(STDDEV):
      return "Stan. Dev.";
    default:
      return "DEFAULT   ";
  };
}


/*!
    Several QualityMetric objects can be added to a single QualityAssessor
    object.  This allows to perform several quality assessments over a
    single mesh sweep.
    \param qm is the QualityMetric that will be used to evaluate the mesh
    quality
    \param func is the wrapper function used over the QualityMetric
    (min, max, etc..)
 */
void QualityAssessor::add_quality_assessment(QualityMetric* qm,
                                             int func,
                                             MsqError &/*err*/)
{ 
  msq_std::list<Assessor>::iterator iter;
  
  iter = find_or_add( qm );
  iter->funcFlags |= func;
  if (func&HISTOGRAM)
    iter->histogram.resize(12);
}

list<QualityAssessor::Assessor>::iterator QualityAssessor::find_or_add( QualityMetric* qm )
{
  list<Assessor>::iterator iter;
  
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


/*!Sets which QualityMetric and QAFunction
combination is used to determine the value return from assess_mesh_quality().
It first ensures that the inputed QAFunction was not HISTOGRAM.  It then
calls add_quality_assessment with the given QualityMetric and QAFunction,
to ensure that this combination will be computed.  Finally, it sets
the stoppingMetric pointer and the stoppingFunction data members.
\param qm Pointer to QualityMetric.     
\param func (QAFUNCTION) Wrapper function for qm (e.g. MINIMUM, MAXIMUM,...).
    */
void QualityAssessor::set_stopping_assessment(QualityMetric* qm,
                                              QAFunction func,
                                              MsqError &err)
{
  if(func==HISTOGRAM){
    MSQ_SETERR(err)("HISTOGRAM DOES NOT GIVE A VALID RETURN VALUE", MsqError::INVALID_ARG);
    return;
  }
  
  stoppingMetric = find_or_add( qm );
  stoppingMetric->funcFlags |= func;
  stoppingFunction = func;
  if (func&HISTOGRAM)
    stoppingMetric->histogram.resize(12);
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
void QualityAssessor::add_histogram_assessment( QualityMetric* qm,
                                                double min_val, 
                                                double max_val,
                                                int intervals,
                                                MsqError &err )
{
  if (min_val >= max_val || intervals < 1) {
    MSQ_SETERR(err)("Invalid histogram range.", MsqError::INVALID_ARG );
    return;
  }
  
  list<Assessor>::iterator assessor = find_or_add( qm );
  assessor->funcFlags |= QualityAssessor::HISTOGRAM;
  assessor->histMin = min_val;
  assessor->histMax = max_val;
  assessor->histogram.resize( intervals + 2 );
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
double QualityAssessor::loop_over_mesh(MeshSet &ms, MsqError& err)
{
    // Clear out any previous data
  reset_data();
  
    // Check for any metrics for which a histogram is to be 
    // calculated and for which the user has not specified 
    // minimum and maximum values.  Enable calculation of
    // the min/max for these metrics in the first pass so
    // we have them for calculating the histograms in the
    // second pass.
    // Element-based metrics are first in list, followed
    // by vertex-based metrics.  Find first vertex-based
    // metric also such that element metrics go from
    // assessList.begin() to elem_end and vertex metrics
    // go from elem_end to assessList.end()
  list<Assessor>::iterator elem_end = assessList.end();
  bool need_second_pass_for_elements = false;
  bool need_second_pass_for_vertices = false;
  list<Assessor>::iterator iter;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
  {
    if (iter->get_metric()->get_metric_type() == QualityMetric::VERTEX_BASED)
      break;

    if (iter->funcFlags&HISTOGRAM && !iter->haveHistRange)
    {
      iter->funcFlags |= MINIMUM|MAXIMUM;
      need_second_pass_for_elements = true;
    }
  }
  elem_end = iter;
  for ( ; iter != assessList.end(); ++iter)
  {
    if (iter->funcFlags&HISTOGRAM && !iter->haveHistRange)
    {
      iter->funcFlags |= MINIMUM|MAXIMUM;
      need_second_pass_for_vertices = true;
    }
  }
  
    // Do element-based metrics
  if (assessList.begin() != elem_end)
  {
    PatchData* pd;
    PatchData local_patch;
    bool more_mesh;
    no_culling_method();
      
    bool first_pass = false;
    do { // might need to loop twice to calculate histograms
      first_pass = !first_pass;
     
      more_mesh = true;
      if (get_global_patch() == 0) {
        pd = &local_patch;
        more_mesh=ms.get_next_patch(*pd, this, err);  MSQ_ERRZERO(err);
      }
      else {
        pd = get_global_patch();
      }
      
        //until there are no more patches
        //there is another get_next_patch at
        //the end of this loop
      while (more_mesh)
      {
        for (unsigned i = 0; i < pd->num_elements(); ++i)
        {
          for (iter = assessList.begin(); iter != elem_end; ++iter)
          {
            if (first_pass)
            {
              double value;
              bool valid = iter->get_metric()->evaluate_element( *pd, 
                                                           &pd->element_by_index(i),
                                                           value, err );
                                                           MSQ_ERRZERO(err);
              
              if (valid)
                iter->add_value(value);
              else 
                iter->add_invalid_value();
            }
            else if (iter->funcFlags&HISTOGRAM && !iter->haveHistRange)
            {
              double value;
              bool valid = iter->get_metric()->evaluate_element( *pd, 
                                                           &pd->element_by_index(i),
                                                           value, err );
                                                           MSQ_ERRZERO(err);
              
              if (valid)
                iter->add_hist_value(value);
            }
          }
        }
        
           // If dealing with local patches, get next element group (PatchData object)
        if (get_patch_type() != PatchData::GLOBAL_PATCH)
          more_mesh = ms.get_next_patch(*pd,this, err); MSQ_ERRZERO(err);
          //Michael:: Since we are doing global right now:
          //Remove this when no longer doing global
        more_mesh=false;
      }
    
    } while (first_pass && need_second_pass_for_elements);
  }
      
    
      // Do vertex-based metrics
  if (assessList.end() != elem_end)
  {
    PatchData* pd;
    PatchData local_patch;
    no_culling_method();
    bool more_mesh;

    bool first_pass = false;
    do { // might need to loop twice to calculate histograms
      first_pass = !first_pass;
     
        //construct the patch we will send to get_next_patch
      more_mesh = true;
      if (get_global_patch() == 0) {
        pd = &local_patch; 
        more_mesh=ms.get_next_patch(*pd, this, err);  MSQ_ERRZERO(err);
      }
      else {
        pd = get_global_patch();
      }
      
        //until there are no more patches
        //there is another get_next_patch at
        //the end of this loop
      while (more_mesh)
      {
        for (unsigned i = 0; i < pd->num_vertices(); ++i)
        {
          for (iter = elem_end; iter != assessList.end(); ++iter)
          {
            if (first_pass)
            {
              double value;
              bool valid = iter->get_metric()->evaluate_vertex( *pd, 
                                                           &pd->vertex_by_index(i),
                                                           value, err );
                                                           MSQ_ERRZERO(err);
              
              if (valid)
                iter->add_value(value);
              else 
                iter->add_invalid_value();
            }
            else if (iter->funcFlags&HISTOGRAM && !iter->haveHistRange)
            {
              double value;
              bool valid = iter->get_metric()->evaluate_vertex( *pd, 
                                                           &pd->vertex_by_index(i),
                                                           value, err );
                                                           MSQ_ERRZERO(err);
              
              if (valid)
                iter->add_hist_value(value);
            }
          }
        }
        
        if (get_patch_type() != PatchData::GLOBAL_PATCH)
          more_mesh = ms.get_next_patch(*pd,this, err); MSQ_ERRZERO(err);
          //Michael:: Since we are doing global right now:
          //Remove this when no longer doing global
        more_mesh=false;
      }
    
    } while (first_pass && need_second_pass_for_vertices);
  }  
  
    // Print results, if requested
  if (printSummary)
    print_summary( this->outputStream );
  
    // If no stopping function, just return zero
  if (!stoppingFunction)
    return 0.0;
  
    // Otherwise return requested value
  if      (stoppingFunction & STDDEV)
    return stoppingMetric->get_stddev();
  else if (stoppingFunction & AVERAGE)
    return stoppingMetric->get_average();
  else if (stoppingFunction & MAXIMUM)
    return stoppingMetric->get_maximum();
  else if (stoppingFunction & MINIMUM)
    return stoppingMetric->get_minimum();
  else if (stoppingFunction & RMS)
    return stoppingMetric->get_rms();
  else 
    MSQ_SETERR(err)("Invalid stopping function for QualityAssessor", 
                    MsqError::INVALID_STATE);
  
  return 0.0;
}

bool QualityAssessor::invalid_elements( ) const
{
  bool result = false;
  msq_std::list<Assessor>::const_iterator iter;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    if (iter->get_invalid_element_count())
      result = true;
  return result;
}

void QualityAssessor::reset_data() 
{
  msq_std::list<Assessor>::iterator iter;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    iter->reset_data();
}

QualityAssessor::Assessor::Assessor( QualityMetric* metric )
  : qualMetric(metric),
    funcFlags(0),
    histMin(1.0),
    histMax(0.0),
    haveHistRange(false)
{
  reset_data();
}

const QualityAssessor::Assessor* QualityAssessor::get_results( QualityMetric* metric ) const
{
  msq_std::list<Assessor>::const_iterator iter;
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
  if ( !(funcFlags & QualityAssessor::HISTOGRAM) )
  {
    MSQ_SETERR(err)("No histogram calculated.", MsqError::INVALID_STATE);
    return;
  }

  if (haveHistRange) {
    lower_bound_out = histMin;
    upper_bound_out = histMax;
  }
  else {
    lower_bound_out = minimum;
    upper_bound_out = maximum;
  }
  
  counts_out = histogram;
}

void QualityAssessor::Assessor::reset_data()
{
  count = 0;
  sum = 0;
  maximum = -DBL_MAX;
  minimum = DBL_MAX;
  sqrSum = 0;
  numInvalid = 0;
  memset( &histogram[0], 0, sizeof(int)*histogram.size() );
}

void QualityAssessor::Assessor::add_value( double metric_value )
{
  if (funcFlags & (QualityAssessor::AVERAGE|QualityAssessor::STDDEV))
    sum += metric_value;
  if (funcFlags & (QualityAssessor::RMS|QualityAssessor::STDDEV))
    sqrSum += metric_value*metric_value;
  if (funcFlags & QualityAssessor::MAXIMUM && metric_value > maximum)
    maximum = metric_value;
  if (funcFlags & QualityAssessor::MINIMUM && metric_value < minimum)
    minimum = metric_value;
  if (funcFlags & QualityAssessor::HISTOGRAM && haveHistRange)
    add_hist_value( metric_value );
  
  ++count;
}

void QualityAssessor::Assessor::add_invalid_value()
{
  ++count;
  ++numInvalid;
}

void QualityAssessor::Assessor::add_hist_value( double metric_value )
{
  double min, max, step;
  if (haveHistRange) {
    min = histMin;
    max = histMax;
  }
  else {
    min = minimum;
    max = maximum;
  }
  
  step = (max - min) / (histogram.size()-2);
  if (metric_value < min)
    ++histogram[0];
  else if (metric_value > max)
    ++histogram[histogram.size()-1];
  else
  {
    unsigned cell = 1+(unsigned)((metric_value - min) / step);
    if (cell + 1 == histogram.size())
      --cell;
    ++histogram[cell];
  }
}

void QualityAssessor::print_summary( msq_stdio::ostream& stream ) const
{
  const int NAMEW = 19;
  const int NUMW = 12;
  
    // Print title
  stream << msq_stdio::endl 
         << "************** " 
         << qualityAssessorName
         << " Summary **************"
         << msq_stdio::endl
         << msq_stdio::endl;
         
    // Get union of function flags
  msq_std::list<Assessor>::const_iterator iter;
  unsigned flags = 0;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    if (iter->get_invalid_element_count())
      stream << ">>>" << msq_stdio::setw(NAMEW) << iter->get_metric()->get_name()
             << " IS INVALID FOR " << iter->get_invalid_element_count()
             << " OF " << iter->get_count() << " VALUES!" << msq_stdio::endl
             << msq_stdio::endl;
    else
      flags |= iter->funcFlags;
  
  if (flags) 
  {
    stream << msq_stdio::setw(NAMEW) << "metric";
    if (flags & AVERAGE)
      stream << msq_stdio::setw(NUMW) << "average";
    if (flags & MINIMUM)
      stream << msq_stdio::setw(NUMW) << "minimum";
    if (flags & MAXIMUM)
      stream << msq_stdio::setw(NUMW) << "maximum";
    if (flags & RMS)
      stream << msq_stdio::setw(NUMW) << "rms";
    if (flags & STDDEV)
      stream << msq_stdio::setw(NUMW) << "std.dev.";
    stream << msq_stdio::endl;

    for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    {
      if (iter->get_invalid_element_count())
        continue;

      stream << msq_stdio::setw(NAMEW) << iter->get_metric()->get_name();

      if (flags & AVERAGE)
      {
        if (iter->funcFlags & AVERAGE)
          stream << msq_stdio::setw(NUMW) << iter->get_average();
        else
          stream << msq_stdio::setw(NUMW) << " ";
      }
      if (flags & MINIMUM)
      {
        if (iter->funcFlags & MINIMUM)
          stream << msq_stdio::setw(NUMW) << iter->get_minimum();
        else
          stream << msq_stdio::setw(NUMW) << " ";
      }
      if (flags & MAXIMUM)
      {
        if (iter->funcFlags & MAXIMUM)
          stream << msq_stdio::setw(NUMW) << iter->get_maximum();
        else
          stream << msq_stdio::setw(NUMW) << " ";
      }
      if (flags & RMS)
      {
        if (iter->funcFlags & RMS)
          stream << msq_stdio::setw(NUMW) << iter->get_rms();
        else
          stream << msq_stdio::setw(NUMW) << " ";
      }
      if (flags & STDDEV)
      {
        if (iter->funcFlags & STDDEV)
          stream << msq_stdio::setw(NUMW) << iter->get_stddev();
        else
          stream << msq_stdio::setw(NUMW) << " ";
      }
      stream << msq_stdio::endl;
    } // for (assessList)
  } // if (flags)
  
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    if (iter->funcFlags & HISTOGRAM && !iter->get_invalid_element_count())
      iter->print_histogram( stream );
}


void QualityAssessor::Assessor::print_histogram( msq_stdio::ostream& stream ) const
{
  const char GRAPH_CHAR = '=';
  const int FLOATW = 12;
  const int GRAPHW = 50;
  
  stream << msq_stdio::endl << "   " << get_metric()->get_name() 
         << " histogram:" << msq_stdio::endl;
  
  double min, max, step;
  if (haveHistRange) {
    min = histMin;
    max = histMax;
  }
  else {
    min = minimum;
    max = maximum;
  }
  step = (max - min) / (histogram.size()-2);
  
  unsigned i;
  int max_interval = 1;
  for (i = 0; i < histogram.size(); ++i)
    if (histogram[i] > max_interval)
      max_interval = histogram[i];
  
  int num_width = 1;
  for (int temp = max_interval; temp > 0; temp /= 10)
    ++num_width;

  char graph_chars[GRAPHW+1];
  memset( graph_chars, GRAPH_CHAR, sizeof(graph_chars) );
  
  for (i = 0; i < histogram.size(); ++i)
  {
    if (0 == i)
    {
      if (0 == histogram[i])
        continue;
      stream << setw(FLOATW) << "under min";
    }
    else if (i+1 == histogram.size())
    {
      if (0 == histogram[i])
        continue;
      stream << setw(FLOATW) << "over max";
    }
    else
    {
      stream << "   " << setw(FLOATW) << min + (i-1)*step;
    }
    
    stream << ": " << setw(num_width) << histogram[i] << ": ";
    
    int num_graph = (GRAPHW * histogram[i]) / max_interval;
    
    graph_chars[num_graph] = '\0';
    stream << graph_chars << msq_stdio::endl;
    graph_chars[num_graph] = GRAPH_CHAR;
  }
  
  stream << msq_stdio::endl;
}
 


} //namespace Mesquite
