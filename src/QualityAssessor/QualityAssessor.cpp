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

#include <math.h>

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <list.h>
#  include <vector.h>
#else
#  include <list>
#  include <vector>
   using namespace std;
#endif

namespace Mesquite {

/*! 
  Constructor creates a QualityAssessor, and adds the given QualityMetric
  pointer and QAFunciton to the list of Assessors (a structure containing
  a QualityMetric pointer and the QAFunctions assossiated with it).  The
  constructor also sets the given QualityMetric and QAFunction as the
  "return value" the combination used to calculate the value returned
  from assess_mesh_quality.  If func == HISTOGRAM, then MAXIMUM
  is used as default QAFunction for purposes of the "return value".
  \param qm Pointer to QualityMetric.
  \param func (QAFUNCTION) Wrapper function for qm (e.g. MINIMUM, MAXIMUM,...).
*/
QualityAssessor::QualityAssessor(QualityMetric* qm, enum QAFunction func,
                                 string name) :
  qualityAssessorName(name)
{
  MsqError err;
  if(func==QualityAssessor::HISTOGRAM || func==QualityAssessor::ALL_MEASURES){
    add_quality_assessment(qm, func, err);  MSQ_ERRRTN(err);
    set_stopping_assessment(qm, QualityAssessor::MAXIMUM, err);  MSQ_ERRRTN(err);
  }
  else{
    set_stopping_assessment(qm,func, err);  MSQ_ERRRTN(err);
  }
  printingTurnedOff=0;
    //sigDig=5;
    //When we are no longer doing a GLOBAL patch, we          
    //must also remove the 'elem_bool=0' below.         
  set_patch_type(PatchData::GLOBAL_PATCH, err, 0); MSQ_ERRRTN(err);
}

QualityAssessor::~QualityAssessor()
{
  list<Assessor*>::iterator pos;
  pos=assessList.begin();
  Assessor* tmp=NULL;
  while ( !(pos==assessList.end())){
    tmp=(*pos);
    ++pos;
    delete tmp;
  }
}

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
                                             enum QAFunction  func,
                                             MsqError &/*err*/)
{ 
  int found=0;
  list<Assessor*>::iterator pos;
  pos=assessList.begin();
  Assessor* assess_ptr=NULL;
    //loop over the assessList (list of Assessor of this QA.)
  while ( !(pos==assessList.end()) && !found){
    assess_ptr=*pos;
      //if this metric (qm) is already in the list
      //then add func to the bit flag (funcFlagBits)
    if(assess_ptr->metric==qm){
      assess_ptr->funcFlagBits |= func;
      found=1;//assessment has been added
    }
    ++pos;
  }
  if(!found){//assessment has not been added (qm was not in list)
      //Michael:  new without associated deleted?
    assess_ptr=new Assessor;//Thus, we need a new Assessor for qm
    assess_ptr->metric=qm;
    assess_ptr->funcFlagBits |= func;
    assess_ptr->minHist = -(MSQ_MAX_CAP+1.0);
    assess_ptr->maxHist = MSQ_MAX_CAP+1.0;    
      //add vertex based metrics to back of list
      //add element based metrics to front of list
    if(assess_ptr->metric->get_metric_type()==QualityMetric::VERTEX_BASED){
      assessList.push_back( assess_ptr );
    }
    else{
      assessList.push_front( assess_ptr );
    }
  }
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
                                              enum QAFunction func,
                                              MsqError &err)
{
  if(func==HISTOGRAM){
    MSQ_SETERR(err)("HISTOGRAM DOES NOT GIVE A VALID RETURN VALUE", MsqError::INVALID_ARG);
  }
  add_quality_assessment(qm,func,err);  MSQ_ERRRTN(err);
  stoppingFunction=func;
  stoppingMetric=qm;
}

  


/*! 
Checks first to see if the QualityMetric, qm, has been added to this
QualityAssessor, and if it has not, adds it.  It then adds HISTOGRAM as a
QAFunciton for that metric.  It then sets the minimum and maximum values
for the histogram.
\param qm Pointer to the QualityMetric to be used in histogram.
\param min_val (double) Minimum range of histogram.
\param max_val (double) Maximum range of histogram.
    */
void QualityAssessor::set_histogram_range(QualityMetric* qm,
                                             double min_val, double max_val,
                                          MsqError &/*err*/)
{
  int found=0;
  list<Assessor*>::iterator pos;
  pos=assessList.begin();
  Assessor* assess_ptr=NULL;
    //loop over the assessList (list of Assessor of this QA.)
  size_t i=0;
  while(i<assessList.size() && !found){
    assess_ptr=*pos;
      //if this metric (qm) is already in the list
      //then add HISTOGRAM to the bit flag (funcFlagBits)
      //and set min and max hist values
    if(assess_ptr->metric==qm){
      found=1;//assessment has been added
    }
    ++pos;
    ++i;
  }
  if(!found){//assessment has not been added (qm was not in list)
    assess_ptr=new Assessor;//Thus, we need a new Assessor for qm
    assess_ptr->metric=qm;    
      //add vertex based metrics to back of list
      //add element based metrics to front of list
    if(assess_ptr->metric->get_metric_type()==QualityMetric::VERTEX_BASED){
      assessList.push_back( assess_ptr );
    }
    else{
      assessList.push_front( assess_ptr );
    }
  }
  assess_ptr->funcFlagBits |= QualityAssessor::HISTOGRAM;
  if(min_val<=max_val){
    assess_ptr->minHist = min_val;
    assess_ptr->maxHist = max_val;
  }
  else{
    assess_ptr->minHist = max_val;
    assess_ptr->maxHist = min_val;
  }
  
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
    //tot_num either equals total_num_elements or total_num_vertices
    //depending on whethere the metric is element or vertex_based.
  int tot_num=0;
    //swich telling us whether two loops are necessary.  Two loops
    //are currently necessary whenever the user specifies that
    //a histogram should be produced, but does not specify
    //a range for that histogram.
  double two_loops=1;
  double return_value=0;
  int total_num_elements=0;
  int total_num_vertices=0;
  list<Assessor*>::iterator pos = assessList.begin();
  int num_elem_based=0;
  int num_metrics=assessList.size();
  //Assessor** assessor_array = new Assessor*[num_metrics];
  vector<Assessor*> assessor_array(num_metrics);
    //array of data holders QAVars.  QAVars is a struct
    //defined in QualityAssessor.hpp which creates a
    //way to store the data we accumulate for the each metric.
    //Thus, we have one QAVars, per metric in our lis.
    //Each QAVar also has an bit flag int, funcFlagBits,
    //which tells us (after we initialize it) what we
    //need to calcuate... max, min, avg, rms, hist, std_dev.
  //QAVars *QAData = new QAVars[num_metrics];
  vector<QAVars> QAData(num_metrics);
  
  int i=0;
    //find num_elem_based (the number of elem based metrics)
    //elem_based metrics are at the beginning of the list
    //vertex based metrics are at the end of the list
    //We loop over the number of metrics (or more correctly
    //but equivalently the number of Assessors) for this
    //QualityAssessor (using its Assessor list).
    //
  while(i<num_metrics){
    assessor_array[i]=(*pos);
    if((*pos)->metric->get_metric_type()==QualityMetric::ELEMENT_BASED)
      num_elem_based=i+1;
      //if histogram must calc min and max
    if(assessor_array[i]->funcFlagBits&QualityAssessor::HISTOGRAM){
      assessor_array[i]->funcFlagBits|=MAXIMUM;
      assessor_array[i]->funcFlagBits|=MINIMUM;
    }
      //Standard Deviation must calc average
    if(assessor_array[i]->funcFlagBits&QualityAssessor::STDDEV){
      assessor_array[i]->funcFlagBits|=AVERAGE;
    }
      //intitialize average
    QAData[i].avgVar=0.0;
      //intitialize each slot in the histogram.
    int ii;
    for(ii=0;ii<MSQ_HIST_SIZE+2;ii++){
      QAData[i].histVar[ii]=0;
    }
      //initialize max, min, rms, stddev., valid
    QAData[i].maxVar=0.0;
    QAData[i].minVar=MSQ_MAX_CAP;
    QAData[i].rmsVar=0.0;
    QAData[i].stdVar=0.0;
    QAData[i].histMax=assessor_array[i]->maxHist;
    QAData[i].histMin=assessor_array[i]->minHist;
    QAData[i].histDelta=(QAData[i].histMax-QAData[i].histMin)/MSQ_HIST_SIZE;
    QAData[i].numInvalid=0;
    ++i;
    ++pos;
  }
    //pointers to elems
  MsqMeshEntity* elems;
  MsqVertex* verts;
  int num_elems=0;
  int num_verts=0;
    //michael temp solution
  int num_pass=0;
    //We loop over the mesh twice here if the user requests
    //a histogram and not a range for that histogram.  In
    //This case we loop over once calculating everything
    //requeseted except the histogram information.  We
    //then use the min and max values to form the range
    //of the histogram.  We then accumulate the histogram
    //data on the second pass over the mesh.
    //
  int metric_counter=0;
  while(num_pass<2 && two_loops){
      //two_loops will be set to one if two loops are necessary
    metric_counter=0;
    two_loops=0;
    double temp_val=0;
    if(num_elem_based){
        //construct the patch we will send to get_next_patch
      PatchData* elem_group;
      PatchData local_elem_group;
        
      no_culling_method();
      
      bool elem_bool=true;
      if (get_global_patch() == 0) {
        elem_group = &local_elem_group;
        elem_bool=ms.get_next_patch(*elem_group, this, err);  MSQ_ERRZERO(err);
      }
      else {
        elem_group = get_global_patch();
      }
      
        //until there are no more patches
        //there is another get_next_patch at
        //the end of this loop
      while(elem_bool){
        
        elems=elem_group->get_element_array(err); MSQ_ERRZERO(err);
        num_elems=elem_group->num_elements();
        
        int element_counter=0;
          //loop over the elems in this group
        while(element_counter<num_elems){
          metric_counter=0;
            //increment the number of elements in mesh (not the group)
            //later, we can delete this, and give a call to meshSet
            //which will return the number of elements.
          ++total_num_elements;
          
            //metrics can be either element_based or vertex based
            //and we must treat them differently.  I didn't realize
            //this originally.  If we get an element group, that
            //ensures that the different groups do not have overlapping
            //elements.  However, they do have overlapping vertices????
            //When we place metrics in a QualityAssessor's list
            //we place element based metrics in the front and
            //vertex_based metrics in the back (push_front or push_back
            //respectively).  So, we can loop over the metrics until
            //we reach num_elem_based, and we know all of these
            //metrics are elment based.
          while(metric_counter<num_elem_based){
              //Message::print_info("\nremove_this_var = %i, num_elems=%i",remove_this_var,num_elems);
              //if first pass or if two passes are required
            if(!num_pass||assessor_array[metric_counter]->maxHist>MSQ_MAX_CAP){
              
              bool b = assessor_array[metric_counter]->metric->
                evaluate_element(*elem_group,&elems[element_counter], 
                temp_val, err);  MSQ_ERRZERO(err);
              if(!b){
                QAData[metric_counter].numInvalid++;
              }
              
                //if we are on the first loop over the mesh, calculate
                //everything we can.  That is accumlate for
                //avergae, max, min, rms, and stddev (stdVar
                //accumalates (temp_val^2) just as rms does)
              if(!num_pass){
                if(assessor_array[metric_counter]->funcFlagBits&AVERAGE)
                  QAData[metric_counter].avgVar+=temp_val;
                if(assessor_array[metric_counter]->funcFlagBits&MAXIMUM)
                  if(temp_val>QAData[metric_counter].maxVar)
                    QAData[metric_counter].maxVar=temp_val;
                if(assessor_array[metric_counter]->funcFlagBits&MINIMUM)
                  if(temp_val<QAData[metric_counter].minVar)
                    QAData[metric_counter].minVar=temp_val;
                if(assessor_array[metric_counter]->funcFlagBits&RMS)
                  QAData[metric_counter].rmsVar+=(temp_val*temp_val);
                  //calculate std_dev    
                if(assessor_array[metric_counter]->funcFlagBits&STDDEV){
                  QAData[metric_counter].stdVar+=(temp_val*temp_val);
                }
              }
                //Calculate histogram if needed
              if(QAData[metric_counter].histMax<MSQ_MAX_CAP &&
                 assessor_array[metric_counter]->funcFlagBits&HISTOGRAM){
                int hist_counter;
                  //FIrst and last bins are for "out of range" values
                  //If temp_val is less that lower histogram bound
                if(temp_val<QAData[metric_counter].histMin-MSQ_MIN){
                  QAData[metric_counter].histVar[0]++;
                }
                  //if temp_val is greater than uper histogram bound
                else if(temp_val>QAData[metric_counter].histMax+MSQ_MIN){
                  QAData[metric_counter].histVar[MSQ_HIST_SIZE+1]++;
                }
                else{
                    //If temp_val is in one of the first MSQ_HIST_SIZE
                    //(the number of divisions of the histogram) then
                    //increment the appropiate slot in the histogram.
                  for(hist_counter=0;hist_counter<MSQ_HIST_SIZE-1;hist_counter++)
                  {
                    if(temp_val<=(QAData[metric_counter].histMin)+
                       (QAData[metric_counter].histDelta*(hist_counter+1))){
                      QAData[metric_counter].histVar[hist_counter+1]++;
                      hist_counter+=MSQ_HIST_SIZE;
                    }
                  }
                    //If temp_val is in the last slot of the histogram
                    //increment the last slot.
                  if(temp_val>(QAData[metric_counter].histMin)+
                     (QAData[metric_counter].histDelta*(MSQ_HIST_SIZE-1))){
                    QAData[metric_counter].histVar[MSQ_HIST_SIZE]++;
                  }
                  
                }//end else (meaning if in range (minHist,maxHist) 
              }//end if histogram
            }//end if first pass or second required for this metric
            metric_counter++;
          }//end while metric_counter is less than num_elem_based
          element_counter++;
        }//end  while element counter < num_elems

          // If dealing with local patches, get next element group (PatchData object)
        if (get_patch_type() != PatchData::GLOBAL_PATCH)
          elem_bool=ms.get_next_patch(*elem_group,this, err); MSQ_ERRZERO(err);
          //Michael:: Since we are doing global right now:
          //Remove this when no longer doing global
        elem_bool=0;
          //Message::print_info("\nInside QA get_next returning %i",elem_bool);        
      }//end  while (elem_bool)

    }//end   if num_elem_based

    if((num_metrics-num_elem_based)>0){
        //construct the patch we will send to get_next_patch
      PatchData* vert_group;
      PatchData local_vert_group;
      no_culling_method();
      metric_counter=0;
      bool vert_bool=true;
      if (get_global_patch() == 0) {
        vert_group = &local_vert_group; 
        vert_bool=ms.get_next_patch(*vert_group, this, err);  MSQ_ERRZERO(err);
      }
      else {
        vert_group = get_global_patch();
      }
      
        //until there are no more patches
        //there is another get_next_patch at
        //the end of this loop
      while(vert_bool){
        
        verts=vert_group->get_vertex_array(err);  MSQ_ERRZERO(err);
        num_verts=vert_group->num_vertices();
        
        int vertex_counter=0;
          //loop over the verts in this group
        while(vertex_counter<num_verts){
            //increment the number of verticess in mesh (not the group)
            //later, we can delete this, and give a call to meshSet
            //which will return the number of verts.
          ++total_num_vertices;
          
          
            //metrics can be either element_based or vertex based
            //and we must treat them differently.  I didn't realize
            //this originally.  If we get an element group, that
            //ensures that the different groups do not have overlapping
            //elements.  However, they do have overlapping vertices????
            //When we place metrics in a QualityAssessor's list
            //we place element based metrics in the front and
            //vertex_based metrics in the back (push_front or push_back
            //respectively).  So, we can loop over the metrics until
            //we reach num_elem_based, and we know all of these
            //metrics are elment based.
          metric_counter=num_elem_based;
          
          while(metric_counter<num_metrics){
           
              //Message::print_info("\nremove_this_var = %i, num_elems=%i",remove_this_var,num_elems);
              //if first pass or if two passes are required
            if(!num_pass||assessor_array[metric_counter]->maxHist>MSQ_MAX_CAP){
              
              if(!assessor_array[metric_counter]->metric->evaluate_vertex(*vert_group,&verts[vertex_counter],temp_val,err)){
                QAData[metric_counter].numInvalid++;
              }
              MSQ_ERRZERO(err);
              
                //if we are on the first loop over the mesh, calculate
                //everything we can.  That is accumlate for
                //avergae, max, min, rms, and stddev (stdVar
                //accumalates (temp_val^2) just as rms does)
              if(!num_pass){
                if(assessor_array[metric_counter]->funcFlagBits&AVERAGE)
                  QAData[metric_counter].avgVar+=temp_val;
                if(assessor_array[metric_counter]->funcFlagBits&MAXIMUM)
                  if(temp_val>QAData[metric_counter].maxVar)
                    QAData[metric_counter].maxVar=temp_val;
                if(assessor_array[metric_counter]->funcFlagBits&MINIMUM)
                  if(temp_val<QAData[metric_counter].minVar)
                    QAData[metric_counter].minVar=temp_val;
                if(assessor_array[metric_counter]->funcFlagBits&RMS)
                  QAData[metric_counter].rmsVar+=(temp_val*temp_val);
                  //calculate std_dev    
                if(assessor_array[metric_counter]->funcFlagBits&STDDEV){
                  QAData[metric_counter].stdVar+=(temp_val*temp_val);
                }
              }
                //Calculate histogram if needed
              if(QAData[metric_counter].histMax<MSQ_MAX_CAP &&
                 assessor_array[metric_counter]->funcFlagBits&HISTOGRAM){
                int hist_counter;
                  //FIrst and last bins are for "out of range" values
                  //If temp_val is less that lower histogram bound
                if(temp_val<QAData[metric_counter].histMin-MSQ_MIN){
                  QAData[metric_counter].histVar[0]++;
                }
                  //if temp_val is greater than uper histogram bound
                else if(temp_val>QAData[metric_counter].histMax+MSQ_MIN){
                  QAData[metric_counter].histVar[MSQ_HIST_SIZE+1]++;
                }
                else{
                    //If temp_val is in one of the first MSQ_HIST_SIZE
                    //(the number of divisions of the histogram) then
                    //increment the appropiate slot in the histogram.
                  for(hist_counter=0;hist_counter<MSQ_HIST_SIZE-1;hist_counter++)
                  {
                    if(temp_val<=(QAData[metric_counter].histMin)+
                       (QAData[metric_counter].histDelta*(hist_counter+1))){
                      QAData[metric_counter].histVar[hist_counter+1]++;
                      hist_counter+=MSQ_HIST_SIZE;
                    }
                  }
                    //If temp_val is in the last slot of the histogram
                    //increment the last slot.
                  if(temp_val>(QAData[metric_counter].histMin)+
                     (QAData[metric_counter].histDelta*(MSQ_HIST_SIZE-1))){
                    QAData[metric_counter].histVar[MSQ_HIST_SIZE]++;
                  }
                  
                }//end else (meaning if in range (minHist,maxHist) 
              }//end if histogram
            }//end if first pass or second required for this metric
            metric_counter++;
          }//end while metric_counter is less than num_metrics (i.e., vertex)
          vertex_counter++;
        }//end  while vertex counter < num_verts
        
          //get next vertex group (PatchData object)
        
        if (get_patch_type() != PatchData::GLOBAL_PATCH)
          vert_bool=ms.get_next_patch(*vert_group,this, err); MSQ_ERRZERO(err);
          //Michael:: Since we are doing global right now:
          //Remove this when no longer doing global
        vert_bool=0;
          //Message::print_info("\nInside QA get_next returning %i",elem_bool);        
      }//end  while (vert_bool)

    }//end   if (num_metrics-num_elem_based)>0


    
    int met_i;
      //Now, we have finished accumulating the data (for one
      //of the passes over the mesh).  We must now:
      //find histDelta (the size of each slot in the histogram)
      //if needed... also finish off avg, rms, stddev
    
    for (met_i=0;met_i<num_metrics;met_i++){
      if(met_i<num_elem_based)
        tot_num=total_num_elements;
      else
        tot_num=total_num_vertices;
        //if first pass
      if(!num_pass){
          //compute histDelta if necessary
        if(assessor_array[met_i]->maxHist>MSQ_MAX_CAP  &&
           assessor_array[met_i]->funcFlagBits&HISTOGRAM){
          QAData[met_i].histMax=QAData[met_i].maxVar;
          QAData[met_i].histMin=QAData[met_i].minVar;
          two_loops=1;
          QAData[met_i].histDelta=(QAData[met_i].maxVar- QAData[met_i].minVar)/
              ( double ) MSQ_HIST_SIZE;
          
        }
          //make sure no div by 0
        if(tot_num==0){
          QAData[met_i].avgVar=0.0;
          QAData[met_i].stdVar=0.0;
          QAData[met_i].rmsVar=0.0;
        }
        else{
            //divide by tot_num to get average
          if(assessor_array[met_i]->funcFlagBits&AVERAGE)
            QAData[met_i].avgVar/=(( double ) tot_num);
            //complete std_dev computation  
          if(assessor_array[met_i]->funcFlagBits&STDDEV){
            QAData[met_i].stdVar=(QAData[met_i].stdVar/ (double) tot_num)
                                       - (QAData[met_i].avgVar*
                                          QAData[met_i].avgVar );
              //make sure to not sqrt a negative number
            if(QAData[met_i].stdVar<MSQ_MIN)
              QAData[met_i].stdVar=0.0;
            else
              QAData[met_i].stdVar=sqrt(QAData[met_i].stdVar);
          }
          
            //complete rms computation
          if(assessor_array[met_i]->funcFlagBits&RMS)
            QAData[met_i].rmsVar=sqrt(QAData[met_i].rmsVar/
                                      (( double ) tot_num));
        }//end else tot_num!=0
      }//end if first pass
      
    }
    num_pass++;
  }//end of num_pass<2
  
    ///////////////////////////////////////////////
    //PRINT TABLE OF VALUES  //////////////////////
    ///////////////////////////////////////////////
    /*!\TODO (Michael) NOTES: delete div 2  when lose pass=2.
       Fix the way and format of display*/
  metric_counter=0;

  int column_counter=0;
  if(!printingTurnedOff){
    
    MSQ_PRINT(1)("\n************  Quality summary of MeshSet  *****************\n ");
    tot_num=0;
      //For each metric, print the data that we have calculated.
    for(metric_counter=0;metric_counter<num_metrics;metric_counter++){
        //if elem based, print element header
      if(metric_counter<num_elem_based){
        if(num_pass==1)
          tot_num=total_num_elements;
        else
          tot_num=total_num_elements/(2);
        MSQ_PRINT(1)("\nELEMENT BASED METRIC :: %s  (%i elements)\n",assessor_array[metric_counter]->metric->get_name().c_str(),tot_num);
      }
        //if metric is vertex_based, print vertex_header
      else{
        if(num_pass==1)
          tot_num=total_num_vertices;
        else
          tot_num=total_num_vertices/(2);
        MSQ_PRINT(1)("\nVERTEX BASED METRIC :: %s (%i vertices)\n",assessor_array[metric_counter]->metric->get_name().c_str(),tot_num);
      }
      column_counter=0;
        //PRINT only a warning if invalid
      if(num_pass==2)
        QAData[metric_counter].numInvalid/=(2);
      else if (num_pass!=1) {
        MSQ_SETERR(err)("Incorrect number of passes over the mesh", MsqError::INTERNAL_ERROR);
        return 0.0;
      }
      
      if(QAData[metric_counter].numInvalid>0){
        MSQ_PRINT(1)("MESH IS INVALID WITH RESPECT TO THIS METRIC (%i values)\n",QAData[metric_counter].numInvalid);
      }
        //mesh was valid with respect to the metric
      else{
        if(assessor_array[metric_counter]->funcFlagBits&MINIMUM)
        {
          MSQ_PRINT(1)("%s = %8.6e   ", get_QAFunction_name(QualityAssessor::MINIMUM).c_str(),QAData[metric_counter].minVar);
          column_counter++;
            //return_value=QAData[metric_counter].minVar;
        }
        if(assessor_array[metric_counter]->funcFlagBits&MAXIMUM)
        {
          MSQ_PRINT(1)("%s = %8.6e   ", get_QAFunction_name(QualityAssessor::MAXIMUM).c_str(),QAData[metric_counter].maxVar);
          column_counter++;
            //return_value=QAData[metric_counter].maxVar;
        }
          //new line every two entries
        if(column_counter==2){
          column_counter=0;
          MSQ_PRINT(1)("\n");
        }
        
        if(assessor_array[metric_counter]->funcFlagBits&AVERAGE)
        {
          MSQ_PRINT(1)("%s = %8.6e   ",get_QAFunction_name(QualityAssessor::AVERAGE).c_str(),QAData[metric_counter].avgVar);
          column_counter++;
            //return_value=QAData[metric_counter].avgVar;
        }
          //new line every two entries
        if(column_counter==2){
          column_counter=0;
          MSQ_PRINT(1)("\n");
        }
        if(assessor_array[metric_counter]->funcFlagBits&RMS)
        {
          MSQ_PRINT(1)("%s = %8.6e   ", get_QAFunction_name(QualityAssessor::RMS).c_str(),QAData[metric_counter].rmsVar);
          column_counter++;
            //return_value=QAData[metric_counter].rmsVar;
        }
          //new line every two entries
        if(column_counter==2){
          column_counter=0;
          MSQ_PRINT(1)("\n");
        }
        if(assessor_array[metric_counter]->funcFlagBits&STDDEV)
        {
          MSQ_PRINT(1)("%s = %8.6e   ", get_QAFunction_name(QualityAssessor::STDDEV).c_str(),QAData[metric_counter].stdVar);
            //return_value=QAData[metric_counter].stdVar;
        }
        if(assessor_array[metric_counter]->funcFlagBits&HISTOGRAM)
        {
          MSQ_PRINT(1)("\n%s  \n", get_QAFunction_name(QualityAssessor::HISTOGRAM).c_str());
          int inner_loop;
          if(QAData[metric_counter].histVar[0])
            MSQ_PRINT(1)("NOTICE:  VALUES BELOW HISTOGRAM RANGE = %d\n",QAData[metric_counter].histVar[0]);
          for (inner_loop=0;inner_loop<MSQ_HIST_SIZE; inner_loop++){
            MSQ_PRINT(1)("Values between %8.6e and %8.6e = %d\n",QAData[metric_counter].histMin+(QAData[metric_counter].histDelta*inner_loop),QAData[metric_counter].histMin+(QAData[metric_counter].histDelta*(inner_loop+1)),QAData[metric_counter].histVar[inner_loop+1]);
          }//end inner_loop       
          if(QAData[metric_counter].histVar[MSQ_HIST_SIZE+1])
            MSQ_PRINT(1)("N       OTICE:  VALUES ABOVE HISTOGRAM RANGE = %d\n",QAData[metric_counter].histVar[MSQ_HIST_SIZE+1]);
        }//end if histo         
      }//end else (else the mesh was valid)
    }//end while loop over metrics
    MSQ_PRINT(1)("\n");
  }

  for(metric_counter=0;metric_counter<num_metrics;metric_counter++){
    if(assessor_array[metric_counter]->metric==stoppingMetric){
      switch(stoppingFunction){
        case(MAXIMUM):
          return_value=QAData[metric_counter].maxVar;
          break;
        case(MINIMUM):
          return_value=QAData[metric_counter].minVar;
          break;
        case(STDDEV):
          return_value=QAData[metric_counter].stdVar;
          break;
        case(RMS):
          return_value=QAData[metric_counter].rmsVar;
          break;
        case(AVERAGE):
          return_value=QAData[metric_counter].avgVar;
          break;
        default:
          MSQ_SETERR(err)("QAFunction used for return value is invalid.",MsqError::INVALID_STATE);
          return 0;
      };
    }
  }
    
  return return_value;
}

} //namespace Mesquite
