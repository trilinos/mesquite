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
#include <math.h>
#include "MsqMessage.hpp"

using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "QualityAssessor::QualityAssessor"
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
                                 std::string name) :
  qualityAssessorName(name)
{
  MsqError err;
  if(func==QualityAssessor::HISTOGRAM || func==QualityAssessor::ALL_MEASURES){
    add_quality_assessment(qm, func, err);
    set_stopping_assessment(qm, QualityAssessor::MAXIMUM, err);
  }
  else{
    set_stopping_assessment(qm,func, err);
  }
  MSQ_CHKERR(err);
  printingTurnedOff=0;
    //sigDig=5;
}

#undef __FUNC__
#define __FUNC__ "QualityAssessor::~QualityAssessor"
QualityAssessor::~QualityAssessor()
{
  Assessor* tmp;
  int i;
  for(i=0;i<assessList.size();i++){
    tmp=assessList.front();
    delete tmp;
    assessList.pop_front();
  }
}

#undef __FUNC__
#define __FUNC__ "QualityAssessor::get_QAFunction_name"
std::string QualityAssessor::get_QAFunction_name(
                              enum QualityAssessor::QAFunction fun)
{
  switch(fun){
    case(AVERAGE):
      return "Average   ";
      break;
    case(HISTOGRAM):
      return "Histogram of metric values: ";
      break;
    case(MAXIMUM):
      return "Maximum   ";
      break;
    case(MINIMUM):
      return "Minimum   ";
      break;
    case(RMS):
      return "RMS       ";
    case(STDDEV):
      return "Stan. Dev.";
      break;
    default:
      return "DEFAULT   ";
  };
}


#undef __FUNC__
#define __FUNC__ "QualityAssessor::add_quality_assessment"
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
                                             MsqError &err)
{ 
  int found=0;
  std::list<Assessor*>::iterator pos;
  pos=assessList.begin();
  Assessor* assess_ptr=NULL;
    //loop over the assessList (list of Assessor of this QA.)
  while (pos!=assessList.end() && !found){
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
    if(assess_ptr->metric->get_evaluation_mode()==QualityMetric::VERTEX){
      assessList.push_back( assess_ptr );
    }
    else{
      assessList.push_front( assess_ptr );
    }
  }
}

#undef __FUNC__
#define __FUNC__ "QualityAssessor::set_stopping_assessment"
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
    err.set_msg("HISTOGRAM DOES NOT GIVE A VALID RETURN VALUE");
  }
  add_quality_assessment(qm,func,err);
  MSQ_CHKERR(err);
  stoppingFunction=func;
  stoppingMetric=qm;
}

  


#undef __FUNC__
#define __FUNC__ "QualityAssessor::set_histogram_range"
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
                                             MsqError &err)
{
  int found=0;
  std::list<Assessor*>::iterator pos;
  pos=assessList.begin();
  Assessor* assess_ptr=NULL;
    //loop over the assessList (list of Assessor of this QA.)
  int i=0;
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
    if(assess_ptr->metric->get_evaluation_mode()==QualityMetric::VERTEX){
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



#undef __FUNC__
#define __FUNC__ "QualityAssessor::assess_mesh_quality"
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
double QualityAssessor::assess_mesh_quality(MeshSet &ms, MsqError &err)
{
    /*! \TODO (Michael) Begin NOTE
      It does not currently work with vertex_based
      metrics.  
    */
    //
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
  std::list<Assessor*>::iterator pos = assessList.begin();
  std::list<QualityAssessor::QAFunction>::iterator func_pos;
  int num_elem_based=0;
  int num_metrics=assessList.size();
  Assessor** assessor_array = new Assessor*[num_metrics];
    //array of data holders QAVars.  QAVars is a struct
    //defined in QualityAssessor.hpp which creates a
    //way to store the data we accumulate for the each metric.
    //Thus, we have one QAVars, per metric in our lis.
    //Each QAVar also has an bit flag int, funcFlagBits,
    //which tells us (after we initialize it) what we
    //need to calcuate... max, min, avg, rms, hist, std_dev.
  QAVars *QAData = new QAVars[num_metrics];
  
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
    if((*pos)->metric->get_evaluation_mode()!=QualityMetric::VERTEX)
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
      //initialize max, min, rms, stddev.
    QAData[i].maxVar=0.0;
    QAData[i].minVar=MSQ_MAX_CAP;
    QAData[i].rmsVar=0.0;
    QAData[i].stdVar=0.0;
    QAData[i].histMax=assessor_array[i]->maxHist;
    QAData[i].histMin=assessor_array[i]->minHist;
    QAData[i].histDelta=(QAData[i].histMax-QAData[i].histMin)/MSQ_HIST_SIZE;
    i++;
    ++pos;
  }
    //pointers to elems
  MsqMeshEntity* elems;
  int num_elems=0;
    //int num_vertices;
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
  int num_nodes_per_elem=0;
  
  while(num_pass<2 && two_loops){
      //two_loops will be set to one if two loops are necessary
    two_loops=0;
    double temp_val=0;
    if(num_elem_based){
        //construct the patch we will send to get_next_element_group
      PatchData elem_group;
        //Michael:: temporary solution to bug
      ms.set_patch_type(MeshSet::ELEMENTS_ON_VERTEX_PATCH, 1);
      ms.copy_culling_method_bits(0);
      
      bool elem_bool;
      elem_bool=ms.get_next_element_group(elem_group, err);
      
        //until there are no more element groups
        //there is another get_next_element_group at
        //the end of this loop
        //int remove_this_var=0;
      while(elem_bool){
          //remove_this_var++;
        
        elems=elem_group.get_element_array(err);
        num_elems=elem_group.num_elements();
        
        int element_counter=0;
          //loop over the elems in this group
        while(element_counter<num_elems){
            //increment the number of elements in mesh (not the group)
            //later, we can delete this, and give a call to meshSet
            //which will return the number of elements.
          ++total_num_elements;
          int metric_counter=0;
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
              //PRINT_INFO("\nremove_this_var = %i, num_elems=%i",remove_this_var,num_elems);
              //if first pass or if two passes are required
            if(!num_pass||assessor_array[metric_counter]->maxHist>MSQ_MAX_CAP){
              
              temp_val=assessor_array[metric_counter]->metric->evaluate_element(elem_group,&elems[element_counter],err);
              num_nodes_per_elem=elems[element_counter].vertex_count();
              
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
        }//end         while element counter < num_elems

          //get next element group (PatchData object)
        elem_bool=ms.get_next_element_group(elem_group, err);
      }//end  while (elem_bool)
      
    }//end   if num_elem_based
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
    /*!\TODO (Michael) NOTES: delete div 2  when lose pass=2
      div 6 is because we use patch instead of group
      fix the way and format of display*/
  int metric_counter=0;
  int column_counter=0;
  if(!printingTurnedOff){
    
    PRINT_INFO("\n************  Quality summary of MeshSet  *****************\n ");
    tot_num=0;
      //Fo  r each metric, print the data that we have calculated.
    for(metric_counter=0;metric_counter<num_metrics;metric_counter++){
        //if elem based, print element header
      if(metric_counter<num_elem_based){
        if(num_pass==1)
          tot_num=total_num_elements/(num_nodes_per_elem);
        else
          tot_num=total_num_elements/(num_nodes_per_elem*2);
        PRINT_INFO("\nELEMENT BASED METRIC :: %s  (%i elements)\n",assessor_array[metric_counter]->metric->get_name().c_str(),tot_num);
      }
        //if metric is vertex_based, print vertex_header
      else{
        tot_num=total_num_vertices/2;
        PRINT_INFO("VERTEX BASED METRIC :: %s (%i vertices)\n",assessor_array[metric_counter]->metric->get_name().c_str(),tot_num);
      }
      column_counter=0;
      if(assessor_array[metric_counter]->funcFlagBits&MINIMUM)
      {
        PRINT_INFO("%s = %8.6e   ", get_QAFunction_name(QualityAssessor::MINIMUM).c_str(),QAData[metric_counter].minVar);
        column_counter++;
          //return_value=QAData[metric_counter].minVar;
      }
      if(assessor_array[metric_counter]->funcFlagBits&MAXIMUM)
      {
        PRINT_INFO("%s = %8.6e   ", get_QAFunction_name(QualityAssessor::MAXIMUM).c_str(),QAData[metric_counter].maxVar);
        column_counter++;
          //return_value=QAData[metric_counter].maxVar;
      }
        //new line every two entries
      if(column_counter==2){
        column_counter=0;
        PRINT_INFO("\n");
      }
      
      if(assessor_array[metric_counter]->funcFlagBits&AVERAGE)
      {
        PRINT_INFO("%s = %8.6e   ",get_QAFunction_name(QualityAssessor::AVERAGE).c_str(),QAData[metric_counter].avgVar);
        column_counter++;
          //return_value=QAData[metric_counter].avgVar;
      }
        //new line every two entries
      if(column_counter==2){
        column_counter=0;
        PRINT_INFO("\n");
      }
      if(assessor_array[metric_counter]->funcFlagBits&RMS)
      {
        PRINT_INFO("%s = %8.6e   ", get_QAFunction_name(QualityAssessor::RMS).c_str(),QAData[metric_counter].rmsVar);
        column_counter++;
          //return_value=QAData[metric_counter].rmsVar;
      }
        //new line every two entries
      if(column_counter==2){
        column_counter=0;
        PRINT_INFO("\n");
      }
      if(assessor_array[metric_counter]->funcFlagBits&STDDEV)
      {
        PRINT_INFO("%s = %8.6e   ", get_QAFunction_name(QualityAssessor::STDDEV).c_str(),QAData[metric_counter].stdVar);
          //return_value=QAData[metric_counter].stdVar;
      }
      if(assessor_array[metric_counter]->funcFlagBits&HISTOGRAM)
      {
         PRINT_INFO("\n%s  \n", get_QAFunction_name(QualityAssessor::HISTOGRAM).c_str());
         int inner_loop;
        if(QAData[metric_counter].histVar[0])
          PRINT_INFO("NOTICE:  VALUES BELOW HISTOGRAM RANGE = %d\n",QAData[metric_counter].histVar[0]/num_nodes_per_elem);
        for (inner_loop=0;inner_loop<MSQ_HIST_SIZE; inner_loop++){
          PRINT_INFO("Values between %8.6e and %8.6e = %d\n",QAData[metric_counter].histMin+(QAData[metric_counter].histDelta*inner_loop),QAData[metric_counter].histMin+(QAData[metric_counter].histDelta*(inner_loop+1)),QAData[metric_counter].histVar[inner_loop+1]/num_nodes_per_elem);
        }//end inner_loop       
        if(QAData[metric_counter].histVar[MSQ_HIST_SIZE+1])
          PRINT_INFO("NOTICE:  VALUES ABOVE HISTOGRAM RANGE = %d\n",QAData[metric_counter].histVar[MSQ_HIST_SIZE+1]/num_nodes_per_elem);
      }//end if histo
    
    }//end while loop over metrics
    PRINT_INFO("\n");
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
          err.set_msg("QAFunction used for return value is invalid.");
      };
    }
  }
    
  delete assessor_array;
  delete QAData;
  return return_value;
}
