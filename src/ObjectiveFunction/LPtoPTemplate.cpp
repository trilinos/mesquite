/*!
  \file   LPtoPTemplate.cpp
  \brief  

  This Objective Function is evaluated using an L P norm to the pth power.
  total=(sum (x_i)^pVal)
  \author Michael Brewer
  \date   2002-01-23
*/
#include <math.h>
#include "LPtoPTemplate.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MsqMessage.hpp"
using  namespace Mesquite;  

using std::cout;
using std::endl;

#undef __FUNC__
#define __FUNC__ "LPtoPTemplate::LPtoPTemplate"

LPtoPTemplate::LPtoPTemplate(QualityMetric *qualitymetric, int Pinput, MsqError &err){
  set_quality_metric(qualitymetric);
  pVal=Pinput;
  if(pVal<1){
    err.set_msg("P_VALUE must be greater than 0.");
  }
  set_feasible(qualitymetric->get_feasible_constraint());
  set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
    //set_use_local_gradient(true);
  set_negate_flag(qualitymetric->get_negate_flag());
}

#undef __FUNC__
#define __FUNC__ "LPtoPTemplate::~LPtoPTemplate"

//Michael:  need to clean up here
LPtoPTemplate::~LPtoPTemplate(){

}

#undef __FUNC__
#define __FUNC__ "LPtoPTemplate::concrete_evaluate"
bool LPtoPTemplate::concrete_evaluate(PatchData &patch, double &fval,
                                      MsqError &err){
  int index=0;
  MsqMeshEntity* elems=patch.get_element_array(err);
  bool obj_bool=true;
    //double check for pVal=0;
  if(pVal==0){
    err.set_msg("pVal equal zero not allowed.  L_0 is not a valid norm.");
    return 0;
  }
  
    //Michael:  this may not do what we want
    //Set currentQM to be the first quality metric* in the list
  QualityMetric* currentQM = get_quality_metric();
  if(currentQM==NULL)
    currentQM=get_quality_metric_list().front();
  if(currentQM==NULL)
    err.set_msg("NULL QualityMetric pointer in LPtoPTemplate");
  int num_elements=patch.num_elements();
  int num_vertices=patch.num_vertices();
  int total_num=0;
  if(currentQM->get_metric_type()==QualityMetric::ELEMENT_BASED)   
    total_num=num_elements;
  else if (currentQM->get_metric_type()==QualityMetric::VERTEX_BASED)
    total_num=num_vertices;
  else
    err.set_msg("Make sure MetricType is initialised in concrete QualityMetric constructor.");
  double *metric_values= new double[total_num];
  if(currentQM->get_metric_type()==QualityMetric::ELEMENT_BASED)
  {
    for (index=0; index<num_elements;index++)
    {
        //if invalid return false after clean-up
      obj_bool=currentQM->evaluate_element(patch, (&elems[index]),
                                           metric_values[index], err);
      MSQ_CHKERR(err);
      if(!obj_bool){
        fval=0.0;
        delete[] metric_values;
        return false;
      }
      
      metric_values[index]=fabs(metric_values[index]);
      MSQ_DEBUG_ACTION(3,{std::cout<< "      o  Quality metric value for element "
                          << index << "\t: " << metric_values[index] << "\n";});
    }
  }
  else if(currentQM->get_metric_type()==QualityMetric::VERTEX_BASED)
  {
    MsqVertex* vertices=patch.get_vertex_array(err);
    for (index=0; index<num_vertices;index++)
    {
        //evaluate metric for this vertex
      obj_bool=currentQM->evaluate_vertex(patch, (&vertices[index]),
                                          metric_values[index], err);
      MSQ_CHKERR(err);
        //if invalid return false after clean-up
      if(!obj_bool){
        fval=0.0;
        delete[] metric_values;
        return false;
      }
      
      metric_values[index]=fabs(metric_values[index]);
    }
  }
  fval=compute_function(metric_values, total_num, err);
  delete[] metric_values;
  return true;
}

#undef __FUNC__
#define __FUNC__ "LPtoPTemplate::compute_analytical_gradient"
/*! \fn LPtoPTemplate::compute_analytical_gradient(PatchData &patch, Vector3D *const &grad, MsqError &err, int array_size)
    \param patch The PatchData object for which the objective function
           gradient is computed.
    \param grad An array of Vector3D, at least the size of the number
           of free vertices in the patch.
    \param array_size is the size of the grad Vector3D[] array and must correspond
           to the number of free vertices in the patch. This argument will be used to
           perform a costly check (counting the number of free vertices) is MSQ_DBG1
           or higher is defined.
*/
bool LPtoPTemplate::compute_analytical_gradient(PatchData &patch,
                                              Vector3D *const &grad,
                                              MsqError &err, int array_size)
{
  int num_elements=patch.num_elements();
  int num_vertices=patch.num_vertices();
  if( num_vertices!=array_size && array_size>=0)
    err.set_msg("Incorrect array size.");
  //Generate vertex to element connectivity if needed
  patch.generate_vertex_to_element_data();
  //vector for storing indices of vertex's connected elems
  std::vector<size_t> elem_on_vert_ind;
  std::vector<size_t> vert_on_vert_ind;
  MsqMeshEntity* elems=patch.get_element_array(err);
  MsqVertex* vertices=patch.get_vertex_array(err);
  bool obj_bool=true;
  // If MSQ_DBG1 is defined, check to make sure that num_vert == array_size.
  MSQ_DEBUG_ACTION(1,{
    int num_vert=patch.num_vertices(); 
    if(num_vert!=array_size && array_size!=-1){
      err.set_msg("Analytical Gradient passed arrays of incorrect size.");
      MSQ_CHKERR(err); cout << num_vert << " instead of " << array_size << endl; }
  });
    
  double big_f=0;
  double temp_value=0;
  int index=0;
  
  //Set currentQM to be the first quality metric* in the list
  QualityMetric* currentQM = get_quality_metric();
  if(currentQM==NULL)
    err.set_msg("LPtoPTemplate has NULL QualityMetric pointer.");
  enum QualityMetric::MetricType qm_type=currentQM->get_metric_type();
  
  int total_num=0;
  if (qm_type==QualityMetric::ELEMENT_BASED)
    total_num=num_elements;
  else if (qm_type==QualityMetric::VERTEX_BASED)
    total_num=num_vertices;
  else
    err.set_msg("Make sure MetricType is initialised in concrete QualityMetric constructor.");

  // in dF = F'F
  double *metric_values=new double[total_num];
  // fill array with element quality metrics
  if(qm_type==QualityMetric::ELEMENT_BASED){
    for (index=0; index<num_elements;++index){
      obj_bool=currentQM->evaluate_element(patch, &elems[index],
                                           metric_values[index], err);
      MSQ_CHKERR(err);
        //if invalid patch
      if(!obj_bool){
        delete[] metric_values;
        return false;
      }
      metric_values[index]=fabs(metric_values[index]);
    }
  }
  // fill array with vertex quality metrics
  else if (qm_type==QualityMetric::VERTEX_BASED) {
    for (index=0; index<num_vertices;++index){
      //evaluate metric for this vertex
      obj_bool=currentQM->evaluate_vertex(patch, &vertices[index],
                                          metric_values[index], err);
      MSQ_CHKERR(err);
        //if invalid patch
      if(!obj_bool){
        delete[] metric_values;
        return false;
      }
      metric_values[index] = fabs( metric_values[index] );
    }
  }
  big_f=get_negate_flag(); //if function negated for minimization, so is gradient.

  MsqFreeVertexIndexIterator free_ind(&patch, err);
  free_ind.reset();
  //position in patch's vertex array
  int vert_count=0;
  //corresponding position in grad array
  double dummy;
  //position in elem array
  size_t elem_pos=0;
  size_t vert_pos=0;
  
  Vector3D grad_vec;
  // Loops over free vertices
  
  for (vert_count=0; vert_count<num_vertices; ++vert_count) {
    grad[vert_count].set(0.,0.,0.);
    if (vertices[vert_count].is_free_vertex()) {
      grad[vert_count].set(0.0,0.0,0.0);
      temp_value=0;
      if(qm_type==QualityMetric::ELEMENT_BASED){
        
        patch.get_vertex_element_indices(vert_count, elem_on_vert_ind,err);
        size_t ele_num_vtces = elem_on_vert_ind.size();
        MsqVertex** ele_free_vtces = new MsqVertex*[ele_num_vtces];
        elem_pos=0;
          //while(el    em_pos<num_elements){
        while(!elem_on_vert_ind.empty()){
          elem_pos=(elem_on_vert_ind.back());
          elem_on_vert_ind.pop_back();
          ele_free_vtces[0] = &vertices[vert_count];
          currentQM->compute_element_gradient(patch, &elems[elem_pos],
                                              ele_free_vtces,
                                              &grad_vec, 1, dummy, err);
          temp_value=pVal;
          for(index=0;index<pVal-1;++index){
            temp_value*=(metric_values[elem_pos]);
          }
            //if pval is odd and met val is negative
          if(metric_values[elem_pos]<0  && pVal%2 ){
            temp_value*=(-1);
          }          
          grad[vert_count] += temp_value*grad_vec;
            //elem_pos++;         
        }
        delete []ele_free_vtces;     
      }
      else{
        
        patch.get_adjacent_vertex_indices(vert_count, vert_on_vert_ind,err);
          //For now we compute the metric for attached vertices and this
          //vertex, the above line gives us the attached vertices.  Now,
          //we must add this vertex.
        vert_on_vert_ind.push_back(vert_count);
        size_t vert_num_vtces = vert_on_vert_ind.size();
        MsqVertex** vert_free_vtces = new MsqVertex*[vert_num_vtces];
        vert_pos=0;
        while(!vert_on_vert_ind.empty()){
          vert_pos=(vert_on_vert_ind.back());
          vert_on_vert_ind.pop_back();
          vert_free_vtces[0] = &vertices[vert_count];
          currentQM->compute_vertex_gradient(patch, vertices[vert_pos],
                                             vert_free_vtces,
                                             &grad_vec, 1, dummy, err);
          temp_value=pVal;
          for(index=0;index<pVal-1;++index){
            temp_value*=(metric_values[vert_pos]);
          }
            //if pval is odd and met val is negative
          if(metric_values[vert_pos]<0  && pVal%2 ){
            temp_value*=(-1);
          }  
          grad[vert_count] += temp_value*grad_vec;
        }
        delete []vert_free_vtces;     
        
      }
      grad[vert_count]*=big_f;
        //PRINT_INFO("  gradx = %f, grady = %f, gradz = %f\n",grad[vert_count][0],grad[vert_count][1],grad[vert_count][2]);
    }//end if free
    
  }
  delete [] metric_values;
  return true;
}
  
	
#undef __FUNC__
#define __FUNC__ "LPtoPTemplate::compute_analytical_hessian"
/*! \fn LPtoPTemplate::compute_analytical_hessian(PatchData &pd, MsqHessian &hessian, MsqError &err)

    For each element, each entry to be accumulated in the Hessian for
    this objective function (\f$ \sum_{e \in E} Q(e)^p \f$ where \f$ E \f$
    is the set of all elements in the patch) has the form:
    \f$ pQ(e)^{p-1} \nabla^2 Q(e) + p(p-1)Q(e)^{p-2} \nabla Q(e) [\nabla Q(e)]^T \f$.

    For \f$ p=2 \f$, this simplifies to
    \f$ 2Q(e) \nabla^2 Q(e) + 2 \nabla Q(e) [\nabla Q(e)]^T \f$.

    For \f$ p=1 \f$, this simplifies to \f$ \nabla^2 Q(e) \f$.

    The \f$ p=1 \f$ simplified version is implemented directly
    to speed up computation. 
    
    \param patch The PatchData object for which the objective function
           hessian is computed.
    \param hessian: this object must have been previously initialized.
*/
bool LPtoPTemplate::compute_analytical_hessian(PatchData &pd,
                                               MsqHessian &hessian,
                                               MsqError &err)
{

  MsqMeshEntity* elements = pd.get_element_array(err); MSQ_CHKERR(err);
  MsqVertex* vertices = pd.get_vertex_array(err); MSQ_CHKERR(err);
  int num_elems = pd.num_elements();
  Matrix3D elem_hessian[MSQ_MAX_NUM_VERT_PER_ENT*(MSQ_MAX_NUM_VERT_PER_ENT+1)/2];
  Matrix3D elem_outer_product[MSQ_MAX_NUM_VERT_PER_ENT*(MSQ_MAX_NUM_VERT_PER_ENT+1)/2];
  Vector3D grad_vec[MSQ_MAX_NUM_VERT_PER_ENT];
  double metric_value;
  double fac1, fac2;
  Matrix3D grad_outprod;
//  Vector3D zero3D(0,0,0);
  QualityMetric* currentQM = get_quality_metric();
  
  MsqVertex* elem_vtx[MSQ_MAX_NUM_VERT_PER_ENT];
  std::vector< size_t > vtx_indices;
  std::vector<size_t>::const_iterator index; 
    
  size_t e;
  int num_vtx;
  int i,j,n;
  
  for (e=0; e<num_elems; ++e) {
    int nve = elements[e].vertex_count();
    
    // Gets a list of free vertices in the element.
    elements[e].get_vertex_indices(vtx_indices);
    num_vtx=0;
    for (index=vtx_indices.begin(); index!=vtx_indices.end(); ++index) {
      if ( vertices[*index].is_free_vertex() ) {
        elem_vtx[num_vtx] = vertices + (*index); ++num_vtx; }
    }
    
    // Computes \nabla^2 Q(e). Only the free vertices will have non-zero entries. 
    currentQM->compute_element_hessian(pd, elements+e, elem_vtx, grad_vec, elem_hessian,
                                       num_vtx, metric_value, err); MSQ_CHKERR(err);
     
    if (pVal == 1) {
      hessian.accumulate_entries(pd, e, elem_hessian, err);
    }
    else if (pVal >= 2) {
      // Computes \nabla Q(e) [\nabla Q(e)]^T 
      n=0;
      for (i=0; i<nve; ++i) {
        for (j=i; j<nve; ++j) {
          if ( vertices[vtx_indices[i]].is_free_vertex() &&
               vertices[vtx_indices[j]].is_free_vertex() ) {
            elem_outer_product[n].outer_product(grad_vec[i], grad_vec[j]);
          } else {
            elem_outer_product[n] = 0.;
          }
          ++n;
        }
      }

      // Computes  pQ(e)^{p-1}
      fac1 = pVal * pow(metric_value, pVal-1);
      // Computes p(p-1)Q(e)^{p-2}
      fac2 = pVal* (pVal-1) * pow(metric_value, pVal-2);

      for (i=0; i<nve*(nve+1)/2; ++i) {
        elem_hessian[i] *= fac1;
        elem_outer_product[i] *= fac2;
      }
      
      hessian.accumulate_entries(pd, e, elem_hessian, err);
      hessian.accumulate_entries(pd, e, elem_outer_product, err);

    } else {
      err.set_msg(" invalid P value.");
      return false;
    }
    
  }
  

  return true;
}
