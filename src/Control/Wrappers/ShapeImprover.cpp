#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ShapeImprover.hpp"
#include "InstructionQueue.hpp"
#include "MeanRatioQualityMetric.hpp"
#include "LPTemplate.hpp"
#include "SteepestDescent.hpp"
#include "TerminationCriterion.hpp"


Mesquite::ShapeImprover::ShapeImprover()
{}

Mesquite::ShapeImprover::~ShapeImprover()
{}

void Mesquite::ShapeImprover::improve_quality(MeshSet &mesh_set,
                                              MsqError &err)
{
    // Create an intruction queue  
  Mesquite::InstructionQueue q;
  
    // Create a mean ratio quality metric
  Mesquite::ShapeQualityMetric* mean_ratio =
    Mesquite::MeanRatioQualityMetric::create_new();
  
    // Build an objective function with it
  Mesquite::LPTemplate* obj_func =
    new Mesquite::LPTemplate(mean_ratio, 2, err);
  
    // Create the steepest descent optimization procedures
  Mesquite::SteepestDescent* pass1 =
    new Mesquite::SteepestDescent( obj_func );
  Mesquite::TerminationCriterion sc1, sc2;
  sc1.add_criterion_type_with_int(Mesquite::TerminationCriterion::NUMBER_OF_ITERATES, 1, err);
  sc2.add_criterion_type_with_int(Mesquite::TerminationCriterion::NUMBER_OF_ITERATES, 1, err);
  pass1->set_outer_termination_criterion(&sc1);
  pass1->set_inner_termination_criterion(&sc2);
  pass1->add_culling_method(Mesquite::PatchData::NO_BOUNDARY_VTX);
  
  q.set_master_quality_improver(pass1, err); MSQ_CHKERR(err);
  q.run_instructions(mesh_set, err); MSQ_CHKERR(err);
}
