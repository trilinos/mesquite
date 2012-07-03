/****************************************************************** 
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
   
  ******************************************************************/
/*!
  \file   NonGradient.cpp
  \brief  
  broken parts marked dmd 
  also exit(-1) ->?, != MSQ errrtn
  implement NonGradient class member functions.
  OF abbreviates ObjectiveFunction
  arrptr extracts the pointer from a std::vector<type>, Mesquite.hpp  
  Amoeba NonGradient Minimization Method, Numerical Recipes in C.
*/

#include "NonGradient.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"
#include <cmath>
#include <iostream>
#include "../../ObjectiveFunction/ObjectiveFunction.hpp"
//#include "../../ObjectiveFunction/ObjectiveFunction.hpp"

namespace MESQUITE_NS {

std::string NonGradient::get_name() const
  { return "NonGradient"; }
  
PatchSet* NonGradient::get_patch_set()
  { return PatchSetUser::get_patch_set(); }

NonGradient::NonGradient(ObjectiveFunction* of) 
  : VertexMover(of),
    PatchSetUser(true),
    projectGradient(false),
    mDimension(0),
    mThreshold(0.0),
    mTolerance(0.0),
    mMaxNumEval(0)
{
         std::cout <<  "  NonGradient constructor"  << std::endl; 
}  

bool
NonGradient::testRowSum( int numRow, int numCol, double* matrix, double* oldRowSum)
{
  bool same = true;
  std::vector<double> rowSum(numRow);
  for (int col=0;col<numCol;col++)    
  {
    for (int row=0;row<numRow;row++)
    {
      rowSum[row] += matrix[row+col*numRow];
    }
  }
  double machEps = 1.e-15;  // better to use system parameters
  for (int row=0;row<numRow;row++)
  {
    if( fabs( rowSum[row] -  oldRowSum[row]) >  machEps * fabs( rowSum[row] ) + fabs(oldRowSum[row]) ) 
    { 
         same = false;
         std::cout << rowSum[row] << "  rowSum[" << row << "]  "  << oldRowSum[row] << std::endl; 
    }
  }
  if( !same )
  {
         std::cout << simplex[0] << " " << simplex[numRow] << " " << simplex[2*numRow] << std::endl; 
         std::cout << simplex[1] << " " << simplex[1+numRow] << " " << simplex[1+ 2*numRow] << std::endl; 
     
  }
  return(same);
}

void 
NonGradient::getRowSum( int numRow, int numCol, std::vector<double>& matrix, std::vector<double>& rowSum)
{
  for (int row=0;row<numRow;row++)
  {
      rowSum[row] = 0.;
  }
  for (int col=0;col<numCol;col++)    
  {
    for (int row=0;row<numRow;row++)
    {
      rowSum[row] += matrix[row+col*numRow];
    }
  }
}

double 
NonGradient::evaluate( double *point,  PatchData &pd, MsqError &err )
{
  return 0.;
  double value;
  // put point in pd. dmd  where 
  assert( pd.num_free_vertices() == 1 );
  size_t vertex_index = 0;
  Vector3D origVec = pd.vertex_by_index(vertex_index);  //[dim]
  Vector3D delta;
  for( int dim = 0; dim<3; dim++)
  {
    delta[dim] = point[dim];
  }
  pd.move_vertex( delta, vertex_index, err );
  OFEvaluator& obj_func = get_objective_function_evaluator();
  //bool feasible = obj_func.evaluate( ObjectiveFunction::CALCULATE, pd, value, OF_FREE_EVALS_ONLY, err ); MSQ_ERRZERO(err);
  bool feasible = obj_func.evaluate( pd, value, err ); MSQ_ERRZERO(err);
  pd.set_vertex_coordinates( origVec, vertex_index, err );
  if( feasible ) 
  {
    return value;
  }
  else
  {
    return -value; // or infinity or zero
  }
}

// Extrapolate by a factor of fac through the face of the simplex
// opposite the high point to a new point.  If the new point is 
// better, swap it with the high point.
double
NonGradient::amotry( std::vector<double>& simplex, 
                 std::vector<double>& height,
                 double psum[], int ihi, double fac, PatchData &pd, MsqError &err)
{
  int numRow = getDimension();
  int numCol = numRow + 1;
  std::vector<double> ptry(numRow); // does this make sense?
  double fac1=(1.0-fac)/static_cast<double>(numRow);
  double fac2=fac1-fac;
  for (int row=0;row<numRow;row++)
  {
    ptry[row]=psum[row]*fac1-simplex[row+ihi*numRow]*fac2;
  }
  std::cout << "Try ";
  double ytry = evaluate(&ptry[0], pd, err); // value at trial point
  if (ytry < height[ihi]) // better than highest (worst)
  {
    height[ihi]=ytry;     // swap ihi and ytry
    for (int row=0;row<numRow;row++)
    {
      psum[row] += (ptry[row]-simplex[row+ihi*numRow]);
      simplex[row+ihi*numRow]=ptry[row];
    }
  }
  return ytry;
}
  
void NonGradient::initialize(PatchData &pd, MsqError &err)
{
  int dimension = 3*pd.num_free_vertices();
  setDimension(dimension);
  int maxNumEval = 1;     // 200,  depends on dimension
  setMaxNumEval(maxNumEval);
  double threshold = 1.e-10; // avoid division by zero
  setThreshold(threshold);
  double minEdgeLen;
  double maxEdgeLen;
  pd.get_minmax_edge_length( minEdgeLen, maxEdgeLen );
  double ftol = minEdgeLen * .2;       // 1.e-4, edge length
  setTolerance(ftol);

  int numRow = dimension;
  int numCol = numRow+1;  
  if( numRow*numCol <= simplex.max_size() )
  { 
    simplex.resize(numRow*numCol); 
    double scale = minEdgeLen;
    for( int col = 1; col < numCol; col++ )
    {
      int row = col-1;
      simplex[ row + col*numRow ] = scale * static_cast<double>(col);
    }
  }
  else
  {
    std::cout<< "memory allocation impossible" << std::endl;
  //  MSQ_ERRRTN(1);  MSQ_ERRRTN(err);  err is a structure or class.
  }
}

void NonGradient::initialize_mesh_iteration(PatchData &/*pd*/, MsqError &/*err*/)
{
}

void NonGradient::optimize_vertex_positions(PatchData &pd, 
                                            MsqError &err)
{
  MSQ_FUNCTION_TIMER( "NonGradient::optimize_vertex_positions" );
  int numRow = getDimension();
  int numCol = numRow+1;  
  std::vector<double> height(numCol); 
  for(int col = 0; col < numCol; col++)
  {
    height[col] =  evaluate(&simplex[col*numRow], pd, err);//  eval patch stuff
  }
  int maxNumEval = getMaxNumEval();
  double threshold = getThreshold();
  double ftol = getTolerance();
  int ilo,inhi,ihi;//height[ilo]<=...<=height[inhi]<=height[ihi] 
  double rtol = 2.*ftol;
  double ysave;
  double ytry;
  std::vector<double> rowSum(numRow);
  getRowSum( numRow, numCol, simplex, rowSum);
  int numEval=0;
  while( rtol >= ftol && numEval < maxNumEval)
  {
    if( numEval > 0 )   // skip before evaluations
    {
      // Reflect highPt through opposite face
      if( height[0] == 0.)
      {
         std::cout << "iterate:reflect pre amotry Zero: fatal error" << std::endl;
         exit(-1);
      }
      if( !testRowSum( numRow, numCol, &simplex[0], &rowSum[0]) )
      {
        std::cout << "prev amo rowsum test failed" << std::endl;
        //MSQ_ERRRTN(-1);
      }
      ytry=amotry(simplex,height,&rowSum[0],ihi,-1.0, pd, err);
      if( !testRowSum( numRow, numCol, &simplex[0], &rowSum[0]) )
      {
         std::cout << "post amo rowsum test failed" << std::endl;
         std::cout << "ytry = " << ytry << std::endl;
         //MSQ_ERRRTN(-2);
      }   
  
      std::cout << " Reflect = "  << ytry << std::endl;
      if( height[0] == 0.)
      {
         std::cout << "iterate:reflect Zero: fatal error" << std::endl;
         exit(-1);
      }
      if (ytry <= height[ilo])   
      {
        // reflect and expand from highPt 
        ytry=amotry(simplex,height,&rowSum[0],ihi,-2.0,pd,err);
        std::cout << " Expand = "  << ytry << std::endl;
      }
      else 
      {
        if (ytry >= height[inhi]) 
        {
          ysave=height[ihi]; // Contract along highPt
          ytry=amotry(simplex,height,&rowSum[0],ihi,0.5,pd,err);
          std::cout << " Contract = "  << ytry << std::endl;
          if (ytry >= ysave)
          { // contract all directions toward lowPt
            for (int col=0;col<numCol;col++)
            {
              if (col != ilo)
              {
                for (int row=0;row<numRow;row++)
                {
                  rowSum[row]=0.5*(simplex[row+col*numRow]+simplex[row+ilo*numRow]);
                  simplex[row+col*numRow]=rowSum[row];
                }
                height[col] = evaluate(&rowSum[0], pd, err); 
                std::cout << "Iter: height[ " << col << "=  "<< height[col] << std::endl;
   
              }
              else
              {
                std::cout << "Iter: height[ " << col << "=  "<< height[col] << std::endl;
              }
            }
            numEval += numRow;
            getRowSum( numRow, numCol, simplex, rowSum);
          }
        }   
        else 
        {
          --numEval;
        }
      } // ytri > h(ilo) 
    } // numEval > 0,  skip evaluation
    ilo=1; // conditional operator or inline if ?:,
    ihi = height[0] > height[1] ? (inhi=1,0) : (inhi=0,1);
    // i.e. (inhi=1,ihi=0)
    for (int col=0;col<numCol;col++)
    {
      if (height[col] <= height[ilo])
      {
        ilo=col;  // ilo := argmin height
      }
      if (height[col] > height[ihi])
      {
        inhi=ihi;
        ihi=col;
      } 
      else  // height[ihi] >= height[col]
        if (col != ihi && height[col] > height[inhi] ) inhi=col;
    }
    rtol=2.0*fabs( height[ihi]-height[ilo] )/
         ( fabs(height[ihi])+fabs(height[ilo])+threshold );
    std::cout << "rtol : " << rtol << std::endl;
    //std::cout << "height: " << ilo << "  " << inhi << "  " << ihi  << std::endl;
    //std::cout << height[ilo] << "  " << height[inhi] << "  " << height[ihi]  << std::endl;
    numEval += 2;
  } //  while( rtol >= ftol && numEval < maxNumEval)
  if (rtol < ftol)
  { 
    if( ilo != 0 )
    {
      double yTemp = height[0];
      height[0] = height[ilo]; // height dimension numCol
      height[ilo] = yTemp;
      for (int row=1;row<numRow;row++)
      { 
          yTemp = simplex[row];
          simplex[row] = simplex[row+ilo*numRow];
          simplex[row+ilo*numRow] = yTemp;
      }
    }
  }
  // back to NonGradient::minimize
  for(int row = 0; row < numRow; row++)
  {
    //solutionPoint[row] = simplex[row];
    //dmd: todo:  set free vertices to this value 
    // no way to return the solution.
  }
  int returnFlag = (numEval >= maxNumEval ? 1:0);
  // minimize diverged
  //  return( returnFlag );  
  // dmd alert user of stagnation
}

void NonGradient::terminate_mesh_iteration(PatchData &/*pd*/, MsqError &/*err*/)
{
    std::cout << "- Executing NonGradient::iteration_complete()\n";
}
  
void NonGradient::cleanup()
{
    std::cout << "- Executing NonGradient::iteration_end()\n";
}
  
} // namespace Mesquite
