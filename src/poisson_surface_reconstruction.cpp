#include "poisson_surface_reconstruction.h"
#include"fd_grad.h"
#include"fd_interpolate.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include<iostream>

void poisson_surface_reconstruction(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & N,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F)
{
  ////////////////////////////////////////////////////////////////////////////
  // Construct FD grid, CONGRATULATIONS! You get this for free!
  ////////////////////////////////////////////////////////////////////////////
  // number of input points
  const int n = P.rows();
  // Grid dimensions
  int nx, ny, nz;
  // Maximum extent (side length of bounding box) of points
  double max_extent =
    (P.colwise().maxCoeff()-P.colwise().minCoeff()).maxCoeff();
  // padding: number of cells beyond bounding box of input points
  const double pad = 8;
  // choose grid spacing (h) so that shortest side gets 30+2*pad samples
  double h  = max_extent/double(30+2*pad);
  // Place bottom-left-front corner of grid at minimum of points minus padding
  Eigen::RowVector3d corner = P.colwise().minCoeff().array()-pad*h;
  // Grid dimensions should be at least 3 
  nx = std::max((P.col(0).maxCoeff()-P.col(0).minCoeff()+(2.*pad)*h)/h,3.);
  ny = std::max((P.col(1).maxCoeff()-P.col(1).minCoeff()+(2.*pad)*h)/h,3.);
  nz = std::max((P.col(2).maxCoeff()-P.col(2).minCoeff()+(2.*pad)*h)/h,3.);
  // Compute positions of grid nodes
  Eigen::MatrixXd x(nx*ny*nz, 3);

      for(int k = 0; k < nz; k++)
      {
    for(int j = 0; j < ny; j++)
    {
    for(int i = 0; i < nx; i++) 
  {
         // Convert subscript to index
         const auto ind = i + nx*(j + k * ny);
         x.row(ind) = corner + h*Eigen::RowVector3d(i,j,k);
      }
    }
  }
  Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);

  ////////////////////////////////////////////////////////////////////////////
   

  Eigen::SparseMatrix<double,Eigen::ColMajor>  W;// n * nx,ny,nz
  fd_interpolate(nx,ny,nz,h,corner,P,W);
   
  Eigen::SparseMatrix<double,Eigen::RowMajor>  G;// XXX * nx,ny,nz
  fd_grad(nx,ny,nz,h,G);
 
 int rows_x = (nx-1)*ny*nz;
 int rows_y = (ny-1)*nx*nz;
 int rows_z = (nz-1)*ny*nx;
 Eigen::VectorXd v  =  Eigen::VectorXd::Zero(rows_x+rows_y+rows_z);
  
   Eigen::VectorXd vx  = Eigen::VectorXd::Zero(rows_x);
   Eigen::VectorXd vy  = Eigen::VectorXd::Zero(rows_y);
   Eigen::VectorXd vz  = Eigen::VectorXd::Zero(rows_z);
 int count_x = -1; int count_y = -1; int count_z = -1;


      for (size_t k = 0; k < nz; k++)
      {
    
    for (size_t j = 0; j < ny; j++)
    {
      for (size_t i = 0; i < nx; i++)
  {

        if(i<nx-1){
        int index_x = i+j*nx+k*nx*ny;
        count_x++;
        vx(count_x) = W.col(index_x).transpose()*N.col(0);
        }

        if(j<ny-1){
        int index_y = i+j*nx+k*nx*ny;
          count_y++;
           vy(count_y) = W.col(index_y).transpose()*N.col(1);
        }

       if(k<nz-1){
        int index_z = i+j*nx+k*nx*ny;
          count_z++;
           vz(count_z) = W.col(index_z).transpose()*N.col(2);
        }
     
      }
       
    }
   
  }
   
  v.head(rows_x)= vx;
  v.segment(rows_x,rows_y) = vy;
  v.tail(rows_z) = vz;
 
  Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b = Eigen::VectorXd::Zero(nx*ny*nz);
 
      A = G.transpose()*G;
      b = G.transpose()*v;
       
      Eigen::BiCGSTAB<  Eigen::SparseMatrix<double> > solver;
       
      solver.compute(A);
       
     g  = solver.solve(b);

     
  //  G^T*G*g = G^T*v;
 



//void fd_grad(
  // const int nx,
  // const int ny,
  // const int nz,
  // const double h,
  // Eigen::SparseMatrix<double,Eigen::RowMajor> & G)
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
