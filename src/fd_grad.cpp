#include "fd_grad.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double,Eigen::RowMajor> & G)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
Eigen::SparseMatrix<double,Eigen::RowMajor> G_x;
Eigen::SparseMatrix<double,Eigen::RowMajor> G_y;
Eigen::SparseMatrix<double,Eigen::RowMajor> G_z;
 
fd_partial_derivative(nx,ny,nz,h,0,G_x);
fd_partial_derivative(nx,ny,nz,h,1,G_y);
fd_partial_derivative(nx,ny,nz,h,2,G_z);

int rows_x = (nx-1)*ny*nz;
int rows_y = (ny-1)*nx*nz;
int rows_z = (nz-1)*ny*nx;

G.resize(rows_x+rows_y+rows_z,nx*ny*nz);

G.topRows(rows_x) = G_x;
G.middleRows(rows_x,rows_y) = G_y;
G.middleRows(rows_x+rows_y,rows_z) = G_z;

}
