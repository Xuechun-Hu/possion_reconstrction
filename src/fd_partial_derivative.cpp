#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double,Eigen::RowMajor> & D)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
  int row = 0;
  int col = nx*ny*nz;
    int count =-1;
 if(dir == 0 ){//x
  row = (nx-1)*ny*nz;
  D.resize(row,col);

        for (size_t k = 0; k < nz; k++){
               for (size_t j = 0; j < ny; j++){
                   for (size_t i = 1; i < nx; i++){
            int index_1 = i + j*nx+k*ny*nx;
            int index = (i-1)+j*nx+k*ny*nx;
            count++;
            D.insert(count,index_1) = h;
            D.insert(count,index) = -h;
        }
     }
  }
  
 }

 if (dir == 1 )
 {//y
   row = (ny-1)*nx*nz;
   D.resize(row,col);
 
     for (size_t k = 0; k < nz; k++){
     for (size_t j = 1; j < ny; j++){ 
      for (size_t i = 0; i < nx; i++)
      {
            int index_1 =  i+j*nx+k*nx*ny ;
            int index =  i+ (j-1)*nx + k*nx*ny ;
              count++;
            D.insert(count,index_1) = h;
            D.insert(count,index) = -h;
        }
     }
  }
 }
 
  if (dir == 2 )
 {//z
    row = (nz-1)*ny*nx;
    D.resize(row,col);
     for (size_t k = 1; k < nz; k++){
    
     for (size_t j = 0; j < ny; j++){
       for (size_t i = 0; i < nx; i++)
  {
            int index_1 = i +j*nx+ k*ny*nx;
            int index = i +j*nx+ (k-1)*ny*nx;
            count++;
            D.insert(count,index_1) = h;
            D.insert(count,index) = -h;
        }
     }
  }
 }

 


}
