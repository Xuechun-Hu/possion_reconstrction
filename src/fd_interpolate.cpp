#include "fd_interpolate.h"
 #include<iostream>

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double,Eigen::ColMajor> & W)
{

  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
    int cols = nx*ny*nz;
    int point_count = P.rows();
   W.resize(point_count,cols);
 
    for (size_t l = 0; l < point_count; l++)
    {
       
      Eigen::RowVector3d  pt = P.row(l);
      Eigen::RowVector3d  diff_pt= pt - corner;
      Eigen::RowVector3i   index_pt = (diff_pt/h).cast<int>();
       
      if(index_pt[0]<0 ){
        index_pt[0] =0;
      }  
    
      if(index_pt[1]<0 ){
        index_pt[1] =0;
      }  
    
      if(index_pt[2]<0 ){
        index_pt[2] =0;
      }  
    
     if(index_pt[0]> nx-2){
       index_pt[0] = nx-2;
     }
      
     if(index_pt[1]> ny-2){
       index_pt[1] = ny-2;
     }

      if(index_pt[2]> nz-2){
       index_pt[2] = nz-2;
     }
        
 
        
      Eigen::RowVector3d interp1 = (index_pt + Eigen::RowVector3i(0,0,0)).cast<double>();
      Eigen::RowVector3d interp2 = (index_pt + Eigen::RowVector3i(1,0,0)).cast<double>();
      Eigen::RowVector3d interp3 = (index_pt + Eigen::RowVector3i(0,1,0)).cast<double>();
      Eigen::RowVector3d interp4 = (index_pt + Eigen::RowVector3i(1,1,0)).cast<double>();
      Eigen::RowVector3d interp5 = (index_pt + Eigen::RowVector3i(0,0,1)).cast<double>();
      Eigen::RowVector3d interp6 = (index_pt + Eigen::RowVector3i(1,0,1)).cast<double>();
      Eigen::RowVector3d interp7 = (index_pt + Eigen::RowVector3i(0,1,1)).cast<double>();
      Eigen::RowVector3d interp8 = (index_pt + Eigen::RowVector3i(1,1,1)).cast<double>();
      Eigen::RowVector3d  pt_interp = corner + interp1*h;
  
      Eigen::RowVector3d txyz = (pt-pt_interp)/h;

      double weight1 = (1-txyz[0])*(1-txyz[1])*(1-txyz[2])  ; 
      double weight2 = (txyz[0])*(1-txyz[1])*(1-txyz[2])  ; 
      double weight3 = (1-txyz[0])*(txyz[1])*(1-txyz[2])  ; 
      double weight4 = (txyz[0])*(txyz[1])*(1-txyz[2])  ; 
      double weight5 = (1-txyz[0])*(1-txyz[1])*(txyz[2])  ; 
      double weight6 = (txyz[0])*(1-txyz[1])*(txyz[2])  ; 
      double weight7 = (1-txyz[0])*(txyz[1])*(txyz[2])  ; 
      double weight8 = (txyz[0])*(txyz[1])*(txyz[2])  ; 
       double weight_sum =weight1+weight2+weight3+weight4+weight5+weight6+weight7+weight8;

     int index1 =   interp1[0]+interp1[1]*nx+interp1[2]*ny*nx;
      int index2 =  interp2[0]+interp2[1]*nx+interp2[2]*ny*nx;
     int index3 =   interp3[0]+interp3[1]*nx+interp3[2]*ny*nx;
      int index4 =  interp4[0]+interp4[1]*nx+interp4[2]*ny*nx;
      int index5 =  interp5[0]+interp5[1]*nx+interp5[2]*ny*nx;
      int index6 =  interp6[0]+interp6[1]*nx+interp6[2]*ny*nx;
      int index7 =  interp7[0]+interp7[1]*nx+interp7[2]*ny*nx;
      int index8 =  interp8[0]+interp8[1]*nx+interp8[2]*ny*nx;
 
      W.insert(l,index1) = weight1/weight_sum;
      W.insert(l,index2) = weight2/weight_sum;
      W.insert(l,index3) = weight3/weight_sum;
      W.insert(l,index4) = weight4/weight_sum;
      W.insert(l,index5) = weight5/weight_sum;
      W.insert(l,index6) = weight6/weight_sum;
      W.insert(l,index7) = weight7/weight_sum;
      W.insert(l,index8) = weight8/weight_sum;
     

    }
  
    //corner, nx,ny,nz;  h  P 



  ////////////////////////////////////////////////////////////////////////////
}
