#include "interF.h"
#include <iostream>

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> shiftMat(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matin,
int value);

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>  interF(
int npts,double alpha,double beta,double gamma)
{
// make Penta diagonal matrix b, one row:

Eigen::Matrix<double, 5, 1> b;
b(0)=beta;
b(1)=-(alpha + 4*beta);
b(2)=(2*alpha + 6 *beta);
b(3)=b(1);
b(4)=b(0);

// Make the penta matrix (for every contour point)

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> B;
A.resize(npts,npts);
B.setIdentity(npts,npts);
//Eigen::Matrix<double,npts,npts> A;
//Eigen::Matrix<double,npts,npts> B;

A=b(0)*shiftMat(Eigen::MatrixXd::Identity(npts,npts),2);
A=A+b(1)*shiftMat(Eigen::MatrixXd::Identity(npts,npts),1);
A=A+b(2)*shiftMat(Eigen::MatrixXd::Identity(npts,npts),0);
A=A+b(3)*shiftMat(Eigen::MatrixXd::Identity(npts,npts),-1);
A=A+b(4)*shiftMat(Eigen::MatrixXd::Identity(npts,npts),-2);

// Calculate the inverse
//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> B;
//A.resize(npts,npts);


//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> C;

//C.resize(npts,npts);
//C = A + gamma* B;

//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> D;

//D.resize(npts,npts);
//D = C.inverse();
//return D;//.inverse();
return (A + gamma* B).inverse();


}

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> shiftMat(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matin,
int value){
int row = matin.rows();
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matout;
matout.resize(row,row);
//if (value<0)
//value = value+row;
int j;
for (int i = 0; i<row;i++){
j = i-value;
if(i-value<0)
j = i-value+row;
if(i-value>=row)
j = i-value-row;

matout.row(i) = matin.row(j);
}
return matout;
}


