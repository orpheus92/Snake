#define _USE_MATH_DEFINES

#include <Eigen/Core>
#include <string>
#include <cmath> 

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>  extForce(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Img,
double wl,
double we,
double wt,
int sigma
//int npts
);
