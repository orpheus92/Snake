#include <Eigen/Core>
#include <Eigen/Dense>
#include <string>
#include <cmath> 

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>  imDev(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Img,
int sigma,
int type
);
