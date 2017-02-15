#include <Eigen/Core>
#include <Eigen/Dense>
#include <string>
#include <cmath>
#include "opencv2/highgui/highgui.hpp"
//#include <iostream>
//#include <stdio.h>
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"

//#include "interpcont.cpp"

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> interpcont(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P,//,
//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P2//,
int npts
);
