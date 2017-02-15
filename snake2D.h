#include <Eigen/Core>
#include <string>
#include <cmath>
#include <iostream>

#include "interpcont.h"
#include "extForce.h"
#include "imDev.h"
#include "interF.h"
#include "GVFimF.h"

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> snake2D(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P, //initial contour
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> input, //input gray-scale image
double gamma, //time step, default 1
int iter, //# of iteration, default 100
int npts, //# of pts to interpolate contours 
int sigma, // sigma to calculate img derivative default 10
double wl, //attraction to lines, < 0 to black line; > 0 to white line; default 0.04
double we, //attraction to edge, default 2.0
double wt, //attraction to end points, default 0.01
double sigma2, //sigma to calculate gradient of edge energy image (give image force), default 20
double alpha, //membrane energy, default 0.2
double beta, //thin plate energy, default 0.2
double delta, //baloon force, default 0.1
double kappa, //weight of external img force, default 2
// the following is used for GVF snake
double mu, //tradeoff between real edge vectors and noise vectors, default 0.2
int Giter, //GVF iteration, default 0
double sigma3 //sigma used to calculate laplacian in GVF, default 1
);

