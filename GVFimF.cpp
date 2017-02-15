#include "GVFimF.h"
#include <iostream>

//Function to calculate Laplacian and update vector field
/*
void Gupdate(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fx,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fy,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &u, 
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &v,
double sigma,
double mu,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> sMag);
*/
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>  GVFimF(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fx,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fy,
double mu,
int Giter,
double sigma3,
int style){
//std::cout<<"in GVF"<<std::endl;
//calculate magnitude,
//std::cout<<"Fx in = "<<Fx.rowwise().sum()<<std::endl;
//std::cout<<"Fy in = "<<Fy.rowwise().sum()<<std::endl;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> sMag;
sMag.resize(Fx.rows(),Fy.cols());
sMag = Fx.array().square()+Fy.array().square();

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> uu;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> vv;
uu.resize(Fx.rows(),Fy.cols());
vv.resize(Fx.rows(),Fy.cols());
uu = Fx;
vv = Fy;

//std::cout<<"=================================================================="<<std::endl;
//std::cout<<"sMag = "<<sMag.colwise().sum()<<std::endl;

//std::cout<<"=================================================================="<<std::endl;
//std::cout<<"vv = "<<vv<<std::endl;
//Function to calculate Laplacian and update vector field
for (int i = 0; i < Giter; i++){
//std::cout<<"vv before = "<<vv<<std::endl;
//Gupdate(Fx, Fy, uu, vv, sigma3, mu, sMag);
//std::cout<<"Fx before = "<<Fx.rowwise().sum()<<std::endl;
//std::cout<<"uu before = "<<uu.rowwise().sum()<<std::endl;
uu = uu + mu*(imDev(uu, sigma3, 3) + imDev(uu, sigma3, 4))- sMag.cwiseProduct(uu-Fx);
//std::cout<<"=================================================================="<<std::endl;
//std::cout<<"Fy before = "<<Fy.rowwise().sum()<<std::endl;
//std::cout<<"vv before = "<<vv.rowwise().sum()<<std::endl;
vv = vv + mu*(imDev(vv, sigma3, 3) + imDev(vv, sigma3, 4))- sMag.cwiseProduct(vv-Fy);
//std::cout<<"vv after= "<<vv<<std::endl;
}

//std::cout<<"Fx = "<<Fx<<std::endl;
 // Uxx=ImageDerivatives2D(u,Sigma,'xx');
 // Uyy=ImageDerivatives2D(u,Sigma,'yy');
  
 // Vxx=ImageDerivatives2D(v,Sigma,'xx');
 // Vyy=ImageDerivatives2D(v,Sigma,'yy');

 // u = u + Mu*(Uxx+Uyy) - sMag.*(u-Fx);
 // v = v + Mu*(Vxx+Vyy) - sMag.*(v-Fy);

if (style == 1)
return uu;
else 
return vv;
}

/*
void Gupdate(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fx,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fy,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &u, 
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &v,
double sigma,
double mu,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> sMag)
{
u = u + mu*(imDev(u, sigma, 3) + imDev(u, sigma, 4))- sMag.cwiseProduct(u-Fx);

v = v + mu*(imDev(v, sigma, 3) + imDev(v, sigma, 4))- sMag.cwiseProduct(v-Fy);

return;
}
*/
