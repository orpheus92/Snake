#include "snake2D.h"
#include <stb_image.h>
#include <iostream>

//snake movement
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> snakeMove(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> intForce,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fx,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fy,
double gamma,
double kappa,
double delta);
//baloon Force
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> baloonF(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P);
//bilinear interp
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> interp2(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> mypts,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P);

//main function for snake
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
double sigma3, //sigma used to calculate laplacian in GVF, default 1
igl::viewer::Viewer& viewer
)
{
// make clockwise contour  (always clockwise due to balloon force)
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> O;
O.resize(2+P.rows(),2);

O << P,
     P.row(0),
     P.row(1);
double area;

// area inside the contour 
area = 0.5*(O.block(1,0,P.rows(),1).cwiseProduct(O.block(2,1,P.rows(),1) - O.block(0,1,P.rows(),1))).sum();

//If the area inside  the contour is positive, change from counter-clockwise to  clockwise
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Pcont;
if (area>0)
{
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Pre;
Pre = P.colwise().reverse();
Pcont = interpcont(Pre, npts);
}
else
{
Pcont = interpcont(P, npts);
}

//Linear interpolation for initial contour
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Eext;

//Calculate external force
//size = image size

//This is correct
Eext = extForce(input,wl,we,wt,sigma);

//Make the external flow field
//size = image size
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fx;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fy;
Fx = -imDev(Eext,sigma2,1)*2*sigma2*sigma2;
Fy = -imDev(Eext,sigma2,2)*2*sigma2*sigma2;

//Calcuate GVF Image Force  Might be needed later 
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fx2;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fy2;
//will be modified
Fx2=GVFimF(Fx,Fy, mu, Giter, sigma3, 1);
Fy2=GVFimF(Fx,Fy, mu, Giter, sigma3, 2);

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> intForce;

//internal Force for snake
//size = npts by npts
//This is correct
intForce = interF(npts,alpha,beta,gamma);
//std::cout<<"internal force = "<<intForce<<std::endl;
Eigen::MatrixXd Hout;
Hout.resize(1,2);
Hout << 0,0;

// Triangulated interior
Eigen::MatrixXd V2;
Eigen::MatrixXi F2;

//igl::viewer::Viewer viewer;
//viewer.launch();
//viewer.launch();
for (int it = 0;it<iter;it++){

Pcont = snakeMove(intForce,Pcont,Fx2,Fy2,gamma,kappa,delta);

//Implement Draw Mesh Here
/*
double rsize =(double)input.rows();
double csize =(double)input.cols();
Eigen::MatrixXd Vout;
Eigen::MatrixXi Eout;

Vout.resize(Pcont.rows()+P.rows()+4,2);
Eout.resize(Pcont.rows()+P.rows()+4,2);
for (int i = 0; i<Pcont.rows();i++){
Vout(i,0) = Pcont(i,0)-csize/2;
Vout(i,1) = Pcont(i,1)-rsize/2;
Eout(i,0) = i;
Eout(i,1) = i+1;
}
Eout(Pcont.rows()-1,1) = 0;
for (int j = 0; j<P.rows();j++){
Vout(j+Pcont.rows(),0) = (double)(P(j,0)-csize/2);
Vout(j+Pcont.rows(),1) = (double)(P(j,1)-rsize/2);
Eout(Pcont.rows()+j,0) = Pcont.rows()+j;
Eout(j+Pcont.rows(),1) = Pcont.rows()+j+1;
}

Eout(Pcont.rows()+P.rows()-1,1) = Eout(Pcont.rows(),0);

//manually set boundary
Vout(Pcont.rows()+P.rows(),0) = -csize/2;
Vout(Pcont.rows()+P.rows(),1) = -rsize/2;
Vout(Pcont.rows()+P.rows()+1,0) = csize/2;
Vout(Pcont.rows()+P.rows()+1,1) = -rsize/2;
Vout(Pcont.rows()+P.rows()+2,0) = csize/2;
Vout(Pcont.rows()+P.rows()+2,1) = rsize/2;
Vout(Pcont.rows()+P.rows()+3,0) = -csize/2;
Vout(Pcont.rows()+P.rows()+3,1) = rsize/2;


Eout(Pcont.rows()+P.rows(),0) = P.rows()+Pcont.rows();
Eout(Pcont.rows()+P.rows(),1) = P.rows()+Pcont.rows()+1;
Eout(Pcont.rows()+P.rows()+1,0) = P.rows()+Pcont.rows()+1;
Eout(Pcont.rows()+P.rows()+1,1) = P.rows()+Pcont.rows()+2;
Eout(Pcont.rows()+P.rows()+2,0) = P.rows()+Pcont.rows()+2;
Eout(Pcont.rows()+P.rows()+2,1) = P.rows()+Pcont.rows()+3;
Eout(Pcont.rows()+P.rows()+3,0) = P.rows()+Pcont.rows()+3;
Eout(Pcont.rows()+P.rows()+3,1) = P.rows()+Pcont.rows();

igl::triangle::triangulate(Vout,Eout,Hout,"a1q",V2,F2);

    viewer.data.clear();
    viewer.data.set_mesh(V2,F2);
    viewer.core.align_camera_center(V2,F2);
    viewer.core.show_texture = false;
   
 viewer.core.animation_max_fps = 30.;

    double tic = get_seconds();
    glfwPollEvents();
        // In microseconds
        double duration = 1000000.*(get_seconds()-tic);
        const double min_duration = 1000000./30;
        if(duration<min_duration)
        {
          std::this_thread::sleep_for(std::chrono::microseconds((int)(min_duration-duration)));
        }


	//for (int xx = 0; xx<1000000;xx++){
	//std::cout<<xx<<std::endl;
	//}
*/
}
double rsize =(double)input.rows();
double csize =(double)input.cols();
Eigen::MatrixXd Vout;
Eigen::MatrixXi Eout;

Vout.resize(Pcont.rows()+P.rows()+4,2);
Eout.resize(Pcont.rows()+P.rows()+4,2);
for (int i = 0; i<Pcont.rows();i++){
Vout(i,0) = Pcont(i,0)-csize/2;
Vout(i,1) = Pcont(i,1)-rsize/2;
Eout(i,0) = i;
Eout(i,1) = i+1;
}
Eout(Pcont.rows()-1,1) = 0;
for (int j = 0; j<P.rows();j++){
Vout(j+Pcont.rows(),0) = (double)(P(j,0)-csize/2);
Vout(j+Pcont.rows(),1) = (double)(P(j,1)-rsize/2);
Eout(Pcont.rows()+j,0) = Pcont.rows()+j;
Eout(j+Pcont.rows(),1) = Pcont.rows()+j+1;
}

Eout(Pcont.rows()+P.rows()-1,1) = Eout(Pcont.rows(),0);

//manually set boundary
Vout(Pcont.rows()+P.rows(),0) = -csize/2;
Vout(Pcont.rows()+P.rows(),1) = -rsize/2;
Vout(Pcont.rows()+P.rows()+1,0) = csize/2;
Vout(Pcont.rows()+P.rows()+1,1) = -rsize/2;
Vout(Pcont.rows()+P.rows()+2,0) = csize/2;
Vout(Pcont.rows()+P.rows()+2,1) = rsize/2;
Vout(Pcont.rows()+P.rows()+3,0) = -csize/2;
Vout(Pcont.rows()+P.rows()+3,1) = rsize/2;


Eout(Pcont.rows()+P.rows(),0) = P.rows()+Pcont.rows();
Eout(Pcont.rows()+P.rows(),1) = P.rows()+Pcont.rows()+1;
Eout(Pcont.rows()+P.rows()+1,0) = P.rows()+Pcont.rows()+1;
Eout(Pcont.rows()+P.rows()+1,1) = P.rows()+Pcont.rows()+2;
Eout(Pcont.rows()+P.rows()+2,0) = P.rows()+Pcont.rows()+2;
Eout(Pcont.rows()+P.rows()+2,1) = P.rows()+Pcont.rows()+3;
Eout(Pcont.rows()+P.rows()+3,0) = P.rows()+Pcont.rows()+3;
Eout(Pcont.rows()+P.rows()+3,1) = P.rows()+Pcont.rows();

igl::triangle::triangulate(Vout,Eout,Hout,"a1q",V2,F2);

    viewer.data.clear();
    viewer.data.set_mesh(V2,F2);
    viewer.core.align_camera_center(V2,F2);
    viewer.core.show_texture = false;
return Pcont; 
}

//function to calculate baloon force
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> baloonF(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P){
//input n by 2
//output n by 2
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> out; 
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> dx; 
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> dy;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> l;
dx.resize(P.rows(),1);
dy.resize(P.rows(),1);
l.resize(P.rows(),1);

Eigen::VectorXd f;
Eigen::VectorXd b;

f.setLinSpaced(P.rows(),1,P.rows()); 
b.setLinSpaced(P.rows(),1,P.rows()); 
f.array() += 4;                
b.array() -= 4; 
for(int i = 0; i<P.rows();i++){

	if (f(i)>P.rows())
	f(i) = f(i)-P.rows();
	if (b(i)<1)
	b(i) = b(i)+P.rows();

	dx(i,0) = P(f(i)-1, 0)-P(b(i)-1, 0);
	dy(i,0) = P(f(i)-1, 1)-P(b(i)-1, 1);
	
}
l = (dx.array().square()+dy.array().square()).cwiseSqrt() ;            

out.resize(P.rows(),2);
out.col(0) = -dy.array() / l.array();  
out.col(1) = dx.array() / l.array();  
return out;
}


// function for linear interpolation
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> interp2(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> mypic,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P)
{
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> out;
out.resize(P.rows(),1);
double x;
double y;
int indx;
int indy;
double q11;
double q12;
double q21;
double q22;

double x1;
double x2;
double y1;
double y2;

//bilinear interpolation
for (int i = 0; i<P.rows();i++){

y = P(i,0);
x = P(i,1);

indx = (int)ceil(x);
indy = (int)ceil(y);

x1 = (double)(indx-1);
x2 = (double)(indx);
y1 = (double)(indy-1);
y2 = (double)(indy);

	if ((indx==0) & (indy==0))	
		out(i,0) = mypic(0,0);
	else if (indx == 0)
		out(i,0) = mypic(indy,0)*(y-y1) + mypic(indy-1,0)*(y2-y);
	else if (indy == 0)
		out(i,0) = mypic(0,indx)*(x-x1) + mypic(0,indx-1)*(x2-x);
	else{
		q22 = mypic(indy,indx);
		q12 = mypic(indy-1,indx);
		q21 = mypic(indy,indx-1);
		q11 = mypic(indy-1,indx-1);
		out(i,0) = q11*(x2-x)*(y2-y) + q21*(x-x1)*(y2-y) + q12*(x2-x)*(y-y1) + q22*(x-x1)*(y-y1);
	}

}

return out;
}

//This function will calculate one iteration of contour Snake movement
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> snakeMove(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> intForce, //internal force
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P,//contour pts
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fx,//external vector field
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fy,//external vec field
double gamma, //time step
double kappa,//external field weight
double delta) //Balloon Force Weight)
{
//Clamp contour to boundary
//rows and cols might be checked later 
//image force
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> imforce;
imforce.resize(P.rows(),2);

//baloon force
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> baforce;

//clamp boundary
P.col(0)=P.col(0).cwiseMax(0).cwiseMin(Fx.cols()-1);
P.col(1)=P.col(1).cwiseMax(0).cwiseMin(Fx.rows()-1);

//Get image force on the contour points
imforce.col(0)=kappa*interp2(Fx,P);
imforce.col(1)=kappa*interp2(Fy,P);

//Get baloon force on the contour points

baforce = delta * baloonF(P);

//Update contour positions

P.col(0) = intForce * (gamma*P.col(0) + imforce.col(0) + baforce.col(0));
P.col(1) = intForce * (gamma*P.col(1) + imforce.col(1) + baforce.col(1));

//Clamp contour to boundary
P.col(0)=P.col(0).cwiseMax(0).cwiseMin(Fx.cols()-1);
P.col(1)=P.col(1).cwiseMax(0).cwiseMin(Fx.rows()-1);


return P;
}


