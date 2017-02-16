#include <iostream>
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include <igl/triangle/triangulate.h>
#include <igl/png/writePNG.h>
#include <igl/png/readPNG.h>
//#include <igl/png/myread.h>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "interpcont.h"
#include "extForce.h"
#include "imDev.h"
#include "snake2D.h"

// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi E;
Eigen::MatrixXd H;


// Input contour Vertices
Eigen::MatrixXd xy;
Eigen::MatrixXi xyedge;
//Eigen::MatrixXd y;

// Triangulated interior
Eigen::MatrixXd V2;
Eigen::MatrixXi F2;

Eigen::MatrixXd V3;
Eigen::MatrixXi F3;

// Function to press key; called every time a keyboard button is pressed 
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
  {
    // Allocate temporary buffers
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(1280,800);
   // Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(1280,800);

    std::cout<<"first if "<<std::endl;
    // Draw the scene in the buffers
    viewer.core.draw_buffer(viewer.data,viewer.opengl,false,R,G,B,A);

    // Save it to a PNG
    igl::png::writePNG(R,G,B,A,"out.png");
  }

  if (key == '2')
  {
    // Allocate temporary buffers
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A,I;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> temp;
    // Read the PNG
    temp = igl::png::readPNG("myim2.png",R,G,B,A,temp);
//igl::png::readPNG("newim.png",R,G,B,A);
   // std::cout<<"R = "<<R<<std::endl;
    // Replace the mesh with a triangulated square
//std::cout<<"temp"<<temp<<std::endl;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P;
	P.resize(7,2);
	P <<
	38,25,
	20,59,
	39,97,
	81,105,
	109,84,
	112,39,
	73,24;

//output contour
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> outcont;	
std::cout<<"before snake"<<std::endl;
outcont = snake2D(P, //initial contour
temp, //input gray-scale image
1, //time step, default 1
400, //# of iteration, default 100
100, //# of pts to interpolate contours 
3, // sigma to calculate img derivative default 10
0, //attraction to lines, < 0 to black line; > 0 to white line; default 0.04 wl
2, //attraction to edge, default 2.0 wedge
0, //attraction to end points, default 0.01 wterm
3, //sigma to calculate gradient of edge energy image (give image force), default 20
0.1, //membrane energy, default 0.2 alpha
0.1, //thin plate energy, default 0.2 beta
-0.1, //baloon force, default 0.1 delta
4, //weight of external img force, default 2 kappa 
// the following is used for GVF snake
0.2, //tradeoff between real edge vectors and noise vectors, default 0.2
100, //GVF iteration, default 0
1 //sigma used to calculate laplacian in GVF, default 1
);
//std::cout<<"after snake"<<std::endl;

std::cout<<"out contour = "<<outcont<<std::endl;

//Eigen::MatrixXd Vout;

   /* Eigen::MatrixXd VV(4,3);
    VV <<
      -0.5,-0.5,0,
       0.5,-0.5,0,
       0.5, 0.5,0,
      -0.5, 0.5,0;
    Eigen::MatrixXi FF(2,3);
    FF <<
      0,1,2,
      2,3,0;
//std::cout<<FF<<std::endl;
    Eigen::MatrixXd UV(4,2);
    UV <<
      0,0,
      1,0,
      1,1,
      0,1;
*/
	/*Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P;
	P.resize(6,6);
	P <<
	1,5,6,4,1,12,
	2,10,4,8,12,12,
	3,13,32,1,7,11,
	5,16,21,1,5,2,
	8,20,12,5,1,8,
	10,24,1,4,11,6;

	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Pp;
	Pp.resize(6,2);
	Pp <<
	1,5,
	2,10,
	3,13,
	5,16,
	8,20,
	10,6;
*/


//	Eigen::VectorXd v;
//v.setLinSpaced(100,1,100); 
//std::cout<<"v = "<<v<<std::endl;
//std::cout<<"power2 = "<<P.array().square()<<"negative = "<<-P<<std::endl;
	//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Pout0, Pout;
//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Pout;
//Pout = interpcont(Pp,20);

//std::cout<<"Pout = "<<Pout<<std::endl;
//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Eext;
//Eext = extForce(P,0.04,2,0.01,10);
//Eext = P.inverse();//imDev(P,10,1);
//std::cout<<"min = "<<P.cwiseMin(18).cwiseMax(5)  <<std::endl;

//std::cout<<"r = "<<Pout.rows()<<"c = "<<Pout.cols()<<"Pout"<<Pout<<std::endl;


  //   Eigen::MatrixXd HH(1,2);

// H << 0,1;
  //  igl::triangle::triangulate(VV,FF,HH,"a0.005q",V3,F3);
  // Plot the generated mesh
  //igl::viewer::Viewer viewer;
  
//viewer.data.set_mesh(V3,F3);
//viewer.launch();
Eigen::MatrixXd Vout;
Eigen::MatrixXi Eout;
Eigen::MatrixXd Hout;
Vout.resize(outcont.rows()+P.rows()+4,2);
Eout.resize(outcont.rows()+P.rows()+4,2);
for (int i = 0; i<outcont.rows();i++){
Vout(i,0) = outcont(i,0)-60;
Vout(i,1) = outcont(i,1)-60;
Eout(i,0) = i;
Eout(i,1) = i+1;
}
Eout(outcont.rows()-1,1) = 0;
for (int j = 0; j<P.rows();j++){
Vout(j+outcont.rows(),0) = (double)(P(j,0)-60);
Vout(j+outcont.rows(),1) = (double)(P(j,1)-60);
Eout(outcont.rows()+j,0) = outcont.rows()+j;
Eout(j+outcont.rows(),1) = outcont.rows()+j+1;
}

Eout(outcont.rows()+P.rows()-1,1) = Eout(outcont.rows(),0);

//manually set boundary
Vout(outcont.rows()+P.rows(),0) = -60;
Vout(outcont.rows()+P.rows(),1) = -60;
Vout(outcont.rows()+P.rows()+1,0) = 60;
Vout(outcont.rows()+P.rows()+1,1) = -60;
Vout(outcont.rows()+P.rows()+2,0) = 60;
Vout(outcont.rows()+P.rows()+2,1) = 60;
Vout(outcont.rows()+P.rows()+3,0) = -60;
Vout(outcont.rows()+P.rows()+3,1) = 60;


Eout(outcont.rows()+P.rows(),0) = 107;
Eout(outcont.rows()+P.rows(),1) = 108;
Eout(outcont.rows()+P.rows()+1,0) = 108;
Eout(outcont.rows()+P.rows()+1,1) = 109;
Eout(outcont.rows()+P.rows()+2,0) = 109;
Eout(outcont.rows()+P.rows()+2,1) = 110;
Eout(outcont.rows()+P.rows()+3,0) = 110;
Eout(outcont.rows()+P.rows()+3,1) = 107;


Hout.resize(1,2);
Hout << 0,0;
igl::triangle::triangulate(Vout,Eout,Hout,"a1q",V2,F2);
//viewer.data.set_mesh(V2,F2);
std::cout<<"Vout = "<< Vout<<std::endl;
    viewer.data.clear();
    viewer.data.set_mesh(V2,F2);
  //  viewer.data.set_uv(UV);
  //  viewer.core.align_camera_center(VV);
    viewer.core.show_texture = true;

    // Use the image as a texture
//std::cout<< temp <<std::endl;
  //  viewer.data.set_texture(R,G,B);

// std::cout<<"G = "<<&G<<std::endl;
// std::cout<<"B = "<<B<<std::endl;
// std::cout<<"R = "<<R<<std::endl;
//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> aaaa;
//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &bbbb;
//aaaa.resize(5,2);
//	aaaa <<
//	96,63,
//	51,147,
//	98,242,
//	280,97,
//	182,59;
//bbbb = aaaa;
//std::cout<<"aaaa = "<< aaaa <<std::endl;
//std::cout<<"bbbb = "<< bbbb <<std::endl;
//aaaa.array() += 1; 
//std::cout<<"aaaa1 = "<< aaaa <<std::endl;
//std::cout<<"bbbb1 = "<< bbbb <<std::endl;
//*aaaa.array() += 1;
//std::cout<<"aaaa2 = "<< aaaa <<std::endl;
//std::cout<<"bbbb2 = "<< bbbb <<std::endl;
  }


  return false;
}




int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  // Create the boundary of a square
  V.resize(12,2);
  E.resize(12,2);
  H.resize(1,2);

xy.resize(10,2);
xyedge.resize(10,2);

 V << 100,100, 200,100, 200,200, 100, 200,
       50,50, 250,50, 250,250, 50, 250,
       0,0, 300,0, 300,300, 0, 300;

/*  V << -1,-1, 1,-1, 1,1, -1, 1,
       -2,-2, 2,-2, 2,2, -2, 2,
       -3,-3, 3,-3, 3,3, -3, 3;
*/
//V = V-150;
  E << 0,1, 1,2, 2,3, 3,0,
       4,5, 5,6, 6,7, 7,4,
       8,9, 9,10, 10,11, 11,8;

  H << 0,1;

  // Triangulate the interior
 // igl::triangle::triangulate(V,E,H,"a0.005q",V2,F2);
igl::triangle::triangulate(V,E,H,"a5q",V2,F2);
  // Plot the generated mesh
  igl::viewer::Viewer viewer;

 // std::cout<<"V2 = "<<V2<<"F2 = "<<F2<<endl;
// Wait for Key? 
  viewer.callback_key_down = &key_down;

  viewer.data.set_mesh(V2,F2);
  viewer.launch();
}
