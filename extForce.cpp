#include "extForce.h"
#include <iostream>
//function to calculate the external force

//function to filter an image with a designed filter 
void imfilter2(Eigen::MatrixXd src_image,Eigen::MatrixXd filter,Eigen::MatrixXd &dst_image);

//IMGAUSSIAN filters a 2D greyscale image with an Gaussian filter.
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> imgaussian(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> I, int Sigma);


Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> extForce(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Img,
double wl,
double we,
double wt,
int sigma
//int npts
){ int total = 2*3*sigma+1;
int row = Img.rows();
int col = Img.cols();
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> x;
x.resize(total,total);
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> y;
y.resize(total,total);
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gx;
Gx.resize(total,total);
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gy;
Gy.resize(total,total);
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gxx;
Gxx.resize(total,total);
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gyy;
Gyy.resize(total,total);
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gxy;
Gxy.resize(total,total);


//Eigen::Matrix<double, total,total> Gx;
//Eigen::Matrix<double, total,total> Gy;
//Eigen::Matrix<double, total,total> Gxx;
//Eigen::Matrix<double, total,total> Gyy;
//Eigen::Matrix<double, total,total> Gxy;





Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Ix;
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Iy;
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Ixx;
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Iyy;
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Ixy;

for (int i = 0; i< total;i++){
	for (int j = 0; j< total;j++){
		x(i,j) = i-sigma*3;
		y(i,j) = j-sigma*3;
//Gaussian Filters
	Gx(i,j)=-x(i,j)/(2*M_PI*pow(sigma,4))*exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));

	Gy(i,j)=-y(i,j)/(2*M_PI*pow(sigma,4))*exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));

        Gxx(i,j)=1/(2*M_PI*pow(sigma,4))*(pow(x(i,j),2)/pow(sigma,2)-1) *exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));

        Gyy(i,j)=1/(2*M_PI*pow(sigma,4))*(pow(y(i,j),2)/pow(sigma,2)-1) *exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));

        Gxy(i,j)=1/(2*M_PI*pow(sigma,6))*(x(i,j)*y(i,j)) *exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));

}//end for j
}//end for i
//Gaussian Based Image Derivatives


//static void 	cv::eigen2cv (const Eigen::Matrix< _Tp, _rows, _cols, _options, _maxRows, _maxCols > &src, Mat &dst)
//cv::eigen2cv (Img, Mat &dst);
imfilter2(Img, Gx, Ix);
imfilter2(Img, Gy, Iy);
imfilter2(Img, Gxx, Ixx);
imfilter2(Img, Gyy, Iyy);
imfilter2(Img, Gxy, Ixy);

//Line Energy
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic>Eline;
Eline = imgaussian(Img,sigma);

//std::cout<<Iyy.cwiseProduct(Ix2) <<std::endl;
//Line Termination Energy
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic>Eterm;
//R = P.cwiseProduct(Q);   

Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic>Ix2 = Ix.array().square();
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic>Iy2 = Iy.array().square();
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic>Ixyx = Ixy.cwiseProduct(Ix);
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic>nom = Iyy.cwiseProduct(Ix2) -2*Ixyx.cwiseProduct(Iy) + Ixx.cwiseProduct(Iy2);
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic>denom1 = Eigen::MatrixXd::Ones(Ix.rows(),Ix.cols())+Ix2 + Iy2;
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic>denom = denom1.array().pow(3/2);


Eterm = nom.cwiseQuotient(denom);

//Edge Energy

Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic>Eedge;
Eedge = (Ix.array().square() + Iy.array().square()).array().sqrt(); 

//Externa =l Energy
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic>Eext;

Eext = (wl*Eline.array() - we*Eedge.array() -wt * Eterm.array()); 

return Eext;
}



void imfilter2(Eigen::MatrixXd src_image, Eigen::MatrixXd filter,Eigen::MatrixXd &dst_image)
{
	dst_image = Eigen::MatrixXd::Zero(src_image.rows(),src_image.cols()) ;	
	
	int start_row = int(filter.rows()/2) ;
	int start_col = int(filter.cols()/2) ;

	Eigen::MatrixXd Mid_Matrix = Eigen::MatrixXd::Zero(src_image.rows()+2*start_row,src_image.cols()+2*start_col) ;

	for (int i = 0; i < src_image.rows(); i++)
	{
		for (int j = 0; j < src_image.cols(); j++)
		{
			Mid_Matrix(i+start_row,j+start_col) = src_image(i,j) ;
		}
	}

	int end_row = Mid_Matrix.rows() -1 - start_row ;
	int end_col = Mid_Matrix.cols() -1 - start_col ;

	int filter_row = filter.rows();
	int filter_col = filter.cols() ;
	
	for (int i = start_row; i <= end_row; i++)
	{
		for (int j = start_col; j <= end_col; j++)
		{			
			int tmp_row = i - start_row  ;
			int tmp_col = j - start_col  ;
			for (int m = 0; m < filter_row; m++)
			{				
				for (int n = 0; n < filter_col; n++)
				{
					dst_image(tmp_row,tmp_col) += Mid_Matrix(tmp_row + m,tmp_col + n)*filter(m,n) ; 
				}
			}
		}
	}
	return ;
}

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> imgaussian(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> I, int Sigma){
int siz = Sigma*6;

if (Sigma>0){
	// Make 1D Gaussian kernel

	Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(2*ceil(siz/2)+1,-ceil(siz/2),ceil(siz/2));     // low:step:hi
	
	Eigen::VectorXd H;
	H = (-(x.array().square()/(2*Sigma*Sigma))).array().exp();
	Eigen::VectorXd HH;

	    HH = H.array()/(H.sum());

	//Eigen::Matrix<double, HH.rows(),1> Hx;
	//Eigen::Matrix<double, 1,HH.rows()> Hy;
	//Hx = HH;
	//Hy = HH.transpose();
	Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Im2;
	Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> out;
	imfilter2(I, HH, Im2);
	imfilter2(Im2, HH.transpose(), out);

	//Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Iyy;
	//Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Ixy;
	return out;
	}
else
return I;

}




