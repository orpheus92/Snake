#include "imDev.h"
#include <iostream>

void imfilter2a(Eigen::MatrixXd src_image,Eigen::MatrixXd filter,Eigen::MatrixXd &dst_image);


Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> imDev(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Img,
int sigma,
int type
// type 1 = x; 2 = y; 3 = xx; 4 = yy; 5 = xy;

){
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> G;
int total = 2*3*sigma+1;

Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> x;
x.resize(total,total);
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> y;
y.resize(total,total);

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> out;
//Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gx;
//Gx.resize(total,total);
//Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gy;
//Gy.resize(total,total);
//Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gxx;
//Gxx.resize(total,total);
//Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gyy;
//Gyy.resize(total,total);
//Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gxy;
G.resize(total,total);

for (int i = 0; i< total;i++){
	for (int j = 0; j< total;j++){
		x(i,j) = i-sigma*3;
		y(i,j) = j-sigma*3;
//Gaussian Filters

switch (type)
	{case 1:
		{G(i,j)=-x(i,j)/(2*M_PI*pow(sigma,4))*exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));
		}
	case 2:
		{G(i,j)=-y(i,j)/(2*M_PI*pow(sigma,4))*exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));
		}
	case 3:
		{G(i,j)=1/(2*M_PI*pow(sigma,4))*(pow(x(i,j),2)/pow(sigma,2)-1) *exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));
		}
	case 4:
		{G(i,j)=1/(2*M_PI*pow(sigma,4))*(pow(y(i,j),2)/pow(sigma,2)-1) *exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));
		}
	case 5:
		{G(i,j)=1/(2*M_PI*pow(sigma,6))*(x(i,j)*y(i,j)) *exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));
		}

	}

}//end for j
}//end for i
imfilter2a(Img, G, out);

return out;

}


void imfilter2a(Eigen::MatrixXd src_image, Eigen::MatrixXd filter,Eigen::MatrixXd &dst_image)
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
