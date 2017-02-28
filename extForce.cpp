#include "extForce.h"
#include <iostream>
//function to calculate the external force

//function to filter an image with a designed filter 
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatrixXd;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> VectorXd;

//Imfilter function that uses correlation
MatrixXd compute2(MatrixXd X, MatrixXd W, bool flip, int stride_x, int stride_y);
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
x.resize(total,1);
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> y;
y.resize(1,total);

Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gx;
Gx.resize(total,1);
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gy;
Gy.resize(1,total);

for (int i = 0; i< total;i++){

		x(i,0) = i-sigma*3;
		y(0,i) = i-sigma*3;
}

Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Ix;
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Iy;
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Ixx;
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Iyy;
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Ixy;

for (int type =0;type<5; type++ )
{if (type == 0){
Gx = -(-x.array().square()/2/sigma/sigma).array().exp()*(x.array())/(2*M_PI*pow(sigma,4));
Gy = (-y.array().square()/2/sigma/sigma).array().exp();
Ix = compute2(Img, Gx, false, 1, 1);
Ix = compute2(Ix, Gy, false, 1, 1);
}
//G(i,j)=-x(i,j)/(2*M_PI*pow(sigma,4))*exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));
else if (type == 1){
Gy = -(-y.array().square()/2/sigma/sigma).array().exp()*(y.array())/(2*M_PI*pow(sigma,4));
Gx = (-x.array().square()/2/sigma/sigma).array().exp();
Iy = compute2(Img, Gx, false, 1, 1);
Iy = compute2(Iy, Gy, false, 1, 1);
}
//G(i,j)=-y(i,j)/(2*M_PI*pow(sigma,4))*exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));
else if (type == 2){
Gx = (-x.array().square()/2/sigma/sigma).array().exp()*(x.array().square()/sigma/sigma-1)/(2*M_PI*pow(sigma,4));
Gy = (-y.array().square()/2/sigma/sigma).array().exp();
Ixx = compute2(Img, Gx, false, 1, 1);
Ixx = compute2(Ixx, Gy, false, 1, 1);
}
//G(i,j)=1/(2*M_PI*pow(sigma,4))*(pow(x(i,j),2)/pow(sigma,2)-1) *exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));
else if (type == 3){
Gy = (-y.array().square()/2/sigma/sigma).array().exp()*(y.array().square()/sigma/sigma-1)/(2*M_PI*pow(sigma,4));
Gx = (-x.array().square()/2/sigma/sigma).array().exp();
Iyy = compute2(Img, Gx, false, 1, 1);
Iyy = compute2(Iyy, Gy, false, 1, 1);
}
//G(i,j)=1/(2*M_PI*pow(sigma,4))*(pow(y(i,j),2)/pow(sigma,2)-1) *exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));
else {
Gx = (-x.array().square()/2/sigma/sigma).array().exp()*(x.array())/(2*M_PI*pow(sigma,6));
Gy = (-y.array().square()/2/sigma/sigma).array().exp()*(y.array());
Ixy = compute2(Img, Gx, false, 1, 1);
Ixy = compute2(Ixy, Gy, false, 1, 1);
}
}
/*
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
}*/
//end for i

//Gaussian Based Image Derivatives

//static void 	cv::eigen2cv (const Eigen::Matrix< _Tp, _rows, _cols, _options, _maxRows, _maxCols > &src, Mat &dst)
//cv::eigen2cv (Img, Mat &dst);

//imfilter2(Img, Gx, Ix);
//imfilter2(Img, Gy, Iy);
//imfilter2(Img, Gxx, Ixx);
//imfilter2(Img, Gyy, Iyy);
//imfilter2(Img, Gxy, Ixy);

//Line Energy
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic>Eline;
Eline = imgaussian(Img,sigma);
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

//Imfilter function that uses correlation 
void imfilter2(Eigen::MatrixXd src_image, Eigen::MatrixXd filter,Eigen::MatrixXd &dst_image)
{
	dst_image = Eigen::MatrixXd::Zero(src_image.rows(),src_image.cols()) ;	
	
	int start_row = int(filter.rows()/2) ;
	int start_col = int(filter.cols()/2) ;
	//std::cout<<"start_row = "<<start_row<<std::endl;
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
					dst_image(tmp_row,tmp_col) -= Mid_Matrix(tmp_row + m,tmp_col + n)*filter(m,n) ; 
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

	Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Im2;
	Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> out;
	imfilter2(I, HH, Im2);
	imfilter2(Im2, HH.transpose(), out);

	return out;
	}
else
return I;

}


MatrixXd compute2(MatrixXd X, MatrixXd W, bool flip, int stride_x, int stride_y)
	{
	int width = X.cols();
	int height = X.rows();
	int kx = W.cols();
	int ky = W.rows();
	int rx = (kx-1)/2;
	int ry = (ky-1)/2;
	MatrixXd out;
	out.resize(height, width);
           for (int x=0; x<width; x+=stride_x)
           {
               int xout = x/stride_x;
   
               for (int y=0; y<height; y+=stride_y)
               {
                   int yout = y/stride_y;
   
                   double sum = 0;

		           for (int x1=x-rx; x1<=x+rx; x1++)
		           {
		               int wx = flip ? x1-x+rx : rx-x1+x;
		               for (int y1=y-ry; y1<=y+ry; y1++)
		               {
		                   if (x1>=0 && y1>=0 && x1<width && y1<height)
		                   {
		                       if (flip)
		                           sum += W(y1-y+ry,wx)*X(y1,x1);
		                       else
		                           sum += W(ry-y1+y,wx)*X(y1,x1);
		                   }
		               }
		           }
                   out(yout,xout) = sum;
               }
           }
	return out;       
	}



