#include "imDev.h"
#include <iostream>

//typedef typename Matrix::Scalar T;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatrixXd;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> VectorXd;


MatrixXd compute(MatrixXd X, MatrixXd W, bool flip, int stride_x, int stride_y);

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
x.resize(total,1);
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> y;
y.resize(1,total);

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> out;
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gx;
Gx.resize(total,1);
Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gy;
Gy.resize(1,total);
//Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gxx;
//Gxx.resize(total,total);
//Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gyy;
//Gyy.resize(total,total);
//Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic> Gxy;
G.resize(total,total);
//std::cout<<"total = "<<total<<std::endl;
for (int i = 0; i< total;i++){
	//for (int j = 0; j< total;j++){
//std::cout<<"i = "<<i<<" i-3sig = "<<i-sigma*3<<std::endl;
		x(i,0) = i-sigma*3;
		y(0,i) = i-sigma*3;
}
//Gaussian Filters
//std::cout<<" x =  "<<x<<std::endl;
//std::cout<<" y =  "<<y<<std::endl;
//std::cout<<"step 1 "<<std::endl;
if (type == 1){
Gx = -(-x.array().square()/2/sigma/sigma).array().exp()*(x.array())/(2*M_PI*pow(sigma,4));
Gy = (-y.array().square()/2/sigma/sigma).array().exp();
}
//G(i,j)=-x(i,j)/(2*M_PI*pow(sigma,4))*exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));
else if (type == 2){
Gy = -(-y.array().square()/2/sigma/sigma).array().exp()*(y.array())/(2*M_PI*pow(sigma,4));
Gx = (-x.array().square()/2/sigma/sigma).array().exp();
}
//G(i,j)=-y(i,j)/(2*M_PI*pow(sigma,4))*exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));
else if (type == 3){
Gx = (-x.array().square()/2/sigma/sigma).array().exp()*(x.array().square()/sigma/sigma-1)/(2*M_PI*pow(sigma,4));
Gy = (-y.array().square()/2/sigma/sigma).array().exp();
}
//G(i,j)=1/(2*M_PI*pow(sigma,4))*(pow(x(i,j),2)/pow(sigma,2)-1) *exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));
else if (type == 4){
Gy = (-y.array().square()/2/sigma/sigma).array().exp()*(y.array().square()/sigma/sigma-1)/(2*M_PI*pow(sigma,4));
Gx = (-x.array().square()/2/sigma/sigma).array().exp();
}
//G(i,j)=1/(2*M_PI*pow(sigma,4))*(pow(y(i,j),2)/pow(sigma,2)-1) *exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));
else {
Gx = (-x.array().square()/2/sigma/sigma).array().exp()*(x.array())/(2*M_PI*pow(sigma,6));
Gy = (-y.array().square()/2/sigma/sigma).array().exp()*(y.array());
}
//G(i,j)=1/(2*M_PI*pow(sigma,6))*(x(i,j)*y(i,j)) *exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));
//std::cout<<"step 2 "<<std::endl;

//}//end for j
//end for i
//std::cout<<"=================================================================="<<std::endl;
//std::cout<<"type = "<<type<<"Gx = "<<Gx<<std::endl;
//std::cout<<"x = "<<x<<std::endl;
//std::cout<<"last vv = "<<Img<<std::endl;
//std::cout<<"before last vout = "<<out.colwise().sum()  <<std::endl;
//std::cout<<"type = "<<type<<"Gy = "<<Gy<<std::endl;
//std::cout<<"type = "<<type<<" imsumc = "<<Img.colwise().sum()<<std::endl;
//std::cout<<"type = "<<type<<" imsumr = "<<Img.rowwise().sum()<<std::endl;
//std::cout<<"=================================================================="<<std::endl;
//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> out;

//use 2D convolution as the filter
//If true, use correlation 
out = compute(Img, Gx, false, 1, 1);
out = compute(out, Gy, false, 1, 1);
//imfilter2a(Img, Gx, out);
//imfilter2a(out, Gy, out);
//std::cout<<"Type = "<< type <<" last out = "<<out.colwise().sum()  <<std::endl;
return out;

}


void imfilter2a(Eigen::MatrixXd src_image, Eigen::MatrixXd filter,Eigen::MatrixXd &dst_image)
{	Eigen::MatrixXd filterconv;
	filterconv = filter.rowwise().reverse().colwise().reverse();
	dst_image = Eigen::MatrixXd::Zero(src_image.rows(),src_image.cols()) ;
//std::cout<<"=================================================================="<<std::endl;	
	//std::cout<<"filter1 = "<<filter<<std::endl;
//std::cout<<"=================================================================="<<std::endl;
	//std::cout<<"filter2 = "<<filterconv<<std::endl;
//std::cout<<"=================================================================="<<std::endl;
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
					dst_image(tmp_row,tmp_col) += Mid_Matrix(tmp_row + m,tmp_col + n)*filterconv(m,n) ; 
				}
			}
		}
	}
//	std::cout<<"fsize = "<<filter.rows()<<" outim r = "<<dst_image.rowwise().sum()<<std::endl;
//	std::cout<<"outim c = "<<dst_image.colwise().sum()<<std::endl;


	return ;
}



MatrixXd compute(MatrixXd X, MatrixXd W, bool flip, int stride_x, int stride_y)
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





