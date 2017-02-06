#include "interpcont.h"
#include <iostream>

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> interpcont(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P,
//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P2
int npts
){//Eigen::Matrix<double,npts,2> out;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> out;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> dis;
//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P2;
out.resize(npts,2);
//P2.resize(10*P.rows(),2);
//upsample by 10 with linear interpolation
//double diffx,diffy;

//std::cout << "I am bad" << P << P2 << "\n";

//for(int i = 0; i < P.rows()-1; i++){

//	diffx = (double)(P(i+1,0)-P(i,0))/10;
//	diffy = (double)(P(i+1,1)-P(i,1))/10;
//std::cout << "diffx = " << diffx << "\n";
//	for(int j = 0; j < 10;j++){
//		P2(10*i+j,0)=(double)P(i,0)+diffx*j;
//		P2(10*i+j,1)=(double)P(i,1)+diffy*j;
		//std::cout<<"tp2 = "<< (double)P(i,1)+diffx*j<<"tp22 = "<<P2(10*i+j,0)<<"i = "<<i<<"j="<<j<<"\n";
//}//end of for j
//}//end of for i
//	diffx = (double)(P(0,0)-P(P.rows()-1,0))/10;
//	diffy = (double)(P(0,1)-P(P.rows()-1,1))/10;
//for (int k = 0;k<10;k++){

//	P2(10*(P.rows()-1)+k,0)=(double)P((P.rows()-1),0)+diffx*k;
//	P2(10*(P.rows()-1)+k,1)=(double)P((P.rows()-1),1)+diffy*k;
//std::cout<<"a="<<(double)P((P.rows()-1),0)+diffx*k<<"b="<<P2(10*(P.rows()-1)+k,0)<<"\n";
//end of for loop
//}

// May be modified later if npts is used as params to set up # of pts on the contour
int row;
row = P.rows();
dis.resize(row,1);
//
for (int ii =0; ii<row-1; ii++){
dis(ii,0) = sqrt(pow(P(ii+1,0)-P(ii,0),2) + pow(P(ii+1,1)-P(ii,1),2));
//std::cout<<"disii = "<<dis(ii,0)<<"xi = "<<pow(P(ii+1,0)-P(ii,0),2)<<"yi = "<<pow(P(ii+1,1)-P(ii,1),2)<<std::endl;
}

for (int x =1; x<row-1; x++){
dis(x,0) = dis(x,0)+dis(x-1,0);
//std::cout<<"disii = "<<dis(ii,0)<<"xi = "<<pow(P(ii+1,0)-P(ii,0),2)<<"yi = "<<pow(P(ii+1,1)-P(ii,1),2)<<std::endl;
}


//dis(0,0)=0;
dis(row-1,0)=dis(row-2,0)+sqrt(pow(P(0,0)-P(row-1,0),2) + pow(P(0,1)-P(row-1,1),2));
//std::cout<<"dis = "<<dis<<std::endl;
double temp;
double tempx;
double tempy;
double totald = dis(row-1,0);
out(0,0)=P(0,0);
out(0,1)=P(0,1);
for (int i = 1; i<npts; i++){

	temp = (double)i/npts*totald;
  //      std::cout<<"temp = "<<temp<<"i = "<<i <<std::endl;
	for (int j = 0; j<row-1; j++){
		if(temp<dis(0,0)){
			tempx = temp/dis(0,0)*(P(1,0)-P(0,0))+P(0,0);
			tempy = temp/dis(0,0)*(P(1,1)-P(0,1))+P(0,1);
		//	std::cout<<"in for 1 "<<"tempx = "<<tempx<<"tempy = "<<tempy<<std::endl;
			out(i,0)=tempx;
			out(i,1)=tempy;
			//break;
		}
		else if(temp>dis(j,0)&(temp<=dis(j+1,0))){
		//	std::cout<<"j = "<< j <<std::endl;
			if(j+2>=row)
			{
			tempx = (temp-dis(j,0))/(dis(j+1,0)-dis(j,0))*(P(0,0)-P(j+1,0))+P(j+1,0);
			tempy = (temp-dis(j,0))/(dis(j+1,0)-dis(j,0))*(P(0,1)-P(j+1,1))+P(j+1,1);
		//	std::cout<<"in for 2 "<<"tempx = "<<tempx<<"tempy = "<<tempy<<std::endl;
			out(i,0)=tempx;
			out(i,1)=tempy;
			}


			else{
			tempx = (temp-dis(j,0))/(dis(j+1,0)-dis(j,0))*(P(j+2,0)-P(j+1,0))+P(j+1,0);
			tempy = (temp-dis(j,0))/(dis(j+1,0)-dis(j,0))*(P(j+2,1)-P(j+1,1))+P(j+1,1);
		//	std::cout<<"in for 2 "<<"tempx = "<<tempx<<"tempy = "<<tempy<<std::endl;
			out(i,0)=tempx;
			out(i,1)=tempy;
			}
			//break;
		}
	}
	//std::cout<<"outi = "<<i<<"" <<std::endl;
	

}


//dis(row-1,1) = sqrt(pow(P2(0,1)-P2(row-1,1),2) + pow(P2(0,2)-P2(row-1,2),2));
//std::cout << "I am good" << P << P2 << "\n";
return out;
}
