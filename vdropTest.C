#include "vdrop.H"
#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Dense>

int main(){

	VectorXd dia(5);
	double n= dia.size();
	VectorXd Vol(dia.size());
	dia<< 0.0005,0.001,0.002,0.003,0.004;
	Vol=1.0/6*M_PI*dia.array().pow(3);
	double dmin=1.0e-6;	
	double sigma=0.0155;
	double rhoc=1000;
	double rhod=700;
	double dynvisc=0.0011;
	double dynvisd=0.01;
	double caliBrk=0.2;
	double OilFrac=0.5;
	VectorXd wr(dia.size());
	wr<< 0.05,0.1,0.2,0.3,0.4;


	vdrop::initialization(dmin, sigma,rhoc,rhod,dynvisc,dynvisd,caliBrk,dia,wr,OilFrac);
	double eps = 0.01;
	
	
	VectorXd tC=VectorXd::Zero(dia.size());
	VectorXd alpha=VectorXd::Zero(dia.size());
	tC(0) = 0;
	tC(1) = 0;
	tC(2) = 2e6;
	tC(3) = 0;
	tC(4) = 0;
	
	
	
	//MatrixXd T=tC.replicate(1, n);
	//cout<<"T is \n"<<T<<endl;
	
	double dt=0.001;
	double Kc=5e-5;
	for (int i = 0; i < 1; i++) 
	{
		vdrop obj;
		double sum_before = 0.0;
		sum_before=Vol.dot(tC);
		cout<<"sum_before \n"<<sum_before<<endl;
		cout<<"tC before is \n"<<tC<<endl;
		obj.vdropCalcCol(dt, eps, tC, Kc);
		double sum_after = 0.0;
		sum_after=Vol.dot(tC);
		cout<<"sum_after \n"<<sum_after<<endl;
		//alpha(1)=(tC(1)*Vol(1))/(tC(1)*Vol(1)+tC(0)*Vol(0));
		//alpha(0)=1-alpha(1);
		cout<<"Time step \n"<<i<<endl;
		cout<<"tC is \n"<<tC<<endl;
		//cout<<"alpha is \n"<<alpha<<endl;
	}
	return 0;
	
	}
