#include "vdrop.H"
#include <cmath>
#include <cstdlib>
#include <iostream>


Eigen::VectorXd vdrop::vol_; 
Eigen::VectorXd vdrop::a_; 
Eigen::VectorXd vdrop::b_; 
Eigen::VectorXd vdrop::c_; 
Eigen::VectorXd vdrop::s_; 


Eigen::VectorXd vdrop::Oh_; 
Eigen::MatrixXd vdrop::Beta_M_; 
Eigen::MatrixXd vdrop::fracMat_; 
Eigen::MatrixXd vdrop::volMat_; 
Eigen::MatrixXd vdrop::indxMat_; 
Eigen::MatrixXd vdrop::Sij_;
Eigen::MatrixXd vdrop::rij_;
Eigen::MatrixXd vdrop::tij_;
Eigen::MatrixXd vdrop::thetaijB_;
Eigen::MatrixXd vdrop::thetaijLS_;
Eigen::MatrixXd vdrop::coaDia_;
Eigen::MatrixXd vdrop::coaVol_;
Eigen::MatrixXd vdrop::coaIndx_; 
Eigen::MatrixXd vdrop::coaFrac_; 
//double vdrop::Tao; 
int vdrop::n_; 
Eigen::VectorXd vdrop::dia_;
Eigen::VectorXd vdrop::wr_;
double vdrop::dmin_; 
double vdrop::sigma_; 
double vdrop::rhoc_; 
double vdrop::rhod_; 
double vdrop::muc_; 
double vdrop::mud_; 
double vdrop::Kb_; 

vdrop::vdrop()
{
//cout << "Object is being created!" << std::endl;
}

vdrop::~vdrop()
{
//cout << "Object is being deleted!" << std::endl;   
}


void vdrop::initialization(double dmin, double sigma,double rhoc,double rhod,double muc,double mud,double Kb, VectorXd dia, VectorXd wr) 
{
	//cout<<"vdrop initialization is starting"<<endl;
	vdrop::dmin_=dmin; 
	vdrop::sigma_=sigma; 
	vdrop::rhoc_=rhoc; 
	vdrop::rhod_=rhod; 
	vdrop::muc_=muc; 
	vdrop::mud_=mud; 
	vdrop::Kb_=Kb;
	vdrop::dia_=dia; 
	vdrop::wr_=wr;
//second droplet matrix****************************************************************
	vdrop::n_=dia_.size();
	vol_=VectorXd(n_);
	vol_=(M_PI/6)*dia_.array().pow(3);
	indxMat_=MatrixXd(n_,n_);
	volMat_=MatrixXd(n_,n_);
	fracMat_=MatrixXd(n_,n_);
	secondDrop();
//Beta matrix**************************************************************************
	Beta_M_=MatrixXd::Zero(n_,n_);
	Beta();
	
//g function vector********************************************************************
	Oh_=VectorXd(n_);
	Oh_ = mud_/sqrt(rhod_*sigma_*dia_.array());
	//double Tao = (mud_/muc_)*sqrt(rhoc_/rhod_);
	
	double a1, a2, a3, a4;
	double b1, b2, b3, b4; 
	double c1, c2, c3, c4,c5;
	double s1, s2, s3, s4;


	a1=-0.2211;  a2=-37.5278; a3=1.9411;  a4=-69.8917;
	b1=207.8005; b2=52.9507; b3=-207.563; b4=52.9613; 
	c1=677920000; c2=-150140000;  c3=11823000; c4=-370610; c5=2812.6;
	s1=5.7626; s2=-90.9438; s3=-4.3021;  s4=0.766;
	
	
	//a1--->e3 value need to change when continuous or dispersed phase properties change

	a_=VectorXd(n_);
	b_=VectorXd(n_);
	c_=VectorXd(n_);
	s_=VectorXd(n_);
	
	
	a_=a1*exp(a2*Oh_.array())+a3*exp(a4*Oh_.array());
	//cout<<"a_ =\n"<<a_<<endl;

	b_=b1*exp(b2*Oh_.array())+b3*exp(b4*Oh_.array());
	
	c_=c1*pow(Oh_.array(),4)+c2*pow(Oh_.array(),3)+
	c3*pow(Oh_.array(),2)+c4*Oh_.array()+c5;
	
	s_=s1*exp(s2*Oh_.array())+s3*exp(s4*Oh_.array());
	//cout<<"c_ =\n"<<c_<<endl;

	
	//cout<<"s =\n"<<s_<<endl;

	
	//cout<<"e_ =\n"<<e_<<endl;
	
//Gamma function vector****************************************************************
	/*
	
	Sij_=MatrixXd::Zero(n_,n_);
	Sij_=(M_PI/4)*pow((dia_.replicate(1, n_)+ dia_.transpose().replicate(n_, 1)).array(),2);
	
	rij_=MatrixXd::Zero(n_,n_);
	for(int i=0;i<n_;i++)
	{
		for(int j=0;j<n_;j++)
		{
			rij_(i,j)=0.5/(2.0/dia_(i)+2.0/dia_(j));
		}
	}
	
	tij_=MatrixXd::Zero(n_,n_);
	tij_=pow((pow(rij_.array(),3)*rhoc_/(16*sigma_)),0.5)*std::log(1e4);
	
	thetaijB_=MatrixXd::Zero(n_,n_);
	thetaijB_=Sij_.array()*abs(wr_.replicate(1, n_).array()-wr_.transpose().replicate(n_, 1).array());
	
	thetaijLS_=MatrixXd::Zero(n_,n_);
	double shearRate=(1-0.75*OilFrac)/(1-OilFrac)*OilFrac*68.9*rhoc_;
	thetaijLS_=1.0/6*pow((dia_.replicate(1, n_).array()+dia_.transpose().replicate(n_, 1).array()).array(),3)*shearRate;
		
	coaDia_=pow(((dia_.replicate(1, n_)).array().pow(3)+ (dia_.transpose().replicate(n_, 1)).array().pow(3)).array(),1.0/3);
	coaVol_=(M_PI/6)*coaDia_.array().pow(3);
	
	coaIndx_=MatrixXd::Zero(n_,n_);
	coaFrac_=MatrixXd::Zero(n_,n_);
	coaDrop();
	//cout<<"coaIndx_ is\n"<<coaIndx_<<endl;
	//cout<<"coaFrac_ is\n"<<coaFrac_<<endl;
	
	*/
	
	
	
	
}
//Preparing work finished**************************************************************


//main function************************************************************************
void vdrop::vdropCalc(double dt, double eps, VectorXd& nDrop,double Kc)
{
	//cout<<"vdropCalc is starting"<<endl;
	VectorXd nSec=VectorXd::Zero(n_);
	VectorXd gd(n_);//brk efficiency
	VectorXd gdn_pre(n_);//brk number pre
	VectorXd gdn(n_); //brk number correct
	VectorXd gene1(n_); 
	G(eps, gd);//calculating breakup effiency value
	//cout<<"gd="<<gd<<endl;
	//cout<<"nDrop at VDROP beginning is =\n"<<nDrop<<"\n"<<endl;
	gdn_pre=gd.array()*nDrop.array();
	
	MatrixXd Fij=MatrixXd::Zero(n_,n_);
	Gamma(eps,Fij,Kc);
	//cout<<"Fij=\n"<<Fij<<"\n"<<endl;
	MatrixXd N_coa=MatrixXd::Zero(n_,n_);
	N_coa=(nDrop.replicate(1, n_).array()* nDrop.transpose().replicate(n_, 1).array()).array()*Fij.array();
	//cout<<"N_coa1=\n"<<N_coa<<"\n"<<endl;
	VectorXd n_coa_pre=VectorXd::Zero(n_);
	n_coa_pre=N_coa.rowwise().sum();
	VectorXd n_coa=VectorXd::Zero(n_);//coa disappeared
	VectorXd coa_n=VectorXd::Zero(n_);//coa generated
	
	//correct coalescence and breakup total number to prevent it>total droplet number
	VectorXd n_react_pre=n_coa_pre+gdn_pre;
	VectorXd gdnCorr=nDrop.array()*gdn_pre.array()/(n_react_pre.array()+1e-30);
	VectorXd n_coaCorr=nDrop.array()*n_coa_pre.array()/(n_react_pre.array()+1e-30);
	gdn=(n_react_pre.array()*dt<=nDrop.array()).select(gdn_pre,gdnCorr/dt);
	n_coa=(n_react_pre.array()*dt<=nDrop.array()).select(n_coa_pre,n_coaCorr/dt);//divided by dt is correct
	//cout<<"n_coa=\n"<< n_coa<<"\n"<<endl;
	VectorXd r=n_coa.array()/(n_coa_pre.array()+1e-30);
	MatrixXd R=MatrixXd::Zero(n_,n_);
	for(int i=0;i<n_;i++)
	{	
		for(int j=0;j<n_;j++)
		{
			R(i,j)=min(r(i),r(j));
		}
	}
	//cout<<"R=\n"<<R<<"\n"<<endl;
	N_coa=N_coa.array() * R.array();
	//cout<<"N_coa2=\n"<<N_coa<<"\n"<<endl;
	n_coa=N_coa.rowwise().sum();
	//cout<<"N_coa2=\n"<<N_coa<<"\n"<<endl;
	//finish correction
	
	gene1=Beta_M_*gdn;//the first breakup droplet
	//cout<<"VDROP line113\n"<<endl;
	//double volBefore = this->TotalVol(nDrop);//cout<<"volBefore is\n"<<volBefore <<endl;

	//for loop to calculate second breakup droplet
	for(int i=1;i<n_;i++)
	{	
		for(int j=0;j<i;j++)
		{
			int indx = indxMat_(j,i);	 
			double v1 = vol_(indx);
			double v2 = vol_(indx+1);
			 
			nSec[indx] += Beta_M_(j,i)*gdn(i)*(1-fracMat_(j,i))*volMat_(j,i)/v1;//interpolate to smaller bin
			nSec[indx+1] += Beta_M_(j,i)*gdn(i)*fracMat_(j,i)*volMat_(j,i)/v2;//interpolate to larger bin
		}
	}
	
	for(int i=0;i<n_;i++)
	{	
		for(int j=0;j<n_;j++)
		{
			int coaidx = coaIndx_(i,j);
			if (coaidx<n_-1)
			{
			double v1 = vol_(coaidx);
			double v2 = vol_(coaidx+1);
			coa_n(coaidx)+=N_coa(i,j)*coaFrac_(i,j)*coaVol_(i,j)/v1;
			coa_n(coaidx+1)+=N_coa(i,j)*(1-coaFrac_(i,j))*coaVol_(i,j)/v2;
			} 
			else
			{
			coa_n(coaidx)+=N_coa(i,j)*coaVol_(i,j)/vol_(coaidx);
			}
		}
	}
	
	coa_n=0.5*coa_n;//calculate twice for di and dj coalescence
	//cout<<"coa_n=\n"<< coa_n<<"\n"<<endl;
	nDrop = nDrop+(-gdn+gene1+nSec-n_coa+coa_n)*dt;
	//cout<<"gdn=\n"<< gdn<<endl;
	//cout<<"gene1=\n"<< gene1<<endl;
	//cout<<"nSec=\n"<< nSec<<endl;
	//cout<<"n_coa=\n"<< n_coa<<endl;
	//cout<<"coa_n=\n"<< coa_n<<endl;
	//cout<<"nDROP after\n"<< nDrop<<endl;
	//Total volume after the breakage

	//double deltaVol = volAfter - volBefore;cout<<"vdrop138"<<endl;
	//double error=deltaVol/volBefore;cout<<"vdrop139"<<endl;
	// to check algorithm
	//cout<<"vol_ matrix is\n"<<volMat_ <<endl;
	//cout<<"indxMat_ matrix is\n"<<indxMat_ <<endl;
	//cout<<"frac matrix is\n"<<fracMat_ <<endl;
	//cout<<"Beta matrix is\n"<<Beta_M_ <<endl;
	//
	//cout<<"break matrix is\n"<<gdn <<endl;
	//cout<<"gene1 matrix is\n"<<Beta_M_*gdn <<endl;
	//cout<<"gene2 matrix is\n"<<nSec<<endl;
	//cout<<"next teim step is\n"<<nDrop <<endl;
	//cout<<"volBefore is\n"<<volBefore <<endl;
	//cout<<"volAfter is\n"<<volAfter <<endl;
	//cout<<"error is "<<error <<endl;

}

//vdrop of Column************************************************************************
void vdrop::vdropCalcBrkOnly(double dt, double eps, VectorXd& nDrop)
{
	//cout<<"vdropCalc is starting"<<endl;
	VectorXd nSec=VectorXd::Zero(n_);
	VectorXd gd(n_);//brk efficiency
	VectorXd gdn_pre(n_);//brk number pre
	VectorXd gdn(n_); //brk number correct
	VectorXd gene1(n_); 
	G(eps, gd);//calculating breakup effiency value
	//cout<<"gd="<<gd<<endl;
	//cout<<"nDrop at VDROP beginning is =\n"<<nDrop<<"\n"<<endl;
	gdn_pre=gd.array()*nDrop.array();
	
	
	//correct coalescence and breakup total number to prevent it>total droplet number
	VectorXd n_react_pre=gdn_pre;
	VectorXd gdnCorr=nDrop.array()*gdn_pre.array()/(n_react_pre.array()+1e-30);
	gdn=(n_react_pre.array()*dt<=nDrop.array()).select(gdn_pre,gdnCorr/dt);
	
	gene1=Beta_M_*gdn;//the first breakup droplet
	//cout<<"VDROP line113\n"<<endl;
	//double volBefore = this->TotalVol(nDrop);//cout<<"volBefore is\n"<<volBefore <<endl;

	//for loop to calculate second breakup droplet
	for(int i=1;i<n_;i++)
	{	
		for(int j=0;j<i;j++)
		{
			int indx = indxMat_(j,i);	 
			double v1 = vol_(indx);
			double v2 = vol_(indx+1);
			 
			nSec[indx] += Beta_M_(j,i)*gdn(i)*(1-fracMat_(j,i))*volMat_(j,i)/v1;//interpolate to smaller bin
			nSec[indx+1] += Beta_M_(j,i)*gdn(i)*fracMat_(j,i)*volMat_(j,i)/v2;//interpolate to larger bin
		}
	}
	
	nDrop = nDrop+(-gdn+gene1+nSec)*dt;
	//cout<<"gdn=\n"<< gdn<<endl;
	//cout<<"gene1=\n"<< gene1<<endl;
	//cout<<"nSec=\n"<< nSec<<endl;
	//cout<<"n_coa=\n"<< n_coa<<endl;
	//cout<<"coa_n=\n"<< coa_n<<endl;
	//cout<<"nDROP after\n"<< nDrop<<endl;
	//Total volume after the breakage

	//double deltaVol = volAfter - volBefore;cout<<"vdrop138"<<endl;
	//double error=deltaVol/volBefore;cout<<"vdrop139"<<endl;
	// to check algorithm
	//cout<<"vol_ matrix is\n"<<volMat_ <<endl;
	//cout<<"indxMat_ matrix is\n"<<indxMat_ <<endl;
	//cout<<"frac matrix is\n"<<fracMat_ <<endl;
	//cout<<"Beta matrix is\n"<<Beta_M_ <<endl;
	//
	//cout<<"break matrix is\n"<<gdn <<endl;
	//cout<<"gene1 matrix is\n"<<Beta_M_*gdn <<endl;
	//cout<<"gene2 matrix is\n"<<nSec<<endl;
	//cout<<"next teim step is\n"<<nDrop <<endl;
	//cout<<"volBefore is\n"<<volBefore <<endl;
	//cout<<"volAfter is\n"<<volAfter <<endl;
	//cout<<"error is "<<error <<endl;

}

//Start calculating Beta*******************************************************************************************************
void vdrop::Beta()
{
	//efMin
	VectorXd efMin=M_PI*sigma_*((dia_.array().pow(3)-pow(dmin_,3)).pow(2.0/3)+pow(dmin_,2)-dia_.array().pow(2));
	
	//efMax
	VectorXd efMax=M_PI*sigma_*(2*(dia_.array()/pow(2,1.0/3)).pow(2)-dia_.array().pow(2));

	
	//ef
	MatrixXd ef=MatrixXd::Zero(n_,n_);
	for (int i = 0; i < n_; i++)
	{
		for (int j = i; j < n_; j++)
		{
		    ef(i, j) = M_PI * sigma_ * (pow((pow(dia_(j), 3.0) - pow(dia_(i), 3.0)), 2.0 / 3) + pow(dia_(i), 2) - pow(dia_(j), 2));

		    ef(i, j) = max(ef(i, j), 0.0);
		}
	}
	
	//cout<<"ef is\n"<<ef <<endl;

// Beta calculation
	for (int i = 1; i < n_; i++)
	{
		double sum = 0.0;
		for (int j = 0; j < i; j++)
		{
		    sum += efMin(i) + efMax(i) - ef(j, i);
		}

		for (int j = 0; j < i; j++)
		{
		    Beta_M_(j, i) = (efMin(i) + efMax(i) - ef(j, i)) / sum;
		}
	}

}
//Finish calculating Beta*******************************************************************************************************




//Start calculating g function**************************************************************************************************
void vdrop::G(double eps,VectorXd& gd)
{	
	//cout<<"In G function gd_initial\n"<<gd<<endl; 
	VectorXd Re = rhoc_*pow(eps,1.0/3.0)*dia_.array().pow(4.0/3.0)/muc_;
	//cout<<"In G function Re:\n"<<Re<<endl;
	VectorXd ans=a_.array()*exp(b_.array()*Oh_.array())+c_.array()*exp(s_.array()*Oh_.array())-a_.array(); 
    //VectorXd ans=a_.array()*pow(log10(Re.array()),
    //b_)+ c_.array()*pow(log10(Re.array()),s_.array()) - e_.array();  
    //cout<<"ans\n"<<ans<<endl;
    
    
    VectorXd gf=pow(10,ans.array());
    //cout<<"gf\n"<<gf<<endl;
        
	gd=(0.638*Kb_/(pow(eps,-1.0/3.0)*dia_.array().pow(2.0/3))*gf.array());
	gd(0)=0;
	//cout<<"In G function gd_final\n"<<gd<<endl;
       
}
//Finish calculating g function*************************************************************************************************



//Start calculating Gamma function**************************************************************************************************
void vdrop::Gamma(double eps,MatrixXd& Fij, double Kc)
{	
       	VectorXd ud2(n_);
       	ud2=pow(1.03*(eps*dia_.array()).pow(1.0/3),2);
       	MatrixXd taoij;
       	taoij=(rij_.array().pow(2.0/3))/pow(eps,1.0/3);
       	MatrixXd Pij;
       	Pij=exp(-tij_.array()/taoij.array());
       	MatrixXd thetaij;
       	thetaij=Sij_.array()*pow((ud2.replicate(1, n_) + ud2.transpose().replicate(n_, 1)).array(),0.5)+thetaijB_.array();
       	//cout<<"thetaij is\n"<<thetaij<<"\n"<<endl;
       	//cout<<"thetaijB is\n"<<thetaijB_<<"\n"<<endl;
       	
       	Fij=Kc*thetaij.array()*Pij.array();
       	
       	//Fij = (coaDia_.array() > dia_(n_ - 1)).select(0, Fij);

		//cout<<"Fij is\n"<<Fij<<"\n"<<endl;
}

void vdrop::GammaCol(double eps,MatrixXd& Fij, double Kc)
{	
       	VectorXd ud2(n_);
       	ud2=pow(1.03*(eps*dia_.array()).pow(1.0/3),2);
       	MatrixXd taoij;
       	taoij=(rij_.array().pow(2.0/3))/pow(eps,1.0/3);
       	MatrixXd Pij;
       	Pij=exp(-tij_.array()/taoij.array());
       	MatrixXd thetaij;
       	thetaij=Sij_.array()*pow((ud2.replicate(1, n_) + ud2.transpose().replicate(n_, 1)).array(),0.5)+thetaijB_.array()+thetaijLS_.array();
       	
       	//cout<<"thetaij is\n"<<thetaij<<"\n"<<endl;
       	//cout<<"thetaijT is\n"<<thetaij-thetaijB_-thetaijLS_<<"\n"<<endl;
       	//cout<<"thetaijB is\n"<<thetaijB_<<"\n"<<endl;
       	//cout<<"thetaijLS is\n"<<thetaijLS_<<"\n"<<endl;
       	Fij=Kc*thetaij.array()*Pij.array();
       	
       	//Fij = (coaDia_.array() > dia_(n_ - 1)).select(0, Fij);

		//cout<<"Fij is\n"<<Fij<<"\n"<<endl;
}
//Finish calculating Gamma function*************************************************************************************************


//Calculating the Matrix of second droplet**************************************************************************************
void vdrop::secondDrop() 
{
	double ds,dl;
	for(int i=1;i<n_;i++)
	{
		dl = dia_[i];
		for(int j=0;j<i;j++)
		{
			ds = dia_[j];
	
			double dnew = 0.0,dvol = 0.0, fractn = 0.0;
			dvol = (M_PI/6)*(pow(dl,3) - pow(ds,3));
			dnew = pow(pow(dl,3) - pow(ds,3),1.0/3);
	
	
			int ismall=0;
			for(int k=0;k<n_;k++)
			{
				if(dnew<=dia_(k))
				{
					ismall=k-1;
					break;
				}
			}

			if(ismall==-1)
			{
				indxMat_(j,i)= 0;
				volMat_(j,i)= dvol;
				fracMat_(j,i)= 1;
			}
	
			else
			{
				fractn=(dnew-dia_[ismall])/(dia_[ismall+1]-dia_[ismall]);
				volMat_(j,i)= dvol;
				fracMat_(j,i)= fractn;
				indxMat_(j,i)= ismall;				
			}
			
		}
	}
	
}
//Finished*************************************************************************************************************************

//Calculating the Matrix of coalescence droplet*************************************************************************
void vdrop::coaDrop() 
{
	for(int i=0;i<n_;i++)
	{
		for(int j=0;j<n_;j++)
		{
			bool found = false;
			for (int k=1;k<n_;k++)
			{
				if (coaDia_(i,j)<=dia_(k))
				{
					coaIndx_(i,j)=k-1;	
					coaFrac_(i,j)=1-(coaDia_(i,j)-dia_(k-1))/(dia_(k)-dia_(k-1));
					found = true;
					break;
				}	
			}
			if (!found)
		    {
		        coaIndx_(i,j) = n_ - 1;
		        coaFrac_(i,j) = 1; 
		    }
		}
	
	}
}
//Finished************************************************************************


//Calculating total volume*********************************************************************************************************
double vdrop::TotalVol(VectorXd& nDrop)
{
	double sum = 0.0;
	sum=vol_.dot(nDrop);
    	return sum;
}
//Finished*************************************************************************************************************************
