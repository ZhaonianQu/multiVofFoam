#ifndef VDROP_H
#define VDROP_H
#include </apps/easybuild/software/Eigen/3.4.0-GCCcore-12.2.0/include/Eigen/Dense>


using namespace std;
using namespace Eigen;

class vdrop
{
    public:
        vdrop();
        virtual ~vdrop();
        //Input
		static double dmin_;	
		static double sigma_;
		static double rhoc_;
		static double rhod_;
		static double muc_;
		static double mud_;
		static double Kb_;
        static VectorXd dia_;
        static VectorXd wr_;
        //using input to initialize
        static int n_;
        static VectorXd vol_;
        static MatrixXd Beta_M_;
        static MatrixXd fracMat_;
        static MatrixXd volMat_;
        static MatrixXd indxMat_;      
        static VectorXd a_;
        static VectorXd b_;
        static VectorXd c_;
        static VectorXd s_;
        

        //static double Tao;
        static VectorXd Oh_;
        
        
        //for coalescence
        static MatrixXd Sij_;
        static MatrixXd rij_;
        static MatrixXd tij_;
        static MatrixXd thetaijB_;
        static MatrixXd thetaijLS_;
        static MatrixXd coaDia_;
        static MatrixXd coaVol_;
        static MatrixXd coaIndx_;
        static MatrixXd coaFrac_;
        //Member functions
		static void initialization(double dmin, double sigma,double rhoc,double rhod,double muc,double mud,double Kb, VectorXd dia,VectorXd wr);
	
        void vdropCalc(double dt, double eps,VectorXd& nDrop,double Kc);
        void vdropCalcBrkOnly(double dt, double eps,VectorXd& nDrop);
        double TotalVol(VectorXd& nDrop);
        double d50Calc(VectorXd& nDrop, double volAfter);
        void G(double eps,VectorXd& gd);  
        
        static void Beta();
        static void secondDrop(); 
        static void coaDrop(); 
        
        void Gamma(double eps,MatrixXd& Fij,double Kc);
        void GammaCol(double eps,MatrixXd& Fij,double Kc);
        static void secondDropCoa(); 
   
};

#endif
