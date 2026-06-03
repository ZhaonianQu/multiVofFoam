IOdictionary tPs
(
    IOobject
    (
        "transportProperties", 
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,  
        IOobject::NO_WRITE     
    )
);


//List<scalar> diameters = tPs.lookup("diameters");

int n= 13;
VectorXd dia(n);
VectorXd wr(n);
//VectorXd wr(n+1);//+1 is the air rising velocity, the first is air according to the phase sort in constant/transportProperties
VectorXd Vol(n);
/*
for (Foam::label i = 0; i < diameters.size(); ++i)
{
    dia(i) = diameters[i];
}
*/
dia<< 
0.0002,
0.00025,
0.00035,
0.000489,
0.000684,
0.00096,
0.00134,
0.00187,
0.0026,
0.0037,
0.0051,
0.0072,
0.01;




/*
dia<<
0.003,  
0.004;


wr<< 
0.1,
0.1,
0.2;
*/
Vol=1.0/6*M_PI*dia.array().pow(3);


int phaseIndexwr = 0;
for (phase& alpha : mixture.phases())
{	
	
	if (phaseIndexwr >=2)
	{
		wr[phaseIndexwr-2]=alpha.wrVector().value()[2];
	}
	++phaseIndexwr;
}

double dmin=1.0e-6;

double sigma=0.0124;


double rhoc=1000;


double rhod=840;


double muc=0.001;


double mud=0.004;


double Kb;
tPs.lookup("Kb") >> Kb;


vdrop::initialization(dmin, sigma,rhoc,rhod,muc,mud,Kb,dia,wr);


