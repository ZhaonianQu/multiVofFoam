double dt=runTime.deltaTValue();

//calculate epsl

//const volScalarField& epsilon=mesh.lookupObject<volScalarField>("epsilon"); 
Foam::tmp<Foam::volTensorField> gradU(fvc::grad(U));
Foam::tmp<Foam::volSymmTensorField> S(symm(fvc::grad(U)));

epsilon =(2*(turbulence->nu()+turbulence->nut())*(S && S));

forAll(mesh.C(),cellI)
{		
	double eps = epsilon[cellI];
	
	if (alphaOil[cellI]>0.001 && eps>0.0001 && alphaOil[cellI]<0.8 && mesh.C()[cellI].z()>=0.6)
	{
		
		vdrop obj;
		VectorXd nDrop = VectorXd::Zero(n);
		int i=0;
		//i is the index for nDrop
		int phaseIndex1 = 0; 
		// Index for all phases
	
		for (phase& alpha : mixture.phases())
		{	
			if (phaseIndex1 >=2)
			{
				double vf=alpha[cellI];
				nDrop(i)=max(vf,0)/Vol(i);
				++i;
			}
			++phaseIndex1;
		}
		//cout<<"SuSp L41"<<endl;
		VectorXd nDrop_old=nDrop;
		
		//double Kc=0;
		//Kc = (alphaOil[cellI] < 0.3) ? 0.0 : caliCoa;//coalescence only occurs when oil holdup > 0.3
		
		obj.vdropCalcBrkOnly(dt, eps, nDrop);//cout<<"SuSp L48"<<endl;
		//cout<<"SuSp L47"<<endl;
		
		VectorXd nDrop_delta(n);
		nDrop_delta=(nDrop-nDrop_old)/dt;
		
		//nDrop_delta is the real source term (the 			 	difference)
		//cout<<"nDrop_delta is \n"<<nDrop_delta<<endl;
		int j=0;
		//j is the index for nDrop
		int phaseIndex2 = 0;
		// Index for all phases
		//cout<<"SuSp L58"<<endl;
		for (phase& alpha : mixture.phases())
		{	
			if (phaseIndex2 >=2)
			{	
				alpha.Su()[cellI]=nDrop_delta(j)*Vol(j);
				++j;
			}
			++phaseIndex2;
		}
	}
	
	else
	{
		int j=0;
		//j is the index for nDrop
		int phaseIndex2 = 0;
		// Index for all phases
		//cout<<"SuSp L58"<<endl;
		for (phase& alpha : mixture.phases())
		{	
			if (phaseIndex2 >=2)
			{	
				alpha.Su()[cellI]=0;
				++j;
			}
			++phaseIndex2;
		}
	}

}
/*
for (phase& alpha : mixture.phases())
{	
		Info<<"phase-"<<alpha.name()<<"'s Su is\n "<<alpha.Su()<<endl;
}
*/










