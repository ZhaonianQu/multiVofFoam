int phaseIndex4=0;
// Index for all phases
alphaOil = dimensionedScalar("zero", dimless, 0.0);

for (phase& alpha : mixture.phases())
{	
	alpha.Su().field() = 0.0;//After mixture.solve(), Su set back to 0
	//Info<<"set to 0\n"<<endl;
	if (phaseIndex4 >=2)
	{	
		alphaOil+=alpha;
		//scalarField alphaIf=alpha.internalField();
		//Info<<alpha.name()<<" is \n"<<alphaIf<<endl;
	}
	++phaseIndex4;
}

forAll(mesh.C(),cellJ)
{	
	int phaseJ=0;
	// Index for dia()
	int phaseIndex5=0;
	// Index for all phases
	double alphaCumu=0;
	
	if (alphaOil[cellJ]>0.0001)
	{
		for (phase& alpha : mixture.phases())
		{	
			if (phaseIndex5 >=2)
			{	
				alphaCumu+=alpha[cellJ]/(alphaOil[cellJ]);
				
				if(alphaCumu>=0.5)
				{
					d50[cellJ]=dia(phaseJ);
					break;
				}
				++phaseJ;	
			}
			++phaseIndex5;
		}
	}
	
	else
	{
		d50[cellJ]=0;
	}

}
