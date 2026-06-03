    MRF.correctBoundaryVelocity(U);
    	

	
	wrMean  = mixture.wrmean();
	
	wr2Mean = mixture.wr2mean();	

	surfaceScalarField phiwrMean
    (
         "phiwrMean",
         fvc::flux(wrMean)
    ); 
    

    /*fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
	  - fvc::div(rho*wr2Mean)
 	  + fvc::div(fvc::interpolate(rho) * phiwrMean, wrMean)
     ==
        fvOptions(rho, U)
        
    );*/
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
      + fvc::div(rho*wr2Mean)
 	  - fvc::div(fvc::interpolate(rho) * phiwrMean, wrMean)
     ==
        fvOptions(rho, U)
        
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
