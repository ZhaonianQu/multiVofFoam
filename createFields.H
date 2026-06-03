Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

volVectorField wrMean
(
	IOobject
	(
		"wrMean",
		runTime.timeName(),
		mesh, 
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedVector ("wrMean", dimVelocity, Foam::vector (0,0,0))
);

volTensorField wr2Mean
(
	IOobject
	(
		"wr2Mean",
		runTime.timeName(),
		mesh, 
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedTensor ("wr2Mean", dimensionSet(0, 2, -2, 0, 0, 0 ,0), Foam::tensor (0,0,0,0,0,0,0,0,0))
);

volScalarField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar ("D", dimensionSet( 0, 2, -1, 0, 0, 0, 0), scalar (0))
);
    


multiphaseMixture mixture(U, phi, D);

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mixture.rho()
);
rho.oldTime();

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);

D = turbulence->nut() + turbulence->nu();

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p_rgh.name());

volScalarField epsilon
(
    IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh, 
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh, 
    dimensionedScalar ("epsilon", dimensionSet( 0, 2, -3, 0, 0, 0, 0), scalar (0))
);

volScalarField alphaOil
(
    IOobject
    (
        "alphaOil",
        runTime.timeName(),
        mesh, 
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh, 
    dimensionedScalar ("alphaOil", dimless, scalar (0))
);

volScalarField d50
(
    IOobject
    (
        "d50",
        runTime.timeName(),
        mesh, 
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh, 
    dimensionedScalar ("d50", dimLength, scalar (0))
);




#include "createMRF.H"
#include "createFvOptions.H"
