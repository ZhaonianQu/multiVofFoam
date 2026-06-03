#include "multiphaseMixture.H"
#include "Time.H"
#include "subCycle.H"
#include "MULES.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "unitConversion.H"
#include "alphaContactAngleFvPatchScalarField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiphaseMixture::calcAlphas()
{
    scalar level = 0.0;
    alphas_ == 0.0;

    for (const phase& ph : phases_)
    {
        alphas_ += level * ph;
        level += 1.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseMixture::multiphaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    phases_(lookup("phases"), phase::iNew(U, phi)),

    mesh_(U.mesh()),
    U_(U),
    phi_(phi),

    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimTime, Zero)
    ),

    alphas_
    (
        IOobject
        (
            "alphas",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),

    nu_
    (
        IOobject
        (
            "nu",
            mesh_.time().timeName(),
            mesh_
        ),
        mu()/rho()
    ),

    sigmas_(lookup("sigmas")),
    dimSigma_(1, 0, -2, 0, 0),
    deltaN_
    (
        "deltaN",
        1e-8/cbrt(average(mesh_.V()))
    )
{
    rhoPhi_.setOriented();

    calcAlphas();
    alphas_.write();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::multiphaseMixture::rho() const
{
    auto iter = phases_.cbegin();

    tmp<volScalarField> trho = iter()*iter().rho();
    volScalarField& rho = trho.ref();

    for (++iter; iter != phases_.cend(); ++iter)
    {
        rho += iter()*iter().rho();
    }

    return trho;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::rho(const label patchi) const
{
    auto iter = phases_.cbegin();

    tmp<scalarField> trho = iter().boundaryField()[patchi]*iter().rho().value();
    scalarField& rho = trho.ref();

    for (++iter; iter != phases_.cend(); ++iter)
    {
        rho += iter().boundaryField()[patchi]*iter().rho().value();
    }

    return trho;
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseMixture::mu() const
{
    auto iter = phases_.cbegin();

    tmp<volScalarField> tmu = iter()*iter().rho()*iter().nu();
    volScalarField& mu = tmu.ref();

    for (++iter; iter != phases_.cend(); ++iter)
    {
        mu += iter()*iter().rho()*iter().nu();
    }

    return tmu;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::mu(const label patchi) const
{
    auto iter = phases_.cbegin();

    tmp<scalarField> tmu =
    (
        iter().boundaryField()[patchi]
       *iter().rho().value()
       *iter().nu(patchi)
    );

    scalarField& mu = tmu.ref();

    for (++iter; iter != phases_.cend(); ++iter)
    {
        mu +=
        (
            iter().boundaryField()[patchi]
           *iter().rho().value()
           *iter().nu(patchi)
        );
    }

    return tmu;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseMixture::muf() const
{
    auto iter = phases_.cbegin();

    tmp<surfaceScalarField> tmuf =
        fvc::interpolate(iter())*iter().rho()*fvc::interpolate(iter().nu());
    surfaceScalarField& muf = tmuf.ref();

    for (++iter; iter != phases_.cend(); ++iter)
    {
        muf +=
            fvc::interpolate(iter())*iter().rho()*fvc::interpolate(iter().nu());
    }

    return tmuf;
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseMixture::nu() const
{
    return nu_;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::nu(const label patchi) const
{
    return nu_.boundaryField()[patchi];
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseMixture::nuf() const
{
    return muf()/fvc::interpolate(rho());
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseMixture::surfaceTensionForce() const
{
    tmp<surfaceScalarField> tstf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "surfaceTensionForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(dimensionSet(1, -2, -2, 0, 0), Zero)
        )
    );

    surfaceScalarField& stf = tstf.ref();
    stf.setOriented();

    forAllConstIters(phases_, iter1)
    {
        const phase& alpha1 = iter1();

        auto iter2 = iter1;

        for (++iter2; iter2 != phases_.cend(); ++iter2)
        {
            const phase& alpha2 = iter2();

            auto sigma = sigmas_.cfind(interfacePair(alpha1, alpha2));

            if (!sigma.good())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << " in list of sigma values"
                    << exit(FatalError);
            }

            stf += dimensionedScalar("sigma", dimSigma_, *sigma)
               *fvc::interpolate(K(alpha1, alpha2))*
                (
                    fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                  - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                );
        }
    }

    return tstf;
}


void Foam::multiphaseMixture::solve()
{
    correct();

    const Time& runTime = mesh_.time();

    volScalarField& alpha = phases_.first();

    const dictionary& alphaControls = mesh_.solverDict("alpha");
    label nAlphaSubCycles(alphaControls.get<label>("nAlphaSubCycles"));
    scalar cAlpha(alphaControls.get<scalar>("cAlpha"));

    solveAlphas(cAlpha);

    // Update the mixture kinematic viscosity
    nu_ = mu()/rho();
}

void Foam::multiphaseMixture::solveAlphas
(
    const scalar cAlpha
)
{
    static label nSolves(-1);
    ++nSolves;

    const word alphaScheme("div(phi,alpha)");
    const word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi_/mesh_.magSf()));
    phic = min(cAlpha*phic, max(phic));

    PtrList<surfaceScalarField> alphaPhiCorrs(phases_.size());
    int phasei = 0;

    for (phase& alpha : phases_)
    {
        alphaPhiCorrs.set
        (
            phasei,
            new surfaceScalarField
            (
                "phi" + alpha.name() + "Corr",
                fvc::flux
                (
                    phi_,
                    alpha,
                    alphaScheme
                )
            )
        );

        surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phasei];

        for (phase& alpha2 : phases_)
        {
            if (&alpha2 == &alpha) continue;

            surfaceScalarField phir(phic*nHatf(alpha, alpha2));

            alphaPhiCorr += fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha,
                alpharScheme
            );
        }

        MULES::limit
        (
            1.0/mesh_.time().deltaT().value(),
            geometricOneField(),
            alpha,
            phi_,
            alphaPhiCorr,
            zeroField(),//Sp
            zeroField(),//Su
            oneField(),//alpha upper limit
            zeroField(),////alpha lower limit
            true
        );

        ++phasei;
    }

    MULES::limitSum(alphaPhiCorrs);

    rhoPhi_ = dimensionedScalar(dimMass/dimTime, Zero);

    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    phasei = 0;

    for (phase& alpha : phases_)
    {
        surfaceScalarField& alphaPhi = alphaPhiCorrs[phasei];
        alphaPhi += upwind<scalar>(mesh_, phi_).flux(alpha);

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha,
            alphaPhi
        );

        rhoPhi_ += alphaPhi*alpha.rho();

        Info<< alpha.name() << " volume fraction, min, max = "
            << alpha.weightedAverage(mesh_.V()).value()
            << ' ' << min(alpha).value()
            << ' ' << max(alpha).value()
            << endl;

        sumAlpha += alpha;

        ++phasei;
    }

    Info<< "Phase-sum volume fraction, min, max = "
        << sumAlpha.weightedAverage(mesh_.V()).value()
        << ' ' << min(sumAlpha).value()
        << ' ' << max(sumAlpha).value()
        << endl;

    // Correct the sum of the phase-fractions to avoid 'drift'
    volScalarField sumCorr(1.0 - sumAlpha);
    for (phase& alpha : phases_)
    {
        alpha += alpha*sumCorr;
    }

    calcAlphas();
}


bool Foam::multiphaseMixture::read()
{
    if (transportModel::read())
    {
        bool readOK = true;

        PtrList<entry> phaseData(lookup("phases"));
        label phasei = 0;

        for (phase& ph : phases_)
        {
            readOK &= ph.read(phaseData[phasei++].dict());
        }

        readEntry("sigmas", sigmas_);

        return readOK;
    }

    return false;
}


// ************************************************************************* //
