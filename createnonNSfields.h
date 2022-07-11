// Nematic tensor
Info<< "Reading nematic tensor Qxx and Qxy\n" << endl;
volScalarField Qxx
(
    IOobject
    (
        "Qxx",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Qxy
(
    IOobject
    (
        "Qxy",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Nematic order parameter
volScalarField S
(
    IOobject
    (
        "S",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    2*sqrt(sqr(Qxx)+sqr(Qxy))
);

// Phase field order parameter

volScalarField Psi
(
    IOobject
    (
        "Psi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField lapPsi
(
    IOobject
    (
        "lapPsi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// ------------------------------------------------------------------------
Info<< "Reading model parameters: \n" << endl;

dimensionedScalar rho
(
    transportProperties.lookup("rho")
);

dimensionedScalar nu
(
    transportProperties.lookup("nu")
);

dimensionedScalar GammaQ
(
    transportProperties.lookup("GammaQ")
);

dimensionedScalar KQ
(
    transportProperties.lookup("KQ")
);

dimensionedScalar A2Q
(
    transportProperties.lookup("A2Q")
);

dimensionedScalar A4Q
(
    transportProperties.lookup("A4Q")
);

dimensionedScalar AcQ
(
    transportProperties.lookup("AcQ")
);

dimensionedScalar Sc
(
    transportProperties.lookup("Sc")
);

dimensionedScalar lambda
(
    transportProperties.lookup("lambda")
);

dimensionedScalar alp
(
    transportProperties.lookup("alp")
);

dimensionedScalar alpGrad
(
    transportProperties.lookup("alpGrad")
);

dimensionedScalar gammav
(
    transportProperties.lookup("gammav")
);

// Turn on elastic stress that drives back flow.
dimensionedScalar isBackFlow
(
    transportProperties.lookup("isBackFlow")
);

// Turn on higher-order flow alignment.
dimensionedScalar isFlowAlignAddi
(
    transportProperties.lookup("isFlowAlignAddi")
);

// Advection of velocity
dimensionedScalar isAdveU
(
    transportProperties.lookup("isAdveU")
);

// Advection of Q
dimensionedScalar isAdveQ
(
    transportProperties.lookup("isAdveQ")
);


// Phase field parameters
dimensionedScalar mPs   // mobility
(
    transportProperties.lookup("mPs")
);
dimensionedScalar aPs   // linear term
(
    transportProperties.lookup("aPs")
);

dimensionedScalar cPs   // square nonlin.
(
    transportProperties.lookup("cPs")
);

dimensionedScalar bPs   // cubic nonlin.
(
    transportProperties.lookup("bPs")
);

dimensionedScalar kPs   // square dif.
(
    transportProperties.lookup("kPs")
);

dimensionedScalar kPsStress   // square dif.
(
    transportProperties.lookup("kPsStress")
);

dimensionedScalar isPsiAdve
(
    transportProperties.lookup("isPsiAdve")
);