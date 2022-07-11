dimensionedScalar PsiPlus = sqrt(-aPs / bPs);

volTensorField gradU(fvc::grad(U));
volScalarField S2(4 * (sqr(Qxx) + sqr(Qxy)));
volScalarField Hxx(fvc::laplacian(KQ, Qxx) - (A2Q * 0.5 * (Psi / PsiPlus + 1) + 0.5 * A4Q * S2) * Qxx);
volScalarField Hxy(fvc::laplacian(KQ, Qxy) - (A2Q * 0.5 * (Psi / PsiPlus + 1) + 0.5 * A4Q * S2) * Qxy);
S = sqrt(S2);

volScalarField Dxx = 0.5 * (Psi / PsiPlus + 1) * gradU.component(tensor::XX);
volScalarField Dxy = 0.5 * (Psi / PsiPlus + 1) * 0.5 * (gradU.component(tensor::XZ) + gradU.component(tensor::ZX));
volScalarField W = 0.5 * (gradU.component(tensor::XZ) - gradU.component(tensor::ZX));
volScalarField FAxx = Dxx - isFlowAlignAddi * 4 * Qxx * (Dxx * Qxx + Dxy * Qxy);
volScalarField FAxy = Dxy - isFlowAlignAddi * 4 * Qxy * (Dxx * Qxx + Dxy * Qxy);
volScalarField AQ = A2Q * 0.5 * (Psi / PsiPlus + 1) + 0.5 * A4Q * S2;

// Intermediate quantities for scalar field
// The first two are not needed any more

volScalarField Psi2_l(sqr(Psi));
volScalarField Psi_l(Psi);
volVectorField gradPsi(fvc::grad(Psi));
volScalarField gradPsiSq(gradPsi & gradPsi);

// Stress tensor originating from nematic fields:
volTensorField qStress(
    IOobject(
        "qStress",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedTensor(
        "qStress",
        dimensionSet(1, -1, -2, 0, 0, 0, 0),
        tensor::zero));

qStress.replace(0, 0.5 * (Psi / PsiPlus + 1) * (alp + alpGrad * gradPsiSq) * Qxx - isBackFlow * lambda * Hxx);
qStress.replace(2, 0.5 * (Psi / PsiPlus + 1) * (alp + alpGrad * gradPsiSq) * Qxy - isBackFlow * (lambda * Hxy + 2 * (Qxx * Hxy - Qxy * Hxx)));
qStress.replace(6, 0.5 * (Psi / PsiPlus + 1) * (alp + alpGrad * gradPsiSq) * Qxy - isBackFlow * (lambda * Hxy - 2 * (Qxx * Hxy - Qxy * Hxx)));
qStress.replace(8,-0.5 * (Psi / PsiPlus + 1) * (alp + alpGrad * gradPsiSq) * Qxx + isBackFlow * lambda * Hxx);

// Has to be changed to simply:
// volScalarField lapPsi(fvc::laplacian(Psi)); 
// Check speed difference

fvScalarMatrix lapPsiEqn(
    fvm::Sp(1, lapPsi) == fvc::laplacian(Psi));
lapPsiEqn.relax();
fvOptions.constrain(lapPsiEqn);
lapPsiEqn.solve();
fvOptions.correct(lapPsi);

// In the case of non-scalar fields, the overridden operators
// '*' and '&' return the outer and inner products of the two 
// operands respectively, so this line implements the stress:
//
//    - k [ (nabla phi) x (nabla phi) - 1/2 delta_ij (nabla phi) · (nabla phi) ]
//
volTensorField pStress(-kPsStress * (gradPsi * gradPsi - 0.5 * (gradPsi & gradPsi) * tensor::I));
