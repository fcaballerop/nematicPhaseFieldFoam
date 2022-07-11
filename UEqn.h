// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
 fvm::ddt(U) + isAdveU*fvm::div(phi, U)
  - nu*fvm::laplacian(U)
  + MRF.DDt(U)
  - 1.0/rho*fvc::div(qStress + pStress)
  + fvm::Sp(gammav/rho,U)
 ==
    fvOptions(U)
);
// + turbulence->divDevSigma(U)
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
