// Equations for Psi:
fvScalarMatrix PsiEqn
(             
 fvm::ddt(Psi) + isPsiAdve*fvm::div(phi, Psi)
 - mPs*aPs*fvc::laplacian(Psi)
 - mPs*bPs*fvc::laplacian(Psi * Psi * Psi)
 + mPs*kPs*fvc::laplacian(lapPsi)
 ==
  fvOptions(Psi)
  //- kPs * fvc::laplacian(lapPsi)
	      );


PsiEqn.relax();
fvOptions.constrain(PsiEqn);
PsiEqn.solve();
fvOptions.correct(Psi);   
