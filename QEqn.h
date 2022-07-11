// Equations for Q:
fvScalarMatrix QxxEqn
(             
 fvm::ddt(Qxx) + isAdveQ*fvm::div(phi, Qxx)
 - GammaQ*KQ*fvm::laplacian(Qxx)
 + GammaQ*fvm::Sp(AQ,Qxx)
 - lambda*FAxx
 + 2*Qxy*W
 ==
 fvOptions(Qxx)
	      );


fvScalarMatrix QxyEqn
(             
 fvm::ddt(Qxy) + isAdveQ*fvm::div(phi, Qxy)
 - GammaQ*KQ*fvm::laplacian(Qxy)
 + GammaQ*fvm::Sp(AQ,Qxy)
 - lambda*FAxy
 - 2*Qxx*W
 ==
 fvOptions(Qxy)
	      );


QxxEqn.relax();
QxyEqn.relax();
fvOptions.constrain(QxxEqn);
fvOptions.constrain(QxyEqn);
QxxEqn.solve();
QxyEqn.solve();
fvOptions.correct(Qxx);
fvOptions.correct(Qxy);	    
