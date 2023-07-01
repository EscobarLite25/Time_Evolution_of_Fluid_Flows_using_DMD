clear all, close all, clc
load CYLINDER_ALL.mat
X = VORTALL(:,1:end-1);
X2 = VORTALL(:,2:end);
[U,S,V] = svd(X,'econ');

%%  Compute DMD (Phi are eigenvectors)
r = 21;  % truncate at 21 modes
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde = U'*X2*V*inv(S);
[W,eigs] = eig(Atilde);
Phi = X2*V*inv(S)*W;

%% Plot DMD modes
for i=10:2:20
    plotCylinder(reshape(real(Phi(:,i)),nx,ny),nx,ny);
    plotCylinder(reshape(imag(Phi(:,i)),nx,ny),nx,ny);
end

