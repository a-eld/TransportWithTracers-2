clc;
clear;
close all;

Nt = 1000;
T = linspace(0,1,Nt);
dt = T(2)-T(1);
n = 2;
g = 30*pi/180;
PosI = [-cos(g);sin(g)];
NI = [sin(g),cos(g)];
a = -120*pi/180;
Nmax = 20000;


SigmaF = (3).*(eye(2));
PhioF = sqrtm(SigmaF);


Phio = @(t) (1-t).*eye(2) + t.*PhioF;
Sigma = @(t) Phio(t)*Phio(t)';
A = @(t)(PhioF-eye(2))*pinv(eye(2)+t.*(PhioF-eye(2)));
dotSigma = @(t) A(t)*Sigma(t)+Sigma(t)*A(t)';

U = @(t)[cos(a*t) sin(a*t);
        -sin(a*t) cos(a*t)];
Pos = @(t) Phio(t)*U(t)*PosI;
dotPos =@(t) (Pos(t+dt)-Pos(t))/dt;
M = @(t)(1./(Pos(t)'*Pos(t))).*dotPos(t)*Pos(t)';


PhioEval = zeros(n,n,Nt);
PosEval = zeros(n,Nt);
dotPosEval = zeros(n,Nt);
N = zeros(Nt,n);
N(1,:) = NI;



for i = 1:Nt
    PhioEval(:,:,i) = Phio(T(i));
    PosEval(:,i) = Pos(T(i));
    dotPosEval(:,i) = dotPos(T(i));
end

for i = 1:Nt-1
 N(i+1,:) = N(i,:)-N(i,:)*M(T(i))*dt;
end

figure,plot(T,squeeze(PhioEval(1,1,:)),T,squeeze(PhioEval(1,2,:)),....
    T,squeeze(PhioEval(2,1,:)),T,squeeze(PhioEval(2,2,:)),LineWidth=2)
legend('$\Phi_{11}^{\rm o}$','$\Phi_{12}^{\rm o}$',...
    '$\Phi_{21}^{\rm o}$','$\Phi_{22}^{\rm o}$','FontSize',19,'interpreter','latex');



