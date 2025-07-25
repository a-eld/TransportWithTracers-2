
function dydt = bvpfuno(yo,n) 
Phi = reshape(yo(1:n*n,1),[n,n]);
P = reshape(yo(n*n+1:2*n*n,1),[n,n]);


Ko = -P*pinv(Phi);
dy1dt = reshape(Ko*Phi,[n*n,1]);
dy2dt = reshape(-Ko'*Ko*Phi-Ko'*P,[n*n,1]);
dydt = [dy1dt;dy2dt];
end

function res = bcfuno(ya,yb,d,sigma_f) 
res = [ya(1:4)-[1;0;0;1];
       yb(1)^2+yb(3)^2-sigma_f(1,1,end);
       yb(2)^2+yb(4)^2-sigma_f(2,2,end);
       yb(1)*yb(2)+yb(3)*yb(4)-sigma_f(1,2,end);
       yb(5:8)-[2*d(1),d(2),0,0;
                d(2),2*d(3),0,0;
                0,0,2*d(1),d(2);
                0,0,d(2),2*d(3)]*yb(1:4)];
end





SigmaF = (3).*(eye(2));
PhioF = sqrtm(SigmaF);

sol_in_omt = bvpinit(T, [1;0;0;1;zeros(4,1)],ones(3,1));
options = bvpset('Nmax',Nmax);
sol_omt = bvp4c(@(t,yo,d)bvpfuno(yo,n),...
    @(ya,yb,d)bcfuno(ya,yb,d,SigmaF), sol_in_omt,options);


figure,plot(sol_omt.x,sol_omt.y(1,:),sol_omt.x,sol_omt.y(2,:),....
    sol_omt.x,sol_omt.y(3,:),sol_omt.x,sol_omt.y(4,:),LineWidth=2)
legend('$\Phi_{11}^{\rm o,n}$','$\Phi_{21}^{\rm o,n}$',...
    '$\Phi_{12}^{\rm o,n}$','$\Phi_{22}^{\rm o,n}$','FontSize',19,'interpreter','latex');
