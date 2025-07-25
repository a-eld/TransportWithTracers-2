function dydt = bvpfun_init(y,n) 
Phi = reshape(y(1:n*n,1),[n,n]);
P = reshape(y(n*n+1:2*n*n,1),[n,n]);


Kd = -P*Phi'*pinv(eye(2)+Phi*Phi');

dy1dt = reshape(Kd*Phi,[n*n,1]);
dy2dt = reshape(-Kd'*Kd*Phi-Kd'*P,[n*n,1]);

dydt = [dy1dt;dy2dt];
end

function g = guessdual(t,vec1,vec2)
  g = interp1(vec1,vec2',t)';
end



function res = bcfun_init(ya,yb,n,sigma1,a) 
U = [cos(a) sin(a);
    -sin(a) cos(a)];
sq_mat = sqrtm(sigma1)*U;

res = [ya(1:4)-[1;0;0;1];
       yb(1:4)-reshape(sq_mat,[n*n,1])];
     
end

vec1 = sol_omt.x;
vec2 = sol_omt.y;


sol_trn = bvpinit(T, @(t)guessdual(t,vec1,vec2));
options = bvpset('Nmax',Nmax,'Stats','on');
sold_trn = bvp5c(@(t,y)bvpfun_init(y,n),...
    @(ya,yb)bcfun_init(ya,yb,n,SigmaF,a), sol_trn,options);


Phi1_trn = [sold_trn.y(1,end),sold_trn.y(2,end);
        sold_trn.y(3,end),sold_trn.y(4,end)];

err_init = SigmaF-Phi1_trn*Phi1_trn';

figure,plot(sold_trn.x,sold_trn.y(1,:),sold_trn.x,sold_trn.y(2,:),...
    sold_trn.x,sold_trn.y(3,:),sold_trn.x,sold_trn.y(4,:),LineWidth=2)
legend('$\Phi_{11}^{\rm trn}$','$\Phi_{21}^{\rm trn}$','$\Phi_{12}^{\rm trn}$',...
    '$\Phi_{22}^{\rm trn}$','FontSize',19,'interpreter','latex');



