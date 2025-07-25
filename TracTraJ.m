function dydt = bvpfunTraJ(t,y,Xt,dotXt,n,N,T) 
N = interp1(T,N,t);
Phi = reshape(y(1:n*n,1),[n,n]);
P = reshape(y(n*n+1:2*n*n,1),[n,n]);


B = eye(2)+Phi*Phi';
M = (1/(Xt'*Xt)).*dotXt*Xt';
R = -((M*B+P*Phi')*N')*pinv(N*B*N');

K = M+R*N;
dy1dt = reshape(K*Phi,[n*n,1]);
dy2dt = reshape(-K'*K*Phi-K'*P,[n*n,1]);

dydt = [dy1dt;dy2dt];
end



function res = bcTraJ(ya,yb,n,sigma1,a) 
U = [cos(a) sin(a);
    -sin(a) cos(a)];
sq_mat = sqrtm(sigma1)*U;

res = [ya(1:4)-[1;0;0;1];
       yb(1:4)-reshape(sq_mat,[n*n,1])];
     
end


function g = guessTraJ(t,vec1,vec2)
  g = interp1(vec1,vec2',t)';
end

vec1 = sol_trn.x;
vec2 = sol_trn.y;


solinit = bvpinit(T, @(t)guessTraJ(t,vec1,vec2));
options = bvpset('Nmax',Nmax,'RelTol',1e-7,'AbsTol',1e-7,'Stats','on');
sold = bvp5c(@(t,y)bvpfunTraJ(t,y,Pos(t),dotPos(t),n,N,T),...
    @(ya,yb)bcTraJ(ya,yb,n,SigmaF,a),solinit,options);


figure,plot(sold.x,sold.y(1,:),sold.x,sold.y(2,:),...
    sold.x,sold.y(3,:),sold.x,sold.y(4,:),LineWidth=2)
ax = gca;
ax.DataAspectRatio=[0.1,1,1];
ax.FontSize = 15;
legend('$\Phi_t^\star(1,1)$','$\Phi_t^\star(2,1)$','$\Phi_t^\star(1,2)$','$\Phi_t^\star(2,2)$',...
                    'FontSize',14,'interpreter','latex','Location','northeastoutside');
xlabel('Time \ t', 'FontSize',19 ,'interpreter','latex');

Phit_mat = zeros(n,n,size(sold.y,2));
sigma_out = zeros(n,n,size(sold.y,2));
pos_out = zeros(2,length(sold.x));
posi = Pos(0);

for k = 1: size(Phit_mat,3)
    Phit_mat(:,:,k) = [sold.y(1,k),sold.y(3,k);
                       sold.y(2,k),sold.y(4,k)];
    sigma_out(:,:,k) =  Phit_mat(:,:,k)*Phit_mat(:,:,k)';
     pos_out(:,k) = Phit_mat(:,:,k)*posi;
end

angle_in = acosd(dot(PosEval(:,1),PosEval(:,end))./(norm(PosEval(:,1))*norm(PosEval(:,end))));
angle_out = acosd(dot(pos_out(:,1),pos_out(:,end))./(norm(pos_out(:,1))*norm(pos_out(:,end))));
err_dual = SigmaF-Phit_mat(:,:,end)*Phit_mat(:,:,end)';


resol = 100;
xo = zeros(2,resol);
yo = zeros(2,resol);
curv = linspace(0,2.*pi,resol);
v = linspace(0,1);
Theta = linspace(0,2*pi,1e2);

sigma__out_(:,:,1) = sigma_out(:,:,1);
sigma__out_(:,:,2) = sigma_out(:,:,end);



L1 =  reshape(0.5.*(sigma__out_(1,1,:)+sigma__out_(2,2,:))+...
           sqrt((0.5.*sigma__out_(1,1,:)-0.5.*sigma__out_(2,2,:)).^2+...
           sigma__out_(1,2,:).^2),1,[]);
L2 = reshape(0.5.*(sigma__out_(1,1,:)+sigma__out_(2,2,:))-...
           sqrt((0.5.*sigma__out_(1,1,:)-0.5.*sigma__out_(2,2,:)).^2+...
           sigma__out_(1,2,:).^2),1,[]);
theta = zeros(1,Nt);


for tt = 1:2
    if (sigma__out_(1,2,tt)==0 && sigma__out_(1,1,tt)>=sigma__out_(2,2,tt))
            theta(tt) = 0;
    elseif (sigma__out_(1,2,tt)==0 && sigma__out_(1,1,tt)<sigma_out(2,2,tt))
         theta(tt) = pi/2;
    else
        theta(tt) = atan2(L1(tt)-sigma__out_(1,1,tt),sigma__out_(1,2,tt));
    end
end




for k = 1:resol
    for i = 1:2
    xo(i,k) = sqrt(L1(i)).*cos(theta(i)).*cos(curv(k))-...
        sqrt(L2(i)).*sin(theta(i)).*sin(curv(k));
    yo(i,k) = sqrt(L1(i)).*sin(theta(i)).*cos(curv(k))+...
        sqrt(L2(i)).*cos(theta(i)).*sin(curv(k));
    end
end



Xo = zeros(length(sold.x),length(Theta));
Yo = zeros(length(sold.x),length(Theta));
Zo = zeros(length(sold.x),length(Theta));

for i =1:length(sold.x)
    for j = 1:length(Theta)
        s = Theta(j);
        Xo(i,j) = sold.x(i);
        Yo(i,j) = [1,0]*(sqrtm(sigma_out(:,:,i)))*[cos(s);sin(s)];
        Zo(i,j) = [0,1]*(sqrtm(sigma_out(:,:,i)))*[cos(s);sin(s)];
    end
end


figure,surf(Xo,Yo,Zo,[0 1 0 2*pi],'LineWidth',0.1,'FaceColor','k',...
    'FaceAlpha',0.1,'MeshStyle','none'), hold on;
plot3(0.*ones(resol,1),xo(1,:),yo(1,:),'k',LineWidth=2); hold on;
fill3(1.*ones(resol,1),xo(end,:),yo(end,:),'k',FaceAlpha=0.125,EdgeColor='none',LineWidth=2);hold on
plot3(1.*ones(resol,1),xo(end,:),yo(end,:),'k',LineWidth=2);hold on
plot3(1.*ones(50,1),xo(end,16:65),yo(end,16:65),'k',LineWidth=2); hold on;
plot3(1.*ones(15,1),xo(end,1:15),yo(end,1:15),'k--',LineWidth=2);hold on;
plot3(1.*ones(35,1),xo(end,66:100),yo(end,66:100),'k--',LineWidth=2);
hold on
    plot3(T,PosEval(1,:),PosEval(2,:),'Color',[0.4660 0.6740 0.1880],LineWidth=3);
   hold on;
   % plot3(sold.x,pos_out(1,:),pos_out(2,:),LineWidth=2.5);
   % hold on;
plot3(0.*ones(resol,1),linspace(0,pos_out(1,1),resol),...
    linspace(0,pos_out(2,1),resol),':','Color', [0.41 0.41 0.41],LineWidth=2); hold on;
plot3(ones(resol,1,1),linspace(0,pos_out(1,end),resol),....
    linspace(0,pos_out(2,end),resol),':','Color', [0.41 0.41 0.41],LineWidth=2); hold on;
plot3(linspace(0,1,resol),zeros(resol,1),zeros(resol,1),...
    ':','Color', [0.41 0.41 0.41],LineWidth=2);hold on;
ax = gca;
ax.DataAspectRatio=[0.125,1,1];
ax.FontSize = 15;
xlabel('Time \ t', 'FontSize',19 ,'interpreter','latex');
% view(-90,0)
