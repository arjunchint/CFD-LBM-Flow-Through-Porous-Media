%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%
clear
Reiter=1;
for Re=[.1 1 10 100]

ny=200;
nx = 4*ny;
Re = 100;
U = .1;
%%%%%%%%%%%%%%%% IMMERSED OBJECT %%%%%%%%%%%%%%%
%Cylinder
r=ny/10+1;
c_x = nx/5+1;   
c_y = ny/2+3;  
%obst = zeros(ny,nx); 
[y,x] = meshgrid(1:ny,1:nx); 
circ = (x-c_x).^2 + (y-c_y).^2 <= r.^2;
obst= circ;
%[rr, cc] = meshgrid(1:ny);
%obst = sqrt((rr-(floor((nx-1)/2))).^2+(cc-(floor((ny-1)/2))).^2)<=r;
obst(:,[1,ny]) = 1;
immersed = find(obst);
circi = find(circ);
%%%%%%%%%%%%%%%% OUTER BOUNDS %%%%%%%%%%%%%%%%%
L = ny-2; 
H = y-nx/ny;
u = 4*U/(L^2)*(H*L-H.^2);
v = zeros(nx,ny);
%%%%%%%%%%%%%%%% INSTANTIATION %%%%%%%%%%%%%%%%%
tau = 0.5+3*U*ny/Re;
%tau = 0.5+6*U*r/Re;
w = [4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36];
cx=[0 1 0 -1 0 1 -1 -1 1];
cy=[0 0 1 0 -1 1 1 -1 -1];
oppdir=[1 4 5 2 3 8 9 6 7];
c=[cx;cy]';
% u = zeros(ny,nx);
% v = zeros(ny,nx);
%rho = ones(ny,nx);
rho = 1;
feq=zeros(9,nx,ny);
f=zeros(9,nx,ny);
fnx=zeros(9,nx,ny);
for i=1:9 
    feq(i,:,:)=w(i)*rho.*(1+3*(c(i,1)*u+c(i,2)*v)+ 9/2*(c(i,1)*u+c(i,2)*v).^2 - 3/2*(u.^2+v.^2)); 
end;
f = feq;

iter=1;
uerror=100;
verror=100;
%%%%%%%%%%%%%%%%%%%%%%%% START ITERATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(iter<1000 || uerror>.001 || verror >.001)
uprev=reshape(u,nx,ny);
vprev=reshape(v,nx,ny);
%%%%%%%%%%%%%%%%%%%%%%%% MACRO VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
rho = sum(f);
% rhou=f(1,:,:)*c(1,1);
% for i=2:9
%     rhou=rhou+f(i,:,:)*c(i,1);
% end;
% rhov=f(1,:,:)*c(1,2);
% for i=2:9
%     rhov=rhov+f(i,:,:)*c(i,2);
% end;
% u=squeeze(rhou./rho);
% v=squeeze(rhov./rho);
% rho=squeeze(rho);    
u=reshape((cx * reshape(f,9,nx*ny)),1,nx,ny)./rho;
v=reshape((cy * reshape(f,9,nx*ny)),1,nx,ny)./rho;

%%%%%% BOUNDARIES
% INLET
H =[2:(ny-1)]-1.5;
u(:,1,2:(ny-1))=4*U/(L^2)*(H*L-H.^2);
v(:,1,[2:(ny-1)])= 0;
rho(:,1,[2:(ny-1)])=1./(1-u(:,1,[2:(ny-1)])).*(sum(f([1,3,5],1,[2:(ny-1)]))+2*sum(f([4,7,8],1,[2:(ny-1)])));
f(2,1,[2:(ny-1)])=f(4,1,[2:(ny-1)])+2/3*rho(:,1,[2:(ny-1)]).*u(:,1,[2:(ny-1)]); 
f(6,1,[2:(ny-1)])=f(8,1,[2:(ny-1)])+1/2*(f(5,1,[2:(ny-1)])-f(3,1,[2:(ny-1)]))+ 1/2*rho(:,1,[2:(ny-1)]).*v(:,1,[2:(ny-1)]) ...
                                + 1/6*rho(:,1,[2:(ny-1)]).*u(:,1,[2:(ny-1)]); 
f(9,1,[2:(ny-1)])=f(7,1,[2:(ny-1)])+1/2*(f(3,1,[2:(ny-1)])-f(5,1,[2:(ny-1)]))-1/2*rho(:,1,[2:(ny-1)]).*v(:,1,[2:(ny-1)]) ...
                                + 1/6*rho(:,1,[2:(ny-1)]).*u(:,1,[2:(ny-1)]); 

% OUTLET
rho(:,nx,[2:(ny-1)])=1;
u(:,nx,[2:(ny-1)])=-1+1./(rho(:,nx,[2:(ny-1)])).*(sum(f([1,3,5],nx,[2:(ny-1)]))+2*sum(f([2,6,9],nx,[2:(ny-1)])));
v(:,nx,[2:(ny-1)])=0;
f(4,nx,[2:(ny-1)])=f(2,nx,[2:(ny-1)])-2/3*rho(:,nx,[2:(ny-1)]).*u(:,nx,[2:(ny-1)]); 
f(8,nx,[2:(ny-1)])=f(6,nx,[2:(ny-1)])+1/2*(f(3,nx,[2:(ny-1)])-f(5,nx,[2:(ny-1)]))- 1/2*rho(:,nx,[2:(ny-1)]).*v(:,nx,[2:(ny-1)]) ...
                                  - 1/6*rho(:,nx,[2:(ny-1)]).*u(:,nx,[2:(ny-1)]); 
f(7,nx,[2:(ny-1)]) = f(9,nx,[2:(ny-1)]) + 1/2*(f(5,nx,[2:(ny-1)])-f(3,nx,[2:(ny-1)]))+ 1/2*rho(:,nx,[2:(ny-1)]).*v(:,nx,[2:(ny-1)]) ...
                                  - 1/6*rho(:,nx,[2:(ny-1)]).*u(:,nx,[2:(ny-1)]); 
%%%%%%%%%%%%%%%%%%%% F EQIILIBRIUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
collOperator=(1-1/tau)*f+1/tau*feq; 
for i=1:9 
feq(i,:,:)=w(i)*rho.*(1+3*(c(i,1)*u+c(i,2)*v)+ 9/2*(c(i,1)*u+c(i,2)*v).^2 - 3/2*(u.^2+v.^2)); 
fnx(i,:,:)=f(i,:,:)-(1/tau).*(f(i,:,:)-feq(i,:,:));
end
%%%%%%%%%%%%%%%%%%%% COLLISION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BOUNCE-BACK
for i=1:9
fnx(i,immersed)=f(oppdir(i),immersed);
end
%%%%%%%%%%%%%%%%% PROPAGATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:9
f(i,:,:)=circshift(fnx(i,:,:),[0,cx(i),cy(i)]);
end
%%%%%%%%%%%%%%%%%% DRAG FORCE %%%%%%%%%%%%%%%%%
if mod(iter,25) == 0
        F = (f(:,circi)+feq(:,circi));
        %         F(:,2:end-1,2:end-1) = 0;
        for ii = 1:9 
            Fx_mat(ii,:,:) = ((F(ii,:,:).*c(ii,1)));
            Fy_mat(ii,:,:) = ((F(ii,:,:).*c(ii,2)));
        end
        Fx(ceil(iter/25)) = sum(sum(sum(Fx_mat)));
        Fy(ceil(iter/25)) = sum(sum(sum(Fy_mat)));
    end
%%%%%%%%%%%%%%%%%%% STABILITIY CRITERIA %%%%%%%%%%%%%%%%%%%%%%%  
uerror=norm((reshape(u,nx,ny)-uprev)/nx/ny,'fro');
verror=norm((reshape(u,nx,ny)-uprev)/nx/ny,'fro');
iter=iter+1
end
%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
finError(Reiter)=uerror;
Reiter=Reiter+1;

x=.5/nx:1/nx:1;
y=.5/ny:1/ny:1;
u(immersed) = nan;
v(immersed) = nan;
UX=reshape(u,nx,ny);
UY=reshape(v,nx,ny);
figure;streamslice(x,y,UX',UY')
title(sprintf('Streams & Re = %g',Re))
xlabel('x')
ylabel('y')
axis tight

[X,Y]=meshgrid(x,y);
[sx,sy] = meshgrid(.05,.05:.05:.95);
figure; 
streamline(stream2(X,Y,UX',UY',sx,sy))
title(sprintf('Streamlines & Re = %g',Re))
xlabel('x')
ylabel('y')

% Streamlines
figure;contourf(x,y,hypot(UX',UY'))
title(sprintf('Contours & Re = %g',Re))
Cd(Reiter)=max(Fx(2:end));
end
