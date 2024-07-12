
% solve 1D Poisson problem by centered finite differences on a 
% cell-centered, non-uniform grid; b.c.s are prescribed by a ghost cell method
% Dirichelet-Neumann boundary conditions are imposed

% A. Nitti, Polytechnic University of Bari (2024)

clc
clear 
close all
kf=1;


%% input parameters
Lx=2;                                   % domain length
nx=40;                                  % number of nodes
pbc=false(1);                           % periodic boundary condition flag
ref='tanh-e';                           % grid refinement law
rhs=@(x) 2.*x./x;                       % right-hand-side function
val=[1,1;0,0;0,4];                      % value at boundaries {\alpha,\beta,\gamma}


%% create grid, assemble FD matrix
gx=getgrid(Lx,nx,ref,pbc);
ndof=nx-1;

A=zeros(ndof);
for i=2:nx-2
    A(i, i)     = -1/(gx.dxn(i)*gx.dxc(i+1)) -1/(gx.dxn(i)*gx.dxc(i));
    A(i, i+1)   =  1/(gx.dxn(i)*gx.dxc(i+1));
    A(i, i-1)   =  1/(gx.dxn(i)*gx.dxc(i));        
end

b=rhs(gx.xp');


%% prescribe boundary conditions
rhse=rhs(gx.xc);
[A,b]=bcs1D(A,b,gx,rhs,pbc,val(:,1),val(:,2));


%% solve linear system
u=A\b;


%% plot result
figure(kf); kf=kf+1;
ua=gx.xp.^2;    % analytical result
plot(gx.xp,ua,'o')
hold on
plot(gx.xp,u,'-x')
xlabel('x');
ylabel('u');
legend('analytical','numerical')


%% compute error 
err=(ua-u')./ua;
disp(['relative error rms: ',num2str(rms(err))])
