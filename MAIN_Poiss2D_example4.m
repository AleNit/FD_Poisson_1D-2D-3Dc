
% solve 2D Poisson problem by centered finite differences on a 
% cell-centered, non-uniform grid; b.c.s are prescribed by a ghost cell method
% Dirichelet-Periodic boundary conditions are imposed.
% Use the DFT in the periodic direction

% A. Nitti, Polytechnic University of Bari (2024)

clc
clear 
close all
kf=1;


%% input parameters
Lx=1;                                           % domain length
Ly=1;                                           % domain length
nx=50;                                          % number of nodes in x dir.
ny=40;                                          % number of nodes in y dir.
pbc=[true(1),false(1)];                         % periodic boundary condition flag [w-e,s-n]
refx='lin';                                     % grid refinement law x dir.
refy='lin';                                  % grid refinement law y dir.
rhs=@(x,y) -x.*sin(5*pi.*y)-exp(-((x-0.5).^2+(y-0.5).^2)./0.02);              % right-hand-side function
valw=@(y) [ones(size(y)), zeros(size(y)), y.*0];            % values at west boundary {\alpha,\beta,\gamma}
vale=@(y) [ones(size(y)), zeros(size(y)), y.*0];        % values at east boundary {\alpha,\beta,\gamma}    
vals=@(x) [ones(size(x)), zeros(size(x)), x.*0];            % values at south boundary {\alpha,\beta,\gamma}    
valn=@(x) [ones(size(x)), zeros(size(x)), x.*0];            % values at north boundary {\alpha,\beta,\gamma}


%% preliminary operations
% create modified wavenumber array
dx=Lx/(nx-1);
kx=0:nx-2;
mwx=2.*(cos(2*pi.*kx./(nx-1))-1)./dx^2;

% create grid
gx=getgrid(Lx,nx,refx,pbc(1));
gy=getgrid(Ly,ny,refy,pbc(2));
ndof=(ny-1);

% transform rhs (this must be done on inner points only!)
[Yp,Xp]=meshgrid(gy.xp,gx.xp);
rhshat=fft(rhs(Xp,Yp),nx-1,1);

% assemble coefficient matrix
A=zeros(ndof);
for j=2:ny-2
    jp=j+1;
    A(j, j)     = -1/(gy.dxn(j)*gy.dxc(j+1)) -1/(gy.dxn(j)*gy.dxc(j));
    A(j, j+1)   =  1/(gy.dxn(j)*gy.dxc(j+1));
    A(j, j-1)   =  1/(gy.dxn(j)*gy.dxc(j));                  
end


%% solve problem
us=zeros(nx-1,ny-1);
for i=1:nx-1

    ip=i+1;

    % isolate rhs of the tridiagonal problem
    b=zeros(ndof,1);
    b(2:ny-2) = rhshat(i,2:ny-2); 

    % apply boundary conditions in the non-periodic direction    
    [AA,b]=bcs1D(A,b,gy,rhshat(i,:),pbc,vals(0)',valn(0)');

    % augment coefficient matrix with modified wavenumber
    AA = AA + mwx(i).*eye(ndof);

    % solve tridiagonal problem
    tmp=AA\b;
    tmp(1)=0;
    us(i,:)=tmp;  

end

% apply inverse DFT
u=real( ifft(us,nx-1,1) );


%% plot result
figure(kf); kf=kf+1;
[Xp,Yp]=meshgrid(gx.xp,gy.xp);
surf(Xp,Yp,u','FaceColor','none')
shading interp
hold on
% axis equal
xlabel('x');
ylabel('y');
zlabel('\phi');
legend('numerical')



