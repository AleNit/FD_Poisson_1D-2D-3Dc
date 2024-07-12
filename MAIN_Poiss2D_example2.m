
% solve 2D Poisson problem by centered finite differences on a 
% cell-centered, non-uniform grid; b.c.s are prescribed by a ghost cell method
% Dirichelet-Neumann boundary conditions are imposed

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
pbc=false(2,1);                                   % periodic boundary condition flag
refx='tanh-i';                                     % grid refinement law x dir.
refy='tanh-i';                                  % grid refinement law y dir.
rhs=@(x,y) -2.*(2.*y.^3-3.*y.^2+1)+6.*(1-x.^2).*(2.*y-1);              % right-hand-side function
valw=@(y) [ones(size(y)), zeros(size(y)), 2.*y.^3-3.*y.^2+1];            % values at west boundary {\alpha,\beta,\gamma}
vale=@(y) [ones(size(y)), zeros(size(y)), y.*0];        % values at east boundary {\alpha,\beta,\gamma}    
vals=@(x) [zeros(size(x)), ones(size(x)), x.*0];            % values at south boundary {\alpha,\beta,\gamma}    
valn=@(x) [zeros(size(x)), ones(size(x)), x.*0];            % values at north boundary {\alpha,\beta,\gamma}



%% create grid, assemble FD matrix
gx=getgrid(Lx,nx,refx,pbc(1));
gy=getgrid(Ly,ny,refy,pbc(2));
ndof=(nx-1)*(ny-1);

A=eye(ndof);
b=zeros(ndof,1);
for j=2:ny-2
    jp=j+1;
    for i=2:nx-2
        ip=i+1;
        c=(j-1)*(nx-1)+i;        
        cjp=j*(nx-1)+i;
        cjm=(j-2)*(nx-1)+i;        
        A(c, c)    = -1/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                     -1/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));
        A(c, c+1)  =  1/(gx.dxn(i)*gx.dxc(ip));
        A(c, c-1)  =  1/(gx.dxn(i)*gx.dxc(i));    
        A(c, cjp)  =  1/(gy.dxn(j)*gy.dxc(jp));
        A(c, cjm)  =  1/(gy.dxn(j)*gy.dxc(j));   
        b(c)       =  rhs(gx.xc(ip),gy.xc(jp));        
    end
end


%% prescribe boundary conditions
[Yc,Xc]=meshgrid(gy.xc,gx.xc);
rhse=rhs(Xc,Yc);
[A,b]=bcs2D(A,b,gx,gy,rhse,pbc,valn(gx.xc'),vale(gy.xc'),vals(gx.xc'),valw(gy.xc'));


%% solve linear system
us=A\b;
u=reshape(us,[nx-1,ny-1]);


%% plot result
figure(kf); kf=kf+1;
[Xp,Yp]=meshgrid(gx.xp,gy.xp);
ua=(1-Xp.^2).*(2.*Yp.^3-3.*Yp.^2+1);    % analytical solution
surf(Xp,Yp,ua)
shading interp
hold on
surf(Xp,Yp,u','FaceColor','none')
xlabel('x');
ylabel('y');
zlabel('\phi');
legend('analytical','numerical')


%% compute error
err=(ua-u')./(ua+1);
disp(['relative error rms: ',num2str(rms(err,'all'))])

