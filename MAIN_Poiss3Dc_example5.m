
% solve 3D Poisson problem on cylindrical coordinates
% by centered finite differences on a cell-centered, non-uniform grid; 
% b.c.s are prescribed by a ghost cell method with generalized {\alpha,\beta,\gamma} coefficients; 
% Boundary conditions: r --> (inner-Dirichelet),  z --> (periodic)
% Solve by DFT in the theta and z directions

% A. Nitti, Polytechnic University of Bari (2024)

clc
clear 
close all
clearAllMemoizedCaches
kf=1;


%% input parameters
Lz=2;                                                                       % axial domain length
Lr=0.5;                                                                     % radial domain length
nz=101;                                                                     % nodes in the axial direction
nr=31;                                                                      % nodes in the radial direction
nt=64;                                                                      % nodes in the tangential direction
refr='tanh-e';                                                              % grid law in dir. z
ua=@(r,t,z) sin(pi.*z).*cos(pi/2.*r);                                       % analytical solution
rhs=@(r,t,z) -(pi.*sin(pi.*z).*(2.*sin((pi.*r)./2) + ...
             5*pi.*r.*cos((pi.*r)./2)))./(4.*r);                            % right-hand-side function
valr1=@(t,z) cat(3,ones(nt-1,nz-1), zeros(nt-1,nz-1), sin(pi.*z).*sin(pi/4));  % values at r=Lr boundary


%% preliminary operations
% create modified wavenumber arrays
dt=2*pi/(nt-1);
kt=0:nt-2;
mwt=2/dt^2.*(cos(2*pi.*kt./(nt-1))-1);

dz=Lz/(nz-1);
kz=0:nz-2;
mwz=2/dz^2.*(cos(2*pi.*kz./(nz-1))-1);

% create and plot grid
gr=getgrid(Lr,nr,refr,false(1));
gt=getgrid(2*pi,nt,'lin',true(1));
gz=getgrid(Lz,nz,'lin',true(1));
ndof=nr-1;

kf=plotgrid(gr.xn,gt.xn,gz.xn,Lr,Lz,kf);

% DFT of rhs and boundary conditions
[rhshat,valr1t]=ftransf2(gr,gt,gz,rhs,valr1);

% assemble coefficient matrix
A=getCoeffMatR(gr,ndof);


%% solve problem
tic

Ier=sparse(1:ndof,1:ndof,1./gr.xp.^2,ndof,ndof);
Iez=speye(nr-1);

us=zeros(nr-1,nt-1,nz-1);
for k=1:nz-1
    for j=1:nt-1
    
        % isolate rhs of the tridiagonal problem    
        b=squeeze(rhshat(:,j,k));
    
        % assign boundary conditions
        [AA,b]=bcs1DR(A,b,gr,b,valr1t(j,k,:));
        
        % augment coefficient matrix with modified wavenumber
        AA = AA + mwt(j).*Ier + mwz(k).*Iez;
    
        % solve the j-th pentadiagonal problem
        sol=AA\b;
    
        % local to global solution array
        us(:,j,k)=sol;
        
    end
end

% apply inverse DFT
u=real( (ifft(ifft(us,nz-1,3),nt-1,2)) );

toc

%% compute error norm and plot results
figure(kf); kf=kf+1;
set(gcf,'Position',[100,100,1000,450])
tiledlayout(1,2);

t1=ceil(nt/8);  % plot solution at theta \approx 45°
nexttile
[Rp,Zp]=meshgrid(gr.xp,gz.xp);
Ua1=ua(Rp,gt.xp(t1),Zp);
surf(Zp,Rp,Ua1)
shading interp
hold on
Un1=squeeze(u(:,t1,:))';
surf(Zp,Rp,Un1,'FaceColor','none')
hold on
axis equal
xlim([0,Lz]); ylim([0,Lr]); 
view(-50,20)
xlabel('z','Interpreter','latex',FontSize=14);
ylabel('r','Interpreter','latex',FontSize=14);
zlabel('$\phi$','Interpreter','latex',FontSize=14);
ang=gt.xp(t1)*180/pi;
tit=strcat('$\theta=',num2str(ang),'^{\circ}$');
title(tit,'Interpreter','latex',FontSize=14)
legend('analytical','numerical')
view(-40,20)

err1=(Ua1-Un1)./(Ua1+1);
disp(['z-r plane 1, relative rmse: ',num2str(rms(err1,'all'))])



t2=ceil(3*nt/8);  % plot solution at theta \approx 135°
nexttile
[Rp,Zp]=meshgrid(gr.xp,gz.xp);
Ua2=ua(Rp,gt.xp(t2),Zp);
surf(Zp,Rp,Ua2)
shading interp
hold on
Un2=squeeze(u(:,t2,:))';
surf(Zp,Rp,Un2,'FaceColor','none')
hold on
axis equal
xlim([0,Lz]); ylim([0,Lr]); 
view(-50,20)
xlabel('z','Interpreter','latex',FontSize=14);
ylabel('r','Interpreter','latex',FontSize=14);
zlabel('$\phi$','Interpreter','latex',FontSize=14);
ang=gt.xp(t2)*180/pi;
tit=strcat('$\theta=',num2str(ang),'^{\circ}$');
title(tit,'Interpreter','latex',FontSize=14)
legend('analytical','numerical')
view(-40,20)

err2=(Ua2-Un2)./(Ua2+1);
disp(['z-r plane 2, relative rmse: ',num2str(rms(err2,'all'))])

