
% solve 3D Poisson problem on cylindrical coordinates
% by centered finite differences on a cell-centered, non-uniform grid; 
% b.c.s are prescribed by a ghost cell method with generalized {\alpha,\beta,\gamma} coefficients; 
% Boundary conditions: r --> (Dirichelet-Neumann),  z --> (periodic)
% Solve by DFT in the theta direction

% A. Nitti, Polytechnic University of Bari (2024)

clc
clear 
close all
clearAllMemoizedCaches
kf=1;


%% input parameters
Lz=1;
Lr=1;
nz=40;
nr=30;
nt=64;
pbz=true(1);                                                               % periodic z boundary
refr='tanh-e';                                                              % grid law in dir. r
refz='lin';                                                                 % grid law in dir. z
rhs=@(r,t,z) pi.*sin(2*pi.*z).*( 2.*cos(pi/2.*r)-17*pi.*r.*sin(pi/2.*r) )./(4.*r);                                               % right-hand-side function
valz0=@(t,r) cat(3,ones(nt-1,nr-1), zeros(nt-1,nr-1), r.*0);                % values at z=0 boundary
valz1=@(t,r) cat(3,ones(nt-1,nr-1), zeros(nt-1,nr-1), r.*0);                % values at z=Lz boundary
valr0=@(t,z) cat(3,ones(nt-1,nz-1), zeros(nt-1,nz-1), z.*0);                % values at r=0 boundary
valr1=@(t,z) cat(3,zeros(nt-1,nz-1), ones(nt-1,nz-1), z.*0);                % values at r=Lr boundary


%% preliminary operations
% create modified wavenumber array
dt=2*pi/(nt-1);
kt=0:nt-2;
mwt=2/dt^2.*(cos(2*pi.*kt./(nt-1))-1);

% create and plot grid
gr=getgrid(Lr,nr,refr,false(1));
gt=getgrid(2*pi,nt,'lin',true(1));
gz=getgrid(Lz,nz,refz,pbz);
ndof=(nr-1)*(nz-1);

kf=plotgrid(gr.xn,gt.xn,gz.xn,Lr,Lz,kf);

% DFT of rhs and boundary conditions
[rhshat,valr0t,valr1t,valz0t,valz1t]=ftransf(gr,gt,gz,rhs,valr0,valr1,valz0,valz1);

% assemble coefficient matrix
A=getCoeffMat(gr,gz,ndof,pbz);


%% solve problem
tic

tmparr=repmat(gr.xp,[nz-1,1]);
tmparr=reshape(tmparr,1,[]);
Ier=sparse(1:ndof,1:ndof,1./tmparr.^2,ndof,ndof);

us=zeros(nr-1,nt-1,nz-1);
for j=1:nt-1

    % isolate rhs of the pentadiagonal problem    
    tmpmat=squeeze(rhshat(:,j,:))';
    b=reshape(tmpmat,[1,ndof])';
    rhse=tmpmat;

    % assign boundary conditions
    [AA,b]=bcs2Dhat(A,b,gz,gr,rhse,pbz, ...
                    squeeze(valr1t(j,:,:)), ...
                    squeeze(valz1t(j,:,:)), ...
                    squeeze(valr0t(j,:,:)), ...
                    squeeze(valz0t(j,:,:)) );
    
    % augment coefficient matrix with modified wavenumber
    AA=AA+mwt(j).*Ier;

    % solve the j-th pentadiagonal problem
    sol=AA\b;
    
    % local to global solution array
    tmparr=reshape(sol,[nz-1,1,nr-1]);
    us(:,j,:)=permute(tmparr,[3,2,1]);

end

% apply inverse DFT
u=real(ifft(us,nt-1,2));
u=u-mean(u(:));        % this is necessary for the chosen boundary conditions

toc

%% compute error norm and plot results
ua=@(r,t,z) sin(2*pi.*z).*sin(pi/2.*r);     % analytical solution

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

