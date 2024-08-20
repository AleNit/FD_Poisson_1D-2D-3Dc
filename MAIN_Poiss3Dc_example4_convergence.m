
% solve 3D Poisson problem on cylindrical coordinates
% by centered finite differences on a cell-centered, non-uniform grid;
% b.c.s are prescribed by a ghost cell method with generalized {\alpha,\beta,\gamma} coefficients;
% non-homogeneous Dirichelet boundary conditions are imposed
% Solve by DFT in the theta direction

% A. Nitti, Polytechnic University of Bari (2024)
% A. Fraccalvieri, Polytechnic University of Bari (2024)


%%  convergence study
clc
clear
close all
clearAllMemoizedCaches
kf=1;

%% input parameters
Lz=2;
Lr=0.5;

%% 3 grids ( fine, medium and coarse grids with constant grid refinement ratio r in each direction )
r = 2;
d = 1/10;
mr = ceil(Lr/d);      % Number of sequencing levels in dir.r
mz = ceil(Lz/d);      % Number of sequencing levels in dir.z
mt = ceil(pi/(asin(d/(2*Lr)))); % Number of sequencing levels in dir.t
% the number of grid points in each coordinate direction
Nz = 2.^[2,1,0].*mz+1;
Nr = 2.^[2,1,0].*mr+1;
Nt = 2.^[2,1,0].*mt+1;

U = {};
RMS1 = zeros(1,length(Nz));
RMS2 = zeros(1,length(Nz));
for i= 1: length(Nz)
    nz = Nz(i);
    nr = Nr(i);
    nt = Nt(i);

    pbz=false(1);                                                               % periodic z boundary
    refr='lin';                                                                 % grid law in dir. r
    refz='lin';                                                                 % grid law in dir. z
    rhs=@(r,t,z)  z.^2.*(pi.*cos(pi.*r)./(2.*r)-pi^2.*sin(pi.*r)./2) + sin(pi.*r);      % right-hand-side function

    valz0=@(t,r) cat(3,ones(nt-1,nr-1), zeros(nt-1,nr-1), r.*0);                % values at z=0 boundary
    valz1=@(t,r) cat(3,ones(nt-1,nr-1), zeros(nt-1,nr-1), (1/2)*Lz^2.*sin(pi.*r));                % values at z=Lz boundary
    valr0=@(t,z) cat(3,ones(nt-1,nz-1), zeros(nt-1,nz-1), z.*0);                % values at r=0 boundary
    valr1=@(t,z) cat(3,ones(nt-1,nz-1), zeros(nt-1,nz-1), (1/2).*z.^2.*sin(pi.*Lr));        % values at r=Lr boundary

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
    index = strcat('num_points_grid',num2str(i));
    U.(index) = u;
    ua=@(r,t,z) (1/2).*z.^(2).*sin(pi.*r);     % analytical solution


    toc

    figure(kf); kf=kf+1;
    nexttile
    t1=ceil(nt/8);  % plot solution at theta \approx 45°
    Un1=squeeze(U.(index)(:,t1,:))';
    [Rp,Zp]=meshgrid(gr.xp,gz.xp);
    Ua1=ua(Rp,gt.xp(t1),Zp);
    surf(Zp,Rp,Ua1);
    shading interp
    hold on
    surf(Zp,Rp,Un1,'FaceColor','none')
    hold on
    axis equal
    xlim([0,Lz]); ylim([0,Lr]);
    view(-50,20)
    xlabel('z','Interpreter','latex',FontSize=14);
    ylabel('r','Interpreter','latex',FontSize=14);
    zlabel('$\phi$','Interpreter','latex',FontSize=14);
    legend('analytical','numerical')
    title(strcat('z-r plane in theta = ', num2str(ceil(nt/t1)*2*pi)),'interpreter','latex','FontSize',14);

    view(-40,20)
    err1=(Ua1-Un1)./(Ua1+1);
    RMS1(i) = rms(err1,'all');
    disp(['z-r plane, relative rmse: ',num2str(RMS1(i))]);

    nexttile
    z1 = ceil(nz/2); % plot solution at  z = Lz/2;
    Un2=squeeze(U.(index)(:,:,z1))';
    [Rp,Tp]=meshgrid(gr.xp,gt.xp);
    Ua2=ua(Rp,Tp,gz.xp(z1));
    %surf(Tp,Rp,Ua2);
    surf(Rp.*cos(Tp), Rp.*sin(Tp),Ua2);
    shading interp
    hold on
    % surf(Tp,Rp,Un2);
    surf(Rp.*cos(Tp), Rp.*sin(Tp),Un2,'FaceColor','none')
    hold on
    axis equal
    zlabel('$\phi$','Interpreter','latex',FontSize=14);
    legend('analytical','numerical')
    title(strcat('$\theta$ -r plane in  z = ', num2str(ceil((z1/nz)*Lz))),'interpreter','latex','FontSize',14);
    err2 =(Ua2-Un2)./(Ua2+1);
    RMS2(i) = rms(err2,'all');
    disp(['theta-r plane, relative rmse: ',num2str( RMS2(i))]);

end


%% theta constant
Normalized_grid_spacing = 2*pi./(Nt-1);
Normalized_grid_spacing = Normalized_grid_spacing/(min(Normalized_grid_spacing));
p =  abs((log(RMS1(2))-log(RMS1(1)))/(log(Normalized_grid_spacing(2))-log(Normalized_grid_spacing(1))));

figure(kf);kf=kf+1;
loglog(Normalized_grid_spacing,RMS1,'-s','color','b');
hold on
y = 10^-4.*(Normalized_grid_spacing(1):1:Normalized_grid_spacing(end)).^2;
loglog(Normalized_grid_spacing(1):1:Normalized_grid_spacing(end),y,'--','color','r');
xlabel('$\Delta \theta / \Delta \theta_{min}$','interpreter','latex','fontsize',14)
ylabel('RMSE','interpreter','latex','fontsize',14);
legend(strcat('p=',num2str(p)),'slope=2', 'location','northwest',...
    'interpreter','latex','fontsize',12)
title('z-r plane','interpreter','latex','fontsize',14);


%% z constant
Normalized_grid_spacing = Lz./(Nz-1);
Normalized_grid_spacing = Normalized_grid_spacing/(min(Normalized_grid_spacing));
figure(kf);kf=kf+1;
p = abs((log(RMS2(2))-log(RMS2(1)))/(log(Normalized_grid_spacing(2))-log(Normalized_grid_spacing(1))));
loglog(Normalized_grid_spacing,RMS2,'-s','color','b');
hold on
y = 10^-4.*(Normalized_grid_spacing(1):1:Normalized_grid_spacing(end)).^2;
loglog(Normalized_grid_spacing(1):1:Normalized_grid_spacing(end),y,'--','color','r');
xlabel('$\Delta \theta / \Delta \theta_{min}$','interpreter','latex','fontsize',14)
ylabel('RMSE','fontsize',14,'interpreter','latex');
legend(strcat('p=',num2str(p)),'slope=2', 'location','northwest',...
    'interpreter','latex','fontsize',12)
title('$\theta$-r plane','interpreter','latex','fontsize',14);




