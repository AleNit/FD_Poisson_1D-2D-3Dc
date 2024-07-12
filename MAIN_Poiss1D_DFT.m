
% solve 1D Poisson problem by centered finite differences on a 
% cell-centered uniform grid with DFT (Discrete Fourier Transform) or DCT
% (Discrete Cosine transform)

% A. Nitti, Polytechnic University of Bari (2024)

clc
clear
close all


% create grid
N=32;
dx=1/N;
x=dx/2:dx:1-dx/2;

% formulate a solution and retrieve corresponding right-hand side
% ua=sin(2*pi.*x);                                % DFT case
% f=-4*pi^2.*sin(2*pi.*x);
ua=cos(pi.*x);                                % DCT case
f=-pi^2.*cos(pi.*x);

% modified wavenumber
kx=0:N-1;
% mx=2*(cos(2*pi.*kx/N)-1)./dx^2;                                % DFT case
mx=2*(cos(pi.*kx/N)-1)./dx^2;                                % DCT case

% apply DFT
% fhat=fft(f);                                % DFT case
fhat=dct(f);                                % DCT case

% solve
uhat=fhat./mx;
uhat(1)=0;

% u=real(ifft(uhat));                                % DFT case
u=idct(uhat);                                % DCT case

figure
plot(x,ua,'x')
hold on
plot(x,u,'o')
xlabel('x');
ylabel('\phi');
legend('analytical','numerical');

