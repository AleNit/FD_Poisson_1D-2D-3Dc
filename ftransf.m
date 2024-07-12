
% apply discrete Fourier transform along the \theta direction over the
% whole r.h.s. and boundary values

function [rhshat,valr0t,valr1t,valz0t,valz1t]=ftransf(gr,gt,gz,rhs,valr0,valr1,valz0,valz1)

nt=length(gt.xn);

% transform rhs (this must be done on inner points only!)
[Tp,Rp,Zp]=meshgrid(gt.xp,gr.xp,gz.xp);
rhshat=fft( rhs(Rp,Tp,Zp),nt-1,2 );

% transform bcs, only the value matrix
[Zp,Tp]=meshgrid(gz.xp,gt.xp);
tmp=valr0(Tp,Zp);
valr0t=tmp;
valr0t(:,:,3)=fft( squeeze(tmp(:,:,3)),nt-1,1 );

tmp=valr1(Tp,Zp);
valr1t=tmp;
valr1t(:,:,3)=fft( squeeze(tmp(:,:,3)),nt-1,1 );

[Rp,Tp]=meshgrid(gr.xp,gt.xp);
tmp=valz0(Tp,Rp);
valz0t=tmp;
valz0t(:,:,3)=fft( squeeze(tmp(:,:,3)),nt-1,1 );

tmp=valz1(Tp,Rp);
valz1t=tmp;
valz1t(:,:,3)=fft( squeeze(tmp(:,:,3)),nt-1,1 );

end

