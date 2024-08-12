
% apply discrete Fourier transform along the \theta and z directions over the
% whole r.h.s. and boundary values

function [rhshat,valr1t]=ftransf2(gr,gt,gz,rhs,valr1)

nt=length(gt.xn);
nz=length(gz.xn);

% transform rhs (this must be done on inner points only!)
[Tp,Rp,Zp]=meshgrid(gt.xp,gr.xp,gz.xp);
rhshat=fft( fft(rhs(Rp,Tp,Zp),nt-1,2),nz-1,3);

% transform bcs, only the value matrix
[Zp,Tp]=meshgrid(gz.xp,gt.xp);
tmp=valr1(Tp,Zp);
valr1t=tmp;
valr1t(:,:,3)=fft( fft(squeeze(tmp(:,:,3)),nt-1,1),nz-1,2 );

end

