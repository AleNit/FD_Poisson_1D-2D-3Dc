
% assign boundary conditions in the form:
%       \alpha*u + \beta*dudx = \gamma.
% This works for the tridiagonal matrices within the cylindrical Poisson solver


function [A,b]=bcs1DR(A,b,gr,rhs,vale)

nx=length(gr.xn);


% initial boundary (the unknown at the inner ghost point cancels out!)
i=1;
ip=i+1;
A(i, i)     = -gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip));
A(i, ip)    =  gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip));
b(i)        =  rhs(i);


% final boundary
i=nx-1;
ip=i+1;
c3=(-2*vale(2)+vale(1)*gr.dxc(i+1))/(-2*vale(2)-vale(1)*gr.dxc(i+1));
c4=-2*vale(3)*gr.dxc(i+1)/(-2*vale(2)-vale(1)*gr.dxc(i+1));
A(i, i)     =  (c3-1)*gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip)) - ...
               gr.xn(i)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(i));
A(i, i-1)   =  gr.xn(i)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(i));
b(i)        =  rhs(i) - c4*gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip));        


end