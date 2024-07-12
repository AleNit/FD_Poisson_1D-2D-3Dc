
% assign boundary conditions in the form:
%       \alpha*u + \beta*dudx = \gamma.
% coefficients are passed by column vectors
% alternatively, periodic boundary conditions can be prescribed


function [A,b]=bcs1D(A,b,gr,rhs,pbc,vali,vale)

nx=length(gr.xn);

if (pbc)    %------------------------------------------------ PERIODIC

    % initial boundary
    i=1;
    im=nx-1;
    A(i, i)     = -1/(gr.dxn(i)*gr.dxc(i+1)) -1/(gr.dxn(i)*gr.dxc(i));
    A(i, i+1)   =  1/(gr.dxn(i)*gr.dxc(i+1));
    A(i, im)    =  1/(gr.dxn(i)*gr.dxc(i));  

    % final boundary
    i=nx-1;
    ip=1;
    A(i, i)     = -1/(gr.dxn(i)*gr.dxc(i+1)) -1/(gr.dxn(i)*gr.dxc(i));
    A(i, ip)    =  1/(gr.dxn(i)*gr.dxc(i+1));
    A(i, i-1)   =  1/(gr.dxn(i)*gr.dxc(i));      

else        %------------------------------------------------ DNR

    % initial boundary
    i=1;
    ip=i+1;
    c1=(2*vali(2)+vali(1)*gr.dxc(i))/(2*vali(2)-vali(1)*gr.dxc(i));
    c2=-2*vali(3)*gr.dxc(i)/(2*vali(2)-vali(1)*gr.dxc(i));
    A(i, i)     = -1/(gr.dxn(i)*gr.dxc(i+1)) + (c1-1)/(gr.dxn(i)*gr.dxc(i));
    A(i, i+1)   =  1/(gr.dxn(i)*gr.dxc(i+1));
    b(i)        =  rhs(i) - c2/(gr.dxn(i)*gr.dxc(i));    

    % final boundary
    i=nx-1;
    ip=i+1;
    c3=(-2*vale(2)+vale(1)*gr.dxc(i+1))/(-2*vale(2)-vale(1)*gr.dxc(i+1));
    c4=-2*vale(3)*gr.dxc(i+1)/(-2*vale(2)-vale(1)*gr.dxc(i+1));
    A(i, i)     =  (c3-1)/(gr.dxn(i)*gr.dxc(i+1)) -1/(gr.dxn(i)*gr.dxc(i));
    A(i, i-1)   =  1/(gr.dxn(i)*gr.dxc(i));
    b(i)        =  rhs(i) - c4/(gr.dxn(i)*gr.dxc(i+1));        

end

end
