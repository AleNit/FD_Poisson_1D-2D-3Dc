
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
    den=(gr.dxc(1)-gr.dxn(1)/2)*gr.dxn(1)/2*gr.dxc(1);
    dp=(gr.dxc(1)-gr.dxn(1)/2)^2/den;
    dm=(gr.dxn(1)/2)^2/den;
    dd=( (gr.dxn(1)/2)^2 - (gr.dxc(1)-gr.dxn(1)/2)^2 )/den;
    c1=( vali(2)*dp + (2*gr.dxc(1)-gr.dxn(1))/2/gr.dxc(1)*( vali(1)+vali(2)*dd ) )/ ...
       ( vali(2)*dm + ( (2*gr.dxc(1)-gr.dxn(1))/2/gr.dxc(1)-1 )*( vali(1)+vali(2)*dd ) );
    c2= - vali(3) / ( vali(2)*dm + ( (2*gr.dxc(1)-gr.dxn(1))/2/gr.dxc(1)-1 )*( vali(1)+vali(2)*dd ) );
    A(i, i)     = -1/(gr.dxn(i)*gr.dxc(i+1)) + (c1-1)/(gr.dxn(i)*gr.dxc(i));
    A(i, i+1)   =  1/(gr.dxn(i)*gr.dxc(i+1));
    b(i)        =  rhs(i) - c2/(gr.dxn(i)*gr.dxc(i));    

    % final boundary
    i=nx-1;
    den=(gr.dxc(nx)-gr.dxn(nx-1)/2)*gr.dxn(nx-1)/2*gr.dxc(nx);
    dp=(gr.dxn(nx-1)/2)^2/den;
    dm=(gr.dxc(nx)-gr.dxn(nx-1)/2)^2/den;
    dd=( (gr.dxc(nx)-gr.dxn(nx-1)/2)^2 - (gr.dxn(nx-1)/2)^2 )/den;
    c3=( vale(2)*dm + ( gr.dxn(nx-1)/2/gr.dxc(nx)-1 )*( vale(1)+vale(2)*dd ) )/ ...
    ( vale(2)*dp + gr.dxn(nx-1)/2/gr.dxc(nx)*( vale(1)+vale(2)*dd ) );
    c4= vale(3) / ( vale(2)*dp + gr.dxn(nx-1)/2/gr.dxc(nx)*( vale(1)+vale(2)*dd ) );
    A(i, i)     =  (c3-1)/(gr.dxn(i)*gr.dxc(i+1)) -1/(gr.dxn(i)*gr.dxc(i));
    A(i, i-1)   =  1/(gr.dxn(i)*gr.dxc(i));
    b(i)        =  rhs(i) - c4/(gr.dxn(i)*gr.dxc(i+1));        

end

end
