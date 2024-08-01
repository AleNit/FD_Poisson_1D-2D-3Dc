
% assemble coefficient matrix in sparse format

function A=getCoeffMatR(gr,ndof)

nr=length(gr.xn);

% store coefficients of the tridiagonal matrix into sparse arrays I,J,VA
I=zeros((nr-3)*3,1);
J=I;
VA=I;

q=1;
for i=2:nr-2

    ip=i+1;
    im=i-1;
    
    I(q)=i;
    J(q)=i;
    VA(q)=-gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip)) ...
          -gr.xn(i)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(i)); ...

    I(q+1)=i;
    J(q+1)=ip;
    VA(q+1)=gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip));

    I(q+2)=i;
    J(q+2)=im;
    VA(q+2)=gr.xn(i)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(i));

    q=q+3;

end

A=sparse(I,J,VA,ndof,ndof);

% assign fictitious boundary conditions, just to add to include elements in the sparse format
[A,~]=bcs1DR(A,ones(ndof,1),gr,ones(ndof,1),[1,0,0],[1,0,0]);

end
