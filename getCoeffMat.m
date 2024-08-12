
% assemble coefficient matrix in sparse format

function [A,b]=getCoeffMat(gr,gt,gz,ndof,bb,cA,cBr,cBt,cBz)

nr=length(gr.xn);
nt=length(gt.xn);
nz=length(gz.xn);

% store coefficients of the pentadiagonal matrix into sparse arrays I,J,VA
I=zeros((nr-3)*(nz-3)*5,1);
J=I;
VA=I;

q=1;
for i=2:nr-2
    ip=i+1;
    for k=2:nz-2
        kp=k+1;
        c=(i-1)*(nz-1)+k;   
        cjp=i*(nz-1)+k;
        cjm=(i-2)*(nz-1)+k;        

        I(q)=c;
        J(q)=c;
        VA(q)=-gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip)) ...
              -gr.xn(i)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(i)) ...
              -1/(gz.dxn(k)*gz.dxc(kp)) -1/(gz.dxn(k)*gz.dxc(k));

        I(q+1)=c;
        J(q+1)=c+1;
        VA(q+1)=1/(gz.dxn(k)*gz.dxc(kp));

        I(q+2)=c;
        J(q+2)=c-1;
        VA(q+2)=1/(gz.dxn(k)*gz.dxc(k));

        I(q+3)=c;
        J(q+3)=cjp;
        VA(q+3)=gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip));

        I(q+4)=c;
        J(q+4)=cjm;
        VA(q+4)=gr.xn(i)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(i));

        q=q+5;

    end
end

A=sparse(I,J,VA,ndof,ndof);

% assign fictitious boundary conditions, just to add to include elements in the sparse format
[A,~]=bcs2Dhat(A,ones(ndof,1),gz,gr,ones(nz-1,nr-1),pbz, ...
                ones(nz-1,3),ones(nr-1,3),ones(nz-1,3),ones(nr-1,3) );

end
