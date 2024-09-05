
% set Dirichelet/Neumann/Periodic boundary conditions for a 2D Poisson problem; 
% corner points are addressed in the west and east boundary sections

function [A,b]=bcs2D(A,b,gx,gy,rhs,pbc,valn,vale,vals,valw)

nx=length(gx.xn);
ny=length(gy.xn);

if (all(pbc))
    error('... bi-periodic case not implemented')
end


%% compute boundary arrays
i=1;
den=(gx.dxc(1)-gx.dxn(1)/2)*gx.dxn(1)/2*gx.dxc(1);
dp=(gx.dxc(1)-gx.dxn(1)/2)^2/den;
dm=(gx.dxn(1)/2)^2/den;
dd=( (gx.dxn(1)/2)^2 - (gx.dxc(1)-gx.dxn(1)/2)^2 )/den;
c1we=( valw(:,2)*dp + (2*gx.dxc(1)-gx.dxn(1))/2/gx.dxc(1)*( valw(:,1)+valw(:,2)*dd ) )/ ...
     ( valw(:,2)*dm + ( (2*gx.dxc(1)-gx.dxn(1))/2/gx.dxc(1)-1 )*( valw(:,1)+valw(:,2)*dd ) );
c2we= - valw(:,3) / ( valw(:,2)*dm + ( (2*gx.dxc(1)-gx.dxn(1))/2/gx.dxc(1)-1 )*( valw(:,1)+valw(:,2)*dd ) );

i=nx-1;
den=(gx.dxc(nx)-gx.dxn(nx-1)/2)*gx.dxn(nx-1)/2*gx.dxc(nx);
dp=(gx.dxn(nx-1)/2)^2/den;
dm=(gx.dxc(nx)-gx.dxn(nx-1)/2)^2/den;
dd=( (gx.dxc(nx)-gx.dxn(nx-1)/2)^2 - (gx.dxn(nx-1)/2)^2 )/den;
c3we=( vale(:,2)*dm + ( gx.dxn(nx-1)/2/gx.dxc(nx)-1 )*( vale(1)+vale(:,2)*dd ) )/ ...
     ( vale(:,2)*dp + gx.dxn(nx-1)/2/gx.dxc(nx)*( vale(:,1)+vale(:,2)*dd ) );
c4we= vale(:,3) / ( vale(:,2)*dp + gx.dxn(nx-1)/2/gx.dxc(nx)*( vale(:,1)+vale(2)*dd ) );

j=1;
den=(gy.dxc(1)-gy.dxn(1)/2)*gy.dxn(1)/2*gy.dxc(1);
dp=(gy.dxc(1)-gy.dxn(1)/2)^2/den;
dm=(gy.dxn(1)/2)^2/den;
dd=( (gy.dxn(1)/2)^2 - (gy.dxc(1)-gy.dxn(1)/2)^2 )/den;
c1sn=( vals(:,2)*dp + (2*gy.dxc(1)-gy.dxn(1))/2/gy.dxc(1)*( vals(:,1)+vals(:,2)*dd ) )/ ...
     ( vals(:,2)*dm + ( (2*gy.dxc(1)-gy.dxn(1))/2/gy.dxc(1)-1 )*( vals(:,1)+vals(:,2)*dd ) );
c2sn= - vals(:,3) / ( vals(:,2)*dm + ( (2*gy.dxc(1)-gy.dxn(1))/2/gy.dxc(1)-1 )*( vals(:,1)+vals(:,2)*dd ) );

j=ny-1;
den=(gy.dxc(ny)-gy.dxn(ny-1)/2)*gy.dxn(ny-1)/2*gy.dxc(ny);
dp=(gy.dxn(ny-1)/2)^2/den;
dm=(gy.dxc(ny)-gy.dxn(ny-1)/2)^2/den;
dd=( (gy.dxc(ny)-gy.dxn(ny-1)/2)^2 - (gy.dxn(ny-1)/2)^2 )/den;
c3sn=( valn(:,2)*dm + ( gy.dxn(ny-1)/2/gy.dxc(ny)-1 )*( valn(:,1)+valn(:,2)*dd ) )/ ...
     ( valn(:,2)*dp + gy.dxn(ny-1)/2/gy.dxc(ny)*( valn(:,1)+valn(:,2)*dd ) );
c4sn= valn(:,3) / ( valn(:,2)*dp + gy.dxn(ny-1)/2/gy.dxc(ny)*( valn(:,1)+valn(:,2)*dd ) );




%% west-east boundaries
if (pbc(1))



    %-------------------------------------------------------- west boundary
    i=1;
    ip=i+1;

    % west-south corner
    j=1;
    jp=j+1;            
    c=(j-1)*(nx-1)+i;        
    cjp=j*(nx-1)+i;   
    cim=j*(nx-1);
    A(c, c)    = -1/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                 -1/(gy.dxn(j)*gy.dxc(jp)) +(c1sn(ip)-1)/(gy.dxn(j)*gy.dxc(j));
    A(c, c+1)    =  1/(gx.dxn(i)*gx.dxc(ip));
    A(c, cim)  =  1/(gx.dxn(i)*gx.dxc(i));    
    A(c, cjp)  =  1/(gy.dxn(j)*gy.dxc(jp));
    b(c)       =  rhs(ip,jp) - c2sn(ip)/(gy.dxn(j)*gy.dxc(j));

    % inner points
    for j=2:ny-2
        jp=j+1;            
        c=(j-1)*(nx-1)+i;        
        cjp=j*(nx-1)+i;
        cjm=(j-2)*(nx-1)+i;        
        cim=j*(nx-1);
        A(c, c)    = -1/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                     -1/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));
        A(c, c+1)    =  1/(gx.dxn(i)*gx.dxc(ip));
        A(c, cim)  =  1/(gx.dxn(i)*gx.dxc(i));    
        A(c, cjp)  =  1/(gy.dxn(j)*gy.dxc(jp));
        A(c, cjm)  =  1/(gy.dxn(j)*gy.dxc(j));   
        b(c)       =  rhs(ip,jp);                
    end    

    % west-north corner
    j=ny-1;     
    jp=j+1;
    c=(j-1)*(nx-1)+i;        
    cjm=(j-2)*(nx-1)+i;        
    cim=j*(nx-1);
    A(c, c)    = -1/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                 +(c3sn(ip)-1)/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));
    A(c, c+1)    =  1/(gx.dxn(i)*gx.dxc(ip));
    A(c, cim)  =  1/(gx.dxn(i)*gx.dxc(i));    
    A(c, cjm)  =  1/(gy.dxn(j)*gy.dxc(j));   
    b(c)       =  rhs(ip,jp) - c4sn(ip)/(gy.dxn(j)*gy.dxc(jp));     



    %-------------------------------------------------------- east boundary
    i=nx-1;
    ip=i+1;

    % east-south corner
    j=1;
    jp=j+1;            
    c=(j-1)*(nx-1)+i;        
    cjp=j*(nx-1)+i;
    cip=(j-1)*(nx-1)+1;
    A(c, c)    = -1/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                 -1/(gy.dxn(j)*gy.dxc(jp)) +(c1sn(ip)-1)/(gy.dxn(j)*gy.dxc(j));
    A(c, cip)    =  1/(gx.dxn(i)*gx.dxc(ip));
    A(c, c-1)  =  1/(gx.dxn(i)*gx.dxc(i));    
    A(c, cjp)  =  1/(gy.dxn(j)*gy.dxc(jp));
    b(c)       =  rhs(ip,jp) - c2sn(ip)/(gy.dxn(j)*gy.dxc(j));  

    % inner points
    for j=2:ny-2
        jp=j+1;            
        c=(j-1)*(nx-1)+i;        
        cjp=j*(nx-1)+i;
        cjm=(j-2)*(nx-1)+i;        
        cip=(j-1)*(nx-1)+1;
        A(c, c)    = -1/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                     -1/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));
        A(c, cip)    =  1/(gx.dxn(i)*gx.dxc(ip));
        A(c, c-1)  =  1/(gx.dxn(i)*gx.dxc(i));    
        A(c, cjp)  =  1/(gy.dxn(j)*gy.dxc(jp));
        A(c, cjm)  =  1/(gy.dxn(j)*gy.dxc(j));   
        b(c)       =  rhs(ip,jp);                
    end  

    % east-north corner
    j=ny-1;     
    jp=j+1;    
    c=(j-1)*(nx-1)+i;        
    cjm=(j-2)*(nx-1)+i;        
    cip=(j-1)*(nx-1)+1;
    A(c, c)    = -1/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                 +(c3sn(ip)-1)/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));
    A(c, cip)    =  1/(gx.dxn(i)*gx.dxc(ip));
    A(c, c-1)  =  1/(gx.dxn(i)*gx.dxc(i));    
    A(c, cjm)  =  1/(gy.dxn(j)*gy.dxc(j));   
    b(c)       =  rhs(ip,jp) - c4sn(ip)/(gy.dxn(j)*gy.dxc(jp));   

    


else




    %-------------------------------------------------------- west boundary
    i=1;            
    ip=i+1;

    % south-west corner
    j=1;
    jp=j+1;
    jm=ny-1;
    c=(j-1)*(nx-1)+i;        
    cjp=j*(nx-1)+i;
    if (pbc(2))
        cjm=(jm-1)*(nx-1)+i;
        A(c, c)    = -1/(gx.dxn(i)*gx.dxc(ip)) +(c1we(jp)-1)/(gx.dxn(i)*gx.dxc(i)) ...
                     -1/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));
        A(c, cjm)  =  1/(gy.dxn(j)*gy.dxc(j));   
        b(c)       =  rhs(ip,jp) - c2we(jp)/(gx.dxn(i)*gx.dxc(i));               
    else
        A(c, c)    =  -1/(gx.dxn(i)*gx.dxc(ip)) +(c1we(jp)-1)/(gx.dxn(i)*gx.dxc(i)) ...
                      -1/(gy.dxn(j)*gy.dxc(jp)) +(c1sn(ip)-1)/(gy.dxn(j)*gy.dxc(j));
        b(c)       = rhs(ip,jp) - c2we(jp)/(gx.dxn(i)*gx.dxc(i)) ...
                    - c2sn(ip)/(gy.dxn(j)*gy.dxc(j));            
    end
    A(c, c+1)  =  1/(gx.dxn(i)*gx.dxc(ip));  
    A(c, cjp)  =  1/(gy.dxn(j)*gy.dxc(jp));

    % west, inner points
    for j=2:ny-2
        jp=j+1;
        c=(j-1)*(nx-1)+i;        
        cjp=j*(nx-1)+i;
        cjm=(j-2)*(nx-1)+i;
        A(c, c)    = -1/(gx.dxn(i)*gx.dxc(ip)) +(c1we(jp)-1)/(gx.dxn(i)*gx.dxc(i)) ...
                     -1/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));
        A(c, c+1)  =  1/(gx.dxn(i)*gx.dxc(ip));  
        A(c, cjp)  =  1/(gy.dxn(j)*gy.dxc(jp));
        A(c, cjm)  =  1/(gy.dxn(j)*gy.dxc(j));   
        b(c)       =  rhs(ip,jp) - c2we(jp)/(gx.dxn(i)*gx.dxc(i));            
    end

    % west-north corner
    j=ny-1;
    jp=j+1;
    c=(j-1)*(nx-1)+i;
    cjm=(j-2)*(nx-1)+i;
    if (pbc(2))
        cjp=(nx-1)+i;
        A(c, c)    = -1/(gx.dxn(i)*gx.dxc(ip)) +(c1we(jp)-1)/(gx.dxn(i)*gx.dxc(i)) ...
                     -1/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));        
        A(c, cjp)  =  1/(gy.dxn(j)*gy.dxc(jp));        
        b(c)       =  rhs(ip,jp) - c2we(jp)/(gx.dxn(i)*gx.dxc(i));  
    else
        A(c, c)    = -1/(gx.dxn(i)*gx.dxc(ip)) +(c1we(jp)-1)/(gx.dxn(i)*gx.dxc(i)) ...
                      +(c3sn(ip)-1)/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));        
        b(c)       =  rhs(ip,jp) - c2we(jp)/(gx.dxn(i)*gx.dxc(i)) ...
                      - c4sn(ip)/(gy.dxn(j)*gy.dxc(jp));            
    end
    A(c, c+1)  =  1/(gx.dxn(i)*gx.dxc(ip));    
    A(c, cjm)  =  1/(gy.dxn(j)*gy.dxc(j));  



    %-------------------------------------------------------- east boundary
    i=nx-1;            
    ip=i+1;

    % south-east corner
    j=1;
    jp=j+1;
    jm=ny-1;
    c=(j-1)*(nx-1)+i;        
    cjp=j*(nx-1)+i;
    if (pbc(2))
        cjm=(jm-1)*(nx-1)+i; 
        A(c, c)    =  (c3we(jp)-1)/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                     -1/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));                
        A(c, cjm)  =  1/(gy.dxn(j)*gy.dxc(j));   
        b(c)       =  rhs(ip,jp) - c4we(jp)/(gx.dxn(i)*gx.dxc(ip));
    else
        A(c, c)    =  (c3we(jp)-1)/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                         -1/(gy.dxn(j)*gy.dxc(jp)) +(c1sn(ip)-1)/(gy.dxn(j)*gy.dxc(j));
        b(c)       =  rhs(ip,jp) - c4we(jp)/(gx.dxn(i)*gx.dxc(ip)) ...
                        - c2sn(ip)/(gy.dxn(j)*gy.dxc(j));            
    end
    A(c, c-1)  =  1/(gx.dxn(i)*gx.dxc(i));    
    A(c, cjp)  =  1/(gy.dxn(j)*gy.dxc(jp));

    % east, inner points
    for j=2:ny-2
        jp=j+1;     % coordinate index y
        c=(j-1)*(nx-1)+i;        
        cjp=j*(nx-1)+i;
        cjm=(j-2)*(nx-1)+i;
        A(c, c)    =  (c3we(jp)-1)/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                     -1/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));        
        A(c, c-1)  =  1/(gx.dxn(i)*gx.dxc(i));    
        A(c, cjp)  =  1/(gy.dxn(j)*gy.dxc(jp));
        A(c, cjm)  =  1/(gy.dxn(j)*gy.dxc(j));   
        b(c)       =  rhs(ip,jp) - c4we(jp)/(gx.dxn(i)*gx.dxc(ip));            
    end    

    % north-east corner
    j=ny-1;
    jp=j+1;     % coordinate index y
    c=(j-1)*(nx-1)+i;
    cjm=(j-2)*(nx-1)+i;
    if (pbc(2))
        cjp=(nx-1)+i;
        A(c, c)    =  (c3we(jp)-1)/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                     -1/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));                
        A(c, cjp)  =  1/(gy.dxn(j)*gy.dxc(jp));
        b(c)       =  rhs(ip,jp) - c4we(jp)/(gx.dxn(i)*gx.dxc(ip));                    
    else
        A(c, c)    =  (c3we(jp)-1)/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                      +(c3sn(ip)-1)/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));        
        b(c)       =  rhs(ip,jp) - c4we(jp)/(gx.dxn(i)*gx.dxc(ip)) ...
                      - c4sn(ip)/(gy.dxn(j)*gy.dxc(jp));
    end
    A(c, c-1)  =  1/(gx.dxn(i)*gx.dxc(i));    
    A(c, cjm)  =  1/(gy.dxn(j)*gy.dxc(j));



end






%% south-north boundaries
if (pbc(2))



    %-------------------------------------------------------- south boundary
    j=1;    
    jm=ny-1;
    for i=2:nx-2
        ip=i+1;
        c=(j-1)*(nx-1)+i;        
        cjp=j*(nx-1)+i;
        cjm=(jm-1)*(nx-1)+i;        
        A(c, c)    = -1/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                     -1/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));
        A(c, c+1)  =  1/(gx.dxn(i)*gx.dxc(ip));
        A(c, c-1)  =  1/(gx.dxn(i)*gx.dxc(i));    
        A(c, cjp)  =  1/(gy.dxn(j)*gy.dxc(jp));
        A(c, cjm)  =  1/(gy.dxn(j)*gy.dxc(j));   
        b(c)       =  rhs(ip,jp);        
    end 


    %-------------------------------------------------------- north boundary
    j=ny-1;
    jp=1;
    for i=2:nx-2
        ip=i+1;
        c=(j-1)*(nx-1)+i;        
        cjp=jp*(nx-1)+i;
        cjm=(j-2)*(nx-1)+i;        
        A(c, c)    = -1/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                     -1/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));
        A(c, c+1)  =  1/(gx.dxn(i)*gx.dxc(ip));
        A(c, c-1)  =  1/(gx.dxn(i)*gx.dxc(i));    
        A(c, cjp)  =  1/(gy.dxn(j)*gy.dxc(jp));
        A(c, cjm)  =  1/(gy.dxn(j)*gy.dxc(j));   
        b(c)       =  rhs(ip,jp);        
    end     



else



    %-------------------------------------------------------- south boundary
    j=1;
    jp=j+1;

    % south, inner points    
    for i=2:nx-2
        ip=i+1;
        c=(j-1)*(nx-1)+i;        
        cjp=j*(nx-1)+i;
        A(c, c)    = -1/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                     -1/(gy.dxn(j)*gy.dxc(jp)) +(c1sn(ip)-1)/(gy.dxn(j)*gy.dxc(j));
        A(c, c+1)  =  1/(gx.dxn(i)*gx.dxc(ip));
        A(c, c-1)  =  1/(gx.dxn(i)*gx.dxc(i));    
        A(c, cjp)  =  1/(gy.dxn(j)*gy.dxc(jp));
        b(c)       =  rhs(ip,jp) - c2sn(ip)/(gy.dxn(j)*gy.dxc(j));
    end     



    %------------------------------------------------------- north boundary
    j=ny-1;
    jp=j+1;

    % north, inner points     
    for i=2:nx-2
        ip=i+1;     % coordinate index x
        c=(j-1)*(nx-1)+i;        
        cjm=(j-2)*(nx-1)+i;
        A(c, c)    = -1/(gx.dxn(i)*gx.dxc(ip)) -1/(gx.dxn(i)*gx.dxc(i)) ...
                     +(c3sn(ip)-1)/(gy.dxn(j)*gy.dxc(jp)) -1/(gy.dxn(j)*gy.dxc(j));
        A(c, c+1)  =  1/(gx.dxn(i)*gx.dxc(ip));
        A(c, c-1)  =  1/(gx.dxn(i)*gx.dxc(i));   
        A(c, cjm)  =  1/(gy.dxn(j)*gy.dxc(j));   
        b(c)       =  rhs(ip,jp) - c4sn(ip)/(gy.dxn(j)*gy.dxc(jp));
    end 



end



end
