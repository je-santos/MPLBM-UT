function [Fy,Fx,Fz]=pointmin(I)
    
Fx=zeros(size(I),class(I));
Fy=zeros(size(I),class(I));
Fz=zeros(size(I),class(I));

J=zeros(size(I)+2,class(I));
J(:,:)=max(I(:));
if(ndims(I)==2)
    J(2:end-1,2:end-1)=I;
    Ne=[-1 -1; -1  0; -1  1; 0 -1; 0  1; 1 -1;  1  0; 1  1];
    for i=1:length(Ne);
       In=J(2+Ne(i,1):end-1+Ne(i,1),2+Ne(i,2):end-1+Ne(i,2));
       check = In<I;
       I(check)= In(check);
       D=Ne(i,:); D=D./sqrt(sum(D.^2));
       Fx(check)= D(1); Fy(check)= D(2);
    end
else
    J(2:end-1,2:end-1,2:end-1)=I;
    Ne=[-1 -1 -1; -1 -1  0; -1 -1  1; -1  0 -1; -1  0  0; -1  0  1; -1  1 -1; -1  1  0; -1  1  1;        
         0 -1 -1;  0 -1  0;  0 -1  1;  0  0 -1;            0  0  1;  0  1 -1;  0  1  0;  0  1  1; 
         1 -1 -1;  1 -1  0;  1 -1  1;  1  0 -1;  1  0  0;  1  0  1;  1  1 -1;  1  1  0;  1  1  1];
    
    for i=1:length(Ne);
       In=J(2+Ne(i,1):end-1+Ne(i,1),2+Ne(i,2):end-1+Ne(i,2),2+Ne(i,3):end-1+Ne(i,3));
       check = In<I;
       I(check)= In(check);
       D=Ne(i,:); D=D./sqrt(sum(D.^2));
       Fx(check)= D(1); Fy(check)= D(2); Fz(check)= D(3);
    end
end

