function [ug,itergr1,residue1]= gauss_seidel_u(u,dx,dy,a1,a2,a3,a4,a5,R1)
%%gauss-seidel sor for u
%%inputs for the code

w=1;
ug=u;
u_old=ug;
%%error and tolerance
residue1=1;
tol=10^-6;
itergr1=0;
 
%%iteration loop
while(residue1>tol)
     for j=2:size(dy,2)-1
         for i=2:size(dx,2)-1
             ug(j,i)= (1-w)*u_old(j,i)+w*((-a1(j,i)*u_old(j,i+1)-a2(j,i)*ug(j,i-1)-a4(j,i)*ug(j-1,i)-a5(j,i)*u_old(j+1,i)+R1(j,i))/a3(j,i));
         end
     end
   
    % ug(:,size(dx,2))=ug(:,size(dx,2)-1);
    % ug(1,:)=ug(2,:);
    % ug(size(dy,2),:)=ug(size(dy,2)-1,:);

     res=zeros([size(dy,2),size(dx,2)]); 
     for j=2:size(dy,2)-1
         for i=2:size(dx,2)-1
             res(j,i)= a1(j,i)*ug(j,i+1)+a2(j,i)*ug(j,i-1)+a3(j,i)*ug(j,i)+a4(j,i)*ug(j-1,i)+a5(j,i)*ug(j+1,i)-R1(j,i);
         end
     end
 residue1= mean(mean(abs(res))); 
 u_old=ug;
 itergr1=itergr1+1;
end
end

