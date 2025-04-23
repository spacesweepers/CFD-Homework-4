function [vg,itergr2,residue2]= gauss_seidel_v(v,dx,dy,a1,a2,a3,a4,a5,R2)
%%gauss-seidel sor for u
%%inputs for the code

w=1;
vg=v;
v_old=vg;
%%error and tolerance
residue2=1;
tol=10^-6;
itergr2=0;
 
%%iteration loop
while(residue2>tol)
     for j=2:size(dy,2)-1
         for i=2:size(dx,2)-1
             vg(j,i)= (1-w)*v_old(j,i)+w*((-a1(j,i)*v_old(j,i+1)-a2(j,i)*vg(j,i-1)-a4(j,i)*vg(j-1,i)-a5(j,i)*v_old(j+1,i)+R2(j,i))/a3(j,i));
         end
     end
    
    % vg(:,size(dx,2))=vg(:,size(dx,2)-1);
    % vg(1,:)=vg(2,:);
    % vg(size(dy,2),:)=vg(size(dy,2)-1,:);

     res=zeros([size(dy,2),size(dx,2)]); 
     for j=2:size(dy,2)-1
         for i=2:size(dx,2)-1
             res(j,i)= a1(j,i)*vg(j,i+1)+a2(j,i)*vg(j,i-1)+a3(j,i)*vg(j,i)+a4(j,i)*vg(j-1,i)+a5(j,i)*vg(j+1,i)-R2(j,i);
         end
     end
 residue2= mean(mean(abs(res))); 
 v_old=vg;
 itergr2=itergr2+1;
end
end


