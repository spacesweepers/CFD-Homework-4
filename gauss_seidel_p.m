function [pg,itergr3,residue3]= gauss_seidel_p(p,dx,dy,c1,c2,c3,c4,c5,R3)
%%gauss-seidel sor for p
%%inputs for the code

w=1.5;
pg=p;
p_old=pg;
%%error and tolerance
residue3=1;
tol=10^-3;
itergr3=0;

%%iteration loop
while(residue3>tol)
%     itergr3
     for j=2:size(dy,2)-1
         for i=2:size(dx,2)-1
             pg(j,i)= (1-w)*p_old(j,i)+w*((-c1(j,i)*p_old(j,i+1)-c2(j,i)*pg(j,i-1)-c4(j,i)*pg(j-1,i)-c5(j,i)*p_old(j+1,i)+R3(j,i))/c3(j,i));
         end
     end
    pg(:,1)=pg(:,2);
    pg(:,size(dx,2))=pg(:,size(dx,2)-1);
    pg(1,:)=pg(2,:);
    pg(size(dy,2),:)=pg(size(dy,2)-1,:);

    res=zeros([size(dy,2),size(dx,2)]);
     for j=2:size(dy,2)-1
         for i=2:size(dx,2)-1
             res(j,i)= c1(j,i)*pg(j,i+1)+c2(j,i)*pg(j,i-1)+c3(j,i)*pg(j,i)+c4(j,i)*pg(j-1,i)+c5(j,i)*pg(j+1,i)-R3(j,i);
         end
     end
 residue3= mean(mean(abs(res)));
 p_old=pg;
 itergr3=itergr3+1;

end




