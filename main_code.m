clc;
clear all;
tic
%% Initializing domain
Lx=1;
Ly=1;
N=129;

% Element 1 is 0, Elements 2 to (N+1) are equally spaced as 1/N, Last
% element is 0
dx=[0 (Lx/N)*ones(1,N) 0];
dy=[0 (Ly/N)*ones(1,N) 0];

delta_t=0.001;  % Time dicretisation
Re=1000;      % Reynold's Number  
Nt=10;
 
%% Initialising u, v and p 
u=zeros(size(dy,2),size(dx,2));
v=zeros(size(dy,2),size(dx,2));
p=zeros(size(dy,2),size(dx,2));

%% Boundary conditions 

% Top left hand corner of the cavity is the origin. Horizontal direction -
% i, Vertical direction is j 

u=zeros(size(dy,2),size(dx,2));
u(:,1)=0; % Left boundary
v(:,1)=0;

u(:,size(dx,2))=0;  % Right boundary
v(:,size(dx,2))=0;

u(1,:)=1;  % Top boundary
v(1,:)=0;

u(size(dy,2),:)=0;  % Bottom boundary
v(size(dy,2),:)=0;

% Initial guesses for u and v
u0=u;
v0=v;
%% Defining the coefficients from the discretisation

% Diffusion terms (Term-1) for calculating the intermdiate velocities u* and v*

a1=zeros(size(dy,2),size(dx,2));
a2=zeros(size(dy,2),size(dx,2));
a3=zeros(size(dy,2),size(dx,2));
a4=zeros(size(dy,2),size(dx,2));
a5=zeros(size(dy,2),size(dx,2));

for j=2:size(dy,2)-1
    for i=2:size(dx,2)-1 
        a1(j,i)= -(delta_t/(2*Re))*(2/((dx(i+1)+dx(i))*dx(i)));
        a2(j,i)= -(delta_t/(2*Re))*(2/((dx(i-1)+dx(i))*dx(i)));
        a4(j,i)= -(delta_t/(2*Re))*(2/((dy(j-1)+dy(j))*dy(j)));
        a5(j,i)= -(delta_t/(2*Re))*(2/((dy(j+1)+dy(j))*dy(j)));
        a3(j,i)= 1+(delta_t/(2*Re))*(((2/(dx(i+1)+dx(i))+(2/(dx(i-1)+dx(i))))*(1/dx(i)))+ ...
            ((2/(dy(j+1)+dy(j))+(2/(dy(j-1)+dy(j))))*(1/dy(j))));

    end
end

% Diffusion terms (Term-2) for calculating the values at time instant n 

b1=-a1;
b2=-a2;
b3= 2-a3;
b4=-a4;
b5=-a5;


% Convective terms (Term-3) for calculating the values at time instant n and (n-1)

for j=2:size(dy,2)-1
for i=2:size(dx,2)-1
    le(j,i)=dx(i+1)/(dx(i)+dx(i+1));
    lw(j,i)=dx(i-1)/(dx(i)+dx(i-1));
    ln(j,i)=dy(j+1)/(dy(j)+dy(j+1));
    ls(j,i)=dy(j-1)/(dy(j)+dy(j-1));
end
end


%% Initialise the RHS of the momentum equations

H1=zeros(size(dy,2),size(dx,2));  % RHS of the u- momentum equation 
H2=zeros(size(dy,2),size(dx,2));  % RHS of the v- momentum equation 

%% Momentum equation solving to get intermdiate u* and v*

k=1;   % Time counter starts 
for t=0:delta_t:Nt
    k 
for j=2:size(dy,2)-1
         for i=2:size(dx,2)-1
             
             co1(j,i) = (b1(j,i)*u(j,i+1)+b2(j,i)*u(j,i-1)+b5(j,i)*u(j+1,i)+b4(j,i)*u(j-1,i)+b3(j,i)*u(j,i));
             co2(j,i) = (1.5/dx(i))*((le(j,i)*u(j,i)+(1-le(j,i))*u(j,i+1))^2-(lw(j,i)*u(j,i)+(1-lw(j,i))*u(j,i-1))^2);
             co3(j,i) = (1.5/dy(j))*((ln(j,i)*u(j,i)+(1-ln(j,i))*u(j+1,i))*(ln(j,i)*v(j,i)+(1-ln(j,i))*v(j+1,i)) ...
                   - (ls(j,i)*u(j,i)+ (1-ls(j,i))*u(j-1,i))*(ls(j,i)*v(j,i)+(1-ls(j,i))*v(j-1,i)));
             co4(j,i) = -(0.5/dx(i))*((le(j,i)*u0(j,i)+(1-le(j,i))*u0(j,i+1))^2-(lw(j,i)*u0(j,i)+(1-lw(j,i))*u0(j,i-1))^2);
             co5(j,i) = -(0.5/dy(j))*((ln(j,i)*u0(j,i)+(1-ln(j,i))*u0(j+1,i))*(ln(j,i)*v0(j,i)+(1-ln(j,i))*v0(j+1,i)) ...
                   - ((ls(j,i)*u0(j,i)+(1-ls(j,i))*u0(j-1,i))*(ls(j,i)*v0(j,i)+(1-ls(j,i))*v0(j-1,i))));  

             H1(j,i)=  co1(j,i) - delta_t*(co2(j,i)+co3(j,i)+co4(j,i)+co5(j,i));
  

             H2(j,i)= (b1(j,i)*v(j,i+1)+b2(j,i)*v(j,i-1)+b3(j,i)*v(j,i)+b4(j,i)*v(j-1,i)+b5(j,i)*v(j+1,i)) ...
         - delta_t*((1.5/dy(j))*((ln(j,i)*v(j,i)+(1-ln(j,i))*v(j+1,i))^2-(ls(j,i)*v(j,i)+(1-ls(j,i))*v(j-1,i))^2)...
         -(0.5/dy(j))*((ln(j,i)*v0(j,i)+(1-ln(j,i))*v0(j+1,i))^2-(ls(j,i)*v0(j,i)+(1-ls(j,i))*v0(j-1,i))^2)...
         + (1.5/dx(i))*((le(j,i)*u(j,i)+(1-le(j,i))*u(j,i+1))*(le(j,i)*v(j,i)+(1-le(j,i))*v(j,i+1))-(lw(j,i)*u(j,i)+(1-lw(j,i))*u(j,i-1))*(lw(j,i)*v(j,i)+(1-lw(j,i))*v(j,i-1)))...
         - (0.5/dx(i))*((le(j,i)*u0(j,i)+(1-le(j,i))*u0(j,i+1))*(le(j,i)*v0(j,i)+(1-le(j,i))*v0(j,i+1))-(lw(j,i)*u0(j,i)+(1-lw(j,i))*u0(j,i-1))*(lw(j,i)*v0(j,i)+(1-lw(j,i))*v0(j,i-1))));      


         end
end

[ug,itergr1,res1]=gauss_seidel_u(u,dx,dy,a1,a2,a3,a4,a5,H1);
[vg,itergr2,res2]=gauss_seidel_v(v,dx,dy,a1,a2,a3,a4,a5,H2);

% Global Mass conservation
  % ucorrect=0; 
  % for j=1:size(dy,2)
  %   ucorrect=ucorrect+(((ug(j,size(dx,2)-1)-ug(j,1))*dy(j))/Ly);
  % end
  % ug(:,size(dx,2)-1)=ug(:,size(dx,2)-1)-ucorrect;
  % ug(:,size(dx,2))=ug(:,size(dx,2)-1);

% Updating the face velocities
for j=2:size(dy,2)-1
         for i=2:size(dx,2)-1

          uface(j,i)= le(j,i)*ug(j,i)+(1-le(j,i))*ug(j,i+1);
          vface(j,i)=ln(j,i)*vg(j,i)+(1-ln(j,i))*vg(j+1,i);   

         end
end

%% Pressure solver 

% Pressure Solver coefficients 

c1=-2*(Re/delta_t)*a1;
c2=-2*(Re/delta_t)*a2;
c4=-2*(Re/delta_t)*a4;
c5=-2*(Re/delta_t)*a5;
c3=-2*(Re/delta_t)*(a3-1);


for j=2:size(dy,2)-1
         for i=2:size(dx,2)-1  
            H3(j,i)= (((le(j,i)*ug(j,i)+(1-le(j,i))*ug(j,i+1)-lw(j,i)*ug(j,i)-(1-lw(j,i))*ug(j,i-1))/dx(i))+...
             (((ln(j,i)*vg(j,i)+(1-ln(j,i))*vg(j+1,i))-ls(j,i)*vg(j,i)-(1-ls(j,i))*vg(j-1,i))/dy(j)))/delta_t;
         end
end

[pg,itergr3,res3]= gauss_seidel_p(p,dx,dy,c1,c2,c3,c4,c5,H3);

%% Velocity updating using the u*, v* and pg

uf=zeros(size(dy,2),size(dx,2));
vf=zeros(size(dy,2),size(dx,2));

uf(:,1)=0; % Left boundary
vf(:,1)=0;
uf(:,size(dx,2))=0; % Right boundary
vf(:,size(dx,2))=0;

uf(1,:)=1;   % Top boundary
vf(1,:)=0;
uf(size(dy,2),:)=0; % Bottom boundary
vf(size(dy,2),:)=0;

for j=2:size(dy,2)-1
         for i=2:size(dx,2)-1

          uf(j,i)=ug(j,i)-(delta_t/dx(i))*(le(j,i)*pg(j,i)+(1-le(j,i))*pg(j,i+1)-lw(j,i)*pg(j,i)-(1-lw(j,i))*pg(j,i-1));
          vf(j,i)=vg(j,i)-(delta_t/dy(j))*(ln(j,i)*pg(j,i)+(1-ln(j,i))*pg(j+1,i)-ls(j,i)*pg(j,i)-(1-ls(j,i))*pg(j-1,i));

         end
end

mean(mean(abs(p-pg)))

u0=u;
u=uf;
v0=v;
v=vf;
p=pg;

k=k+1;

end

toc

%% Plotting

% Lx=1;
% Ly=1;
% N=size(dx,2)-2;
% x=zeros(1,N+2);
% y=zeros(1,N+2);
% for i=2:size(dx,2)-1
%     x(i)=x(i-1)+(dx(i-1)+dx(i))/2;
%     y(i)=y(i-1)+(dy(i-1)+dy(i))/2;
% end
% 
% figure 
%  [X,Y] = meshgrid(x,y);
%  colormap('jet');
% contourf(X,Y,uf)
%  axis ij
%  xlabel('X','fontsize',14);
%  ylabel('Y','fontsize',14);
%  pbaspect ([1 1 1])
%  colorbar

%%

%vorticity
dx1=0.0078;
dy1=0.0078;
 xx=0:dx1:Lx;
 yy=0:dy1:Ly;
 [X1,Y1] = meshgrid(xx,yy);
for j=2:size(y,2)-1
    for i=2:size(x,2)-1
    omega(j,i)=((vf(j,i+1)-vf(j,i-1))/(2*dx1))-((uf(j+1,i)-uf(j-1,i))/(2*dy1));
    end
end
contour(X1,Y1,flip(omega),[-4,-3,-2,-1,0,1,2,3,4,5],'ShowText','on');
 pbaspect ([1 1 1])


%%  Exporting for plotting in Tecplot
% I=16;
% J=16;
% Xp = (reshape(X,[1 , size(X,1)*size(Y,2)]))';
% Yp = (reshape(Y,[1 , size(X,1)*size(Y,2)]))';
% Uvel = (reshape(uf,[1 , size(X,1)*size(Y,2)]))';
% Vvel = (reshape(vf,[1 , size(X,1)*size(Y,2)]))';
% 
% Uvel(isnan(Uvel))=0;
% Vvel(isnan(Vvel))=0;
% 
% OutputFolder='C:\Users\Chintan\Desktop\OneDrive - Johns Hopkins\Springs_2024\CFD\Assignment-4\';
% plotlineu_1='TITLE = "INITIAL"';
% plotlineu_2='VARIABLES = "x","y","u","v"';
% plotlineu_3=['ZONE T="zeros" I=',num2str(I,'%d'),', J=', num2str(J,'%d'),', F=POINT'];
% vec=[Xp,Yp,Uvel,Vvel];
% fidv = fopen([OutputFolder,'vel.dat'], 'wt+');
%     fprintf(fidv, '%s\n', plotlineu_1);
%     fprintf(fidv, '%s\n', plotlineu_2);
%     fprintf(fidv, '%s\n', plotlineu_3);
%     fprintf(fidv, '%f %f %f %f\n',vec');
%     fclose(fidv);
save('Re1000_N128.mat', 'uf', 'vf', 'dx', 'dy');