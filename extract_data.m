%% plotting
Lx=1;
Ly=1;
N=size(dx,2)-2;
x=zeros(1,N+2);
y=zeros(1,N+2);
for i=2:size(dx,2)
    x(i)=x(i-1)+(dx(i-1)+dx(i))/2;
    y(i)=y(i-1)+(dy(i-1)+dy(i))/2;
end
    %%
figure
%  
 [X,Y] = meshgrid(x,y);
 colormap('jet');
 axis ij
 xlabel('X','fontsize',14);
 ylabel('Y','fontsize',14);
 pbaspect ([1 1 1])
 colorbar
 % contour(X,Y,flip(qg),100)
%  contour(X,Y,flip(qg),'LevelListMode','manual','LevelList',[0.12:-0.002:0.001,10^-3,10^-4,10^-7,10^-6,-9*10^-6,-10^-1,-0.1*10^-2:-0.1*10^-2:-0.9*10^-2,-0.1*10^-3:-0.1*10^-3:-0.9*10^-3,-0.1*10^-4:-0.1*10^-4:-0.9*10^-4,-10^-5,-0.1*10^-6:-0.1*10^-6:-0.9*10^-6,-0.1*10^-7:-0.1*10^-7:-0.9*10^-7,-10^-8]);

%%
dyy=[0 (1/127)*ones(1,127) 0];
yy=zeros(1,129);
for i=2:129
     yy(i)=yy(i-1)+(dyy(i-1)+dyy(i))/2;
end
ughiaRe100=[0 -0.03717 -0.04192 -0.04775 -0.06434 -.10150 -0.15662 -0.2 -0.020581 0.15641  0.68717 0.73722 1];
ughiaRe1000=[0 -0.18109 -0.20196 -0.22220 -0.29730 -0.38289 -0.27805 -0.10648 -0.06080 0.05702 0.18719 0.33304 0.46604 0.51117 0.57492 0.65928 1];
vghiaRe100=[0 0.09233 0.10091 0.10890 0.12317 0.16077 0.17507  0.05454  -0.22445 -0.239 -0.10313 -0.08864 -0.07391 -0.05906 0 ];
vghiaRe1000=[0 0.27485 0.29012 0.30353 0.32627 0.37095 0.33075 0.32235 0.02526 -0.31966 -0.42665 -0.51550 -0.39188 -0.33714 -0.27669 -0.21388 0];
% yghiaRe100=[yy(1) yy(8) yy(9) yy(10) yy(14) yy(23) yy(37) yy(59) yy(65) yy(80) yy(95) yy(110) yy(123) yy(124) yy(125) yy(126) yy(129)];
yghiaRe100=[yy(1) yy(8) yy(9) yy(10) yy(14) yy(23) yy(37)  yy(65)  yy(95) yy(110) yy(123) yy(124) yy(125) yy(126) yy(129)];
yghiaRe1000=[yy(1) yy(9) yy(10) yy(11) yy(13) yy(21) yy(30) yy(31) yy(65) yy(104) yy(111) yy(117) yy(122) yy(123) yy(124) yy(125) yy(129)];
%%
a=load('Re100_N16.mat','uf','vf','dx','dy');
N=size(a.dx,2)-2;
x1=zeros(1,N+2);
for i=2:size(a.dx,2)
    x1(i)=x1(i-1)+(a.dx(i-1)+a.dx(i))/2;
end
b=load('Re100_N32.mat','uf','vf','dx','dy');
N=size(b.dx,2)-2;
x2=zeros(1,N+2);
for i=2:size(b.dx,2)
    x2(i)=x2(i-1)+(b.dx(i-1)+b.dx(i))/2;
end
c=load('Re100_N64.mat','uf','vf','dx','dy');
N=size(c.dx,2)-2;
x3=zeros(1,N+2);
for i=2:size(c.dx,2)
    x3(i)=x3(i-1)+(c.dx(i-1)+c.dx(i))/2;
end
d=load('Re100_N128.mat','uf','vf','dx','dy');
N=size(d.dx,2)-2;
x4=zeros(1,N+2);
for i=2:size(d.dx,2)
    x4(i)=x4(i-1)+(d.dx(i-1)+d.dx(i))/2;
end
e=load('Re1000_N256.mat','uf','vf','dx','dy');
N=size(e.dx,2)-2;
x5=zeros(1,N+2);
for i=2:size(e.dx,2)
    x5(i)=x5(i-1)+(e.dx(i-1)+e.dx(i))/2;
end
%%
 figure
 plot(0.5*(a.uf(:,size(a.uf,2)*0.5)+a.uf(:,size(a.uf,2)*0.5-1)),fliplr(x1),'k--', 'LineWidth', 2)
 hold on
 plot(0.5*(b.uf(:,size(b.uf,2)*0.5)+b.uf(:,size(b.uf,2)*0.5-1)),fliplr(x2),'g--', 'LineWidth', 2)
 plot(0.5*(c.uf(:,size(c.uf,2)*0.5)+c.uf(:,size(c.uf,2)*0.5-1)),fliplr(x3),'b--.', 'LineWidth', 2)
 plot(0.5*(d.uf(:,size(d.uf,2)*0.5)+d.uf(:,size(d.uf,2)*0.5-1)),fliplr(x4),'m--.', 'LineWidth', 2)
% plot(0.5*(e.uf(:,size(e.uf,2)*0.5)+e.uf(:,size(e.uf,2)*0.5-1)),fliplr(x5),'c-.')
 plot( ughiaRe100,yghiaRe100(1:13), 'k*', 'LineWidth', 2);

 xlabel('U','fontsize',14);
 ylabel('Y','fontsize',14);
 pbaspect ([1 1 1])
 legend('N=16','N=32','N=64','N=128','Ghia et al,1982','location','southeast');
 
 %%
 
 figure
plot(x1, -0.5*(a.vf(size(a.vf,2)*0.5,:) + a.vf(size(a.vf,2)*0.5-1,:)), 'r--', 'LineWidth', 2)
hold on
plot(x2, -0.5*(b.vf(size(b.vf,2)*0.5,:) + b.vf(size(b.vf,2)*0.5-1,:)), 'g-', 'LineWidth', 2)
plot(x3, -0.5*(c.vf(size(c.vf,2)*0.5,:) + c.vf(size(c.vf,2)*0.5-1,:)), 'k-.', 'LineWidth', 2)
plot(x4, -0.5*(d.vf(size(d.vf,2)*0.5,:) + d.vf(size(d.vf,2)*0.5-1,:)), 'm-.', 'LineWidth', 2)
% plot(x5, -0.5*(e.vf(size(e.vf,2)*0.5,:) + e.vf(size(e.vf,2)*0.5-1,:)), 'c-.', 'LineWidth', 2)
plot(yghiaRe100, vghiaRe100, 'k*', 'LineWidth', 1.5);  % you can increase this too if needed

xlabel('X','fontsize',14);
ylabel('V','fontsize',14);
pbaspect([1 1 1])
legend('N=16','N=32','N=64','N=128','Ghia et al, 1982','location','southeast');

