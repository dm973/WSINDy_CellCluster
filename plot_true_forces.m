clear all; clc; close all;

exper = '111';

data_dr = ['~/Desktop/data_dr/',exper,'/'];
save_dr = data_dr;
input_data = findfilestrloc(save_dr,'sim',1);
load([save_dr,input_data],'f_xv_true_cell','h_xv_true_cell','d_xv_true_cell','inds_cell_true');

n = 200;
r = 0.2;
v =[0 1];
[rr,th,xx,yy] = build_polar_grid(n,r,v);
x = xx(1,:);
y = yy(:,1);

figure(1);
c_range = [-15 15]; force_ind=1;
F_dat = f_xv_true_cell{1}(rr,th,0); surfplot;
legend({'$f^\star_{a-r}$'},'interpreter','latex','fontsize',12)
saveas(gcf,'~/Desktop/ftrue.png')

figure(2); 
force_ind=2;
F_dat = h_xv_true_cell{1}(rr,th,0); surfplot;
legend({'$f^\star_{align}$'},'interpreter','latex','fontsize',12)
saveas(gcf,'~/Desktop/htrue.png')

figure(3); 
c_range = [-inf inf]; force_ind=3;
r = 0.2;
[rr,th] = meshgrid(linspace(0,r,n),linspace(0,pi,n));
xx=rr;
yy=th;
x = xx(1,:);
y = yy(:,1);
F_dat = d_xv_true_cell{1}(rr,th,0); surfplot;
legend({'$f^\star_{drag}$'},'interpreter','latex','fontsize',12)
saveas(gcf,'~/Desktop/dtrue.png')

