clear all;clc

close all
clearvars -except inputyear mode    
% h5_files=dir(['c:\Users\Yifan\Documents\GitHub\Coupling-Model\PICmodel\*.h5']);
h5_files=dir(['c:\Users\Yifan\Documents\Github\Coupling-Model\model-parallel-v.2\GridsData.h5']);
h5_files=struct2cell(h5_files);
h5_files=h5_files(1,:)';

% important: the face index is from 1-6, not like 0-5 in C code

for roll=1:numel(h5_files)
close all
clear global
clearvars -except roll h5_files 
tic;
filesinfo=h5info(h5_files{roll});
% h1=h5disp(h5_files{roll});
% gpsinfo=hdf5info('D:\DATA_calculation\TEC_mat\2015\gps150314g.001.hdf5');
data_const = h5read(h5_files{roll},'/ArrayOfGrids_const');
data=h5read(h5_files{roll},'/ArrayOfGrids_200');
toc;

% control panel
gridsize = 17; % fsize+1
showsize = 3; % shell
templevel = 1; % level overlap with other possible region
%%%%%%%%%%%%%%

k_level = [(templevel+2):(gridsize-1-templevel)];


face=[1,2,3,4,5,6];
posx=data_const.pos3.x(:,:,:,face)/1e3/(6371);
posy=data_const.pos3.y(:,:,:,face)/1e3/(6371);
posz=data_const.pos3.z(:,:,:,face)/1e3/(6371);
Bx_const=data_const.b3.x(:,:,:,face);
By_const=data_const.b3.y(:,:,:,face);
Bz_const=data_const.b3.z(:,:,:,face);

posx_const=data_const.pos3.x(showsize,:,:,face)/1e3/(6371);
posy_const=data_const.pos3.y(showsize,:,:,face)/1e3/(6371);
posz_const=data_const.pos3.z(showsize,:,:,face)/1e3/(6371);

Vx_const=data_const.v3.x(showsize,:,:,face);
Vy_const=data_const.v3.y(showsize,:,:,face);
Vz_const=data_const.v3.z(showsize,:,:,face);
Ex_const=data_const.e3.x(showsize,:,:,face);
Ey_const=data_const.e3.y(showsize,:,:,face);
Ez_const=data_const.e3.z(showsize,:,:,face);


denH_const=data_const.densityH(showsize,:,:,face);
denH_showsize = data.densityH(showsize,:,:,face);

Vx_showsize=data.v3.x(showsize,:,:,face);
Vy_showsize=data.v3.y(showsize,:,:,face);
Vz_showsize=data.v3.z(showsize,:,:,face);
Ex_showsize=data.e3.x(showsize,:,:,face);
Ey_showsize=data.e3.y(showsize,:,:,face);
Ez_showsize=data.e3.z(showsize,:,:,face);



figure(1)
plot3(reshape(posx,[],1),reshape(posy,[],1),reshape(posz,[],1),'o');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('pos')

figure(2)
quiver3(reshape(posx,[],1),reshape(posy,[],1),reshape(posz,[],1),reshape(Bx_const,[],1),reshape(By_const,[],1),reshape(Bz_const,[],1),3,'b');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('B-const')

figure(3)
quiver3(reshape(posx_const,[],1),reshape(posy_const,[],1),reshape(posz_const,[],1),reshape(Ex_const,[],1),reshape(Ey_const,[],1),reshape(Ez_const,[],1),3,'b');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
view(-55,45)
title('E-const')

figure(4)
quiver3(reshape(posx_const,[],1),reshape(posy_const,[],1),reshape(posz_const,[],1),reshape(Vx_const,[],1),reshape(Vy_const,[],1),reshape(Vz_const,[],1),3,'b');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
view(-55,45)
title('v-const')

figure(5)
scatter3(reshape(posx_const,[],1),reshape(posy_const,[],1),reshape(posz_const,[],1),10,reshape(denH_const,[],1),'filled');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
view(-55,45)
title('denH-const')
colorbar()
set(gca,'Xlim',[-5,5],'Ylim',[-5,5],'Zlim',[-5,5])


figure(6)
scatter3(reshape(posx_const,[],1),reshape(posy_const,[],1),reshape(posz_const,[],1),10,reshape(denH_showsize,[],1),'filled');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
view(-55,45)
title('denH-showsize')
colorbar
set(gca,'Xlim',[-5,5],'Ylim',[-5,5],'Zlim',[-5,5])


figure(7)
quiver3(reshape(posx_const,[],1),reshape(posy_const,[],1),reshape(posz_const,[],1),reshape(Vx_showsize,[],1),reshape(Vy_showsize,[],1),reshape(Vz_showsize,[],1),3,'b');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('v-showsize')

figure(8)
quiver3(reshape(posx_const,[],1),reshape(posy_const,[],1),reshape(posz_const,[],1),reshape(Ex_showsize,[],1),reshape(Ey_showsize,[],1),reshape(Ez_showsize,[],1),3,'b');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('E-showsize')



%equatorial plane velocity
face=[1,2,4,5];
posx_equa=data_const.pos3.x(k_level,(gridsize+1)/2,:,face)/1e3/(6371);
posy_equa=data_const.pos3.y(k_level,(gridsize+1)/2,:,face)/1e3/(6371);
posz_equa=data_const.pos3.z(k_level,(gridsize+1)/2,:,face)/1e3/(6371);

Vx=data.v3.x(k_level,(gridsize+1)/2,:,face);
Vy=data.v3.y(k_level,(gridsize+1)/2,:,face);
Vz=data.v3.z(k_level,(gridsize+1)/2,:,face);

figure(10)
quiver3(reshape(posx_equa,[],1),reshape(posy_equa,[],1),reshape(posz_equa,[],1),reshape(Vx,[],1),reshape(Vy,[],1),reshape(Vz,[],1),1,'b');hold on

grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')    
title('v-ion')
set(gca,'Xlim',[-5,5],'Ylim',[-5,5],'Zlim',[-5,5])

%equatorial plane density 
N_equa = data.densityO(k_level,(gridsize+1)/2,:,face);

figure(11)
scatter3(reshape(posx_equa,[],1),reshape(posy_equa,[],1),reshape(posz_equa,[],1),10,reshape(N_equa,[],1),'filled');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('den_O')
colorbar
set(gca,'Xlim',[-5,5],'Ylim',[-5,5],'Zlim',[-5,5])

%equatorial plane density 
N_equa = data.densityH(k_level,(gridsize+1)/2,:,face);

figure(12)
scatter3(reshape(posx_equa,[],1),reshape(posy_equa,[],1),reshape(posz_equa,[],1),10,reshape(N_equa,[],1),'filled');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('den_H')
colorbar
set(gca,'Xlim',[-5,5],'Ylim',[-5,5],'Zlim',[-5,5])
        

% vertical plane velocity noon/midnight
face=[1,3,4,6];
posx_vert=data_const.pos3.x(k_level,:,(gridsize+1)/2,face)/1e3/(6371);
posy_vert=data_const.pos3.y(k_level,:,(gridsize+1)/2,face)/1e3/(6371);
posz_vert=data_const.pos3.z(k_level,:,(gridsize+1)/2,face)/1e3/(6371);
Vx=data.vO3.x(k_level,:,(gridsize+1)/2,face);
Vy=data.vO3.y(k_level,:,(gridsize+1)/2,face);
Vz=data.vO3.z(k_level,:,(gridsize+1)/2,face);

figure(20)
quiver3(reshape(posx_vert,[],1),reshape(posy_vert,[],1),reshape(posz_vert,[],1),reshape(Vx,[],1),reshape(Vy,[],1),reshape(Vz,[],1),1.5,'b');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('vO')
set(gca,'Xlim',[-5,5],'Ylim',[-0.5,0.5],'Zlim',[-5,5])
    

% vertical plane density noon/midnight
N_vert = data.densityO(k_level,:,(gridsize+1)/2,face);

figure(21)
scatter3(reshape(posx_vert,[],1),reshape(posy_vert,[],1),reshape(posz_vert,[],1),10,reshape(N_vert,[],1),'filled');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('den_O')
colorbar
set(gca,'Xlim',[-5,5],'Ylim',[-5,5],'Zlim',[-5,5])
   

% vertical plane velocity noon/midnight
Vx=data.vH3.x(k_level,:,(gridsize+1)/2,face);
Vy=data.vH3.y(k_level,:,(gridsize+1)/2,face);
Vz=data.vH3.z(k_level,:,(gridsize+1)/2,face);

figure(22)
quiver3(reshape(posx_vert,[],1),reshape(posy_vert,[],1),reshape(posz_vert,[],1),reshape(Vx,[],1),reshape(Vy,[],1),reshape(Vz,[],1),1.5,'b');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('vH')
set(gca,'Xlim',[-5,5],'Ylim',[-0.5,0.5],'Zlim',[-5,5])
    

% vertical plane density noon/midnight
N_vert = data.densityH(k_level,:,(gridsize+1)/2,face);

figure(23)
scatter3(reshape(posx_vert,[],1),reshape(posy_vert,[],1),reshape(posz_vert,[],1),10,reshape(N_vert,[],1),'filled');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('den_H')
colorbar
set(gca,'Xlim',[-5,5],'Ylim',[-5,5],'Zlim',[-5,5])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%       dawn/dust  O
face=[3,6];
posx_vert=data_const.pos3.x(k_level,(gridsize+1)/2,:,face)/1e3/(6371);
posy_vert=data_const.pos3.y(k_level,(gridsize+1)/2,:,face)/1e3/(6371);
posz_vert=data_const.pos3.z(k_level,(gridsize+1)/2,:,face)/1e3/(6371);
Vx=data.vO3.x(k_level,(gridsize+1)/2,:,face);
Vy=data.vO3.y(k_level,(gridsize+1)/2,:,face);
Vz=data.vO3.z(k_level,(gridsize+1)/2,:,face);
N_vert = data.densityO(k_level,(gridsize+1)/2,:,face);

face1=[2,5];
posx_vert1=data_const.pos3.x(k_level,:,(gridsize+1)/2,face1)/1e3/(6371);
posy_vert1=data_const.pos3.y(k_level,:,(gridsize+1)/2,face1)/1e3/(6371);
posz_vert1=data_const.pos3.z(k_level,:,(gridsize+1)/2,face1)/1e3/(6371);
Vx1=data.vO3.x(k_level,:,(gridsize+1)/2,face1);
Vy1=data.vO3.y(k_level,:,(gridsize+1)/2,face1);
Vz1=data.vO3.z(k_level,:,(gridsize+1)/2,face1);
N_vert1 = data.densityO(k_level,:,(gridsize+1)/2,face1);

x= [reshape(posx_vert,[],1); reshape(posx_vert1,[],1)];
y= [reshape(posy_vert,[],1); reshape(posy_vert1,[],1)];
z= [reshape(posz_vert,[],1); reshape(posz_vert1,[],1)];
vx= [reshape(Vx,[],1); reshape(Vx1,[],1)];
vy= [reshape(Vy,[],1); reshape(Vy1,[],1)];
vz= [reshape(Vz,[],1); reshape(Vz1,[],1)];
n= [reshape(N_vert,[],1); reshape(N_vert1,[],1)];

figure(30)
quiver3(x,y,z,vx,vy,vz,2,'b');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('v_O')
set(gca,'Xlim',[-5,5],'Ylim',[-5,5],'Zlim',[-5,5])

% vertical plane density dawn/dust O

figure(31)
scatter3(x,y,z,10,n,'filled');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('den_O')
colorbar


% vertical plane velocity dawn/dust  H
face=[3,6];
Vx=data.vH3.x(k_level,(gridsize+1)/2,:,face);
Vy=data.vH3.y(k_level,(gridsize+1)/2,:,face);
Vz=data.vH3.z(k_level,(gridsize+1)/2,:,face);
N_vert = data.densityH(k_level,(gridsize+1)/2,:,face);

facd=[2,5];
Vx1=data.vH3.x(k_level,:,(gridsize+1)/2,face1);
Vy1=data.vH3.y(k_level,:,(gridsize+1)/2,face1);
Vz1=data.vH3.z(k_level,:,(gridsize+1)/2,face1);
N_vert1 = data.densityH(k_level,:,(gridsize+1)/2,face1);

vx= [reshape(Vx,[],1); reshape(Vx1,[],1)];
vy= [reshape(Vy,[],1); reshape(Vy1,[],1)];
vz= [reshape(Vz,[],1); reshape(Vz1,[],1)];
n= [reshape(N_vert,[],1); reshape(N_vert1,[],1)];


figure(32)
quiver3(x,y,z,vx,vy,vz,2,'b');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('v_H')
set(gca,'Xlim',[-5,5],'Ylim',[-5,5],'Zlim',[-5,5])

% vertical plane density dawn/dust H

figure(33)
scatter3(x,y,z,10,n,'filled');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('den_H')
colorbar

%%%

end
