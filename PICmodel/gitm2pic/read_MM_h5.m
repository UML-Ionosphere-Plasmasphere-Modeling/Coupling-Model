clear all;clc

close all
clearvars -except inputyear mode    
% h5_files=dir(['c:\Users\Yifan\Documents\GitHub\Coupling-Model\PICmodel\*.h5']);
h5_files=dir(['c:\Users\Yifan\Documents\GITM-M-Modeling\PICmodel\*.h5']);
h5_files=struct2cell(h5_files);
h5_files=h5_files(1,:)';


for roll=1:numel(h5_files)
close all
clear global
clearvars -except roll h5_files 
tic;
filesinfo=h5info(h5_files{roll});
% h1=h5disp(h5_files{roll});
% gpsinfo=hdf5info('D:\DATA_calculation\TEC_mat\2015\gps150314g.001.hdf5');
data_const = h5read(h5_files{roll},'/ArrayOfGrids_const');
data=h5read(h5_files{roll},'/ArrayOfGrids_1');
toc;

% control panel
gridsize = 17;
showsize = 1;
%%%%%%%%%%%%%%

for face=[3,6]
posx=data_const.pos3.x(:,:,:,face)/1e3/(6371);
posy=data_const.pos3.y(:,:,:,face)/1e3/(6371);
posz=data_const.pos3.z(:,:,:,face)/1e3/(6371);

posx_const=data_const.pos3.x(showsize,:,:,face)/1e3/(6371);
posy_const=data_const.pos3.y(showsize,:,:,face)/1e3/(6371);
posz_const=data_const.pos3.z(showsize,:,:,face)/1e3/(6371);
Vx_const=data_const.v3.x(showsize,:,:,face);
Vy_const=data_const.v3.y(showsize,:,:,face);
Vz_const=data_const.v3.z(showsize,:,:,face);
Ex_const=data_const.e3.x(showsize,:,:,face);
Ey_const=data_const.e3.y(showsize,:,:,face);
Ez_const=data_const.e3.z(showsize,:,:,face);

Ex=data.e3.x(:,:,:,face);
Ey=data.e3.y(:,:,:,face);
Ez=data.e3.z(:,:,:,face);
Bx=data_const.b3.x(:,:,:,face);
By=data_const.b3.y(:,:,:,face);
Bz=data_const.b3.z(:,:,:,face);
Vx=data.v3.x(:,:,:,face);
Vy=data.v3.y(:,:,:,face);
Vz=data.v3.z(:,:,:,face);
N=data.density(:,:,:,face);
figure(1)
plot3(reshape(posx,gridsize*gridsize*gridsize,1),reshape(posy,gridsize*gridsize*gridsize,1),reshape(posz,gridsize*gridsize*gridsize,1),'*');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
figure(2)
quiver3(reshape(posx,gridsize*gridsize*gridsize,1),reshape(posy,gridsize*gridsize*gridsize,1),reshape(posz,gridsize*gridsize*gridsize,1),reshape(Bx,gridsize*gridsize*gridsize,1),reshape(By,gridsize*gridsize*gridsize,1),reshape(Bz,gridsize*gridsize*gridsize,1),3,'b');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
figure(3)
quiver3(reshape(posx_const,1,gridsize*gridsize,1),reshape(posy_const,1,gridsize*gridsize,1),reshape(posz_const,1,gridsize*gridsize,1),reshape(Ex_const,1,gridsize*gridsize,1),reshape(Ey_const,1,gridsize*gridsize,1),reshape(Ez_const,1,gridsize*gridsize,1),3,'b');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
figure(4)
quiver3(reshape(posx_const,1,gridsize*gridsize,1),reshape(posy_const,1,gridsize*gridsize,1),reshape(posz_const,1,gridsize*gridsize,1),reshape(Vx_const,1,gridsize*gridsize,1),reshape(Vy_const,1,gridsize*gridsize,1),reshape(Vz_const,1,gridsize*gridsize,1),3,'b');hold on
grid on
box on
xlabel('X')
ylabel('Y')
zlabel('Z')
end
end
