% %  need work with Irf and geodetic tools
% % which could be download:
% https://github.com/irfu/irfu-matlab
% https://ww2.mathworks.cn/matlabcentral/fileexchange/15285-geodetic-toolbox
clear
% load 3DALL_t110331_120402
load_data=textread('trans_t110331_120402.txt');
alt_profile=unique(load_data(:,3));
set_alt=alt_profile(49);
trans_data=load_data(find(load_data(:,3)==set_alt),:);
h5_files=dir(['D:\DATA_calculation\read_MM_H5\*.h5']);
h5_files=struct2cell(h5_files);
h5_files=h5_files(1,:)';
filesinfo=h5info('GridsData.h5');

data_const = h5read('GridsData.h5','/ArrayOfGrids_const');
% control panel
gridsize = 17;
%%%%%%%%%%%%%%
pos3x=data_const.pos3.x/1e3/(6371);
pos3y=data_const.pos3.y/1e3/(6371);
pos3z=data_const.pos3.z/1e3/(6371);
% % % % % now xx = toward sun thus "pos3y" used. remember change it when rotate x and y dirction 
xx=reshape(pos3y,[17^3*6 1]);
yy=reshape(pos3x,[17^3*6 1]);
zz=reshape(pos3z,[17^3*6 1]);
px=reshape(yy,[17,17,17,6]);
 epoch = EpochUnix(repmat(iso2epoch('2011-03-31T12:04:02Z'),[17^3*6 1]));
transfer_martix(:,1)=xx;
transfer_martix(:,2)=yy;
transfer_martix(:,3)=zz;
    TsVec3DRTP = TSeries(epoch,transfer_martix,'vec_xyz');
outputgeo=irf.geocentric_coordinate_transformation(TsVec3DRTP,'sm>geo');
[x1,y1,z1]=ell2xyz(trans_data(:,2),trans_data(:,1),trans_data(:,3));
x1=x1/1e3/6371;
y1=y1/1e3/6371;
z1=z1/1e3/6371;
for i=1:numel(outputgeo.data(:,1))
distance=sqrt((x1-outputgeo.data(i,1)).^2+(y1-outputgeo.data(i,2)).^2+(z1-outputgeo.data(i,3)).^2);
[~,loc]=min(distance);
i2p_no(i,:)=trans_data(loc,4);%% O+ density
i2p_nh(i,:)=trans_data(loc,11);%% H+ density
i2p_te(i,:)=trans_data(loc,14);%% electron temperature
i2p_ti(i,:)=trans_data(loc,15);%% ion tempreature
i2p_ve(i,:)=trans_data(loc,16);%% electron temperature
i2p_vn(i,:)=trans_data(loc,17);%% ion tempreature
i2p_vup(i,:)=trans_data(loc,18);%% electron temperature

end
i2p_no=rotate(i2p_no,[17,17,17,6]);
i2p_nh=rotate(i2p_nh,[17,17,17,6]);
i2p_te=rotate(i2p_te,[17,17,17,6]);
i2p_ti=rotate(i2p_ti,[17,17,17,6]);
i2p_ve=rotate(i2p_ve,[17,17,17,6]);
i2p_vn=rotate(i2p_vn,[17,17,17,6]);
i2p_vup=rotate(i2p_vup,[17,17,17,6]);
% % i2p are the parameters going to be transfered to PIC model in
% % 17*17*17*6, you can decide the output type of these variables   
% % 1 longitude in rad
% % 2 latitude  in rad
% % 3 altitude  in m
% % 4 O(4sp)+ density
% % 5 O2+ density
% % 6 N2+ density
% % 7 N+ density
% % 8 NO+ density
% % 9 O(2D)+ density
% % 10 O(2P)+ density
% % 11 H+ density
% % 12 He+ density
% % 13 electron density
% % 14 electron temperature
% % 15 ion temperature
% % 16 eastward ion velocity
% % 17 northward ion velocity
% % 18 upward velocity

