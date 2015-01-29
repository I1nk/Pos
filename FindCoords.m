close all
clear
clc

an = 0;
an1 = 0;

x = -cosd(35);
z = sind(35);
y = 0;

x1 = cosd(36);
z1 = sind(36);
y1 = 0;

xx = [x;y;z];
yy = [x1;y1;z1];

rots = [cosd(an), 0, sind(an);
         0, 1, 0;
         -sind(an), 0, cosd(an)];
     
xx = rots * xx
yy = rots * yy

rots = [1, 0, 0;
         0, cosd(an1), -sind(an1);
         sind(an1), 0, cosd(an1)];

xx = rots * xx
yy = rots * yy

v1 1.24 1
