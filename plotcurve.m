figure(1)
tmp=dlmread('band1chrwr1.dat');
plot(tmp(:,1),tmp(:,2)),hold on
tmp=[dlmread('band2chrwr1L.dat');dlmread('band2chrwr1c.dat');flipud(dlmread('band2chrwr1R.dat'))];
plot(tmp(:,1),tmp(:,2)),hold on
tmp=[dlmread('band3chrwr1L.dat');dlmread('band3chrwr1c.dat');flipud(dlmread('band3chrwr1R.dat'))];
plot(tmp(:,1),tmp(:,2)),hold on
tmp=dlmread('band4chrwr1L.dat');
plot(tmp(:,1),tmp(:,2)),hold on
tmp=dlmread('band4chrwr1R.dat');
plot(tmp(:,1),tmp(:,2)),hold on

tmp=dlmread('band1rwar1.dat');
plot(tmp(:,1),tmp(:,2)),hold on
tmp=[dlmread('band2rwar1L.dat');dlmread('band2rwar1c.dat');flipud(dlmread('band2rwar1R.dat'))];
plot(tmp(:,1),tmp(:,2)),hold on
tmp=[dlmread('band3rwar1L.dat');dlmread('band3rwar1c.dat');flipud(dlmread('band3rwar1R.dat'))];
plot(tmp(:,1),tmp(:,2)),hold on
tmp=dlmread('band4rwar1L.dat');
plot(tmp(:,1),tmp(:,2)),hold on
tmp=dlmread('band4rwar1R.dat');
plot(tmp(:,1),tmp(:,2)),hold on

figure(2)
tmp=dlmread('band1chrwr0p5.dat');
plot(tmp(:,1),tmp(:,2)),hold on
tmp=[dlmread('band2chrwr0p5L.dat');dlmread('band2chrwr0p5c.dat');flipud(dlmread('band2chrwr0p5R.dat'))];
plot(tmp(:,1),tmp(:,2)),hold on
tmp=dlmread('band3chrwr0p5L.dat');
plot(tmp(:,1),tmp(:,2)),hold on
tmp=flipud(dlmread('band3chrwr0p5R.dat'));
plot(tmp(:,1),tmp(:,2)),hold on
tmp=dlmread('band4chrwr0p5L.dat');
plot(tmp(:,1),tmp(:,2)),hold on
tmp=dlmread('band4chrwr0p5R.dat');
plot(tmp(:,1),tmp(:,2)),hold on


tmp=dlmread('band1rwar0p5c.dat');
plot(tmp(:,1),tmp(:,2)),hold on
tmp=[dlmread('band2rwar0p5L.dat');dlmread('band2rwar0p5c.dat');flipud(dlmread('band2rwar0p5R.dat'))];
plot(tmp(:,1),tmp(:,2)),hold on
tmp=[dlmread('band3rwar0p5L.dat');dlmread('band3rwar0p5c.dat');flipud(dlmread('band3rwar0p5R.dat'))];
plot(tmp(:,1),tmp(:,2)),hold on
tmp=dlmread('band4rwar0p5L.dat');
plot(tmp(:,1),tmp(:,2)),hold on
tmp=dlmread('band4rwar0p5R.dat');
plot(tmp(:,1),tmp(:,2)),hold on