% plot_hermes2.m
% Extremely simple script to plot x, y, and altitude from Hermes2.csv
% Expected columns in the CSV: x_meters, y_meters, altitude_meters

T = readtable('Hermes2_lat.csv');

% 3D trajectory (x, y, altitud_late)
figure;
plot3(T.latitude, T.longitude, T.gpsAltitude, '-');
grid on;
xlabel('lat');
ylabel('long');
zlabel('altitude [m]');
title('Trajectory (lat, long, altitude)');


