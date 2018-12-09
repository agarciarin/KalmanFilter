clear all;
close all;
format long;

diftime = xlsread('vuelta1_data.xlsx', 'data_IMU', 'C6:C72108'); %OJO adelantar una fila la difftime
quat_imu = xlsread('vuelta1_data.xlsx', 'data_IMU', 'D5:G72107');
accel_imu = xlsread('vuelta1_data.xlsx', 'data_IMU', 'H5:J72107');
gyro_imu = xlsread('vuelta1_data.xlsx', 'data_IMU', 'K5:M72107');
euler_imu = xlsread('vuelta1_data.xlsx', 'data_IMU', 'N5:P72107');    
head = xlsread('vuelta1_data.xlsx', 'data_IMU', 'Q5:Q72107');
geod_gnss = xlsread('vuelta1_data.xlsx', 'data_GNSS', 'F11:H371');
speed_gnss = xlsread('vuelta1_data.xlsx', 'data_GNSS', 'N11:N371');
% Qk = xlsread('vuelta1_data.xlsx', 'Qk', 'A1:J10');
% Rk = xlsread('vuelta1_data.xlsx', 'Rk', 'A1:C3');

save vuelta1_data.mat diftime quat_imu accel_imu gyro_imu euler_imu head geod_gnss speed_gnss



