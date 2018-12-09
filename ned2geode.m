clear all;
close all;
format long;


r_kalmanNED = xlsread('vuelta1_resultados.xlsx', 'estimacion_Kalman', 'B3:D72106');
geod_gnss = xlsread('vuelta1_data.xlsx', 'data_GNSS', 'F11:H11'); %posicion inicial GNSS

N1 = length(r_kalmanNED);

% Coordenadas estimacion Kalman en geodesicas
for n=1:N1
    n
    [r_kalmanGEO(n,1),r_kalmanGEO(n,2),r_kalmanGEO(n,3)] = ned2geodetic(r_kalmanNED(n,1), r_kalmanNED(n,2), r_kalmanNED(n,3), geod_gnss(1),geod_gnss(2), geod_gnss(3), wgs84Ellipsoid);
end



% %****************Guardar datos*****************
xlswrite('vuelta1_resultados.xlsx', r_kalmanGEO, 'estimacion_kalman', 'L3:N100000');