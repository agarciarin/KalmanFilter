clear all;
close all;
format long;

%Cargar datos de IMU y GNSS
load('vuelta1_data');

%Matrices de covarianza y valor de la gravedad
% Qk = xlsread('vuelta1_data.xlsx', 'Qk', 'A1:J10');
% Rk = xlsread('vuelta1_data.xlsx', 'Rk', 'A1:C3');
cov_accel = 0.02;
cov_gyro = 0.002;
cov_quat = 0.00000002;
cov_gnssX = 0.01;
cov_gnssY = 0.01;
cov_gnssZ = 0.01;
Qk = [cov_accel 0 0 0 0 0 0 0 0 0;
      0 cov_accel 0 0 0 0 0 0 0 0;
      0 0 cov_accel 0 0 0 0 0 0 0;
      0 0 0 cov_gyro 0 0 0 0 0 0;
      0 0 0 0 cov_gyro 0 0 0 0 0;
      0 0 0 0 0 cov_gyro 0 0 0 0;
      0 0 0 0 0 0 cov_quat 0 0 0;
      0 0 0 0 0 0 0 cov_quat 0 0;
      0 0 0 0 0 0 0 0 cov_quat 0;
      0 0 0 0 0 0 0 0 0 cov_quat];
Rk = [cov_gnssX 0 0; 0 cov_gnssY 0; 0 0 cov_gnssZ];
g = [0; 0; 9.81];


% pasar a rotacion 'ZYX' (Yaw-Pitch-Roll)
% Ejes X-Y-Z Body
quat(:,1) = quat_imu(:,1);
quat(:,2) = quat_imu(:,3);
quat(:,3) = quat_imu(:,2);
quat(:,4) = -quat_imu(:,4);

accel(:,1) = accel_imu(:,2); % eje X_b
accel(:,2) = accel_imu(:,1); % eje Y_b
accel(:,3) = accel_imu(:,3); % eje Z_b

gyro(:,1) = gyro_imu(:,2);  % roll (sobre X_b)
gyro(:,2) = gyro_imu(:,1);  % pitch (sobre Y_b)
gyro(:,3) = -gyro_imu(:,3); % yaw (sobre Z_b)

euler(:,1) = euler_imu(:,2);    % roll (sobre X_b)
euler(:,2) = euler_imu(:,1);    % pitch (sobre Y_b)
euler(:,3) = -euler_imu(:,3);   % yaw (sobre Z_b)


% numero de ciclos
N1 = length(diftime);
N2 = length(geod_gnss);
dtGNSS = 1;

% Escala de tiempo datos IMU
s=0;
for n=1:N1
    T1(1,n) = s;
    s = s + diftime(n);
end

% Escala de tiempo posicion y velocidad IMU, estimacion Kalman
s=0;
for n=1:N1+1
    if(n <= N1)   
        T2(1,n) = s;
        s = s + diftime(n);
    else
        T2(1,n) = s;
    end
end

% Escala de tiempo datos GNSS
T3 = 0:dtGNSS:(N2-1);

time_imu = T1;
time_gnss = T3;


% Cordenadas GNSS en horizonte local(NED)
for n=1:N2
    [ned_gnss(n,1),ned_gnss(n,2),ned_gnss(n,3)] = geodetic2ned(geod_gnss(n,1),geod_gnss(n,2),geod_gnss(n,3), geod_gnss(1,1),geod_gnss(1,2),geod_gnss(1,3), wgs84Ellipsoid);
end
ned_gnss = ned_gnss + [0, 0, 0.8]; %centro de gravedad del vehiculo

%Eliminar la gravedad de las aceleraciones en ejes body
for n=1:N1
    C_bi = quat2dcm(quat(n,:));
    g_b(n,:) = (C_bi*g)';   
    accel(n,1) = accel(n,1) + g_b(n,1);
    accel(n,2) = accel(n,2) + g_b(n,2);
    accel(n,3) = accel(n,3) - g_b(n,3);
end



% Matriz "C" lineal y constante, obtenemos las variables de posicion
C = [1 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0];

% variables
syms x y z vx vy vz q0 q1 q2 q3 fx fy fz wx wy wz dt
r_s = [x; y; z];
v_s = [vx; vy; vz];
q_s = [q0; q1; q2; q3];
f_s = [fx; fy; fz];
w_s = [wx; wy; wz];
% vector de estado
X_s = [r_s; v_s; q_s];

C_ib=[1-2*(q2^2+q3^2) 2*(-q0*q3+q1*q2) 2*(q0*q2+q1*q3);
      2*(q0*q3+q1*q2) 1-2*(q1^2+q3^2) 2*(-q0*q1+q2*q3);
      2*(-q0*q2+q1*q3) 2*(q0*q1+q2*q3) 1-2*(q1^2+q2^2)];
Q=[-q1 -q2 -q3;
    q0 q3 -q2;
    -q3 q0 q1;
    q2 -q1 q0];

% funciones de las variables de estado
rr = v_s*dt + r_s;
vv =(C_ib*f_s+g)*dt + v_s;
qq = (0.5*Q*C_ib*w_s)*dt + q_s;

% Matriz "A" simbolica
for i=1:10;
   A_s(:,i) = [diff(rr, X_s(i)); diff(vv, X_s(i)); diff(qq, X_s(i))];
end
 


% valores iniciales
r = ned_gnss(1,:);
v = zeros(1,3);
% Valores IMU
for n=1:N1;
    n
    C_bi = quat2dcm(quat(n,:));
    C_ib = C_bi';
      
    r(n+1,:) = v(n,:)*diftime(n) + r(n,:);
    v(n+1,:) = ((C_ib*accel(n,:)')*diftime(n))'+ v(n,:);
%     v(n+1,:) = ((C_ib*accel(n,:)'-g)*diftime(n))'+ v(n,:);
end

% % Grafica quaternions IMU
% figure
% plot(T1, quat(:,1), 'r', T1, quat(:,2), 'g', T1, quat(:,3), 'b', T1, quat(:,4), 'k'); grid;
% xlabel('T (s)');
% ylabel('quaternions');
% title('Quat IMU');
% legend('qW', 'qX', 'qY','qZ');
% 
% % Grafica aceleraciones IMU
% figure
% plot(T1, accel(:,1), 'r', T1, accel(:,2), 'g', T1, accel(:,3), 'b'); grid;
% xlabel('T (s)');
% ylabel('aceleracion (m/s^2)');
% title('Aceleraciones IMU');
% legend('aX', 'aY', 'aZ');
% 
% % Grafica velocidades angulares IMU
% figure
% plot(T1, gyro(:,3), 'r', T1, gyro(:,2), 'g', T1, gyro(:,1), 'b'); grid;
% xlabel('T (s)');
% ylabel('vel angular (rad/s)');
% title('Vel angular IMU');
% legend('ejeYaw', 'ejePitch', 'ejeRoll');
% 
% Grafica angulos euler IMU
figure
plot(T1, euler(:,3), 'r', T1, euler(:,2), 'g', T1, euler(:,1), 'b'); grid;
xlabel('T (s)');
ylabel('angulos euler (º)');
title('Ang Euler IMU');
legend('Yaw', 'Pitch', 'Roll');

% Grafica heading(N-E) IMU
figure
plot(T1, head(:,1), 'k'); grid;
xlabel('T (s)');
ylabel('heading (º)');
title('Rumbo magnetico IMU');
legend('azimut(º)');
% 
% % Grafica posicion IMU
% figure
% plot(T2, r(:,1), 'r', T2, r(:,2), 'g', T2, r(:,3), 'b'); grid;
% xlabel('T (s)');
% ylabel('posicion (m)');
% title('Posicion en ejes NED obtenidos por la IMU');
% legend('rX_imu', 'rY_imu', 'rZ_imu');
% 
% % Grafica velocidad IMU
% figure
% plot(T2, v(:,1), 'r', T2, v(:,2), 'g', T2, v(:,3), 'b'); grid;
% xlabel('T (s)');
% ylabel('velocidad (m/s)');
% title('Velocidad en ejes NED obtenidos por la IMU');
% legend('vX_imu', 'vY_imu', 'vZ_imu');
% 
% %--------------Comparacion posicion IMU-GNSS en ejesNED-------------
% figure
% hold on;
% plot(T2, r(:,1), 'r', T2, r(:,2), 'g', T2, r(:,3), 'b'); grid;
% plot(T3, ned_gnss(:,1), 'r+', T3, ned_gnss(:,2), 'g+', T3, ned_gnss(:,3), 'b+');
% xlabel('T (s)');
% ylabel('posicion (m)');
% title('Posicion NED IMU vs GNSS');
% legend('rX_imu', 'rY_imu', 'rZ_imu', 'rX_gnss', 'rY_gnss', 'rZ_gnss');





% valores iniciales kalman
r_kalman = ned_gnss(1,:);
v_kalman = zeros(1,3);
angRad_kalman(1,1) = deg2rad(euler(1,3)); % yaw
angRad_kalman(1,2) = deg2rad(euler(1,2)); % pitch
angRad_kalman(1,3) = deg2rad(euler(1,1)); % roll

% Matriz de estado y "P" iniciales
X = [r(1,:), v(1,:), quat(1,:)]';
P = Qk;

% Filtro de KALMAN
marg = 0.007; %margen de error en sincronizacion IMU y GNSS 
for n=1:N1;
    n
    %A = double(subs(A_s, [X_s; f_s; w_s], [r(n,:)'; v(n,:)'; q(n,:)'; f_real(n,:)'; w_real(n,:)']));
    A = double(subs(A_s, {x, y, z, vx, vy, vz, q0, q1, q2, q3, fx, fy, fz, wx, wy, wz, dt}, {r(n,1), r(n,2), r(n,3), v(n,1), v(n,2), v(n,3), quat(n,1), quat(n,2), quat(n,3), quat(n,4), accel(n,1), accel(n,2), accel(n,3), gyro(n,1), gyro(n,2), gyro(n,3), diftime(n)}));
    % 1) Prediccion
    Xmenos = A*X;
    Pmenos = A*P*A'+Qk;
    
    su=0;
    for j=1:N2;
        if (((time_gnss(j)-marg/2) <= time_imu(n)) && (time_imu(n) <= (time_gnss(j)+marg/2)))
            j
            su=1;
            % 2) Correccion
            K = Pmenos*C'*(inv(C*Pmenos*C'+ Rk));
            P = (eye(10)-K*C)*Pmenos;
            X = Xmenos+K*(ned_gnss(j,:)'-C*Xmenos);
        elseif (su ~= 1)
            X = Xmenos;
        end
    end
    r_kalman(n+1,:) = X(1:3)';
    v_kalman(n+1,:) = X(4:6)';
    [angRad_kalman(n+1,1), angRad_kalman(n+1,2), angRad_kalman(n+1,3)] = quat2angle(X(7:10)');
end
ang_kalman = rad2deg(angRad_kalman); %Yaw-Pitch-Roll


% % Grafica posicion estimada Kalman
% figure
% plot(T2, r_kalman(:,1), 'r', T2, r_kalman(:,2), 'g', T2, r_kalman(:,3), 'b'); grid;
% xlabel('T (s)');
% ylabel('posicion Kalman (m)');
% title('Posicion en funcion del tiempo estimado');
% legend('rX_{Kalman}', 'rY_{Kalman}', 'rZ_{Kalman}');
% 
% % Grafica velocidad estimada Kalman
% figure
% plot(T2, v_kalman(:,1), 'r', T2, v_kalman(:,2), 'g', T2, v_kalman(:,3), 'b'); grid;
% xlabel('T (s)');
% ylabel('velocidad Kalman (m/s)');
% title('Velocidad en funcion del tiempo estimado');
% legend('vX_{Kalman}', 'vY_{Kalman}', 'vZ_{Kalman}');
% 
% % Graficas angulos de euler estimados Kalman
% figure
% plot(T2, ang_kalman(:,1), 'r', T2, ang_kalman(:,2), 'g', T2, ang_kalman(:,3), 'b'); grid;
% xlabel('T (s)');
% ylabel('Angulos Kalman (º)');
% title('Angulos en funcion del tiempo estimado');
% legend('Yaw_{Kalman}', 'Pitch_{Kalman}', 'Roll_{Kalman}');



%-------------Graficas comparacion Posicion IMU-GNSS-Kalman---------------
% Comparacion posicion ejeX IMU-GNSS-Kalman
figure
plot(T2, r(:,1), 'r', T3, ned_gnss(:,1), 'g+', T2, r_kalman(:,1), 'b'); grid;
xlabel('T (s)');
ylabel('posicion (m)');
title('Posicion eje X: IMU-GNSS-Kalman');
legend('rX_{IMU}', 'rX_{GNSS}', 'rX_{Kalman}');

% Comparacion posicion ejeY IMU-GNSS-Kalman
figure
plot(T2, r(:,2), 'r', T3, ned_gnss(:,2), 'g+', T2, r_kalman(:,2), 'b'); grid;
xlabel('T (s)');
ylabel('posicion (m)');
title('Posicion eje Y: IMU-GNSS-Kalman');
legend('rY_{IMU}', 'rY_{GNSS}', 'rY_{Kalman}');

% Comparacion posicion ejeZ IMU-GNSS-Kalman
figure
plot(T2, r(:,3), 'r', T3, ned_gnss(:,3), 'g+', T2, r_kalman(:,3), 'b'); grid;
xlabel('T (s)');
ylabel('posicion (m)');
title('Posicion eje Z: IMU-GNSS-Kalman');
legend('rZ_{IMU}', 'rZ_{GNSS}', 'rZ_{Kalman}');


%-------------Graficas comparacion Posicion GNSS-Kalman---------------
% Comparacion posicion ejeX GNSS-Kalman
figure
plot(T3, ned_gnss(:,1), 'g+', T2, r_kalman(:,1), 'b'); grid;
xlabel('T (s)');
ylabel('posicion (m)');
title('Posicion eje X: GNSS-Kalman');
legend('rX_{GNSS}', 'rX_{Kalman}');

% Comparacion posicion ejeY GNSS-Kalman
figure
plot(T3, ned_gnss(:,2), 'g+', T2, r_kalman(:,2), 'b'); grid;
xlabel('T (s)');
ylabel('posicion (m)');
title('Posicion eje Y: GNSS-Kalman');
legend('rY_{GNSS}', 'rY_{Kalman}');

% Comparacion posicion ejeZ GNSS-Kalman
figure
plot(T3, ned_gnss(:,3), 'g+', T2, r_kalman(:,3), 'b'); grid;
xlabel('T (s)');
ylabel('posicion (m)');
title('Posicion eje Z: GNSS-Kalman');
legend('rZ_{GNSS}', 'rZ_{Kalman}');


% %-------------Graficas comparacion Velocidad IMU-Kalman---------------
% % Comparacion velocidad ejeX IMU-Kalman
% figure
% plot(T2, v(:,1), 'r', T2, v_kalman(:,1), 'b'); grid;
% xlabel('T (s)');
% ylabel('velocidad (m/s)');
% title('Velocidad eje X: IMU-Kalman');
% legend('vX_IMU', 'vX_Kalman');
% 
% % Comparacion velocidad ejeY IMU-Kalman
% figure
% plot(T2, v(:,2), 'r', T2, v_kalman(:,2), 'b'); grid;
% xlabel('T (s)');
% ylabel('velocidad (m/s)');
% title('Velocidad eje Y: IMU-Kalman');
% legend('vY_IMU', 'vY_Kalman');
% 
% % Comparacion velocidad ejeZ IMU-Kalman
% figure
% plot(T2, v(:,3), 'r', T2, v_kalman(:,3), 'b'); grid;
% xlabel('T (s)');
% ylabel('velocidad (m/s)');
% title('Velocidad eje Z: IMU-Kalman');
% legend('vZ_IMU', 'vZ_Kalman');


%-------------Graficas comparacion Angulos IMU-Kalman---------------
% Comparacion angulos ejeZ IMU-Kalman
figure
plot(T1, euler(:,3), 'r', T2, ang_kalman(:,1), 'b'); grid;
xlabel('T (s)');
ylabel('angulos (º)');
title('Angulos euler eje Z: IMU-Kalman');
legend('Yaw', 'Yaw_{Kalman}');

% Comparacion angulos ejeY IMU-Kalman
figure
plot(T1, euler(:,2), 'r', T2, ang_kalman(:,2), 'b'); grid;
xlabel('T (s)');
ylabel('angulos (º)');
title('Angulos euler eje Y: IMU-Kalman');
legend('Pitch', 'Pitch_{Kalman}');

% Comparacion angulos ejeX IMU-Kalman
figure
plot(T1, euler(:,1), 'r', T2, ang_kalman(:,3), 'b'); grid;
xlabel('T (s)');
ylabel('angulos (º)');
title('Angulos euler eje X: IMU-Kalman');
legend('Roll', 'Roll_{Kalman}');





%-------------Graficas comparacion recorrido horizontal-------------
% Comparacion posicion horizontal GNSS-Kalman
figure
hold on;
plot(ned_gnss(:,2), ned_gnss(:,1), 'c.'); grid;
plot(r_kalman(:,2), r_kalman(:,1), 'm');
axis equal
xlabel('posY_{NED} (m)');
ylabel('posX_{NED} (m)');
title('Posicion ejes X-Y: NED');
legend('Recorrido_{GNSS}', 'Recorrido_{Kalman}');



% Calcular la velocidad horizontal(abs) de Kalman
for n=1:N1+1;
    speed_kalman(n) = sqrt(v_kalman(n,1)^2 + v_kalman(n,2)^2);
end

%------------Grafica velocidad horizontal(abs) GNSS-Kalman-----------
figure
hold on;
plot(T3, speed_gnss(:,1), 'c'); grid;
plot(T2, speed_kalman, 'm-.');
xlabel('T (s)');
ylabel('velocidad (m/s)');
title('Velocidad horizontal(abs) GNSS-Kalman');
legend('groundSpeed_{gnss}', 'groundSpeed_{kalman}');



% Calcular la distancia GNSS
dist_GNSS=0;
for n=1:N2-1;
    dist_GNSS = dist_GNSS + sqrt((ned_gnss(n+1,1)-ned_gnss(n,1))^2 + (ned_gnss(n+1,2)-ned_gnss(n,2))^2 + (ned_gnss(n+1,3)-ned_gnss(n,3))^2);
end

% Calcular la distancia estimada Kalman
dist_Kalman=0;
for n=1:N1;
    dist_Kalman = dist_Kalman + sqrt((r_kalman(n+1,1)-r_kalman(n,1))^2 + (r_kalman(n+1,2)-r_kalman(n,2))^2 + (r_kalman(n+1,3)-r_kalman(n,3))^2);
end


%**********************Datos obtenemos************************
% Datos IMU
xlswrite('vuelta1_resultados.xlsx', T2', 'estimacion_IMU', 'A3:A100000');
xlswrite('vuelta1_resultados.xlsx', r, 'estimacion_IMU', 'B3:D100000');
xlswrite('vuelta1_resultados.xlsx', v, 'estimacion_IMU', 'E3:G100000');

% Datos Kalman
xlswrite('vuelta1_resultados.xlsx', T2', 'estimacion_Kalman', 'A3:A100000');
xlswrite('vuelta1_resultados.xlsx', r_kalman, 'estimacion_Kalman', 'B3:D100000');
xlswrite('vuelta1_resultados.xlsx', v_kalman, 'estimacion_Kalman', 'E3:G100000');
xlswrite('vuelta1_resultados.xlsx', ang_kalman, 'estimacion_Kalman', 'H3:J100000');

% Distancia total recorrida
xlswrite('vuelta1_resultados.xlsx', dist_GNSS, 'dist_vuelta', 'B2');
xlswrite('vuelta1_resultados.xlsx', dist_Kalman, 'dist_vuelta', 'C2');







