clear all;
close all;
format long;

% valores que nos dan
diftime = xlsread('Prueba1/Prueba0.ods', 'Data', 'C3:C1934'); %OJO adelantar una fila la difftime
quat = xlsread('Prueba1/Prueba0.ods', 'Data', 'D2:G1933');
accel = xlsread('Prueba1/Prueba0.ods', 'Data', 'H2:J1933');
gyro = xlsread('Prueba1/Prueba0.ods', 'Data', 'K2:M1933');
euler = xlsread('Prueba1/Prueba0.ods', 'Data', 'N2:P1933');    
head = xlsread('Prueba1/Prueba0.ods', 'Data', 'Q2:Q1933');
% r_gps = xlsread('Cl2.xlsx', 'REAL', 'H4:J504');
% Qk = xlsread('Cl2.xlsx', 'Qk', 'A1:J10');
% Rk = xlsread('Cl2.xlsx', 'Rk', 'A1:C3');
g = [0; 0; 9.79991];

% numero de ciclos
N = length(diftime);
% % dt entre medidas
% dtMedidasGPS = 0.1;
% escalas de tiempos
s=0;
for n=1:N
    T1(1,n) = s;
    s = s + diftime(n);
end

s=0;
for n=1:N+1
    if(n <= N)   
        T2(1,n) = s;
        s = s + diftime(n);
    else
        T2(1,n) = s;
    end
end


% Matriz "C" lineal y constante, obtenemos las variables de posicion
C = [1 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0];

% variables
syms x y z vx vy vz q0 q1 q2 q3 fx fy fz wx wy wz
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
syms dt %diftime
rr = v_s*dt + r_s;
vv =((C_ib*f_s+g)*dt) + v_s;
qq = ((0.5*Q*C_ib*w_s)*dt) + q_s;

% Matriz "A" simbólica
for i=1:10;
   A_s(:,i) = [diff(rr, X_s(i)); diff(vv, X_s(i)); diff(qq, X_s(i))];
end
 
% valores iniciales
r = zeros(1,3);
v = zeros(1,3);
% ang(1,:) = [0.428116 1.00897 49.1488];
% ang_rad = ang*pi/180;
% q = angle2quat(ang_rad(1,3), ang_rad(1,2), ang_rad(1,1)); %OJO a la rotacion
% r_estim = zeros(1,3);
% v_estim = zeros(1,3);
% ang_rad_estim = zeros(1,3);


% % Valores IMU(Comprobacion con euler y quat de la IMU)
% for n=1:N;
%     n
%     C_bi = quat2dcm(q(n,:));
%     C_ib = C_bi';
%     Q = [-q(n,2) -q(n,3) -q(n,4);
%           q(n,1) q(n,4) -q(n,3);
%           -q(n,4) q(n,1) q(n,2);
%           q(n,3) -q(n,2) q(n,1)];
%       
%     dift = diftime(n);
%     r(n+1,:) = v(n,:)*dift + r(n,:);
% %     v(n+1,:) = ((C_ib*accel(n,:)')*dift)'+ v(n,:);
%     v(n+1,:) = ((C_ib*accel(n,:)'-g)*dift)'+ v(n,:);
%     q(n+1,:) = ((0.5*Q*C_ib*gyro(n,:)')*dift)'+ q(n,:);
%     [ang_rad(n+1,3), ang_rad(n+1,2), ang_rad(n+1,1)] = quat2angle(q(n+1,:));
% end
% ang = ang_rad*180/pi;


% Valores IMU
for n=1:N;
    n
    C_bi = quat2dcm(quat(n,:));
    C_ib = C_bi';
      
    dift = diftime(n);
    r(n+1,:) = v(n,:)*dift + r(n,:);
%     v(n+1,:) = ((C_ib*accel(n,:)')*dift)'+ v(n,:);
    v(n+1,:) = ((C_ib*accel(n,:)'-g)*dift)'+ v(n,:);
end


% Gráfica angulos euler IMU
figure
plot(T1, euler(:,3), 'r', T1, euler(:,2), 'g', T1, euler(:,1), 'b'); grid;
xlabel('T (s)');
ylabel('angulos euler (º)');
title('Ang Euler IMU');
legend('Yaw', 'Pitch', 'Roll');

% Gráfica heading IMU
figure
plot(T1, head(:,1), 'k'); grid;
xlabel('T (s)');
ylabel('heading (º)');
title('Rumbo magnetico IMU');
legend('psi(º)');

% Gráfica posicion IMU
figure
plot(T2, r(:,1), 'r', T2, r(:,2), 'g', T2, r(:,3), 'b'); grid;
xlabel('T (s)');
ylabel('posicion (m)');
title('Posicion en ejes inerciales obtenidos por la IMU');
legend('r_X', 'r_Y', 'r_Z');

% Gráfica velocidad IMU
figure
plot(T2, v(:,1), 'r', T2, v(:,2), 'g', T2, v(:,3), 'b'); grid;
xlabel('T (s)');
ylabel('velocidad (m/s)');
title('Velocidad en ejes inerciales obtenidos por la IMU');
legend('v_X', 'v_Y', 'v_Z');




% % Matriz de estado y "P" iniciales
% X = [r(1,:), v(1,:), q(1,:)]';
% P = Qk;
% 
% % Filtro de KALMAN
% j=0;
% for n=1:N;
%     %A = double(subs(A_s, [X_s; f_s; w_s], [r(n,:)'; v(n,:)'; q(n,:)'; f_real(n,:)'; w_real(n,:)']));
%     A = double(subs(A_s, {x, y, z, vx, vy, vz, q0, q1, q2, q3, fx, fy, fz, wx, wy, wz}, {r(n,1), r(n,2), r(n,3), v(n,1), v(n,2), v(n,3), q(n,1), q(n,2), q(n,3), q(n,4), f_real(n,1), f_real(n,2), f_real(n,3), w_real(n,1), w_real(n,2), w_real(n,3)}));
%     % 1) Predicción
%     Xmenos = A*X;
%     Pmenos = A*P*A'+Qk;
%     j=j+1;
%     if j == dtMedidas/dt;
%         j=0;
%         % 2) Corrección
%         K = Pmenos*C'*(inv(C*Pmenos*C'+ Rk));
%         P = (eye(10)-K*C)*Pmenos;
%         X = Xmenos+K*(r_gps(n,:)'-C*Xmenos);
%     else
%         X = Xmenos;
%     end
%     r_estim(n,:) = X(1:3)';
%     v_estim(n,:) = X(4:6)';
%     [ang_rad_estim(n,1), ang_rad_estim(n,2), ang_rad_estim(n,3)] = quat2angle(X(7:10)');
% end
% 
% ang = ang_rad*180/pi;
% ang_estim = ang_rad_estim*180/pi;
% 
% % % Gráficos
% % % Gráfica 3D posición estimada
% % figure
% % comet3(r_estim); grid;
% % xlabel('Rx estim(m)');
% % ylabel('Ry estim(m)');
% % zlabel('Rz estim(m)');
% % title('Trayectoria de la aeronave ESTIMADO');
% % 
% % % Gráfica posicion estimada
% % figure
% % plot(T1, r_estim(:,1), 'r', T1, r_estim(:,2), 'g', T1, r_estim(:,3), 'b'); grid;
% % xlim([0 N/100]);
% % xlabel('T (s)');
% % ylabel('R (m)');
% % title('Posicion en funcion del tiempo estimado');
% % legend('Rx_est', 'Ry_est', 'Rz_est');
% % 
% % % Gráficas angulos de euler estimados
% % figure
% % plot(T1, ang_estim(:,1), 'r', T1, ang_estim(:,2), 'g', T1, ang_estim(:,3), 'b'); grid;
% % xlim([0 N/100]);
% % xlabel('T (s)');
% % ylabel('Angulos (º)');
% % title('Angulos en funcion del tiempo estimado');
% % legend('Yaw_est', 'Pitch_est', 'Roll_est');
% % 
% % % Comparacion con datos de la imu (posicion)
% % figure
% % hold on;
% % plot(T1, r_estim(:,1), 'r+-', T1, r_estim(:,2), 'g+-', T1, r_estim(:,3), 'b+-'); grid;
% % plot(T2, r_gps(:,1), 'ro-', T2, r_gps(:,2), 'go-', T2, r_gps(:,3), 'bo-');
% % plot(T2, r(:,1), 'r', T2, r(:,2), 'g', T2, r(:,3), 'b');
% % xlim([0 N/100]);
% % xlabel('T (s)');
% % ylabel('Re Real/Medido/Estimado (m)');
% % title('Comparacion posicion: Real/Medida GPS/Estimada');
% % legend('Rx_estimado', 'Ry_estimado', 'Rz_estimado', 'Rx_gps', 'Ry_gps', 'Rz_gps', 'Rx', 'Ry', 'Rz');
% % 
% % % Comparacion con datos de la imu (angulos)
% % figure
% % hold on;
% % plot(T1, ang_estim(:,1), 'r+-', T1, ang_estim(:,2), 'g+-', T1, ang_estim(:,3), 'b+-'); grid;
% % plot(T2, ang(:,1), 'r', T2, ang(:,2), 'g', T2, ang(:,3), 'b');
% % xlim([0 N/100]);
% % xlabel('T (s)');
% % ylabel('Ang Real/Estimado (º)');
% % title('Comparacion angulos: Real/Estimados');
% % legend('Yaw_estimado', 'Pitch_estimado', 'Roll_estimado', 'Yaw', 'Pitch', 'Roll');
% % 
% % % Comparacion con datos de la imu (velocidad)
% % figure
% % hold on;
% % plot(T1, v_estim(:,1), 'r+-', T1, v_estim(:,2), 'g+-', T1, v_estim(:,3), 'b+-'); grid;
% % plot(T2, v(:,1), 'r', T2, v(:,2), 'g', T2, v(:,3), 'b');
% % xlim([0 N/100]);
% % xlabel('T (s)');
% % ylabel('Velocidad Real/Estimado (º)');
% % title('Comparacion velocidad: Real/Estimada');
% % legend('Vx_estimado', 'Vy_estimado', 'Vz_estimado', 'Vx', 'Vy', 'Vz');
% % 
% % % % Datos obtenemos
% % % xlswrite('datos estimados KALMAN.xls', r_estim, 'B4:D504');
% % % xlswrite('datos estimados KALMAN.xls', ang_estim, 'E4:G504');
% % % xlswrite('datos estimados KALMAN.xls', v_estim, 'H4:J504');