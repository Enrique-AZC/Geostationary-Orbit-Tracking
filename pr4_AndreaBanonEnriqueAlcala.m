%%  PRÁCTICA 4: Filtro Extendido de Kalman para seguimiento de una órbita geostacionaria
% --------------------------------------------------------------
% Autores: Enrique Alcalá-Zamora Castro
%          Andrea Bañón Triguero
% --------------------------------------------------------------
% Descripción:
% En esta práctica se implementa un Filtro Extendido de Kalman (EKF)
% para estimar la trayectoria de un satélite en órbita geostacionaria
% usando mediciones ruidosas del radio.
%
% - Predicción del estado con dinámica orbital.
% - Corrección con mediciones del radio.
% - Uso del Jacobiano para linealizar el sistema.
% - Comparación con la órbita real y análisis del error.
% --------------------------------------------------------------

% Cargar las variables del archivo Variables_Orbita.mat
load('Variables_Orbita.mat');

% Parámetros iniciales
GM = GM; % Constante gravitacional
T = 0.0001; % Intervalo de predicción
Tmed = 0.1; % Intervalo de medición
N = length(tt); % Número de mediciones
steps = Tmed / T; % Número de predicciones entre mediciones

% Inicialización del estado
x0 = 42164; % km
y0 = 0; % km
vx0 = 0; % km/h
vy0 = 11068; % km/h
x_est = [x0; y0; vx0; vy0];

% Matriz de covarianza inicial
P = diag([50, 50, 1, 1]);

% Matrices de ruido
Q = [0 0 0 0; 0 0 0 0; 0 0 sig2w 0; 0 0 0 sig2w]*T; % Ruido de proceso
R = sig2v; % Ruido de medición

% Inicialización de almacenamiento de resultados
x_estimado = zeros(4, N);
radio_estimado = zeros(1, N);
error_estimado_k = zeros(1, N);

med = 1; % Contador de mediciones

% Bucle de Filtro de Kalman
for t = T:T:24
    r = norm(x_est(1:2));
    a = -GM / r^3 * x_est(1:2);
    x_est = x_est + T * [x_est(3); x_est(4); a]; % Integración por Euler

    % Jacobiano del sistema
    Gx = [0 0 1 0;
          0 0 0 1;
          ((-GM/r^3 + 3*GM*x_est(1)^2/r^5) * T), ((3*GM*x_est(1)*x_est(2)/r^5) * T), 0, 0;
          ((3*GM*x_est(1)*x_est(2)/r^5) * T), ((-GM/r^3 + 3*GM*x_est(2)^2/r^5) * T), 0, 0];
    Fx = eye(4) + Gx;

    % Predicción de la covarianza
    P = Fx * P * Fx' + Q;

    if mod(t, Tmed)==0
        % Medición estimada
        z_pred = norm(x_est(1:2));

        % Jacobiano de la medición
        H = [x_est(1)/z_pred, x_est(2)/z_pred, 0, 0];

        % Ganancia de Kalman
        K = P * H' / (H * P * H' + R);

        % Corrección del estado
        x_est = x_est + K * (zmed(med) - z_pred);

        % Actualización de la covarianza
        P = P - K * H * P;

        % Almacenamiento de resultados
        x_estimado(:, med) = x_est;
        radio_estimado(med) = norm(x_est(1:2));
        error_estimado_k(med) = sqrt(P(1,1) + P(2,2));

        % Incrementar contador de mediciones
        med = med + 1;
    end
end

% Gráfica: Órbita estimada vs. verdadera
figure;
plot(Cxy_true(:,1), Cxy_true(:,2), 'b', 'DisplayName', 'Órbita Verdadera');
hold on;
plot(x_estimado(1, :), x_estimado(2, :), 'r', 'DisplayName', 'Órbita Estimada');
legend;
xlabel('X (km)'); ylabel('Y (km)');
title('Órbita estimada vs. verdadera');
grid on;

% Gráfica: Radio estimado vs. verdadero vs. medido
figure;
plot(tt, sqrt(sum(Cxy_true.^2, 2)), 'b', 'DisplayName', 'Radio Verdadero');
hold on;
plot(tt, radio_estimado, 'r', 'DisplayName', 'Radio Estimado');
plot(tt, zmed, 'g--', 'DisplayName', 'Radio Medido');
legend;
xlabel('Tiempo (h)'); ylabel('Radio (km)');
title('Evolución del radio: Estimado, Verdadero y Medido');
grid on;

% Calcular error de posición real
error_posicion = sqrt(sum((Cxy_true - x_estimado(1:2, :)').^2, 2));

% Gráfica: Error de posición vs. error estimado
figure;
plot(tt, error_posicion, 'b', 'DisplayName', 'Error Real');
hold on;
plot(tt, error_estimado_k, 'r--', 'DisplayName', 'Error Estimado');
legend;
xlabel('Tiempo (h)'); ylabel('Error (km)');
title('Error de posición vs. estimado');
grid on;

% Comparación de desviaciones estándar
desv_std_ruido = sqrt(sig2v);
desv_std_error_medida = std(zmed - radio_estimado);

