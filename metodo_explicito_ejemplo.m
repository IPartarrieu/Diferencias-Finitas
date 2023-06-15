%% Tarea 1: Diferencias finitas - MÉTODO EXPLÍCITO
clear all
close all
%
% Coeficiente de difusión
k = 1e-6; % [m^2s^-1]

% Definimos las vairables de ingreso
L = input('\nIngrese distancia de la roca en metros: ');  % Distancia de la roca
disp([num2str(L),' [m]']);
%
R = input('\nIngrese temperatura de la roca en celcius: ');  % Temperatura roca
disp([num2str(R),' [ºC]']);
%
W = input('\nIngrese ancho del dique en metros: ');  % Distancia de dique
disp([num2str(W),' [m]']);
%
D = input('\nIngrese temperatura del dique en celcius: ');  % Temperatura dique
disp([num2str(D),' [ºC]']);
%
T = input('\nIngrese el tiempo de simulación en días: ');  % Tiempo de simulación
T = T*86400; % pasamos los días a segundos
disp([num2str(T),' [s]']);

%%
% % % Caso de ejemplo
% % L = 100;
% % R = 300;
% % W = 5;
% % D = 1200;
% % T = 2*24*3600;
% % disp('VARIABLES DE PRUEBA LISTAS');
%%
nx = input('\nIngrese los puntos de grillas espaciales a simular: ');
dx = L/(nx-1);
%
NT = condicion_cfl(nx,L,T,k);
disp(['Para que se satisfaga la condición de CFL nt debe ser mayor que ', num2str(NT,'%0.2f')])
%
nt = input('\nIngrese el paso de tiempo para simular: ');
dt = T/(nt-1);
%
CFL = dx^2/(2*k);
if dt < CFL
    disp('Se cumple la condición de CFL.');
elseif dt >= CFL
    disp('No se cumple la condición de CFL.');
end

%% Tiempo 0
tic
% Creamos la grilla
x = -L/2:dx:L/2;
% Configuramos el perfil de temperatura inicial (u^l=0 para todo i)
u = ones(size(x))*R;
u(find(abs(x)<=W/2)) = D;
t = 0; % partimos desde tiempo 0
%
%% Solución explícita
for n = 1:nt
    u_new = zeros(1,nx);
    for i = 2:nx-1
        % ecuación explícita de calor discretizada (ecuación de calor transitoria dependiente del tiempo)
        u_new(i) = u(i) + ((k*dt)/dx^2)*(u(i+1)-2*u(i)+u(i-1));
    end
    % Establecemos nuevamente las condiciones de borde
    u_new(1) = u(1);
    u_new(nx) = u(nx);
    % Restablecemos temperatura y tiempo
    u = u_new;
    t = t+dt;
    % Ploteamos resultado
    figure(1),clf
    plot(x,u_new,'LineWidth',1.5)
    xlabel('x [m]','FontSize',16)
    ylabel('Temperatura [ºC]','FontSize',16)
    title(['dx: ',num2str(dx,'%0.2f'),' [m] - dt: ',num2str(dt/3600,'%0.2f'),' [hrs]'],'FontSize',20)
    axis([min(x) max(x) R D])
    grid minor
    drawnow
end
TT = toc;
disp(['Tiempo total de ejecución: ',num2str(TT)]);
%
%