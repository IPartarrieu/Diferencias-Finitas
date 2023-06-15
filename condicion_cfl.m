function NT = condicion_cfl(nx,L,T,k)
% caso específico para la tarea
%k = 1e-6;
%L = 100;
%T = 2*86400
dx = L/(nx-1);
CFL = (dx^2)/(2*k);
NT = T/CFL + 1;
end