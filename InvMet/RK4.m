
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge-Kutta 4th order solver %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y_0 unitial condition
% f the function of the system (u' = f(u))
% dt time step
% tmax stop time
function Y_t = RK4(Y_0, f, dt, tmax)
    
    Y = Y_0;
    t = 0.0;
    
    Y_t = [t;Y_0];
    n = ceil(tmax/dt);
    
     for k=1:n
        t = t + dt;
        k1 = f(Y);
        k2 = f(Y+dt/2*k1);
        k3 = f(Y+dt/2*k2);
        k4 = f(Y+dt*k3);
        Y = Y + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
        Y_t = [Y_t, [t;Y]];
    end
end


