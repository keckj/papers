
% Main function
function question8
    
    clear all;

    % Solver configuration
    dt   = 0.001; %s
    tmax = 0.0  ; %s

    % Unknown true state
    u_t = [-4.62, -6.61, 17.94]';        %true initial state of the system
    Y_t = RK4(u_t, @Lorenz, dt, tmax);  %true state of the system (true trajectory, reference state)
    
    % Initial guess
    u_b = [-5   , -7   , 17   ]';        % a priori initial state of the system
    Y_b = RK4(u_b, @Lorenz, dt, tmax);  % a priori state of the system (background trajectory)
    
    % Generate observations
    Y_obs = Y_t;                        % observations are perfect
    
    %Kalman filtering
    A = [1,2;3,4]
    B = [1,0;0,1]
    
end

function [x_a, x_f, P_a, P_f] = KalmanFilter(x_b, P_b, n)
    
    % Initialization
    x_a = x_b;
    P_a = P_b;
    
    % Forecast - Correction
    for k=1:n
        Rk = R(k);
        Qk = Q(k);
        Hk = H(k);
        y_obs = Y(k);
        
        % Forecast
        x_f = Mk*x_a;
        P_f = Mk*P_a*Mk' + Qk;
        
        % Analysis
        K = P_f * Hk' \ (Hk * P_f * Hk' + Rk);
        x_a = x_f + K * (y_obs - Hk * x_f);
        P_a = P_f - K * Hk * P_f;
    end
end