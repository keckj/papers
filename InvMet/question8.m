
% Main function
function question8
    
    clear all;
    digits(32);
    
    % Solver configuration
    dt   = 0.01; %s
    Tmax = 3;  %s
    n = ceil(Tmax/dt);

    % Unknown true state
    u_t = [-4.62, -6.61, 17.94]';        %true initial state of the system
    Y_t = RK4(u_t, @Lorenz, dt, Tmax);   %true state of the system (true trajectory, reference state)
    
    % Initial guess
    u_b = [-5   , -7   , 17   ]';        % a priori initial state of the system
    Y_b = RK4(u_b, @Lorenz, dt, Tmax);   % a priori state of the system (background trajectory)     
    P_b = 1.0*eye(3);
    
    % Generate observations
    Yobs = Y_t;                         % observations are perfect
    
    % Kalman filtering
    [Xa, Xf, Xobs] = KalmanFilter(u_b, P_b, Yobs, dt, Tmax, ...
        @GetObservations, @ModelMatrix, @ModelErrorCovariance, ...
        @ObservationOperator, @ObservationErrorCovariance);
    
    % Print last value
    custom_plot(Y_t, Y_b, Xa);   
end    

function Yk = GetObservations(k, Y_obs) 
    Yk = [Y_obs(2,k),Y_obs(3,k),Y_obs(4,k)]';
end


function Fk = F(u_k)
    xk = u_k(1);
    yk = u_k(2); 
    zk = u_k(3);
    Fk = [-10       ,+10  , +0   ; ...
         +28-zk/2.0 ,-1   , -xk/2; ...
         +yk/2      ,+xk/2, -8/3];
end

function Mk = ModelMatrix(k,x_a,dt)
    u_k = x_a;

    %F defined such that f(x) = Fx     
    %RK4 Finding Runge Kutta 4 Linear Model
    % xk+1 = xk + dt/6(k1+ 2*k2 + 2*k3 + k4)
    % k1 = f(xk)         <=> k1 = F1*xk
    % k2 = f(xk+dt/2*k1) <=> k2 = F2*(xk+dt/2*k1)
    % k3 = f(xk+dt/2*k2) <=> k3 = F3*(xk+dt/2*k2);
    % k2 = f(xk+dt  *k3) <=> k4 = F4*(xk+dt  *k3)
    % k1 + 2*k2 + 2*k3 + k4
    % = 6*F*xk + dt*F*  [k1 + k2 + k3]
    % = 6*F*xk + dt*F*  [F*xk + F*(xk+dt/2*k1) + F*(xk+dt/2*k2)]
    % = 6*F*xk + 3*dt*F*F*xk + 1/2*dt*dt*F*F* [k1+k2]
    % = 6*F*xk + 3*dt*F*F*xk + 1/2*dt*dt*F*F* [F*xk + F*(xk+dt/2*k1)]
    % = 6*F*xk + 3*dt*F*F*xk + 1*dt*dt*F*F*F*xk + 1/4*dt*dt*dt*F*F*F*F*xk
    % => dt/6*(k1+2*k2+2*k3+k4) = Ak*xk
    % with Ak = dt*F + 1/2*dt*dt*F*F + 1/6*dt*dt*dt*F*F*F + 1/24*dt*dt*dt*dt*F*F*F*F

    % xk+1 = (Ak+I)*xk
    %Ak = A1 + 1/2*A2 + 1/6*A3 + 1/24*A4;    
    
    F1 =   F(u_k);
    k1 = F1*(u_k);
    u1 = u_k + dt/2*k1;
    
    F2 =   F(u1);
    k2 = F2*(u1);
    u2 = u_k + dt/2*k2;
    
    F3 =   F(u2);
    k3 = F3*(u2);
    u3 = u_k + dt  *k3;
    
    F4 =   F(u3);
    k3 = F4*(u3);
    
    I = eye(3); 
    A1 = I + dt*F1; A2 = dt*F2; A3 = dt*F3; A4 = dt*F4;    
    Mk = A1 + 1/2.0*A2*A1 + 1/3.0*A3*A2*A1 + 1/4.0*A4*A3*A2*A1;
end

function Qk = ModelErrorCovariance(k)
    Qk = 1.0*eye(3);
end

function Hk = ObservationOperator(k)
    Hk = eye(3);
end

function Rk = ObservationErrorCovariance(k)
    Rk = 0.1*eye(3);
end


function [Xa, Xf, Xobs] = KalmanFilter(x_b, P_b, Yobs, dt, Tmax, ...
    GetObservations, ...
    ModelMatrixes, ModelErrorCovarianceMatrixes, ...
    ObservationOperatorMatrixes, ObservationErrorCovarianceMatrixes)
    
    % Initialization
    n = ceil(Tmax/dt);
    x_a = x_b;
    P_a = P_b;
    
    Xa = [];
    Xf = [];
    Xobs = [];
    t = 0;
 
    % Forecast - Correction
    for k=1:n

        Mk = ModelMatrixes(k,x_a,dt);
        Qk = ModelErrorCovarianceMatrixes(k);
        Hk = ObservationOperatorMatrixes(k);
        Rk = ObservationErrorCovarianceMatrixes(k);
        Yk = GetObservations(k+1, Yobs);
        
        % Forecast
        x_f = Mk*x_a;
        P_f = Mk*P_a*Mk' + Qk;
                           
        % Analysis X_f
        %Kk = P_f * Hk' \ (Hk * P_f * Hk' + Rk);
        Kk = 0.0;
        x_a = x_f + Kk * (Yk - Hk * x_f);
        P_a = P_f - Kk * Hk * P_f;
        
        % Keep values          
        Xobs = [Xobs, [t; Hk*Yk]];
        Xf =   [Xf  , [t; x_f  ]];        
        Xa =   [Xa  , [t; x_a  ]];         
        t = t+dt;   
    end
end

% Plot solution 
function custom_plot(u_t, u_b, u_a) 
    tt = u_t(1,:);      
    xt = u_t(2,:);
    yt = u_t(3,:);
    zt = u_t(4,:);        
    tb = u_b(1,:);
    xb = u_b(2,:);
    yb = u_b(3,:);
    zb = u_b(4,:);
    ta = u_a(1,:);
    xa = u_a(2,:);
    ya = u_a(3,:);
    za = u_a(4,:);

    subplot(3,5,[1,2]);
    plot(tt,xt,'blue',tb,xb,'green',ta,xa,'red'); 
    legend('xt','xb','xa');
    title('x(t)');
    xlabel('t');
    ylabel('x');
    
    subplot(3,5,[6,7]);
    plot(tt,yt,'blue',tb,yb,'green',ta,ya,'red'); 
    legend('yt','yb','ya');
    title('y(t)');
    xlabel('t');
    ylabel('y');
    
    subplot(3,5,[11,12]);
    plot(tt,zt,'blue', tb,zb,'green', ta,za,'red');   
    legend('zt','zb','za');
    title('z(t)');
    xlabel('t');
    ylabel('z');
    
    subplot(3,5,[3,4,5, 8,9,10, 13,14,15]);
    plot(xt,zt,'blue', xb,zb,'green', xa,za,'red'); 
    legend('zt(xt)','zb(xb)','za(xa)');
    title('Butterfly z(x)');
    xlabel('x');
    ylabel('z');
end