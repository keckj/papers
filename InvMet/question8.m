
% Main function
function question8
    
    clear all;
    digits(32);
    
    % Solver configuration
    dt   = 0.001; %s
    Tmax = 30;  %s
    n = ceil(Tmax/dt);

    % Unknown true state
    u_t = [-4.62, -6.61, 17.94]';        %true initial state of the system
    Y_t = RK4(u_t, @Lorenz, dt, Tmax);   %true state of the system (true trajectory, reference state)
    
    % Initial guess
    u_b = [-5   , -7   , 17   ]';        % a priori initial state of the system
    Y_b = RK4(u_b, @Lorenz, dt, Tmax);   % a priori state of the system (background trajectory)     
    P_b = 1.0*eye(3);
    
    % Generate observations
    Yobs = Y_t;   % these observations are perfect

    % Kalman filtering
    [Xa, Xf, Xobs] = KalmanFilter(u_b, P_b, Yobs, dt, Tmax, ...
       @GetObservations, @ModelMatrix, @ModelErrorCovariance, ...
        @ObservationOperator, @ObservationErrorCovariance);

    % Print last value
    custom_plot(Y_t, Xobs, Y_b, Xf, Xa);   
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

    %Matrix Fk is defined such that f(xk) = Fk.xk     
    %RK4 Finding Runge Kutta 4 Model Matrix
    % xk+1 = xk + dt/6(k1+ 2*k2 + 2*k3 + k4)
    % k1 = f(xk)         <=> k1 = F0*xk
    % k2 = f(xk+dt/2*k1) <=> k2 = F1*(xk+dt/2*k1)
    % k3 = f(xk+dt/2*k2) <=> k3 = F2*(xk+dt/2*k2);
    % k2 = f(xk+dt  *k3) <=> k4 = F3*(xk+dt  *k3)
    % => dt/6*(k1+2*k2+2*k3+k4) = Ak*xk (1)
    % => xk+1 = (Ak+I)*xk 
    
    F0 = F(u_k);
    k1 = F0*u_k;
    u1 = u_k + dt/2*k1;
    
    F1 = F(u1);
    k2 = F1*u1;
    u2 = u_k + dt/2*k2;
    
    F2 = F(u2);
    k3 = F2*u2;
    u3 = u_k + dt  *k3;
    
    F3 = F(u3);
    k4 = F3*(u3);
    
    I = eye(3); %identity

    % "Just" expand all the terms in (1)
    % 3 hours later you should obtain :
    Ak = 1/24 * (                                 ...
         1*dt*dt*dt*dt*(F3*F2*F1*F0)            + ...
         2*   dt*dt*dt*(F2*F1*F0 + F3*F2*F1)    + ...
         4*      dt*dt*(F1*F0 + F2*F1 + F3*F2)  + ...
         4*         dt*(F0+F3)                  + ...
         8*         dt*(F1+F2)                    ...
        );
    
    Mk = I + Ak;
end

function Qk = ModelErrorCovariance(k)
    Qk = 0.5*eye(3);
end

function Hk = ObservationOperator(k)
    Hk = eye(3);
end

function Rk = ObservationErrorCovariance(k)
    Rk = 1.0*eye(3);
end


function [Xa, Xf, Xobs] = KalmanFilter(x_b, P_b, Yobs, dt, Tmax, ...
    GetObservations, ...
    ModelMatrixes, ModelErrorCovarianceMatrixes, ...
    ObservationOperatorMatrixes, ObservationErrorCovarianceMatrixes)
    
    % Initialization
    n = ceil(Tmax/dt);
    x_a = x_b;
    P_a = P_b;
    
    Xa = [0;x_b];
    Xf = [];
    Xobs = [];
    t = 0;
 
    % Forecast - Correction
    for k=1:n+1

        Mk = ModelMatrixes(k,x_a,dt);
        Qk = ModelErrorCovarianceMatrixes(k);
        Hk = ObservationOperatorMatrixes(k);
        Rk = ObservationErrorCovarianceMatrixes(k);
        Yk = GetObservations(k, Yobs);
        
        % Forecast
        x_f = Mk*x_a;
        P_f = Mk*P_a*Mk' + Qk;
                           
        % Analysis X_f
       
        HRi = Hk'/ Rk;
        Kk = (inv(P_f) + HRi * Hk) \ HRi;
        %Kk = P_f * Hk' \ (Hk * P_f * Hk' + Rk);
        %Kk = 0.0; (to check if Mk is really a RK4 step ...)
        
        x_a = x_f + Kk * (Yk - Hk * x_f);
        P_a = P_f - Kk * Hk * P_f;
        
        % Keep values
        Xobs = [Xobs, [t; Hk*Yk]];
        Xf =   [Xf  , [t; x_f  ]];    
        
        if(k == n+1)
            continue;
        end;
        
        t = t+dt;   
        Xa = [Xa  , [t; x_a  ]];  
    end
end

% Plot solution 
function custom_plot(u_t, u_obs, u_b, u_f, u_a) 
    tt = u_t(1,:);      
    xt = u_t(2,:);
    yt = u_t(3,:);
    zt = u_t(4,:);
    tobs = u_obs(1,:);      
    xobs = u_obs(2,:);
    yobs = u_obs(3,:);
    zobs = u_obs(4,:);        
    tb = u_b(1,:);
    xb = u_b(2,:);
    yb = u_b(3,:);
    zb = u_b(4,:);
    tf = u_f(1,:);
    xf = u_f(2,:);
    yf = u_f(3,:);
    zf = u_f(4,:);
    ta = u_a(1,:);
    xa = u_a(2,:);
    ya = u_a(3,:);
    za = u_a(4,:);

    subplot(4,7,[1,2]);
    plot(tt,xt,'blue',tobs, xobs, 'black', tb,xb,'green',ta,xa,'red'); 
    legend('xt','xobs','xb','xa');
    title('x(t)');
    xlabel('t');
    ylabel('x');
    
    subplot(4,7,[8,9]);
    plot(tt,yt,'blue',tobs, yobs, 'black', tb,yb,'green',ta,ya,'red'); 
    legend('yt','yobs','yb','ya');
    title('y(t)');
    xlabel('t');
    ylabel('y');
    
    subplot(4,7,[15,16]);
    plot(tt,zt,'blue',tobs,zobs,'black',tb,zb,'green', ta,za,'red');   
    legend('zt','zobs','zb','za');
    title('z(t)');
    xlabel('t');
    ylabel('z');
    
    subplot(4,7,[3,4,5, 10,11,12, 17,18,19]);
    plot(xb,zb,'green', xa,za,'red',xt,zt,'blue'); 
    legend('zb(xb)','za(xa)','zt(xt)');
    title('Butterfly z(x)');
    xlabel('x');
    ylabel('z');
    
    subplot(4,7,[6,7]);
    Etb = (xt-xb).*(xt-xb) + (yt-yb).*(yt-yb) + (zt-zb).*(zt-zb);
    plot(tt,Etb,'blue'); 
    legend('||Ut - Ub||^2');
    title('Background vs True state');
    xlabel('t');
    ylabel('Etb');

    subplot(4,7,[13,14]);
    Etobs = (xt-xobs).*(xt-xobs) + (yt-yobs).*(yt-yobs) + (zt-zobs).*(zt-zobs);
    plot(tt,Etobs,'red'); 
    legend('||Ut - Uobs||^2');
    title('Observations vs True state');
    xlabel('t');
    ylabel('Etobs');
    
    subplot(4,7,[20,21]);
    Etf = (xt-xf).*(xt-xf) + (yt-yf).*(yt-yf) + (zt-zf).*(zt-zf);
    plot(tf,Etf,'green'); 
    legend('||Ut - Uf||^2');
    title('True State vs Forecast');
    xlabel('t');
    ylabel('Etf');
    
    subplot(4,7,[22,23,24,25,26,27,28]);
    Eta = (xt-xa).*(xt-xa) + (yt-ya).*(yt-ya) + (zt-za).*(zt-za);
    plot(tt,Eta,'magenta'); 
    legend('||Ut - Ua||^2');
    title('True State vs Analysis');
    xlabel('t');
    ylabel('Eta');
    
end