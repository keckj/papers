
% Main function
function question6
    
    % initial condition & configuration
    u_0  = [-4.62,-6.61,17.94]';
    tmax = 30.0;   %s
    dt   = 0.001;  %s
    
    % solving the system with Runge-Kutta of 4th order
    u_t = RK4(u_0, @Lorenz, dt, tmax);
    
    % plot the solution
    custom_plot(u_t);
end

% Plot solution 
function custom_plot(u_t) 
    
    t = u_t(1,:);
    x = u_t(2,:);
    y = u_t(3,:);
    z = u_t(4,:);

    subplot(3,5,[1,2]);
    plot(t,x,'blue');
    title('x(t)');
    xlabel('t');
    ylabel('x');
    
    subplot(3,5,[6,7]);
    plot(t,y,'green');
    title('y(t)');
    xlabel('t');
    ylabel('y');
    
    subplot(3,5,[11,12]);
    plot(t,z,'red');
    title('z(t)');
    xlabel('t');
    ylabel('z');
    
    subplot(3,5,[3,4,5, 8,9,10, 13,14,15]);
    plot(x,z,'cyan');
    title('Butterfly z(x)');
    xlabel('x');
    ylabel('z');
    
end