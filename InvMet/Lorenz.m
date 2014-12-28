
% Lorenz function describing the system 
% u' = f(u)
% with u = [x,y,z]';
function Y = Lorenz(X) 
    echo off;
    x = X(1);
    y = X(2);
    z = X(3);
    Y = [ 10*(y-x), ...
          28*x - y - x*z, ...
          x*y - 8/3*z]';   
end  

