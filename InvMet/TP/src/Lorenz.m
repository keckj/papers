
% Lorenz function describing the system 
% u' = f(u)
% with u = [x,y,z]';
function Y = Lorenz(X) 
    x = X(1);
    y = X(2);
    z = X(3);
    Y = [ 10.0*(y-x), ...
          28.0*x - y - x*z, ...
          x*y - 8.0/3.0*z]';   
end  

