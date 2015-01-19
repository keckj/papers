
function main() 
    N = 64;
    xmin = -3.14/2;
    xmax = 3.14/2;
    
    [x,f] = sample(N, xmin, xmax, @sin);
    %plot(x,f,'x');
    DB(0,0,0,1,10);
end


function [x,f] = sample(n, xmin, xmax, func)   
    x = xmin + rand(1,n)*(xmax - xmin);
    x = sort(x);
    f = arrayfun(func, x);
end

function y = phi(x) 
    if(abs(x) > 1)
        y = 0;
    else 
        y = 1;
    end
end

function phi_jk = DB(j,k,xmin,xmax,pts)
    dx = (xmax - xmin)/(pts-1);
    for x = xmin:dx:xmax
        phi(2^-j * x + k);
    end
    
end

