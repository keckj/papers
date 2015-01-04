
function main() 
    N = 1024;c
    xmin = -3.14/2;
    xmax = 3.14/2;
    
    [x,f] = sample(N, xmin, xmax, @sin);
    plot(x,f);
end


function [x,f] = sample(n, xmin, xmax, func)   
    x = xmin + rand(1,n)*(xmax - xmin);
    x = sort(x);
    f = arrayfun(func, x);
end

