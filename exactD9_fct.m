function u = exactD9_fct(x,y,example)

% returns a vector

switch example
    case 1
        u = sin(pi*x).*sin(pi*y);
    otherwise
        u = x.*(1-x).*y.^2.*(1-y);
end