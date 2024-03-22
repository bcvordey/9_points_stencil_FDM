function f = rhsD9_fct(x,y,example,h,fourth_order)

% Return the values of the source term
% 
% Input
% x and y are vectors
% example: integer, to select the correct example

%
% Output:
% f: vector, the source term

% if the number of iunput arguments is only 2, set h to be 1
if nargin == 3
    example = 1;
end

switch example
    case 1
        f = 2*pi^2*sin(pi*x).*sin(pi*y) - fourth_order * (1/12)*h^2*4*pi^4*sin(pi*x)*sin(pi*y);
    otherwise
        f = -2*(x - 3*x.*y + (-1+y).*y.^2 + x.^2.*(-1+3*y)) - fourth_order * (1/12)*h^2*(24*y-8); 
end