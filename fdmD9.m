function fdmD9(l,example,fourth_order)

% FDl approxilation to the solution of
% -Laplacian(u(x,y)) = f(x,y) in (0,1)^2
% equipped with Dirichlet/Neumann boundary conditions
%
% Input:
% l, l: integer...
% example: integer, selection of the test problel
%
% Output
% a figure of the numerical solution
%
% also prints out the infinity norl of the error
%
%
% Possible ilprovelnts:
%
% add the pcg parameters tolerance and laxiter as input parameters

% General 9 points stencil 

%step size
h = 1/l;


%%%%liddle diagonals
e1 = ones((l-1)*(l-1),1);e1(1:l-1:end) = 0;
e = ones((l-1)*(l-1),1);
e_1 = ones((l-1)*(l-1),1);e_1(l-1:l-1:end) = 0; 



%%%%Lower Diagonals that vary
e_l_2 = ones((l-1)*(l-1),1);e_l_2(1:l-1:end) = 0;
e_l = ones((l-1)*(l-1),1);e_l(l-1:l-1:end) = 0;

%%%%% Upper Diagonals that very
el_2 = ones((l-1)*(l-1),1);el_2(l-1:l-1:end) = 0;
el = ones((l-1)*(l-1),1);el(1:l-1:end) = 0;



%%%%Creating A9 using spdiags
A = spdiags([-e_l -4*e -e_l_2 -4*e_1 20*e -4*e1 -el_2 -4*e -el],[-l,-(l-1),-(l-2),-1,0,1,(l-2),(l-1),l],(l-1)*(l-1),(l-1)*(l-1));
A = (1/(6*h^2))*A;
 
% full(A)



% grid points
x = linspace(0,1,l+1)';
y = linspace(0,1,l+1);


% evaluate the right hand side function on the interior grid points
f = rhsD9_fct(x(2:l),y(2:l),example,h,fourth_order);



% evaluate the exact solution function on the interior grid points
u = exactD9_fct(x,y,example);

%%%lodify RHS before solving  

%corners
f(1,1) = f(1,1) + (1/(6*h^2))*sum(u(1:3,1)) + (1/(6*h^2))*sum(u(1,2:3)); %5 bdry points
f(1,l-1) = f(1,l-1)+ (1/(6*h^2))*sum(u(1:3,l+1)) + (1/(6*h^2))*sum(u(1,(l-1):l)); %5 bdry points
f(l-1,1) = f(l-1,1) + (1/(6*h^2))*sum(u((l-1):l+1,1)) + (1/(6*h^2))*sum(u(l+1,2:3)); 
f(l-1,l-1) = f(l-1,l-1) + (1/(6*h^2))*sum(u(l+1,(l-1):l+1)) + (1/(6*h^2))*sum(u((l-1):(l),l+1));




%top bot left and right bdrys
for i = 2:l-2
    f(1,i) = f(1,i) + (1/(6*h^2))*sum(u(1,i:i+2));
    f(l-1,i) = f(l-1,i) + (1/(6*h^2))*sum(u(l+1,i:i+2));
    f(i,1) = f(i,1) + (1/(6*h^2))*sum(u(i:i+2,1));
    f(i,l-1) = f(i,l-1) + (1/(6*h^2))*sum(u(i:i+2,l+1));
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = reshape(f,(l-1)*(l-1),1);

tic
u_in = pcg(A,f,1e-10,10000); % avoid using a direct solver whenever it is possible
toc


%add the BC to the u and Reshpe
u_in = reshape(u_in,l-1,l-1);
u(2:l,2:l) = u_in; 
u = reshape(u,l+1,l+1);

exact = exactD9_fct(x,y,example);
 

 norm = max(max(abs(u - exact)))

 figure('Name','FDM D9 Approxilation')
 surf(x,y,u')
 xlabel("x")
 ylabel("y")

 
 exact = reshape(exact',l+1,l+1);
 figure('Name','Exact Solution')
 surf(x,y,exact)
 xlabel("x")
 ylabel("y")
end
