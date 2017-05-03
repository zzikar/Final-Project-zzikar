%Matlab Code to solve Poisson's equation with Gauss Seidel Method with the following conditions in the problem statement. 
% Jose Chavez  1161146 
clear all; clc; 

%% Given Conditions 
ax = 0;
ay = 0;
bx = 2*pi;
by = 2*pi; 
MI=input('Value of X Intenal Nodes='); % Number of points on the internal nodes for N and M%
NI=input('Value of Y Internal Nodes='); 
tic
M=NI+2; %Number of points including exterior boundary points for Ne and Me%
N=MI+2; 
% this generates the x and y values that will be used to calculate 
x = linspace(0,2*pi,N); 
y = linspace(0,2*pi,M); 
%% 

U = ones(N,M); %U initial guess %

% For loop solving for right hand side with F equation with i,j indices% 
for i=1:length(x); 
    for j=1:length(y); 
F(i,j) = cos ( (0.5*pi)* (2*((x(i)-ax) / (bx - ax))+1 )).*sin( pi*((y(j)-ay) / (by -ay))); 
    end 
end 

%% Boundary Conditions for "top" and "bottom" 

% Bottom boundary values 
U(1,:) = ((x - ax).^2 ) .* sin( (pi *(x- ax)) / (2*(bx-ax)) ) ; 

% Top boundary values 
U(MI+2,:)= cos (pi*(x-ax)).*cosh(bx-x); 

% place these known values in the solution grid 

%% Left and Right Boundary points 
%   Using the given neumann condition yields special cases of the Gauss-siedel iteration that can be used along entire "side" boundaries. 
%   F and U is computed in solution grid 
% Multipliers that are used in the iterations. 
dx = 2*pi/(MI+1); 
B = 1/dx.^2;
dy = 2*pi/(NI+1); 
C = 1/dy.^2;
den= -2*(B+C); 

% Normalize Multipliers%
B = B/den; 
C = C/den; 
F = F/den; 
den = 1; 
error=10; 
error_iterations=0;
% check for diagonal dominance of elements 
abs(den) >= abs(2*B+2*C)
while error>10^-10; 
    W=U; 

for i = 2:MI+1; 
     
    % Left boundary 
    U(i,1) = den*(  F(i,1) - (2*B)*U(i,2) - C*U(i-1,1) - C*U(i+1,1) );
    
    % Right Boundary 
    U(i,MI+2) = den*(  F(i,MI+2) - (2*B)*U(i,MI+1) - C*U(i-1,MI+2) - C*U(i+1,MI+2) ); 
end 

%% Gauss-Siedel iterating the general U equation%

for i = 2:MI+1; 
    for j = 2:NI+1; 
        U(i,j) = den*(  F(i,j) - C*U(i+1,j) - C*U(i-1,j)- B*U(i,j+1) - B*U(i,j-1) ); 
    end 
end 
error=abs(max(max(((W-U)./W)))); 
error_iterations=error_iterations+1;
end 
toc 
error_iterations
figure 
subplot(1,2,1),surf(U),xlabel('x axis'),ylabel('y axis'),title('F=cos(x)sin(y)');

subplot(1,2,2),contour(U),xlabel('x axis'),ylabel('y axis'),title('F=cos(x)sin(y)');
