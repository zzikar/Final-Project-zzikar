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
xvalues = linspace(0,2*pi,M);
yvalues = linspace(0,2*pi,N);

F = Functionzz(xvalues,yvalues);
% F=zeros(numXp2,numYp2);
   tic;% starting the timer
% Bounds for the Equation 

Er=10^-10 ; %Value of error for system convergence
U=zeros(M,N);
W=zeros(M,N);
%% Defining Boundary Conditions for "top" and "bottom"

% Bottom boundary condition

U(1,:) = ((xvalues - ax).^2 ) .* sin( pi *(xvalues - ax) / (2*(bx-ax)) ) ;
W(1,:)=U(1,:)

% Top boundary condition
U(N,:) = cos(pi*(xvalues-ax)).*cosh(bx-xvalues);
W(N,:)=U(N,:)
Error = zeros(N,M-2);

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
