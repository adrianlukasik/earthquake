

%earthquake model
%Mx''=Kx'+F 
%where M i diagona matrix; K 3diag matrix symmetric
%k=[k_0,...,k_{n-1}] - coefficient 
%   such that on i-th floor acts -k(i-1)(x(i)-x(i-1)) + k(i)(x(i+1)-x(i))
%F=[F_1,0,..0] with Gsin(gamma t) gamma frequency of an earthquake, G magnitude
n=10;
m=10000*ones(n,1);
m(1) = 30000;
k=5000*ones(n,1);
k(1) = 10000;
gamma=3;
G=5000;
x0=v0=zeros(n,1);


function y=pw(x,t)
  %vector field for lsode() for the system x''=Ax+F
  %with F =[GG*sin(gam*t);zeros(2:end,1)]
  global AA gam GG
  n=length(AA);
  y=zeros(size(x));
  y(1:n)=x(n+1:end);
  y(n+1:end)=AA*x(1:n);
  y(n+1)=y(n+1)+GG*sin(gam*t);
end



function [X,tt,f,p]=earthquake(k,m,gamma,G,T=15,x0=zeros(size(m)),v0=zeros(size(x0)))
%function solving the system Mx''=Ax''+F(t) x(0)=x0 x'(0)=v0
%modelling an earthquake of a building: M diagonal matrix with masses of the storeys
%with the diagonal vector m; K stiffness matrix - 3diag symmetric 
%matrix having [k(2),..,k(n)] on sub/super diagonals and -k(i)-k(i+1) on main diagonal (k(n+1)=0)
%F(t)=[G*sin(gamma*t);zeros(2:end,1)] - force applied to the 1st floor - G magnitude
% gamma - the frequency of the earthquake we have 2*pi/gamma \in [2,3]  
%INPUT
%k, m vectors wit hmasses and forces coefficient
%G,gamma - magnitude and frequncy of athe earthquake
%T max time for whch we model EQ
%x0,v0 - (default zero vectors) initial position and velocities of the floors
%OUTPUT
%X discrete solution for 100 points on [0,T]
%tt - vector of discrete times 
%f,p - vectors with frequancies and periods  of the building
global AA gam GG
gam = gamma;
GG=G/m(1);  
n=length(m);
k(n+1)=0; %artifial extension of the vector k by zero
M=diag(m);
K=-diag(k(1:n)+k(2:end))+diag(k(2:n),1)+diag(k(2:n),-1);
AA=M\K;
d=eig(AA);
f=sqrt(-d);
p=(2*pi)./f;

x=[x0;v0];

tt=linspace(0,T);  
X=lsode(@pw,x,tt);

end  

 [X,t,f,p]=earthquake(k,m,gamma,G);
 
 f
 p
 
 fok=!sum((p>2)&(p<3));
 if(fok)
  disp("Building periods not in [2,3] - OK");
 end
 
 plot(t,X(:,1),";1st floor;",t,X(:,length(m)),";last floor;");
 pause(2);
 
 aplitudes = zeros(1, length(m));
 for i = 1 : length(m)
   aplitudes(i) = max(abs(X(:, i)));
 endfor
 aplitudes
 
 F = zeros(1, length(t));
 for i = 1 : length(t)
   F(i) = max(abs(X(i, 1 : length(m))));
 endfor
 plot(t, F, ";F(t);");
 
 
