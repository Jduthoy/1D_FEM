function [ e ] = FEM1D( h1,g,ne,c,f )
%%%%%%% Julius Duthoy%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Math 520     %%%%%%%%%%%%%%%%%%%%%%%%%%
%Test2: FEM1D(6,50,4,17/4,@(x) 1-x)
%Future code let -u'(x1)=h1 and u(x2)=g with x1 and x2 inputs
%Problem: u"+u'+c*u=f
%        -u'(0)=h1 and u(1)=g

%%%%%%%%%%List of Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ln = #Local Nodes
ln=2;
% gn = #Global Nodes (1 <= A <= gn)
gn=ne+1;
% ne = #elements- input into the function

% n_eq = #Nonzero Global Basis Fncts
n_eq=0; %Determined from ID Matrix

%%%%%%%%%% Set up Global Nodes & Constant Interval h %%%%%%%%%%%%%%%%%%
x=[0:4/ne:4];
h=x(2:ne+1)-x(1:ne);

%%%%%%%%%% ID(A)=P Matrix%%%%%%%%%%%%%%%%%%%%%%%%%% 
for A=1:gn
    if A==1
        ID(A)=0;
    end
    if A~=1
        n_eq=n_eq+1;
        ID(A)=n_eq;
    end
end
ID

%%%%%% IEN(a,e)=A Matrix %%%%%%%%%%%%%%%%%%%%%%%%
for a=1:ln
    for e=1:ne
        IEN(a,e)=e+(a-1);
    end
end
IEN

%%%%%%%%% LM(a,e)=ID(IEN) Matrix %%%%%%%%%%%%%%%%%%%%%%%%
LM = ID(IEN)


%%%%%%% Initialize Basis fncts%%%%%%%%%%%%%%%%%%%%%%%%% 
N1=@(z) 1-z;
N2=@(z) z;

% Create K and F 
K=zeros(n_eq,n_eq); 
F=zeros(n_eq,1);
% Constructing stiffness matrix
%K=zeros(neq,neq);
%for e=1:ne
%    for a=1:nln
%        if LM(a,e)~=0
%            if a==1
%            F(LM(a,e))=F(LM(a,e))+quad(@(z)f(h*z+x(e)).*n1,0,1);
%            end
%        for b=1:nln
%            if LM(b,e)~=0
%                K(LM(a,e),LM(b,e))=K(LM(a,e),LM(b,e))+k(a,b);
%            end
%        end
%        end
%    end
%end
%K
%F 
%%%%%%%%% Taken from Class %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e=1:ne
       k(1,1)=-(1/h(e))*(-1)*(-1)+h(e).*quad(@(z) N1(z).*N1(z),0,1);%mini k11 
       k(1,2)=-(1/h(e))*(-1)*(1)+h(e).*quad(@(z) N1(z).*N2(z),0,1);%mini k12 etc
       k(2,1)=-(1/h(e))*(1)*(-1)+h(e).*quad(@(z) N2(z).*N1(z),0,1);
       k(2,2)=-(1/h(e))*(1)*(1)+h(e).*quad(@(z) N2(z).*N2(z),0,1);
%%%%%%%% Constructing Stiffness Matrix and F %%%%%%%%%%%%%%%%%%%%%%%
    if LM(1,e)~=0 && LM(2,e)~=0
        K(LM(1,e):LM(2,e),LM(1,e):LM(2,e))=K(LM(1,e):LM(2,e),LM(1,e):LM(2,e))+k;
        F(LM(1,e))=F(LM(1,e))+h(e).*quad(@(z) f(x(e)+h(e).*z).*(1-z),0,1);
        F(LM(2,e))=F(LM(2,e))+h(e).*quad(@(z) f(x(e)+h(e).*z).*z,0,1);
    elseif LM(1,e)~=0
        K(LM(1,e),LM(1,e))=K(LM(1,e),LM(1,e))+ k(1,1);
        F(LM(1,e))=F(LM(1,e))+h(e).*quad(@(z) f(x(e)+h(e).*z).*(1-z),0,1);
    elseif LM(2,e)~=0
        K(LM(2,e),LM(2,e))=K(LM(2,e),LM(2,e))+k(2,2);
        F(LM(2,e))=F(LM(2,e))+h(e).*quad(@(z) f(x(e)+h(e).*z).*z,0,1);
    end
    K
%{    
%%%%%%%% Adding in BC's and IC's %%%%%%%%%%%%%%%%%%%%%%%%%%
%   if x1==0 
    if e==ne
        F(LM(1,e))=F(LM(1,e))-k(1,2)*g; % since u(1)=g is at last element 
    end
%   if x2==1
    if e==1
        F(LM(1,e))=F(LM(1,e))-h1;% u(0)=h1 is at first element
    end
%%%% Future: Let user state the BC's and check for each condition
    
end
 
%%%%%%%%% Solve for d in U^h=d1N1+...%%%%%%%%%%%%%%%%%%%%
d=K\F; % May use Cholesky since K is SPD...

%%%% Convert from local to Global using ID%%%%%%%%%%%%%%%%%%%%
for A=1:gn
    if ID(A)~=0
        u(A)=d(ID(A));
    else
        u(A)=g;
    end
end

 

%%%%%%%%%%Comparison to Exact Soln%%%%%%%%%%%%%%%%%%%%%%%%%
sol=@(x) x.^2-2*sin(x)+2*cos(x)*(tan(1)+sec(1))-2;
UU=arrayfun(sol,x); %UU=Actual u calculated by hand
%arrayfun(fun,inputs) uses x array as inputs to sol to form the actual

%%%%%%%%%% Plots the Actual U and my FEM approximation
plot(x, u, 'm*', x, UU, 'k')
legend('Approx u', 'Actual u')
%title('u"+u=x^2, u(1)=g, -u_x(0)=h' )
str=sprintf('u"+u=x^2, u(1)=%f, -u_x(0)=%f',g,h1);
% Couldn't get f to change to whatever the user inputs
title(str);
xlabel('x')
ylabel('u')

format long
e=max(abs(UU-u)); %Calculates the error

%}
end

