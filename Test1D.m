function [ e ] = Test1D( H,g,ne,ln,f )
%%%%%%% Julius Duthoy%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Math 520     %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Linear Basis Elements 1D %%%%%%%%%%%%%%
%Test1D(0,0,4,2,@(x) 1-x)
% u"+u'-(15/4)u= 1-x
%BC's u(0)=5=g u'(4)=-12=H
% Test1D(-12,5,4,2,@(x) 1-x)
%%%%%%%%%%List of Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%ln  =# Local Nodes
%gn  =# Global Nodes 
%ne  =# elements- input into the function
%n_eq=# Nonzero Global Basis Fncts

gn=ne+1;

%%%%%%%%%% Set up Global Nodes & Constant Interval h %%%%%%%%%%%%%%%%%%
x = [0:4/ne:4];
h=x(2:ne+1)-x(1:ne);

%%%%%%%%%% ID(A)=P Matrix%%%%%%%%%%%%%%%%%%%%%%%%%% 

n_eq=0;%Determined from ID Matrix below

for A=1:gn
    if A==1
        ID(A)=0; %only works for u(1)=g not u(0)=g
    else
        n_eq=n_eq+1;
        ID(A)=n_eq;
    end
end
%ID

%%%%%% IEN(a,e)=A Matrix %%%%%%%%%%%%%%%%%%%%%%%%
for e=1:ne 
    for a=1:ln
        IEN(a,e)=e+(a-1);
    end
end
%IEN

%%%%%%%%% LM(a,e)=ID(IEN) Matrix %%%%%%%%%%%%%%%%%%%%%%%%
LM=ID(IEN);


%%%%%%% Initialize Local Basis fncts %%%%%%%%%%%%%%%%%%%%%%%%% 

N1=@(z) 1-z;
dN1=@(z) -1+z-z;

N2=@(z) z;
dN2=@(z) 1+z-z;


N={'1-z','z'};
dN={'-1+z-z','1+z-z'};


%%%%%%%%%Stiffness Matrix and F Vector%%%%%%%%%%%%%%%%%%%%%%%%%
K=zeros(n_eq,n_eq); 
F=zeros(n_eq,1);

%fnct=inline(N{1,2});
%F(LM(2,1))=F(LM(2,1))+h(1).*quad(@(z) f(x(LM(2,1))+h(1).*z).*fnct(z), 0, 1)


for e=1:ne
 
       %Set Up Big K
       for a=1:ln
            if LM(a,e)~=0
                if LM(2,e)~=0
                    fnct=inline(N{1,a});
                    F(LM(a,e))=F(LM(a,e))+h(e).*quad(@(z) f(x(LM(2,e))+h(e).*z).*fnct(z), 0, 1);
           
                end
                for b=1:ln
                        fnct_a=inline(N{1,a});
                        fnct_b=inline(N{1,b});
                        der_a=inline(dN{1,a});
                        der_b=inline(dN{1,b});
                        %k(a,b)=-(1/h(e)).*quad(@(z) der_a(z).*der_b(z), 0, 1)-(15/4)*h(e).*quad(@(z) fnct_a(z).*fnct_b(z), 0, 1)+quad(@(z) der_a(z).*fnct_b(z),0,1);
                        k(b,a)=-(1/h(e)).*quad(@(z) der_b(z).*der_a(z), 0, 1)-(15/4)*h(e).*quad(@(z) fnct_b(z).*fnct_a(z), 0, 1)+quad(@(z) der_b(z).*fnct_a(z),0,1);
                    if LM(b,e)~=0
                        K(LM(a,e),LM(b,e))=K(LM(a,e),LM(b,e))+k(b,a);
                    end
                end
            end
       end

% Placing B.C.s in proper entries
    if e==1
        F(LM(2,e)) = F(LM(2,e)) - k(1,2)*g;
    end
  

    if e==ne
      F(LM(2,e)) = F(LM(2,e)) - H; %*N2(1)
    end
       
end


% Solve for d
d = K\F;

% Get u from d
for A=1:gn
    if ID(A)~=0
        u(A)=d(ID(A));
    else
        u(A)=g;
    end
end
      
%%%%%%%%%%Comparison to Exact Soln%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% zero BC's Solution%%%%%%%%%%%%%%%%%%%%%%%
%sol= @(z)-(1/(255*(5+3*exp(16))))*4*exp((-5*z)/2)*(exp((5*z)/2)*...
%    (55-75*z)-55*exp(4*z)+30*exp(4*z+10)+exp(((5*z)/2)+16)*...
%    (33-45*z)-33*exp(16)-30*exp(10));
%%%%%%%% Actual BC's of the problem%%%%%%%%%%%%%%%%%%%%%%%%
sol=@(z)1/(225*(5+3*exp(16)))*exp(-5*z/2)*(20*exp(5*z/2)*(15*z-11) + 12*exp(5*z/2 + 16)*(15*z-11)+ 5845*exp(4*z) -...
  5520*exp(4*z+10) + 3507*exp(16) + 5520*exp(10));
UU=arrayfun(sol, x);%UU=Actual u calculated by hand
%arrayfun(function,inputs) uses x array as inputs to sol to form the actual

%%%%%%%%%% Plots the Actual U and my FEM approximation
plot(x, u, 'm*', x, UU, 'k')
legend('Approx u', 'Actual u')
str=sprintf('u"+u`-(15/4)u=1-x, u(1)=%f, -u_x(0)=%f',g,H);
% Couldn't get f to change to whatever the user inputs
title(str);
xlabel('x')
ylabel('u')

%format long
e=max(abs(UU-u)); %Calculates the error

%}

end

