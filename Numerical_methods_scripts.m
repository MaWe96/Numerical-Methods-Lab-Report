% Solution 1.2: absolute errors in aprroximating f'(1.2) for f(x)=x^x^x,
% using second and forth central order difference formula.
k=linspace(1,15,15);
h=2.^(-k);
CDF2=(((1.2+h).^((1.2+h).^(1.2+h)))-((1.2-h).^((1.2-h).^(1.2-h))))./(2*h);
CDF4=(-(1.2+2.*h).^((1.2+2.*h).^(1.2+2.*h))+8*(1.2+h).^((1.2+h).^(1.2+h))-8*(1.2-h).^((1.2-h).^(1.2-h))+(1.2-2.*h).^((1.2-2.*h).^(1.2-2.*h)))./(12*h);
Correctvalue=1.637932983825062173566;
error1=abs(Correctvalue-CDF2);
error2=abs(Correctvalue-CDF4);
loglog(h,error1,'b',h,h.^2,'b--',h,error2,'g',h,h.^4,'g--')
title('Errors for second and fourth order central difference approximations')
xlabel('Step size h')
ylabel('Error')

% Solution 2.1: Bisection, Newton-Raphson, Secant.
clc; clear
x=linspace(0,6,100);
R=3; L=8; V=100;
f=(R.^2*acos((R-x)./R)-(R-x).*sqrt(2.*R.*x-x.^2)).*L-V;
y=0*x;
plot(x,f,x,y)
grid on
title('Interval [0,6]')
xlabel('h-axis')
ylabel('f(h)')

% Algorithm - Bisection method.
function t = bisection(f,a,b,TOL,MaxIter)
tic
f(a);
f(b);
for n=1:MaxIter
    p=(a+b)/2;
    if f(p)==0
        break
    elseif f(a)*f(p)<0
        b=p;
    elseif f(b)*f(p)<0
        a=p;
    end
    if (b-a)/(2^n) < =TOL
        break
    end
end
n
p
toc

% Algorithm - Newton-Raphson method.
function r = newtonraphson(f,fprime,p,TOL,MaxNiter)
tic
for n=1:MaxNiter
    f(p);
    fprime(p);
    if fprime(p)==0
        disp('method failed')
    break
end
pn = p-(f(p)/fprime(p));
if abs((pn-p)/pn)<=TOL
    break
end
p=pn;
end
p
n
toc

% Algorithm - Secant method.
function s = secant(f,p0,p1,TOL,MaxIter)
tic
f(p0);
for n=1:MaxIter
    p2=p1-(((p1-p0)*f(p1))/(f(p1)-f(p0)));
    if abs((p2-p1)/p2)<=TOL
        break
    else
        p0=p1;
        p1=p2;
    end
end
p1
n
toc

% Algorithm - Newton interpolating polynomial.
function pvalue = newtoninterpolationfunc(year)
x=[1955 1960 1965 1970 1975 1980 1985 1990];
y=[2.7730 3.0349 3.3396 3.7004 4.0795 4.4580 4.8709 5.3272];
n=7; [m s]=size(x); M=[x' y'];
for k=3:s+1
    for i=k-1:s
        M(i,k)=(M(i,k-1)-M(i-1,k-1))/(M(i,1)-M(i-(k-2),1));
    end
end
terms=ones(1,n+n);
for k=1:n
    for i=1:n
        terms(i+k)=terms(i+k)*(year-x(k));
    end
end
for k=1:n+1
    terms(k)=terms(k)*M(k,k+1);
end
pvalue=0;
for k=1:n+1
    pvalue=pvalue+terms(k);
end

% Plotting script
clear; clc
for k=1955:1990
    polyvalues(k-1954)=newtoninterpolationfunc(k);
end
xvalues=1955:1990;
plot(xvalues,polyvalues)
hold on
plot ([1955 1960 1965 1970 1975 1980 1985 1990],[2.7730 3.0349 3.3396 3.7004 4.0795 4.4580
4.8709 5.3272],'o')
title('World population')
xlabel('Year')
ylabel('Population (Billion)')

% Algorithm - Least squares line (curve fitting).
clear; clc
x=[1955 1960 1965 1970 1975 1980 1985 1990];
y=[2.7730 3.0349 3.3396 3.7004 4.0795 4.4580 4.8709 5.3272];
[m n]=size(x);
sumx=sum(x);
sumy=sum(y);
sum2x=sum(x.^2);
sumxy=sum(x.*y);
M=[sum2x sumx; sumx n];
b=[sumxy; sumy];
p=linsolve(M,b);
p=p';
xx=linspace(1950,2000,100);
pvalues=polyval(p,xx);
plot(x,y,'o',xx,pvalues,'g')
title('World population')
xlabel('Year')
ylabel('Population (Billion)')
polyval(p,1976)

% Algorithm - Trapezoidal rule.
function integral = trapezoidalsum(f,a,b,N)
integral=0;
xk=linspace(a,b,N+1);
h=(b-a)/N;
for k=0:N
    if k==0||k==N
        integral=integral+f(xk(k+1));
    else
        integral=integral+2*f(xk(k+1));
    end
end
integral=(h/2)*integral;

% Simpsons function.
function integral = simpsons(f,a,b,N)
integral=0;
xk=linspace(a,b,N+1);
h=(b-a)/N;
t=0;
for k=0:N
    if k==0||k==N
        integral=integral+f(xk(k+1));
        t=t+1;
    elseif t==1
        integral=integral+4*f(xk(k+1));
        t=t-1;
    else
        integral=integral+2*f(xk(k+1));
        t=t+1;
    end
end
integral=(h/3)*integral;

% Analysis of results script.
I=0.746824132812427;
a=0;
b=1;
k=1:10;
h=2.^(-k);
f=@(x) exp(-x.^2);
m=1;
for n=1:10
    errort(m)=abs(I-trapezoidalsum(f,a,b,2^(n)));
    m=m+1;
end
m=1;
for n=1:10
    errors(m)=abs(I-simpsons(f,a,b,2^(n)));
    m=m+1;
end
loglog(h,errort,'b')
hold on
loglog(h,h.^2,'b--')
hold on
loglog(h,errors,'g')
hold on
loglog(h,h.^4,'g--')

% Matlab ODE45 solver
f=@(t,y) [((-0.4/1000)*y(2)*y(1));((0.4/1000)*y(2)*y(1)-0.04*y(2));0.04*y(2)];
y0=[997 3 0];
tspan=[0 150];
[t y]=ode45(f,tspan,y0);
plot(t,y(:,1),'--b',t,y(:,2),'--r',t,y(:,3),'--g')
title('Solving SIR model using ode45')
xlabel('Time')
ylabel('Number of people')

% Algorithm - 4th order Runge-Kutta function.
function r = rk4(f,a,b,y0,N)
q=(b-a)/N;
t=linspace(a,b,N+1)';
yn=y0;
for j=1:N
    k1=f(t(j),yn(j,:))';
    k2=f(t(j)+(q/2),yn(j,:)+(q/2)*k1)';
    k3=f(t(j)+(q/2),yn(j,:)+(q/2)*k2)';
    k4=f(t(j)+q,yn(j,:)+q*k3)';
    yn(j+1,:)=yn(j,:)+((1/6)*q*(k1+2*k2+2*k3+k4));
end
[t yn];
plot(t,yn(:,1),'b',t,yn(:,2),'r',t,yn(:,3),'g')