clear all;
clc;
alpha=0.5;
L=1;
T=1;
M=3;
N=M;
h=L/M;
tau=T/N;
x=0:h:L;
t=0:tau:T;

%f=@(i,j) 4/gamma(2.5)*(j*tau)^(2-alpha)*sin(pi*i*h)*cos(pi*i*h)+4*(j*tau)^2*sin(2*pi*i*h)*pi^2;
f=@(i,j) 3.009011112*(j*tau)^(2-alpha)*sin(pi*i*h)*cos(pi*i*h)+4*(j*tau)^2*sin(2*pi*i*h)*pi^2;
for i=0:M
    for j=0:N
     Exact(i+1,j+1)=(j*tau)^2*sin(2*pi*i*h);
    end
end

rau=@(i) -sin(2*pi*i*h);
u=zeros(M+1,N+1);

A=zeros(N+1,N+1);

for i=2:N+1
        A(i,i)=-1/h^2;
end

B=zeros(N+1,N+1);

for i=1:N+1
    for j=1:N+1
        if i==1
          B(i,1)=1;
          B(i,N+1)=-1;
        elseif i==j
          B(i,i)=(w(0,alpha)/tau^alpha+2/h^2);
        elseif  j<i && j>1 &&  i>2
          B(i,j)=(w(i-j,alpha)/tau^alpha);
        elseif j==1
            B(i,j)=-sum_w(i-2,alpha)/tau^alpha;
        end
        
    end
end
Alpp{1}=zeros(N+1,N+1);
Betaa{1}=zeros(N+1,1);

for n=1:M-1
    phi(1)=rau(n);
    for j=1:N
        phi(j+1,1)=f(n,j);
    end
    Alpp{n+1}=-(B+A*Alpp{n})\A;
    Betaa{n+1}=(B+A*Alpp{n})\(phi-A*Betaa{n});
end

for n=M-1:-1:1
u(n+1,:)=Alpp{n+1}*u(n+2,:)'+Betaa{n+1};
end

error=max(max(abs(u(:,:)-Exact(:,:))))

function total = sum_w(k,alpha)
    total = 0;
    for m = 0:k
        total = total + w(m, alpha);
    end
end

function value=w(m,alpha)
if m==0
    value=1;
else
    value=(1-(alpha+1)/m)*w(m-1,alpha);
end
end

