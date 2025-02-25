tempi=1:5;
esponente=1:5;
epsilon=1e-6;
quantita=1:5;
lambda=6.8e-4;

for numero=1:5
    n=2^(numero+4);
    esponente(numero)=n;
    pi=3.14159265358979323846264338327950288419716939937510582097494459230781640628; 
a=-6;
b=6;
h=(b-a)/n;
x_lambda=zeros(n,1);


K= @(x,s) (1+cos(pi*(s-x)/3)).*(abs(s-x)<=3);
g= @(s) (6-abs(s)).*(1+0.5*cos(pi*s/3))+9/(2*pi)*sin(pi*abs(s)/3);

A=zeros(n,n);
b_vect=zeros(n,1);
B=zeros(n,n);
for i=1:n
    
   b_vect(i,1)=(integral(g,a+(i-1)*h,a+i*h))/sqrt(h);

    for j=1:n
        A(i,j)=(integral2(K,a+(j-1)*h,a+j*h,a+(i-1)*h,a+i*h))/h;
        
    end
end
tic
[U,S,V]=svd(A);
s=svd(A);
for it=1:n
x_lambda=x_lambda+s(it)/(s(it)^2+lambda^2)*(((U(:,it))')*b_vect)*V(:,it);
end
tempi(numero)=toc;
quantita(numero)=sum(s/s(1)>epsilon);
s(1) %serve a vedere come si comporta il valore singolare più alto-
end
figure(1)
loglog(esponente,tempi,"r-");

figure(2)
plot(esponente,quantita,"k-");
figure(3)
plot(esponente,quantita./esponente,"b-");

