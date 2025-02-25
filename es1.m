%Per prima cosa genero la matrice A e il vettore dei termini noti b
pi=3.14159265358979323846264338327950288419716939937510582097494459230781640628; 
n=64;
a=-6;
b=6;
h=(b-a)/n;

K= @(x,s) (1+cos(pi*(s-x)/3)).*(abs(s-x)<=3);
g= @(s) (6-abs(s)).*(1+0.5*cos(pi*s/3))+9/(2*pi)*sin(pi*abs(s)/3);

A=zeros(n,n);
b_vect=zeros(n,1);

for i=1:n
    
    b_vect(i,1)=(integral(g,a+(i-1)*h,a+i*h))/sqrt(h);

    for j=1:n
        A(i,j)=(integral2(K,a+(j-1)*h,a+j*h,a+(i-1)*h,a+i*h))/h;
       
    end
end
%faccio SVD di A
s=svd(A);
[U,S,V]=svd(A);

%perturbo b con un vettore le cui componenti sono generate con una
%distribuzione Gaussiana di media 0 e std^2=1e-8
mu=0;
sigma=1e-8;
rng("default")
R=mvnrnd(mu,sigma,n);
b_vect=b_vect+R;

%discretizzo per creare il grafico
x=1:10000;
y=1:10000;

for k=1:10000
    lambda=10^(-6+k*6.5*10^(-4));
    [x(k),y(k)]=norme_lambda(U,s,V,b_vect,n,lambda);

end
%calcolo i valori U^T*b_vect e li metto in z per tracciare il grafico
z=zeros(n,1);
for it=1:n
    z(it)=abs((U(:,it))'*b_vect);
 
   
end
%dal grafico sembra che il valore di lambda ottimale è circa 1e-3, lo
%confermo con un cerchietto nero che effettivamente si trova nei pressi
%dello spigolo della L
lambda=6.8*10^(-4);
[lambda2,lambda1]=norme_lambda(U,s,V,b_vect,n,lambda);
    
figure(1)
semilogy(1:n,s,"r.",1:n,z,"bo",1:n,z./s,"kx");

figure(2)
semilogx(y,x,"g.",lambda1,lambda2,"ko");




