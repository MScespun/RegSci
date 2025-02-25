M1=zeros(64);
M2=zeros(128);
M3=zeros(256);
M4=zeros(512);
M5=zeros(1024);
%matrice 1
    n=2^(6)
   
    pi=3.14159265358979323846264338327950288419716939937510582097494459230781640628; 
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
M1=A;
%matrice 2
    n=2^(7)
   

h=(b-a)/n;




A=zeros(n,n);
b_vect=zeros(n,1);

for i=1:n
    
    b_vect(i,1)=(integral(g,a+(i-1)*h,a+i*h))/sqrt(h);

    for j=1:n
        A(i,j)=(integral2(K,a+(j-1)*h,a+j*h,a+(i-1)*h,a+i*h))/h;
        
    end
end
M2=A;
%matrice 3
    n=2^(8)
   

h=(b-a)/n;




A=zeros(n,n);
b_vect=zeros(n,1);

for i=1:n
    
    b_vect(i,1)=(integral(g,a+(i-1)*h,a+i*h))/sqrt(h);

    for j=1:n
        A(i,j)=(integral2(K,a+(j-1)*h,a+j*h,a+(i-1)*h,a+i*h))/h;
        
    end
end
M3=A;
%matrice 4
    n=2^(9)
   

h=(b-a)/n;




A=zeros(n,n);
b_vect=zeros(n,1);

for i=1:n
    
    b_vect(i,1)=(integral(g,a+(i-1)*h,a+i*h))/sqrt(h);

    for j=1:n
        A(i,j)=(integral2(K,a+(j-1)*h,a+j*h,a+(i-1)*h,a+i*h))/h;
        
    end
end
M4=A;
%matrice 5
    n=2^(10)
   

h=(b-a)/n;

A=zeros(n,n);
b_vect=zeros(n,1);
B=zeros(n,n);
for i=1:n
    
    b_vect(i,1)=(integral(g,a+(i-1)*h,a+i*h))/sqrt(h);

    for j=1:n
        A(i,j)=(integral2(K,a+(j-1)*h,a+j*h,a+(i-1)*h,a+i*h))/h;
        
    end
end
M5=A;
