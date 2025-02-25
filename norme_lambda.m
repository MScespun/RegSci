function [norma_x,norma_r]= norme_lambda(U,s,V,b_vect,n,lambda)
epsilon=10^(-4);
x_lambda_err=zeros(n,1);
r_lambda_err=zeros(n,1);
x_lambda=zeros(n,1);
r_lambda=zeros(n,1);
for it=1:n
    x_lambda_err=x_lambda_err+s(it)/(s(it)^2+lambda^2)*V(:,it);
    r_lambda_err=r_lambda_err+lambda^2/(s(it)^2+lambda^2)*U(:,it);
    x_lambda=x_lambda+s(it)/(s(it)^2+lambda^2)*(((U(:,it))')*b_vect)*V(:,it);
    r_lambda=r_lambda+lambda^2/(s(it)^2+lambda^2)*(((U(:,it))')*b_vect)*U(:,it);
end
x_lambda_err=x_lambda_err*epsilon;
r_lambda_err=r_lambda_err*epsilon;
norma_x=norm(x_lambda+x_lambda_err);
norma_r=norm(r_lambda+r_lambda_err);