function func = GaussSiedelFunction(A,b,init,tol,maxit)
%% Gauss Seidel Method
%% Solution of x in Ax=b using Gauss Seidel Method
if nargin<3,error('at least 3 input arguments required'),end
if nargin<5 | isempty(maxit),maxit=50;end
if nargin<4 | isempty(tol),tol=1e-3;end
[m,n]=size(A)
if m~=n, error('Matrix A must be square.');end
C=A

%% Algorithm: Gauss Seidel Method
x=init
for i=1:n
  C(i,i)=0
  x(i)=0
end
for i=1:n
  C(i,1:n)=C(i,1:n)/A(i,i)
end
for i=1:n
  d(i)=b(i)/A(i,i) 
end

iter=0
while (1)
    x_old=x;
    
    for i=1:n
      sum=0
      for j=1:n
          sum=sum+C(i,j)*x(j);
      endfor
      x(i)=d(i)-sum
      if x(i)~=0
        eps(i)=abs((x(i)-x_old(i))/x(i))*100
      end
    endfor
    
    iter=iter+1
    if max(eps)<=tol | iter>=maxit, break, end
end

fprintf('The solution of the system is:') 
for i=1:length(x)
  fprintf('\n%d%f\n',x(i))
endfor
fprintf('\nIt was obtained in %d iterations.',iter)
fprintf('\n')
