clear;clc;



A=[1,2,-2; 1,1,1; 2,2,1];
b=[7;2;5];
x0=zeros(3,1);
tol = 1e-5;

[x_j, iterations_j, error_j, maxeig_j] = jacobi2(A, b, x0, tol);
[x_gs, iterations_gs, error_gs, maxeig_gs] = gs2(A, b, x0, tol);




function [x, i, err, spectralj] = jacobi2(A, b, x0, tol)
[n, m]=size(A);
U = zeros(n);
L = zeros(n);
D = zeros(n);

i = 1;
    if m ~=n 
        fprintf('matrix not square')
        return
    end
    for i=1:n
        for j=1:n
            if j > i
                U(i,j) = A(i,j);
            elseif j < i
                L(i,j) = A(i,j);
            else
                D(i,j) = A(i,j);
            end 
        end 
    end
    x = D\(b - U*x0 - L*x0);
    x-x0;
    err = norm(x-x0, inf);
   
while err > tol
    x0 = x;
    x = D\(b - U*x0 - L*x0);
    i = i+1;
    err = norm(x-x0, inf);
end 

 T = inv(D)* (U+L)
 e = eig(T)
 spectralj = norm(e, inf)


end
    
                

function [x, j, err, spectralj] = gs2(A, b, x0, tol)
[n, m]=size(A);
U = zeros(n);
L = zeros(n);
D = zeros(n);
i = 1;

    if m ~=n 
        fprintf('matrix not square')
        return
    end
    for i=1:n
        for j=1:n
            if j > i
                U(i,j) = A(i,j);
            elseif j < i
                L(i,j) = A(i,j);
            else
                D(i,j) = A(i,j);
            end 
        end 
    end
    x = inv(D)\(b - U*x0 - L*x0);
    j=1;
    err (j)= norm(x-x0, inf);
    A;
    for i = 1:10
    %while err(j) > tol
    x0 = x;
    x = (D+L)\(b - U*x0);
    j=j+1;
    err (j) = norm(x-x0, inf);
    end 
    
    T = inv(D-L)* (U)
    e = eig(T)
    spectralj = norm(e, inf)
% figure
% semilogy(1:j, err, 'bo-');
% title('Error vs. Iterations: Gauss-Seidel')
% xlabel('iterations')
% ylabel ('error')

end

