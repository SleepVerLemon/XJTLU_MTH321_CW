function [Tk,k] = dshiftqr(A, tolerance)
%
%

Tk = hess(A); % Tk is upper hessenberg here

for k=1:1000
    n = length(Tk(1));
    [mu1, mu2] = eig(Tk(n-1:n,n-1:n));
    H = Tk - (mu1 + mu2)*eye(n) + mu1*mu2*eye(n);
    [Q, R] = qr(H);
    Tk = R * Q;

    lowertri = tril(Tk, -1);
    if norm(lowertri, 'fro') < tolerance
        break;
    end
end
end

