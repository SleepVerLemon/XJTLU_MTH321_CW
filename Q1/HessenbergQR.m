function [Tk, k] = HessenbergQR(A, tolerance)
%
%

Tk = hess(A); % Tk is upper hessenberg here

for k=1:1000
    [Q, R] = qr(Tk);
    Tk = R * Q;

    lowertri = tril(Tk, -1);
    if norm(lowertri, 'fro') < tolerance
        break;
    end
end

end

