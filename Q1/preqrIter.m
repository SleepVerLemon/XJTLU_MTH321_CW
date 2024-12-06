function [Tk,k] = preqrIter(Q0, A, tolerance)
%
%

Tk = Q0' * A * Q0;
for k=1:1000
    [Q, R] = qr(Tk);
    Tk = R * Q;

    lowertri = tril(Tk, -1);
    if norm(lowertri, 'fro') < tolerance
        break;
    end
end

end

