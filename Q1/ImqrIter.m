function [Tk, k] = ImqrIter(A, tolerance)
%
%

Tk = hess(A);
n = length(A(1));

for k = 1:1000
    mu = Tk(end,end);
    [Q, R] = qr(Tk - mu*eye(n));
    Tk = R + Q + mu * eye(n);

    % Apply deflation
    for i = n:-1:2
        if abs(Tk(i,i-1)) < tolerance
            Tk = Tk(1:i-1, 1:i-1);
            n = i-1;
            break;
        end
    end
    
    % inside tolerance ending 
    H = tril(Tk, -1);
    if norm(H, 'fro') < tolerance
        break;
    end
end

end

