function [lambda1,k,q] = powerIter(A, q, epsilon)
% 
% 

% n = 0;
lambda0 = 0;
for k = 1:1000
    z = A * q;
    q = z / norm(z);
    lambda1 = q' * A * q;
    if abs(lambda1 - lambda0) < epsilon
        break
    end
    lambda0 = lambda1;
end

end

