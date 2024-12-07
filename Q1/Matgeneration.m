function A = Matgeneration(n, seed)
    % SymMatgeneration generates a random symmetric matrix of size n x n
    % Input parameters:
    %   n    - Size of matrix
    %   seed - Seed for random number generation (optional)

    if nargin > 1
        rng(seed); % If a seed is provided, set the seed for the random number generator
    end

    A = randn(n, n); % Generate a random matrix
    % symmetricMatrix = 0.5 * (A + A'); % Generating symmetric matrix
end