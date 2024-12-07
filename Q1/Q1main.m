% Parameters Initialization
matrix_size = [20, 100, 500, 1000, 1500];
num_size = length(matrix_size);
tolerance = 1e-6;

methods = {'Power Iteration', 'Preliminary QR', 'Hessenberg QR', 'Implicit QR'};
num_methods = length(methods);
times = zeros(num_size, num_methods);
iterations = zeros(num_size, num_methods);
errors = zeros(num_size, num_methods);

% Set the random number seed
seed = 2142482;
rng(seed);

% Generating the matrix
A20 = Matgeneration(matrix_size(1), seed);
A100 = Matgeneration(matrix_size(2), seed);
A500 = Matgeneration(matrix_size(3), seed);
A1000 = Matgeneration(matrix_size(4), seed);
A1500 = Matgeneration(matrix_size(5), seed);
A = {A20, A100, A500, A1000, A1500};

%% 

% Power Iteration

for s = 1:num_size
q0 = randn(matrix_size(s),1);
tic;
[lambda_power, k, q] = powerIter(cell2mat(A(s)), q0, tolerance);
times(s, 1) = toc;
iterations(s, 1) = k;
errors(s, 1) = norm(cell2mat(A(s))*q - lambda_power*q);
end

disp ('--- Power Iteration ---') ;

%% 

% Preliminary QR Iteration

for s = 1:num_size
Q0 = randn(matrix_size(s));
tic;
[Tk, k] = preqrIter(cell2mat(A(s)), Q0, tolerance);
times(s, 2) = toc;
iterations(s, 2) = k;
errors(s, 2) = norm(tril(Tk, -1), 'fro');
end

disp ('--- Preliminary QR Iteration ---') ;

%% 

% Hessenberg QR Iteration

for s = 1:num_size
tic;
[Tk, k] = HessenbergQR(cell2mat(A(s)), tolerance);
times(s, 3) = toc;
iterations(s, 3) = k;
errors(s, 3) = norm(tril(Tk, -1), 'fro');
end

disp ('--- Hessenberg QR Iteration ---') ;

%% 
% 
% % Double-Shift QR Iteration
% 
% for s = 1:num_size
% tic;
% [Tk, k] = dshiftqr(cell2mat(A(s)), tolerance);
% times(s, 4) = toc;
% iterations(s, 4) = k;
% errors(s, 4) = norm(tril(Tk, -1), 'fro');
% end
% 
% disp ('--- Double-Shift QR Iteration ---') ;
% 
%% 

% Implicit QR Iteration

for s = 1:num_size
tic;
[Tk, k] = ImqrIter(cell2mat(A(s)), tolerance);
times(s, 4) = toc;
iterations(s, 4) = k;
errors(s, 4) = norm(tril(Tk, -1), 'fro');
end

disp ('--- Implicit QR Iteration ---') ;

%%

% Plot

% Plot results for execution time comparison
figure;
for m = 1:num_methods
    plot(matrix_size, times(:, m), '-o', 'DisplayName', methods{m});
    hold on;
end
title('Execution Time Comparison');
xlabel('Matrix Size (n x n)');
ylabel('Time (s)');
legend;
grid on;

% Plot results for iteration count comparison
figure;
for m = 1:num_methods
    plot(matrix_size, iterations(:, m), '-o', 'DisplayName', methods{m});
    hold on;
end
title('Iteration Count Comparison');
xlabel('Matrix Size (n x n)');
ylabel('Number of Iterations');
legend;
grid on;

% Plot results for convergence error comparison
figure;
for m = 1:num_methods
    plot(matrix_size, errors(:, m), '-o', 'DisplayName', methods{m});
    hold on;
end
title('Convergence Error Comparison');
xlabel('Matrix Size (n x n)');
ylabel('Error');
legend;
grid on;

% Display metrics
for m = 1:num_methods
    disp(['--- ', methods{m}, ' ---']);
    disp(table(matrix_size', times(:, m), iterations(:, m), errors(:, m), ...
        'VariableNames', {'Matrix Size', 'Time (s)', 'Iterations', 'Error'}));
end
