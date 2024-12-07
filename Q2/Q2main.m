% Initialization and Parameters Setup

N_values = [10, 20, 40, 80, 100]; % Steps
minposi_eigen = zeros(length(N_values), 1);
eigen = cell(length(N_values), 1);

for N = 1:length(N_values)
    n = N_values(N);
    dx = 1/n; % pacing
    x = linspace(0, 1, n+1); % Grid

    % a(x) = 10, 1/4 < x < 3/4
    %        20, else
    a = zeros(n-1, 1);
    for i = 1:n-1
        if x(i+1) < 0.25 || x(i+1)>0.75
            a(i) = 20;
        else
            a(i) = 10;
        end
    end

    % A
    main_diag = 2 + dx^2 *a;
    off_diag = -1 * ones(n-2,1);
    A = dx^2 * (diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1)); % (N-1*1)

    [V, D] = eig(A);
    eigenvalues = diag(D);
    eigen(N) = {eigenvalues};

    neg_eigenvalues = eigenvalues(eigenvalues < 0);
    posi_eigenvalues = eigenvalues(eigenvalues > 0);
    smallest_positive = min(posi_eigenvalues);
    minposi_eigen(N) = smallest_positive;

    % plot
    figure;
    plot(x(2:end - 1), V(:, 1), 'r-', 'LineWidth', 1.5, 'DisplayName', 'First Eigenfunction');
    hold on;
    plot(x(2:end - 1), V(:, 2), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Second Eigenfunction');
    legend show;
    title(sprintf('Eigenfunctions for Grid Size N = %d', N));
    xlabel('x');
    ylabel('u(x)');
    grid on;

    saveas(gcf, sprintf('eigenfunctions_N%d.png', N));
    hold off;
end


%%

% Display Results
disp('Mesh Refinement Analysis Results:');
for N = 1:length(N_values)
    fprintf('Grid Size (N = %d): Smallest Positive Eigenvalue = %.6f\n', N_values(N), minposi_eigen(N));
end

% Mesh Refinement Analysis
figure;
plot(N_values, minposi_eigen, '-o', 'LineWidth', 2);
title('Convergence Analysis of Smallest Positive Eigenvalue');
xlabel('Number of Grid Points (N)');
ylabel('Smallest Positive Eigenvalue');
grid on;

% Eigenfunction Visualization for Largest Grid Size
N = N_values(end); % Use the largest grid size for detailed visualization
dx = 1 / N;
x = linspace(0, 1, N + 1);

% Reconstruct matrix A for largest N
a = zeros(N - 1, 1);
for i = 1:N - 1
    if x(i + 1) < 0.25 || x(i + 1) > 0.75
        a(i) = 20;
    else
        a(i) = 10;
    end
end
main_diag = 2 + dx^2 * a;
off_diag = -1 * ones(N - 2, 1);
A = dx^2*(diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1));

% Solve again for the largest N
[V, D] = eig(A);
eigenvalues = diag(D);

% Visualize eigenfunctions for the first two eigenvalues
figure;
plot(x(2:end - 1), V(:, 1), 'r-', 'LineWidth', 1.5, 'DisplayName', 'First Eigenfunction');
hold on;
plot(x(2:end - 1), V(:, 2), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Second Eigenfunction');
legend show;
title('Eigenfunctions for Largest Grid Size');
xlabel('x');
ylabel('u(x)');
grid on;
hold off;

set(gca, 'FontSize', 12);
legend('FontSize', 12);
