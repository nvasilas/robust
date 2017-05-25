function [p] = pp_data()
    A = [0 1 0; 0 0 1; 1 2 -1];
    B = [0; 0; 1];
    Q = eye(length(A));
    Q_0 = inv(Q);
    R = 1;

    d_t1 = [0; 0; 1];
    e_t1 = [0; 1; 0];
    d_t2 = [0; 0; 1];
    e_t2 = [0; 0; 1];

    n = length(A);
    k = 2; % # of uncertainties

    D = [d_t1, d_t2];
    E = [e_t1, e_t2]';

    % relative stability degree
    a = 1;

    p = struct('A', A, 'B', B, 'Q', Q, ...
        'Q_0', Q_0, 'R', R, 'n', n, 'k', k, ...
        'D', D, 'E', E, 'a', a, 'sort_array', @sort_array);
end

function [V] = sort_array(V_, D, n)
    % this function sorts the array
    % of the generalized eigenvectors
    % according to the form required by
    % the step 2 of the algorithm
    if n < 1
        fprintf('Error n must >= 1\n');
        return
    end
    [~, idx] = sort(diag(D));
    V = V_(:, idx(1));
    for i = 2 : 2*n
        V = horzcat(V, V_(:, idx(i)));
    end
end
