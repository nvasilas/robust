function pp_algorithm()
    % get data
    [p] = pp_data();

    % check if uncertainties exists
    if p.k == 0
        fprintf('Solving the nominal system for a = %f\n', p.a);
        fprintf('No uncertainties are present, k = 0\n');
        [Gnom, Pnom, error_flag] = pp_nominal(p);
        if error_flag
            fprintf('Error H, no distinct eigenvalues\n');
            return
        end
        return
    end

    % step 1 of algorithm
    fprintf('Solving the LMI of step 1\n');
    [S] = pp_step1(p);
    % based on S set T, U for k = 2
    T = p.D(:,1)*p.D(:,1)'*S(1, 1)^2 ...
        + p.D(:,2)*p.D(:,2)'*S(2, 2)^2;
    U = p.E(1,:)'*p.E(1,:)*(1/S(1, 1)^2) ...
        + p.E(2,:)'*p.E(2,:)*(1/S(2, 2)^2)';
    p.('T') = T;
    p.('U') = U;

    % step 2 of algorithm
    [P_u, error_flag] = pp_step2(p);
    if error_flag
        fprintf('Error H, no distinct eigenvalues\n');
        return
    end
    if any(eig(P_u) > eps)
        fprintf('Error, P_u is not negative definite\n');
        return
    end
    DQ = -2*p.a*P_u;
    if all(abs(imag(DQ)) > eps)
        fprintf('Error DQ is complex, LMI does');
        fprintf('not accept complex numbers\n');
        return
    else
        % fixing numerical zeros
        DQ = real(DQ);
        p.('DQ') = DQ;
    end

    % step 3 of algorithm
    fprintf('Solving the LMI of step 3\n');
    [P_f] = pp_step3(p);
    if all(eig(P_f) > eps)
        fprintf('Optimal solution found\n');
        p.('P_f') = P_f;
        % plot data, use 'TEX' for latex plots
        pp_makeplot(p, '');
    end
end
