function pp_algorithm()
    % get data
    [p] = pp_data();

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
    [P_u] = pp_step2(p);
    if any(eig(P_u) > eps)
        fprintf('Error, P_u is not negative definite\n');
        return
    end
    DQ = -2*p.a*P_u;
    if all(abs(imag(DQ)) > eps)
        fprintf('Error DQ is complex, LMI does not accept complex numbers\n');
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
        pp_makeplot(p);
    end
end
