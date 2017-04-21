function [err, tr_p] = gcc_linear_bound_all()
    A = [-1, 0; 1, -2];
    B = [0; 1];
    Q = eye(size(A, 1));
    R = 10;

    A_1 = @(r1) r1*[0, 0; 3.0, 0];
    A_2 = @(r2) r2*[0, 0; 0, 1.5];
    DA = @(r1, r2) A_1(r1) + A_2(r2);
    U = @(P, r1, r2, e) e*P + ...
	(1/e)*(DA(r1, r2)'*P*DA(r1, r2));

    options = optimset('Display','off');
    err = [];
    Pnot = care(A, B, Q, R);
    cnt = 0;
    for e = 0.1:0.1:1
	[P, fval] = fsolve(@(p) gare(p, A, B, Q, R, ...
	    U(p, 1, 1, e)), Pnot, options);
	if all(eig(P) > eps)
	    eig(P)
	    cnt = cnt + 1;
	    err(cnt) = e;
	    tr_p(cnt) = trace(P);
	end
    end
end

function f = gare(p, A, B, Q, R, U)
    P = reshape(p(1:length(A)^2), size(A));
    F = P*A + A'*P - P*B*(R\(B'*P)) + Q + U;
    f = F(:);
end
