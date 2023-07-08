function DeltaZ = iteracaoMetodoNewton(
    wcc, Alpha, Pd0,
    PgMin, PgMax, Tmin, Tmax,
    barraVTheta,
    DeltaPd, Pg, Theta,
    mu, Lambdas, pi1, pi2, pi3, pi4, pi5, s1, s2, s3, s4, s5,
    Um, Ag, A, B, Xinv,
    Winv
)
    vetorPi = [
        pi1; pi2; pi3; pi4; pi5;
    ];
    vetorS = [
        s1; s2; s3; s4; s5;
    ];
    Bred = removeColuna(barraVTheta, B);
    Ared = removeLinha(barraVTheta, A);

    ndes = size(vetorS)(1); % quantidade de restrições de desigualdade
    e    = ones(ndes, 1);

    L_u = [
        wcc*Alpha + Um'*Lambdas - pi5;
        Ag'*Lambdas + pi3 - pi1;
        -B*Lambdas + A*Xinv*(pi4 - pi2);
    ];

    G_u = Ag*Pg - (Pd0 - Um*DeltaPd) - Bred*Theta;
    G_u = removeLinha(barraVTheta, G_u);

    H_u = [
        s1 - Pg + PgMin;
        s2 - Xinv*Ared'*Theta + Tmin;
        s3 + Pg - PgMax;
        s4 + Xinv*Ared'*Theta - Tmax;
        s5 - DeltaPd;
    ];

    S_u = diag(vetorS)*vetorPi - mu*e;

    % p(z) = dL/dx = 0
    Pz = [
        L_u;
        G_u;
        H_u;
        S_u;
    ];

    DeltaZ = -Winv*Pz;
end
