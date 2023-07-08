function DeltaZ = iteracaoMetodoNewton(
    wcc, Alpha, Pd0,
    PgMin, PgMax, Tmin, Tmax,
    barraVTheta,
    DeltaPd, Pg, Theta,
    mu, Lambda, pi1, pi2, pi3, pi4, pi5, s1, s2, s3, s4, s5,
    Um, Ag, A, B, Xinv,
    Winv
)
    % reconstrução de variáveis auxiliares
    vetorPi = [pi1; pi2; pi3; pi4; pi5];
    vetorS  = [s1; s2; s3; s4; s5];
    Bred    = removeColuna(barraVTheta, B);
    Ared    = removeLinha(barraVTheta, A);
    ndes    = size(vetorS)(1); % quantidade de restrições de desigualdade
    e       = ones(ndes, 1);

    % matriz hessiana
    W    = montarMatrizW(Um, Ag, Ared, Xinv, Bred, vetorS, vetorPi);
    Winv = inv(W);

    L_u = evaluarGradienteLagrangeano(wcc, Alpha, Um, Lambda, Ag, B, A, Xinv, pi1, pi2, pi3, pi4, pi5);
    G_u = evaluarRestricoesIgualdade(Ag, Pg, Pd0, Um, DeltaPd, Bred, Theta);
    H_u = evaluarRestricoesDesigualdade(Pg, PgMin, Xinv, Ared, Theta, Tmin, PgMax, Tmax, DeltaPd, s1, s2, s3, s4, s5);
    S_u = diag(vetorS)*vetorPi - mu*e;

    % p(z) + W*dz = 0
    Pz = [
        L_u;
        removeLinha(barraVTheta, G_u);
        H_u;
        S_u;
    ];

    DeltaZ = -Winv*Pz;
end
