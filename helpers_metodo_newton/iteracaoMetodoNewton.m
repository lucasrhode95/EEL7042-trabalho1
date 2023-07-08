function DeltaZ = iteracaoMetodoNewton(
    wcc, Alpha, Pd0,
    PgMin, PgMax, Tmin, Tmax,
    barraVTheta,
    DeltaPd, Pg, Theta,
    mu, Lambda, pi1, pi2, pi3, pi4, pi5, s1, s2, s3, s4, s5,
    Um, Ag, A, B, Xinv
)
    % reconstrução de variáveis auxiliares
    vetorPi = [pi1; pi2; pi3; pi4; pi5];
    vetorS  = [s1; s2; s3; s4; s5];
    Bred    = removeColuna(barraVTheta, B);
    Ared    = removeLinha(barraVTheta, A);

    % matriz hessiana
    W    = montarMatrizW(Um, Ag, Ared, Xinv, Bred, vetorS, vetorPi);
    Winv = inv(W);

    % p(z) + W*dz = 0
    Pz = evaluarPz(
        wcc, Alpha, Um, Lambda, Ag, B, A, Xinv, pi1, pi2, pi3, pi4, pi5,
        Pg, Pd0, DeltaPd, Theta,
        PgMin, Tmin, PgMax, Tmax, s1, s2, s3, s4, s5,
        mu,
        barraVTheta
    );

    DeltaZ = -Winv*Pz;
end
