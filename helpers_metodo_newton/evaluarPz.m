function Pz = evaluarPz(
    wcc, Alpha, Um, Lambda, Ag, B, A, Xinv, pi1, pi2, pi3, pi4, pi5,
    Pg, Pd0, DeltaPd, Theta,
    PgMin, Tmin, PgMax, Tmax, s1, s2, s3, s4, s5,
    mu,
    barraVTheta
)
    vetorPi = [pi1; pi2; pi3; pi4; pi5];
    vetorS  = [s1; s2; s3; s4; s5];
    Bred    = removeColuna(barraVTheta, B);
    Ared    = removeLinha(barraVTheta, A);

    L_u = evaluarGradienteLagrangeano(wcc, Alpha, Um, Lambda, Ag, Bred, Ared, Xinv, pi1, pi2, pi3, pi4, pi5);
    G_u = evaluarRestricoesIgualdade(Ag, Pg, Pd0, Um, DeltaPd, Bred, Theta);
    H_u = evaluarRestricoesDesigualdade(Pg, PgMin, Xinv, Ared, Theta, Tmin, PgMax, Tmax, DeltaPd, s1, s2, s3, s4, s5);
    S_u = evaluarGradienteLagrangeanoS(vetorS, vetorPi, mu);

    Pz = [
        L_u;
        G_u;
        H_u;
        S_u;
    ];
end
