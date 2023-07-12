function Pz = evaluarPz(
    Q, Pg, b, Iaux, Mat_hor, Pgsolar, Pd, PgMin, Cap, PgMax,
    mu, Lambda, pi1, pi2, pi3, pi4, s1, s2, s3, s4
)
    vetorPi = [pi1; pi2; pi3; pi4];
    vetorS  = [s1; s2; s3; s4];

    L_u = evaluarGradienteLagrangeano(Q, Pg, b, Iaux, Mat_hor, Lambda, pi1, pi2, pi3, pi4);
    G_u = evaluarRestricoesIgualdade(Iaux, Pg, Pgsolar, Pd);
    H_u = evaluarRestricoesDesigualdade(Pg, PgMin, Mat_hor, Cap, PgMax, s1, s2, s3, s4);
    S_u = evaluarGradienteLagrangeanoS(vetorS, vetorPi, mu);

    Pz = [
        L_u;
        G_u;
        H_u;
        S_u;
    ];
end
