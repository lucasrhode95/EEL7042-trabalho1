function DeltaZ = iteracaoMetodoNewton(
    Q, Pg, b, Iaux, Mat_hor, Pgsolar, Pd, PgMin, Cap, PgMax,
    mu, Lambda, pi1, pi2, pi3, pi4, s1, s2, s3, s4
)
    % reconstrução de variáveis auxiliares
    vetorPi = [pi1; pi2; pi3; pi4];
    vetorS  = [s1; s2; s3; s4];

    % matriz hessiana
    W = montarMatrizW(Q, Iaux, Mat_hor, vetorS, vetorPi);

    % p(z) + W*dz = 0
    Pz = evaluarPz(
        Q, Pg, b, Iaux, Mat_hor, Pgsolar, Pd, PgMin, Cap, PgMax,
        mu, Lambda, pi1, pi2, pi3, pi4, s1, s2, s3, s4
    );

    DeltaZ = -W\Pz;
end
