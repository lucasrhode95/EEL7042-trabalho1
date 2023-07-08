function Lu = evaluarGradienteLagrangeano(
    wcc, Alpha, Um, Lambda,
    Ag, B, A, Xinv,
    pi1, pi2, pi3, pi4, pi5
)
    Lu = [
        wcc*Alpha + Um'*Lambda - pi5;
        Ag'*Lambda + pi3 - pi1;
        -B*Lambda + A*Xinv*(pi4 - pi2);
    ];
end
