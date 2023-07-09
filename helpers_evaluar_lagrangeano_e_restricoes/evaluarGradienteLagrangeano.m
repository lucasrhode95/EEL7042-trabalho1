% gradiente de L com relação às variáveis de otimização
% ∇L(u)
function Lu = evaluarGradienteLagrangeano(
    wcc, Alpha, Um, Lambda,
    Ag, Bred, Ared, Xinv,
    pi1, pi2, pi3, pi4, pi5
)
    Lu = [
        wcc*Alpha - Um'*Lambda - pi5;         % deltaPd
        -Ag'*Lambda + pi3 - pi1;              % Pg
        Bred'*Lambda + Ared*Xinv*(pi4 - pi2); % Theta
    ];
end
