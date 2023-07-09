% gradiente de L com relação às variáveis de otimização
% ∇L(u)
%
% π1 = π_Pgmin
% π2 = π_Tmin
% π3 = π_Pgmax
% π4 = π_Tmax
% π5 = π_Pdmin
function Lu = evaluarGradienteLagrangeano(
    wcc, Alpha, Um, Lambda,
    Ag, Bred, Ared, Xinv,
    pi1, pi2, pi3, pi4, pi5
)
    Lu = [
        wcc*Alpha + Um'*Lambda - pi5;         % deltaPd
        Ag'*Lambda + pi3 - pi1;              % Pg
        -Bred'*Lambda + Ared*Xinv*(pi4 - pi2); % Theta
    ];
end
