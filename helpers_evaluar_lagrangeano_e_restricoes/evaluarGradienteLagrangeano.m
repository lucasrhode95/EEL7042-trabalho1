function Lu = evaluarGradienteLagrangeano(Q, Pg, b, Iaux, Mat_hor, Lambda, pi1, pi2, pi3, pi4)
    Lu = 2*Q*Pg + b + Iaux'*Lambda + pi3 - pi1 + Mat_hor'*(pi2 - pi4);
end
