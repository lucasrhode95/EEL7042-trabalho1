function G = evaluarRestricoesIgualdade(Ag, Pg, Pd0, Um, DeltaPd, Bred, Theta)
    G = Ag*Pg - Pd0 + Um*DeltaPd - Bred*Theta;
end
