function G = evaluarRestricoesIgualdade(Iaux, Pg, Pgsolar, Pd)
    G = Iaux*Pg + Pgsolar - Pd;
end
