function mu = calcularMu(vetorS, vetorPi, dimVar, betaAceleracao)
    mu = vetorS'*vetorPi/(2*dimVar*betaAceleracao);
end
