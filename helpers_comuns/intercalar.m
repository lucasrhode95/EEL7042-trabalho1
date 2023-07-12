function vetorIntercalado = intercalar(vetor1, vetor2)
    vetor1 = vetor1(:)'; % transforma em vetor-linha
    vetor2 = vetor2(:)';
    vetorIntercalado = [vetor1; vetor2];
    vetorIntercalado = vetorIntercalado(:);
end
