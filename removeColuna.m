function novaMatriz = removeColuna(linha, matrizOriginal)
    novaMatriz = matrizOriginal;
    novaMatriz(:, linha) = [];
end
