function novaMatriz = removeLinha(linha, matrizOriginal)
    novaMatriz = matrizOriginal;
    novaMatriz(linha, :) = [];
end
