function L_s = evaluarGradienteLagrangeanoS(vetorS, vetorPi, mu)
    L_s  = vetorS.*vetorPi - mu; % não precisa do vetor de `e` em MATLAB, nem criar matriz diagonal
end
