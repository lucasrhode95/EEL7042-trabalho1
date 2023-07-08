% retorna lista de barras conectadas à barra X
%
% acabou não sendo utilizada no cálculo do FPO, mas pode ser útil no futuro
function barrasAdjacentes = getBarrasAdjacentes(barraX, matrizIncidenciaBarraLinha)
    barrasAdjacentes = [];
    nb = size(matrizIncidenciaBarraLinha)(1);
    nl = size(matrizIncidenciaBarraLinha)(2);

    linhasConectadasABarra = getLinhasConectadasABarra(barraX, matrizIncidenciaBarraLinha);
    for (linha = linhasConectadasABarra)
        for (barra = 1:nb)
            incidencia = matrizIncidenciaBarraLinha(barra, linha);
            if (barra == barraX || incidencia == 0)
                continue;
            endif

            barrasAdjacentes(end+1) = barra;
        end
    endfor
end
