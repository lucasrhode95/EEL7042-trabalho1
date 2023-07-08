% retorna lista de quais linhas estão conectadas à barra X
function linhasConectadas = getLinhasConectadasABarra(barraX, matrizIncidenciaBarraLinha)
    linhasConectadas = [];
    nl = size(matrizIncidenciaBarraLinha)(2);

    for (linhaX = 1:nl)
        incidencia = matrizIncidenciaBarraLinha(barraX, linhaX);
        if (incidencia ~= 0)
            linhasConectadas(end+1) = linhaX;
        endif
    endfor
end
