function B = montarMatrizB(impedanciasDeLinha, matrizIncidenciaBarraLinha)
    nb = size(matrizIncidenciaBarraLinha)(1);
    B = zeros(nb, nb);

    admitancias = 1./impedanciasDeLinha;
    for (i = 1:nb)
        for (j = 1:nb)
            % a função `getAdmitanciaEntreBarras` retorna o somatório das
            % impedâncias das linhas conectadas à barra X, caso i=j.
            B(i, j) = getAdmitanciaEntreBarras(i, j, admitancias, matrizIncidenciaBarraLinha);

            % se não for diagonal, multiplica por -1
            if (i ~= j)
                B(i, j) *= -1;
            endif
        endfor
    endfor
end
