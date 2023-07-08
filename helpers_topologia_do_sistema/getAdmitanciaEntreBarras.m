% retorna a admitância que conecta duas barras. Caso barras não se conectem,
% retorna 0. Caso se conectem por mais de uma linha em paralelo, retorna o
% somatório das admitâncias (i.e. linhas em paralelo)
%
% caso seja passado argumentos tal que barra1=barra2, a função retorna o
% somatório das admitâncias de todas as linhas que se conectam à barra. O que
% condiz diretamente com a diagonal da matriz B
function admitancia = getAdmitanciaEntreBarras(barra1, barra2, admitancias, matrizIncidenciaBarraLinha)
    linhasConectadasABarra1 = getLinhasConectadasABarra(barra1, matrizIncidenciaBarraLinha);
    linhasConectadasABarra2 = getLinhasConectadasABarra(barra2, matrizIncidenciaBarraLinha);

    linhaComumAsDuasBarras = intersect(linhasConectadasABarra1, linhasConectadasABarra2);
    if (isempty(linhaComumAsDuasBarras))
        admitancia = 0;
    else
        admitancia = sum(admitancias(linhaComumAsDuasBarras));
    endif
end
