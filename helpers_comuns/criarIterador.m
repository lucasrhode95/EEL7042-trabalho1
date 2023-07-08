function iterador = criarIterador(vetor)
    cursor   = 0;

    function subvetor = proximo(passo)
        subvetor = vetor((1:passo) + cursor);
        cursor += passo;
    end

    iterador = @proximo;
end
