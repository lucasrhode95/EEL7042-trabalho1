function [vetor1, vetor2] = desintercalar(vetor)
    vetor = vetor(:)'; % transforma em vetor-linha

    vetor1 = vetor(1:2:end);
    vetor2 = vetor(2:2:end);
end
