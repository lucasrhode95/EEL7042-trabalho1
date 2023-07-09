% Essa função de montagem da matriz W não tem como foco a performance, mas sim
% o entendimento do procedimento de montagem da matriz Hessiana, ou matriz W.
%
% Uma grande melhoria de performance seria possível mantendo um cache da matriz
% W gerada, e somente atualizar as células referentes ao vetorS e vetorPi.
function W = montarMatrizW(Um, Ag, Ared, Xinv, Bred, vetorS, vetorPi)
    % para não precisar passar todos argumentos, inferimos usando as dimensões
    % dos argumentos recebidos
    nb     = size(Ag)(1);
    ng     = size(Ag)(2);
    nc     = size(Um)(2);
    nl     = size(Ared)(2);
    nr     = nb - 1;
    nigual = nb;
    nvar   = nc + ng + nr;    % quantidade de variáveis otimizadas (ou dimensão do vetor u)
    ndes   = size(vetorS)(1); % quantidade de restrições de desigualdade


    L_u_u = [
    %    deltaPd           Pg           Theta
        zeros(nc,nc), zeros(nc,ng), zeros(nc,nr); % deltaPd
        zeros(ng,nc), zeros(ng,ng), zeros(ng,nr); % Pg
        zeros(nr,nc), zeros(nr,ng), zeros(nr,nr); % Theta
    ];

    L_u_y = [
    %    lambda
        -Um';  % deltaPd
        -Ag';  % Pg
        Bred';  % Theta
    ];

    % Obs: transformamos o problema para que só houvessem PI máx e nenhum PI min
    L_u_pi = [
    %       pi1           pi2           pi3           pi4           pi5
        zeros(nc,ng), zeros(nc,nl), zeros(nc,ng), zeros(nc,nl),  -eye(nc,nc); % deltaPd
         -eye(ng,ng), zeros(ng,nl),   eye(ng,ng), zeros(ng,nl), zeros(ng,nc); % Pg
        zeros(nr,ng),   -Ared*Xinv, zeros(nr,ng),    Ared*Xinv, zeros(nr,nc); % Theta
    ];

    L_u_s   = zeros(nvar, ndes);
    L_y_y   = zeros(nigual, nigual);
    L_y_pi  = zeros(nigual, ndes);
    L_y_s   = zeros(nigual, ndes);
    L_pi_pi = zeros(ndes, ndes);
    L_pi_s  = eye(ndes, ndes);
    L_s_pi  = diag(vetorS);
    L_s_s   = diag(vetorPi);

    W = [
        L_u_u    L_u_y   L_u_pi  L_u_s;
        L_u_y'   L_y_y   L_y_pi  L_y_s;
        L_u_pi'  L_y_pi' L_pi_pi L_pi_s;
        L_u_s'   L_y_s'  L_s_pi  L_s_s;
    ];
end
