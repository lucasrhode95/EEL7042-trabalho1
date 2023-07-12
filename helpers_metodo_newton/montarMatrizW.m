% Essa função de montagem da matriz W não tem como foco a performance, mas sim
% o entendimento do procedimento de montagem da matriz Hessiana, ou matriz W.
%
% Uma grande melhoria de performance seria possível mantendo um cache da matriz
% W gerada, e somente atualizar as células referentes ao vetorS e vetorPi.
function W = montarMatrizW(Q, Iaux, Mat_hor, vetorS, vetorPi)
    nt     = size(Q)(1);
    nigual = size(Mat_hor)(1);
    nvar   = size(Mat_hor)(2);
    ndes   = size(vetorS)(1);

    L_u_u   = 2*Q;
    L_u_y   = Iaux';
    L_u_pi  = [-eye(nvar, nt), Mat_hor', eye(nvar, nt), -Mat_hor'];
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
