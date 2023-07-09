% contagem das variáveis primais e duais
dimDeltaPd = size(DeltaPd)(1);
dimPg      = size(Pg)(1);
dimTheta   = size(Theta)(1);
dimU       = dimDeltaPd + dimPg + dimTheta;
dimLambda  = size(Lambda)(1);
dimPi      = size(vetorPi)(1);
dimS       = size(vetorS)(1);
dimVar     = dimU + dimLambda + dimPi + dimS;

% deltaZ = [
%     deltaU;  % dimU
%     deltaY;  % dimLambda
%     deltaPi; % dimPi
%     deltaS;  % dimS
% ];
deltaZ = iteracaoMetodoNewton(
    wcc, Alpha, Pd0,
    PgMin, PgMax, Tmin, Tmax,
    barraVTheta,
    DeltaPd, Pg, Theta,
    mu, Lambda, pi1, pi2, pi3, pi4, pi5, s1, s2, s3, s4, s5,
    Um, Ag, A, B, Xinv
);
deltaZ = sigmaInterioridade*deltaZ; % para garantir interioridade, conforme apostila

% separar variáveis MPI
next        = criarIterador(deltaZ);
deltaU      = next(dimU);
deltaLambda = next(dimLambda);
deltaPi     = next(dimPi);
deltaS      = next(dimS);

% calcular sigmas e mu
alfaPrimal = calcularAlfa(vetorS, deltaS);
alfaDual   = calcularAlfa(vetorPi, deltaPi);

% atualizar variáveis
u       += alfaPrimal*deltaU;
Lambda  += alfaDual*deltaLambda;
vetorPi += alfaDual*deltaPi;
vetorS  += alfaPrimal*deltaS;

% atualizar mi
mu = calcularMu(vetorS, vetorPi, dimVar, betaAceleracao);

% separar/atualizar variáveis originais
next    = criarIterador(u);
DeltaPd = next(dimDeltaPd);
Pg      = next(dimPg);
Theta   = next(dimTheta);

next = criarIterador(vetorS);
s1 = next(size(s1)(1));
s2 = next(size(s2)(1));
s3 = next(size(s3)(1));
s4 = next(size(s4)(1));
s5 = next(size(s5)(1));

next = criarIterador(vetorPi);
pi1 = next(size(pi1)(1));
pi2 = next(size(pi2)(1));
pi3 = next(size(pi3)(1));
pi4 = next(size(pi4)(1));
pi5 = next(size(pi5)(1));

% condição de tolerância
Pz  = evaluarPz(
    wcc, Alpha, Um, Lambda, Ag, B, A, Xinv, pi1, pi2, pi3, pi4, pi5,
    Pg, Pd0, DeltaPd, Theta,
    PgMin, Tmin, PgMax, Tmax, s1, s2, s3, s4, s5,
    mu,
    barraVTheta
);
