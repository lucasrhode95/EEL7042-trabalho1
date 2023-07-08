
% variáveis de folga S
s1     = Pg - PgMin;
s2     = -Tmin;
s3     = -Pg + PgMax;
s4     = Tmax;
s5     = DeltaPd;
vetorS = [s1; s2; s3; s4; s5];
% variáveis de folga PI
pi1     = s1./mu;
pi2     = s2./mu;
pi3     = s3./mu;
pi4     = s4./mu;
pi5     = s5./mu;
vetorPi = [pi1; pi2; pi3; pi4; pi5];

% contagem das variáveis primais e duais
dimDeltaPd = size(DeltaPd)(1);
dimPg      = size(Pg)(1);
dimTheta   = size(Theta)(1);
dimU       = dimDeltaPd + dimPg + dimTheta;
dimY       = size(Lambda)(1);
dimPi      = size(vetorPi)(1);
dimS       = size(vetorS)(1);
dimVar     = dimU + dimY + dimPi + dimS;

% método de Newton
u = [
    DeltaPd; % dimDeltaPd
    Pg;      % dimPg
    Theta;   % dimTheta
];
% deltaZ = [
%     deltaU;  % dimU
%     deltaY;  % dimY
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
deltaS      = next(dimS);
deltaLambda = next(dimY);
deltaPi     = next(dimPi);

% calcular sigmas e mu
alfaPrimal = calcularAlfa(vetorS, deltaS);
alfaDual   = calcularAlfa(vetorPi, deltaPi);
mu         = calcularMu(vetorS, vetorPi, dimVar, betaAceleracao);

% atualizar variáveis
u       =       u + alfaPrimal*deltaU;
vetorS  =  vetorS + alfaPrimal*deltaS;
Lambda  =  Lambda + alfaDual*deltaLambda;
vetorPi = vetorPi + alfaDual*deltaPi;

% condição de tolerância
Pz  = evaluarPz(
    wcc, Alpha, Um, Lambda, Ag, B, A, Xinv, pi1, pi2, pi3, pi4, pi5,
    Pg, Pd0, DeltaPd, Theta,
    PgMin, Tmin, PgMax, Tmax, s1, s2, s3, s4, s5,
    mu,
    barraVTheta
);

% separar/atualizar variáveis originais
next    = criarIterador(u);
DeltaPd = next(dimDeltaPd);
Pg      = next(dimPg);
Theta   = next(dimTheta);
