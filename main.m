clear();
clc();
addpath(genpath("helpers_comuns"));
addpath(genpath("helpers_topologia_do_sistema"));
addpath(genpath("helpers_metodo_newton"));
addpath(genpath("helpers_evaluar_lagrangeano_e_restricoes"));
addpath(genpath("helpers_MPI"));

%% QUANTIDADE DE PATAMARES / DISCRETIZAÇÃO TEMPORAL
nT = 3;
T  = ones(nT, 1);
%% LIMITE CAPACIDADE BATERIA
Cap = 500*T;
%% PATAMARES DE CARGA
Pd = [1.188; 2.960; 4.118];
% GERAÇÃO FOTOVOLTAICA
Pgsolar = [0; 3.5; 0];
%% LIMITE DE POTÊNCIA MÁXIMA
maxTermica = 7*T;
maxBateria = 3*T;
PgMax      = intercalar(maxTermica, maxBateria);
%% LIMITE DE POTÊNCIA MÍNIMA
minTermica = 0*T;
minBateria = -3*T;
PgMin      = intercalar(minTermica, minBateria);
%% CUSTOS DE GERAÇÃO
aTermica = 5*T;
bTermica = 2*T;
aBateria = 1e-9*T;
bBateria = 0*T;
Q        = diag(intercalar(aTermica, aBateria));
b        = intercalar(bTermica, bBateria);

%% RESTRIÇÃO DE IGUALDADE
% Iaux*Pg + Pgsolar - Pd = 0
Iaux = [
    1  -1   0   0  0   0;
    0   0   1  -1  0   0;
    0   0   0   0  1  -1;
];
%% RESTRIÇÃO DESIGUALDADE
%  -Pg + Pgmin <= 0
%  Pg  - Pgmax <= 0
%  -Mat_hor*Pg <= Cap
Mat_hor = [
    0   148   0     0   0   0;
    0   148   0   428   0   0;
    0   148   0   428   0 168;
];


%% MPI
% fatores de aceleração
sigmaInterioridade = 0.9995;
betaAceleracao     = 10;

% definição de chute inicial
mu      = 0.1;
Pg      = (PgMin + PgMax)/2;
Lambda  = ones(nT, 1); % lambdas = 1;
% variáveis de folga S
s1     = Pg - PgMin;
s2     = Cap - Mat_hor*Pg;
s3     = -Pg + PgMax;
s4     = Mat_hor*Pg;
vetorS = [s1; s2; s3; s4];
% variáveis de folga PI
pi1     = mu./s1;
pi2     = mu./s2;
pi3     = mu./s3;
pi4     = mu./s4;
vetorPi = [pi1; pi2; pi3; pi4];

% condições de parada
limiteIteracoes = 100;
toleranciaPz    = 1e-6;
toleranciaMu    = 1e-6;

u = Pg;
iteracoes        = 0;
alcancouPrecisao = false;
while (iteracoes < limiteIteracoes) && ~alcancouPrecisao
    iteracoes++;

    ITERACAO_MPI

    alcancouPrecisao = norm(Pz) <= toleranciaPz && mu <= toleranciaMu; % atingiu precisão desejada?
end

if alcancouPrecisao
    fprintf("Solução em %d iterações:\n", iteracoes);
else
    fprintf("[AVISO] PROBLEMA NÃO CONVERGIU DEPOIS DE %d ITERAÇÕES\n\n", limiteIteracoes);
endif

printVetor(" Pg", u);

restoredefaultpath();
