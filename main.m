clear();
clc();
addpath(genpath("helpers_topologia_do_sistema"));
addpath(genpath("helpers_metodo_newton"));
addpath(genpath("helpers_evaluar_lagrangeano_e_restricoes"));

nb     = 5; % número de barras
nc     = 3; % número de cargas
ng     = 3; % número de geradores
wcc    = 1; % índice de ponderação
nl     = 6; % número de linhas
nigual = 5; % número de restrições de igualdade (se adicionar mais restrições, é preciso atualizar a matriz W)

% peso dos cortes de carga ($/MWpu)
% (problema não deu valores, assumindo peso uniforme = 1)
Alpha = ones(nc, 1);

%% DEMANDA
% matriz de demanda inicial
Pd0 = [
    1;   % carga 1 na barra 1
    2.5; % carga 2 na barra 2
    1;   % carga 3 na barra 3
    0;   % nenhuma carga na barra 4
    0;   % nenhuma carga na barra 5
];
% matriz de incidência de cargas cortáveis
Um = zeros(nb, nc);
Um(1,1) = 1; % carga 1 na barra 1
Um(2,2) = 1; % carga 2 na barra 2
Um(3,3) = 1; % carga 3 na barra 3

%% GERADORES
% limites geradores
PgMax = [
    5;
    3;
    5;
];
PgMin = [
    0;
    0;
    0;
];
% matriz de incidência de geradores
Ag = zeros(nb, ng);
Ag(1,1) = 1; % gerador 1 na barra 1
Ag(3,2) = 1; % gerador 2 na barra 3
Ag(4,3) = 1; % gerador 3 na barra 4

%% LINHAS
% impedância das linhas
impedanciasDeLinha = [
    0.1680; % de 1 para 2
    0.1260; % de 2 para 3
    0.2100; % de 3 para 5
    0.3360; % de 3 para 4
    0.2520; % de 5 para 4
    0.1260; % de 5 para 1
];
X    = diag(impedanciasDeLinha);
Xinv = inv(X);
% limite de fluxo das linhas
Tmax = [
    1.5; % de 1 para 2
    1.5; % de 2 para 3
    1.5; % de 3 para 5
    1.5; % de 3 para 4
    1.5; % de 5 para 4
    0.6; % de 5 para 1
];
Tmin = -Tmax;

% matrix incidência barra-ramo A(barra, linha)
A = zeros(nb, nl);
% linha 1 (b1->b2)
A(1, 1) =  1; % barra 1
A(2, 1) = -1; % barra 2
% linha 2 (b2->b3)
A(2, 2) =  1; % barra 2
A(3, 2) = -1; % barra 3
% linha 3 (b3->b5)
A(3, 3) =  1; % barra 3
A(5, 3) = -1; % barra 5
% linha 4 (b3->b4)
A(3, 4) =  1; % barra 3
A(4, 4) = -1; % barra 4
%linha 5 (b4->b5)
A(5, 5) = -1; % barra 5
A(4, 5) =  1; % barra 4
% linha 6 (b1->b5)
A(1, 6) =  1; % barra 1
A(5, 6) = -1; % barra 5

% matriz B
B = montarMatrizB(impedanciasDeLinha, A);

%% Definição da Barra Vθ
nr          = nb - 1; % quantidade de barras menos barra Vθ
barraVTheta = 1;


%% MPI
% definição de chute inicial
mu      = 0.1;
Pg      = (PgMax+PgMin)/2; % geração a 50% da capacidade
DeltaPd = 0.05*Um'*Pd0;    % corte de carga inicial = 5%
Theta   = zeros(nr, 1);    % ângulos = 0 rad
Lambda  = ones(nigual, 1); % lambdas = 1;

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

restoredefaultpath();
