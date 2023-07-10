clear();
clc();
addpath(genpath("helpers_comuns"));
addpath(genpath("helpers_topologia_do_sistema"));
addpath(genpath("helpers_metodo_newton"));
addpath(genpath("helpers_evaluar_lagrangeano_e_restricoes"));
addpath(genpath("helpers_MPI"));

nb     = 5; % número de barras
nc     = 3; % número de cargas
ng     = 3; % número de geradores
nl     = 6; % número de linhas
nigual = 5; % número de restrições de igualdade (se adicionar mais restrições, é preciso atualizar a matriz W)

% peso dos cortes de carga ($/MWpu)
% (problema não deu valores, assumindo peso uniforme = 1)
wcc   = 1; % índice de ponderação
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
% fatores de aceleração
sigmaInterioridade = 0.9995;
betaAceleracao     = 10;

% definição de chute inicial
mu      = 0.1;
Pg      = (PgMin + PgMax)/2;
DeltaPd = [1 1 1]';
Theta   = zeros(nr, 1);
Lambda  = ones(nigual, 1); % lambdas = 1;
% variáveis de folga S
s1     = Pg - PgMin;
s2     = -Tmin;
s3     = -Pg + PgMax;
s4     = Tmax;
s5     = DeltaPd;
vetorS = [s1; s2; s3; s4; s5];
% variáveis de folga PI
pi1     = mu./s1;
pi2     = mu./s2;
pi3     = mu./s3;
pi4     = mu./s4;
pi5     = mu./s5;
vetorPi = [pi1; pi2; pi3; pi4; pi5];

% condições de parada
limiteIteracoes = 100;
toleranciaPz    = 1e-6;
toleranciaMu    = 1e-6;

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

printVetor("ΔPd", DeltaPd);
printVetor(" Pg", Pg);
printVetor("  θ", Theta);
sum(DeltaPd)
sum(Pg)

restoredefaultpath();
