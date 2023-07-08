clear();
clc();

nb     = 5; % número de barras
nc     = 3; % número de cargas
ng     = 3; % número de geradores
wcc    = 1; % índice de ponderação
nl     = 6; % número de linhas
nigual = nb; % número de restrições de igualdade;
ndes   = 2*ng + 2*nl + nc; % número de restrições de desigualdade

% peso dos cortes de carga ($/MWpu)
% (problema não deu valores, assumindo peso uniforme = 1)
alpha = ones(nc, 1);

%% DEMANDA
% matriz de demanda inicial
Pd0 = zeros(nb, 1);
Pd0(1) = 1;   % carga 1 na barra 1
Pd0(2) = 2.5; % carga 2 na barra 2
Pd0(3) = 1;   % carga 3 na barra 3
% matriz de incidência de cargas cortáveis
Um = zeros(nb, nc);
Um(1,1) = 1; % carga 1 na barra 1
Um(2,2) = 1; % carga 2 na barra 2
Um(3,3) = 1; % carga 3 na barra 3

%% GERADORES
% limites geradores
PgMax = [5, 3, 5]';
PgMin = [0, 0, 0]';
% matriz de incidência de geradores
Ag = zeros(nb, ng);
Ag(1,1) = 1; % gerador 1 na barra 1
Ag(3,2) = 1; % gerador 2 na barra 3
Ag(4,3) = 1; % gerador 3 na barra 4

%% LINHAS
% impedância das linhas
X = zeros(nl, nl);
X(1, 1) = 0.1680; % de 1 para 2
X(2, 2) = 0.1260; % de 2 para 3
X(3, 3) = 0.2100; % de 3 para 5
X(4, 4) = 0.3360; % de 3 para 4
X(5, 5) = 0.2520; % de 5 para 4
X(6, 6) = 0.1260; % de 5 para 1
Xinv = inv(X);
% limite de fluxo das linhas
Tmax = zeros(nl, 1);
Tmax(1) = 1.5; % de 1 para 2
Tmax(2) = 1.5; % de 2 para 3
Tmax(3) = 1.5; % de 3 para 5
Tmax(4) = 1.5; % de 3 para 4
Tmax(5) = 1.5; % de 5 para 4
Tmax(6) = 0.6; % de 5 para 1
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
B = montarMatrizB(diag(X), A);

%% Barra Vθ
nr = nb - 1; % quantidade de barras menos barra Vθ
barraVTheta = 1;
Bred = removeColuna(barraVTheta, B);
Ared = removeLinha(barraVTheta, A);
Agred = removeLinha(barraVTheta, Ag);



%% MPI
mu      = 0.1;
Pg      = (PgMax+PgMin)/2;
DeltaPd = 0.05*Um'*Pd0; % corte de carga inicial = 5%
Theta   = zeros(nr, 1);
Lambdas = ones(nr, 1);


% variáveis de folga S
s1 = Pg - PgMin;
s2 = -Tmin;
s3 = -Pg + PgMax;
s4 = Tmax;
s5 = DeltaPd;
vetorS = [s1; s2; s3; s4; s5];
% variáveis de folga PI
pi1 = s1./mu;
pi2 = s2./mu;
pi3 = s3./mu;
pi4 = s4./mu;
pi5 = s5./mu;
vetorPi = [pi1; pi2; pi3; pi4; pi5];


%% MATRIZ HESSIANA
W = montarMatrizW(Um, Ag, Ared, Xinv, Bred, vetorS, vetorPi);
