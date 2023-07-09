% vetorS = [s1; s2; s3; s4; s5]
%
% s1 = s_Pgmin
% s2 = s_Tmin
% s3 = s_Pgmax
% s4 = s_Tmax
% s5 = s_Pdmin
%
% π1 = π_Pgmin
% π2 = π_Tmin
% π3 = π_Pgmax
% π4 = π_Tmax
% π5 = π_Pdmin
function L_s = evaluarGradienteLagrangeanoS(vetorS, vetorPi, mu)
    % L_s  = -mu./vetorS + vetorPi;
    L_s  = vetorS.*vetorPi - mu; % não precisa do vetor `e` em MATLAB, nem criar matriz diagonal
end
