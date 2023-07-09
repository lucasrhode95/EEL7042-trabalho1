% s1 = s_Pgmin
% s2 = s_Tmin
% s3 = s_Pgmax
% s4 = s_Tmax
% s5 = s_Pdmin
function H = evaluarRestricoesDesigualdade(
    Pg, PgMin, Xinv, Ared, Theta, Tmin, PgMax, Tmax, DeltaPd,
    s1, s2, s3, s4, s5
)
    H = [
        s1 - Pg + PgMin;
        s2 - Xinv*Ared'*Theta + Tmin;
        s3 + Pg - PgMax;
        s4 + Xinv*Ared'*Theta - Tmax;
        s5 - DeltaPd;
    ];
end
