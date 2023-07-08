function alfa = calcularAlfa(z, deltaZ)
    negativos = z < 0;
    zi        = z(negativos);
    dzi       = deltaZ(negativos);
    dummy     = min(-zi./dzi);
    alfa      = min([dummy, 1]);
end
