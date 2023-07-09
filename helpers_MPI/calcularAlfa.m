function alfa = calcularAlfa(z, deltaZ)
    negativos = z < 0;
    zi        = z(negativos);
    dzi       = deltaZ(negativos);
    dummy     = -zi./dzi;
    alfa      = min([dummy; 1]);
end
