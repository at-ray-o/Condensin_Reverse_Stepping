function dwellProb = dwFwd2(k1p,t)
    dwellProb = k1p.*exp(-k1p.*t);
end