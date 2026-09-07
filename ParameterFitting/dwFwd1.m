function dwellProb = dwFwd1(k1p,k1m,k2p,t)
    dwellProb = k1p.*(k1m+k2p).*(exp(-k1p.*t)-exp(-(k1m+k2p).*t))./(k1m-k1p+k2p);
end