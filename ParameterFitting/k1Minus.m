function rate = k1Minus(prefactor,theta,force,deltaR,kBT)
 rate = prefactor*exp((1-theta)*force*deltaR/kBT);
end

