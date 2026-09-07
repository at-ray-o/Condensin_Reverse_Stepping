function rate = k1Plus(prefactor,theta,force,deltaR,kBT)
 rate = prefactor*exp(-theta*force*deltaR/kBT);
end

