function rate = k2Plus(prefactor,force,forceLengthData,delF)
    length = SampleLofF(forceLengthData,force,delF);
    rate = prefactor*exp((length-52)*0.58);
end

