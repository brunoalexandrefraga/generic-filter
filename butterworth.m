d = fdesign.lowpass('Fp,Fst,Ap,Ast',2000,3200,0.8,50,48000);
Hd = design(d,'butter');
fvtool(Hd);