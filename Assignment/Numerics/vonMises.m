function esv = vonMises(sig)
sigkk = sig(1) + sig(2) + sig(3);
s11 = sig(1) - sigkk/3;
s22 = sig(2) - sigkk/3;
s33 = sig(3) - sigkk/3;
s12 = sig(4);
esv = sqrt(3/2*(s11.^2 + s22.^2 + s33.^2 + 2*s12.^2));
end