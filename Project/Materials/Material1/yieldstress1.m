function sy = yieldstress1(peteff, dmat)
% m1yieldstress computes the yield stress for the given state
%
% INPUT
%       peteff:     the plastic effective strain
%       dmat:       matieral parameters,
%                       - kinf
%                       - h
%                       - sy0:  Initial yield stress
%
% OUTPUT
%       sy:         the new yield stress
kinf = dmat(1);
h = dmat(2);
sy0 = dmat(3);

sy = sy0 + kinf*(1 - exp(-h/kinf*peteff));
end