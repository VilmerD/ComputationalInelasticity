EPS = eps;
load('checkroutines_Lab.mat');

% Test stress
my_es = my_stress(eps, mp);
rel_des = (my_es - es)./es;
ndes = norm(rel_des);
grade = ndes < sqrt(EPS);
if grade, es_grade  = 'Passed'; else, es_grade = 'Failed'; end
clear grade rel_des ndes

% Test material tangent
my_Dt = my_tangent(eps, mp);
rel_dDt = (my_Dt - Dt)./Dt;
ndDt = norm(rel_dDt, 'fro');
grade = ndDt < sqrt(EPS);
if grade, d_grade = 'Passed'; else, d_grade = 'Failed'; end
clear grade rel_dD ndD

% Print results
fprintf(['my_stress: ' es_grade '\n']);
fprintf(['my_tangent: ' d_grade '\n']);
clear es_grade d_grade EPS