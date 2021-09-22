%%% Combines the valid parameter values acquired in the script
%%% Parameter_Sweep_Transmission_Coeffs

%%% When saving the combined variables use same rho value as considered in
%%% the original Parameter_Sweep_Transmission_Coeffs script

load('var0year')
var1 = var;
load('var03rdyear')
var2 = var;
load('var06thyear')
var3 = var;

var4 = intersect(var1, var2, 'rows');   %intersection of var1 and var2

varintersect = intersect(var3, var4, 'rows');   %intersection of valid
                                    %parameters for all time starts given

save('varallrho083')