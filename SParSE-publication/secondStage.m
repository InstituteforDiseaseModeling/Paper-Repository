function [num, denom, counter] = secondStage(ic, k, ies, Kce, ofile_prefix)

x0 = ic.x0;
t0 = ic.t0;
tf = ic.tf;

M       = length(k);
ie_len  = length(ies);
num     = zeros(ie_len,M);
denom   = zeros(ie_len,M);
counter = zeros(1,ie_len);

switch ofile_prefix
    case 'revIsom'
        parfor i=1:Kce
            [num_out, denom_out, counter_out] = RI_solveOnce_gamma(x0, t0, tf, k, ies);
            counter = counter + counter_out;
            num     = num + num_out;
            denom   = denom + denom_out;
        end
        
    case 'birthDeath'
        parfor i=1:Kce
            [num_out, denom_out, counter_out] = BD_solveOnce_gamma(x0, t0, tf, k, ies);
            counter = counter + counter_out;
            num     = num + num_out;
            denom   = denom + denom_out;
        end
        
    case 'SIRS'
        parfor i=1:Kce
            [num_out, denom_out, counter_out] = SIRS_solveOnce_gamma(x0, t0, tf, k, ies);
            counter = counter + counter_out;
            num     = num + num_out;
            denom   = denom + denom_out;
        end
end
