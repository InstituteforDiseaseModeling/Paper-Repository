function [counter, localExts] = firstStage(ic, k, gammas, target_e, event_type, Kce, ofile_prefix)

x0 = ic.x0;
t0 = ic.t0;
tf = ic.tf;
counter   = 0;
localExts = zeros(1,Kce);

switch ofile_prefix
    case 'revIsom'
        parfor i=1:Kce
            [counter_out, localExt]  = RI_solveOnce_ie(x0, t0, tf, k, gammas, target_e, event_type);
            if counter_out
                counter = counter + counter_out;
            end
            localExts(i) = localExt;
        end
        
    case 'birthDeath'
        parfor i=1:Kce
            [counter_out, localExt]  = BD_solveOnce_ie(x0, t0, tf, k, gammas, target_e, event_type);
            if counter_out
                counter = counter + counter_out;
            end
            localExts(i) = localExt;
        end
        
    case 'SIRS'
        parfor i=1:Kce
            [counter_out, localExt]  = SIRS_solveOnce_ie(x0, t0, tf, k, gammas, target_e, event_type);
            if counter_out
                counter = counter + counter_out;
            end
            localExts(i) = localExt;
        end
end


