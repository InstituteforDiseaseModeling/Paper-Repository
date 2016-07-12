function localExts = opStage(ic, k, event_type, Kce, ofile_prefix)

x0 = ic.x0;
t0 = ic.t0;
tf = ic.tf;
localExts = zeros(1,Kce);

switch ofile_prefix
    case 'revIsom'
        parfor gd=1:Kce
            localExts(gd) = RI_solveOnce_ie_op(x0, t0, tf, k, event_type);
        end
        
    case 'birthDeath'
        parfor gd=1:Kce
            localExts(gd) = BD_solveOnce_ie_op(x0, t0, tf, k, event_type);
        end
        
    case 'SIRS'
        parfor gd=1:Kce
            localExts(gd) = SIRS_solveOnce_ie_op(x0, t0, tf, k, event_type);
        end        
end