function localExt  = BD_solveOnce_ie_op(x0, t0, tf, k, event_type)

x        = x0;
t        = t0;
localExt = x;

while(t < tf)
    % compute propensity
    b  = [k(1) k(2)*x];
    b0 = b(1) + b(2);
    
    % time to the next reaction
    tau = log(1/rand) / b0;
    t   = t + tau;
    
    if(t <= tf)
        if (b(1)/b0) > rand
            x = x+1;
        else
            x = x-1;
        end
        
        currentEevntVal = x;
        % record extreme values
        if (event_type == 1)
            localExt = max(localExt, currentEevntVal);
        else
            localExt = min(localExt, currentEevntVal);
        end
    end
end