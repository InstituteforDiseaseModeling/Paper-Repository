function  [num_out, denom_out, counter_out] = SIRS_solveOnce_gamma(x0, t0, tf, k, ies)

ielen        = length(ies);
counter_out  = zeros(1,ielen);
num_out      = zeros(ielen, 3);
denom_out    = zeros(ielen, 3);
stoch_matrix = [-1 1 0; 0 -1 1; 1 0 -1];

x      = x0;
t      = t0;
a      = zeros(1,3);
n      = zeros(1,3);
lambda = zeros(1,3);

% ies has the rarest one in the first element
ie_ind  = ielen;
while(t < tf)
    % check for a rare event
    if (x(2) == ies(ie_ind))
        counter_out(ie_ind) = 1;
        
        num_out(ie_ind,:)   = n;
        denom_out(ie_ind,:) = lambda;
        
        if ie_ind == 1
            break;
        else
            ie_ind = ie_ind - 1;
        end
    end
    
    % compute propensity
    a    = k.*x;
    a(1) = a(1) * x(2);
    a0   = sum(a);
    
    % time to the next reaction
    tau = log(1/rand) / a0;
    t   = t + tau;
    if t < tf
        j = find(cumsum(a)/a0 > rand);
        j = j(1);
        x      = x + stoch_matrix(j,:);
        n(j)   = n(j) + 1;
        lambda = lambda + a*tau;
    end
end
