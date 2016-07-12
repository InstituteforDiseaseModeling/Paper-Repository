function ks = generateKs(kslen, k1min, k1max, k2min, k2max, k3min, k3max)

if nargin == 5
    ks = [k1min + (k1max-k1min)*rand(kslen,1) k2min + (k2max-k2min)*rand(kslen,1)];
elseif nargin == 7
    ks = [k1min + (k1max-k1min)*rand(kslen,1), k2min + (k2max-k2min)*rand(kslen,1), k3min + (k3max-k3min)*rand(kslen,1)];
else
    error('unknown input type')
end


