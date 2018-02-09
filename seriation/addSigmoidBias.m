function [fBias, gBias] = addSigmoidBias(funh, x, w, mu, lambda)
% given function handle funh, add mu*sigmoid(lambda*w'*x) to the function

[f, g] = funh(x);
% wnorm = norm(w);
% xnorm = norm(x);
sigm = 1 / (1 + exp(-lambda * w'*x));
fBias = f - mu * sigm;
gBias = g - (mu * lambda * sigm^2) * w;

end