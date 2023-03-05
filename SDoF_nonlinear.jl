using PyPlot;

# simulation interval
dt = 1.0e-4;
t_span = [0., 5.];
t = t_span[1]:dt:t_span[2]

# define damping constants and time relation
m = 2; Dₐ = 0.004; kᵧ = 60
damp(x) = Dₐ*(1 - (x/kᵧ)^m);

ζ = damp.(t);

# define frequenci constants and time relation
kf = -0.7; shift = 50;
freq(x) = exp(kf*x) + shift;

f = freq.(t);

ω = 2π*f;
ϕ = cumsum(ω)*dt;

A = 100 .* exp.(-ζ .* ω .* t);

y = A .* sin.(ϕ)

plot(t,y, linewidth = 1.0); plot(t,A, "--k", linewidth = 2.0 )
