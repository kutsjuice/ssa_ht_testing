using PyPlot;

# simulation interval
dt = 1.0e-4;
t_span = [0., 5.];
t = t_span[1]:dt:t_span[2]

# define damping constants and time relation
m = 0.2; Dₐ = 0.004; kᵧ = 0.8
damp(x) = Dₐ*(exp(-x/kᵧ)+1);

ζ = damp.(t);

# define frequenci constants and time relation
kf = 1.7;  f_amp = 4.5; shift = 50 - f_amp;
freq(x) = shift + f_amp*(1 - exp(-kf*x));

f = freq.(t);

ω = 2π*f;
ϕ = cumsum(ω)*dt;

A = 100 .* exp.(-ζ .* ω .* t);

y = A .* sin.(ϕ)

# plot(t,y, linewidth = 1.0); plot(t,A, "--k", linewidth = 2.0 )
subplot(3,2,1)
plot(A, ζ)
# ylim([0.003,0.005])
xscale("log")
yscale("log")
subplot(3,2,2)
plot(t,ζ)
yscale("log")

subplot(3,2,3)
plot(A, f)
xscale("log")

subplot(3,2,4)
plot(t,f)

subplot(3,1,3)
plot(t,y)
gcf().set_size_inches(12,8)
tight_layout()
