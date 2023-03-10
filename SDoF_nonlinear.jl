using PyPlot;


function diff14(vec::AbstractVector{<: Number}, time_step::Number)
    der = similar(vec);
    s_coefs = [1/12, -2/3, 2/3, -1/12];
    f_coefs = [-25/12, 4, -3, 4/3, -1/4];
    b_coefs = [1/4, -4/3, 3, -4, 25/12];
    N = size(vec, 1);
    M = size(vec, 2);
    # almost all points calculates using midpoint scheme
    for i = 3:(length(vec)-2)
        der[i, :] = (vec[i-2,:].*s_coefs[1] + vec[i-1,:].*s_coefs[2] + vec[i+1,:].*s_coefs[3] + vec[i+2,:].*s_coefs[4]) /time_step;
    end
    # first points are calculated using forward scheme of the same order
    for i = 1:2
        der[i,:] = sum(hcat([vec[i+j-1,:].*f_coefs[j] for j in eachindex(f_coefs)]...), dims =2)./time_step;
    end
    # last points are calculated using bacckward scheme of the same order
    for j in 1:M
        for i = N-1:N
            der[i,j] = (vec[i-4:i, j]'*b_coefs)/time_step;
        end
    end
    return der;
end

function frequency(z, dt)
    y = real(z);
    ỹ = imag(z);
    ẏ = diff14(y, dt);
    ỹ̇ = diff14(ỹ, dt);

    A_sq = [(y[i]^2 + ỹ[i]^2) for i in eachindex(z)];

    ω = similar(y);
    for i in eachindex(ω)
        ω[i] = (y[i]*ỹ̇[i] - ẏ[i]*ỹ[i])/A_sq[i];
    end
    return ω
end

function undampedfrequency1(ω, A, dt)
    Ȧ=diff14(A,dt);
    Ä=diff14(Ȧ,dt);
    dω=diff14(ω,dt);
    ω₀ = similar(ω);
    for i in eachindex(ω)
        ω₀[i]=sqrt(ω[i]^2-(Ä[i]/A[i])+(2Ȧ[i]^2)/(A[i]^2)+(Ȧ[i]*dω[i] )/(A[i]*ω[i]))
    end
    return ω₀
end

function undampedfrequency2(ω, A, dt)
    Ȧ=diff14(A,dt);
    Ä=diff14(Ȧ,dt);
    ω₀ = similar(ω);
    for i in eachindex(ω);
    ω₀[i]=sqrt(ω[i]^2-(Ä[i]/A[i]));
    end
    return ω₀
end

function instdamping(ω, A, dt)
    Ȧ=diff14(A,dt);
    dω=diff14(ω,dt);
    β = similar(ω);
    for i in eachindex(ω);
        β[i]=-Ȧ[i]/A[i]-(dω[i])/(2*ω[i]);
    end
    return β
end

function decrement(ω, A, dt)
    Ȧ=diff14(A,dt);
    Ä=diff14(Ȧ,dt);
    dω=diff14(ω,dt);
    δ = similar(ω);
    for i in eachindex(ω);
        δ[i]=-(2*π*Ȧ[i]*ω[i])/(A[i]*(ω[i]^2-(Ä[i]/A[i])))-(π*dω[i])/(ω[i]^2-(Ä[i]/A[i]));
    end
    return δ
end


function damping(env::AbstractVector{<:Number}, time::AbstractVector{<:Number}, step::Number)
    Ȧ = diff14(env, dt);
    g = - Ȧ ./ env;
    βₐₚ = Vector{Float64}(undef, length(time));
    for i in 2:length(time);
        βₐₚ[i] = (step * time[i] / ( step + time[i])) * ( βₐₚ[i-1] / step + g[i] / time[i]);
    end
    return βₐₚ
 end

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
ϕ2 = cumsum(ω.*ζ)*dt;

######A = 100 .* exp.(-ζ .* ω .* t);
######A2 = 100 .* exp.(-ζ .* ϕ);
A=100 .* exp.(-ϕ2);
plot(t,A)
plot(t,A2) 
plot(t,A3)

y = A .* sin.(ϕ);
ỹ = A .* cos.(ϕ); # to obtain imagine part of the signal we compute it directly to avoid using of HT 

Ã = [sqrt(y[i]^2 + ỹ[i]^2) for i in eachindex(y)];
plot(t,y);plot(t,ỹ);plot(t, A); plot(t, Ã, "k--")
ẏ = diff14(y, dt);


plot(t,y, linewidth = 1.0); plot(t,A, "--k", linewidth = 2.0 )
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


full = y + ỹ * im

ŵ = frequency(full, dt);
plot(t,ω);plot(t, -ŵ, "--")


ω₀₁=undampedfrequency1(ω,A,dt);
ω₀₂=undampedfrequency2(ω,A,dt);

plot(t, ω₀₁,"-"); plot(t, ω₀₂,"--"); plot(t, -ŵ,"--");
plot(t, ω₀₁ - ω₀₂)
plot(t, ω₀₁ - -ŵ)
ylim([0,0.02])


β = instdamping(ω, A, dt);
plot(t, β);plot(t, ζ.*ω)
δ = decrement(ω, A, dt);

subplot(1,2,1)
plot(t, ζ); plot(t, δ/2π);
yscale("log");

subplot(1,2,2)
plot(A, ζ); plot(A, δ/2π)
yscale("log");xscale("log")