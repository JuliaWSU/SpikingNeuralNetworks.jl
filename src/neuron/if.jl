
abstract type AbstractIFParameter end
struct IFParameter# <: AbstractIFParameter
    τm::Float32# = 20ms
    τe::Float32# = 5ms
    τi::Float32# = 10ms
    Vt::Float32# = -50mV
    Vr::Float32# = -60mV
    El::Float32# = Vr
    function IFParameter()
        new(20.0,5.0,10.0,-50.0,-60.0,-35.0)
    end

end

abstract type AbstractIF end

struct IF <: AbstractIF
    param::IFParameter
    N::Int32
    v::Vector{Float32} 
    ge::Vector{Float32}
    gi::Vector{Float32}
    fire::Vector{Bool}
    I::Vector{Float32}
    records::Dict
    function IF(;N::Int64,param::SpikingNeuralNetworks.IFParameter)
        VFT = Vector{Float32}
        v::VFT = param.Vr .+ rand(N) .* (param.Vt - param.Vr)
        ge::VFT = zeros(N)
        gi::VFT = zeros(N)
        fire::Vector{Bool} = zeros(Bool, N)
        I::VFT = zeros(N)
        records::Dict = Dict()
        new(param,N,v,ge,gi,fire,I,records)
    end
end


spike!(neuron::IF, t::Integer; dt::Real = 1.0) = (t > 0) && (neuron.voltage = neuron.vreset)
function spike!(neurons::T, spikes; dt::Real = 1.0) where T<:AbstractArray{<:IF}
    map!((x, y, s) -> (s > 0) ? y : x, neurons.v, neurons.v, neurons.Vr, spikes)
end
"""
    [Integrate-And-Fire Neuron](https://neuronaldynamics.epfl.ch/online/Ch1.S3.html)
"""
IF

#function lif!(t::CuVector{<:Real}, I::CuVector{<:Real}, V::CuVector{<:Real}; vrest::CuVector{<:Real}, R::CuVector{<:Real}, tau::CuVector{<:Real})
function integrate!(p::IF, param::IFParameter, dt::Float32)
    @unpack N, v, ge, gi, fire, I = p
    τm, τe, τi, Vt, Vr, El = (100.0,5.0,10.0,0.1,0.0,-35.0)

    # account for leakage
    @. V = V * exp(-dt / τm) + Vr

    # apply update step
    @. V += dt * (ge[i] + gi[i] -I[i]) * (R / τm)
    @. ge += dt * -ge / τe
    @. gi += dt * -gi / τi
    @. fire = v > Vt
end
#=
function integrate!(p::IF, param::IFParameter, dt::Float32)
    @unpack N, v, ge, gi, fire, I = p
    τm, τe, τi, Vt, Vr, El = (100.0,5.0,10.0,0.1,0.0,-35.0)

    #τm = 100
    #vreset = 0.0
    #vth = 0.1
    R = 1.75
    #@unpack τm, τe, τi, Vt, Vr, El = param
    @inbounds for i = 1:N

        v[i] += v[i] * exp(-dt / τm) + Vr
        v[i] += dt * (ge[i] + gi[i] -I[i]) * (R / τm)

        #v[i] += dt * (ge[i] + gi[i] - (v[i] - El) + I[i]) / τm
        ge[i] += dt * -ge[i] / τe
        gi[i] += dt * -gi[i] / τi
        fire[i] = v[i] > Vt
    end
    @inbounds for i = 1:N
        v[i] = ifelse(fire[i], Vr, v[i])
    end
end
=#
#IF(τm::Real, vreset::Real, R::Real = 1.0) = IF{Float32, Int}(vreset, 0, 0, τm, vreset, R)

isactive(neuron::IF, t::Integer; dt::Real = 1.0) = (neuron.v > 0)
getvoltage(neuron::IF) = neuron.v