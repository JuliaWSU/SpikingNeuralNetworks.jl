
abstract type AbstractIFParameter end

struct IFParameter <: AbstractIFParameter 
    τm::Float32 
    τe::Float32 
    τi::Float32 
    Vt::Float32 
    Vr::Float32 
    El::Float32 
end

function IFParameter()
    return IFParameter(20,5.0,10.0,-50.0,-60.0,-49)
end    
function IFParameter(;τm::Float32=20ms,τe::Float32=5.0,τi::Float32=10.0,Vt::Float32=-50.0,Vr::Float32=-60.0,El::Float32=-49mV)
    return IFParameter(τm,τe,τi,Vt,Vr,El)
end
abstract type AbstractIF end

@snn_kw mutable struct IF{VFT=Vector{Float32},VBT=Vector{Bool}} <: AbstractIF
    param::IFParameter = IFParameter()
    N::Int32 = 10
    v::Vector{Float32} = param.Vr .+ rand(N) .* (param.Vt - param.Vr)
    ge::Vector{Float32} = zeros(N)
    gi::Vector{Float32} = zeros(N)
    fire::VBT = zeros(Bool, N)
    I::Vector{Float32} = zeros(N)
    records::Dict = Dict()
end

"""
    [Integrate-And-Fire Neuron](https://neuronaldynamics.epfl.ch/online/Ch1.S3.html)
"""
IF

#integrate!(::SpikingNeuralNetworks.IF{Vector{Float32}, Vector{Bool}}, ::SpikingNeuralNetworks.IFParameter, ::Float16)

function integrate!(p::IF, param::IFParameter, dt::Float64)
    """
    Integrating over cells, not integrating over time. 
    An ideal place for parallelisation.
    """

    @unpack N, v, ge, gi, fire, I = p
    @unpack τm, τe, τi, Vt, Vr, El = param
    @inbounds for i = 1:N
        v[i] += dt * (ge[i] + gi[i] - (v[i] - El) + I[i]) / τm
        ge[i] += dt * -ge[i] / τe
        gi[i] += dt * -gi[i] / τi
    end
    @inbounds for i = 1:N
        fire[i] = v[i] > Vt
        v[i] = ifelse(fire[i], Vr, v[i])
    end
end
#integrate!(p::IF, param::IFParameter, dt::Float64) = integrate!(p::IF, param::IFParameter, dt::Float32)
function integrate!(p::IF, param::IFParameter, dt::Float32)
    @unpack N, v, ge, gi, fire, I = p
    @unpack τm, τe, τi, Vt, Vr, El = param
    @inbounds for i = 1:N
        v[i] += dt * (ge[i] + gi[i] - (v[i] - El) + I[i]) / τm
        ge[i] += dt * -ge[i] / τe
        gi[i] += dt * -gi[i] / τi
    end
    @inbounds for i = 1:N
        fire[i] = v[i] > Vt
        v[i] = ifelse(fire[i], Vr, v[i])
    end
end
