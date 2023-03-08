#valgrind --smc-check=all-non-file julia if_net.jl

#@snn_kw 
abstract type AbstractIFParameter end

struct IFParameter <: AbstractIFParameter #<: AbstractIFParameter
    τm::Float32 # = 20#ms
    τe::Float32 # = 5#ms
    τi::Float32 # = 10#ms
    Vt::Float32 # = -50#mV
    Vr::Float32 # = -60#mV
    El::Float32 # = -65
end

    #function IFParameter{Vector{Float32}}(τm::Float32,τe::Float32,τi::Float32,Vt::Float32,Vr::Float32,El::Float32)# where T<:Integer
    #    new(τm, τe, τi, Vt, Vr, El)
    #end
    #function IFParameter()
        #(τm::Float32,τe::Float32,τi::Float32,Vt::Float32,Vr::Float32,El::Float32)# where T<:Integer
    #    new(20.0, 5.0, 10.0, -50.0, -60.0, -49.0)
    #end

    #IFParameter{T<:Float32}  = new(20.0,5.0,10.0,-50.0,-60.0,-65.0)
function IFParameter()#;τm::Float32,τe::Float32,τi::Float32,Vt::Float32,Vr::Float32,El::Float32)#20.0,5.0,10.0,-50.0,-60.0,-65.0)
    FT = Float32
    τm::FT = 20ms
    τe::FT = 5ms
    τi::FT = 10ms
    Vt::FT = -50mV
    Vr::FT = -60mV
    El::FT = -65.0
    return IFParameter(20ms,5.0,10.0,-50.0,-60.0,-65.0)
end    


#IFParameter = IFParameter()

#=
function IFParameter()
    FT=Float32
    τm::FT = 20ms
    τe::FT = 5ms
    τi::FT = 10ms
    Vt::FT = -50mV
    Vr::FT = -60mV
    El::FT = Vr
end    
IFParameter = IFParameter()
=#
abstract type AbstractIF end

@snn_kw mutable struct IF{VFT=Vector{Float32},VBT=Vector{Bool}} <: AbstractIF
    param::IFParameter = IFParameter()
    #@show(param)
    N::Int32 = 100
    v::VFT = ones(N).*param.Vr #.+ rand(N) .* (param.Vt - param.Vr)
    ge::VFT = zeros(N)
    gi::VFT = zeros(N)
    fire::VBT = zeros(Bool, N)
    I::VFT = zeros(N)
    records::Dict = Dict()
end

"""
    [Integrate-And-Fire Neuron](https://neuronaldynamics.epfl.ch/online/Ch1.S3.html)
"""
IF

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