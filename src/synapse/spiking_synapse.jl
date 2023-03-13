#using Revise
#using SpikingNeuralNetworks
@snn_kw struct SpikingSynapseParameter16{FT=Float16}
    τpre::FT = 20ms
    τpost::FT = 20ms
    Wmax::FT = 0.01
    ΔApre::FT = 0.01 * Wmax
    ΔApost::FT = -ΔApre * τpre / τpost * 1.05
end

@snn_kw struct SpikingSynapseParameter{FT=Float32}
    τpre::FT = 20ms
    τpost::FT = 20ms
    Wmax::FT = 0.01
    ΔApre::FT = 0.01 * Wmax
    ΔApost::FT = -ΔApre * τpre / τpost * 1.05
end

@snn_kw mutable struct SpikingSynapse{VIT=Vector{Int32},VFT=Vector{Float32},VBT=Vector{Bool}}
    param::SpikingSynapseParameter = SpikingSynapseParameter()
    rowptr::VIT # row pointer of sparse W
    colptr::VIT # column pointer of sparse W
    I::VIT      # postsynaptic index of W
    J::VIT      # presynaptic index of W
    index::VIT  # index mapping: W[index[i]] = Wt[i], Wt = sparse(dense(W)')
    W::VFT  # synaptic weight
    tpre::VFT = zero(W) # presynaptic spiking time
    tpost::VFT = zero(W) # postsynaptic spiking time
    Apre::VFT = zero(W) # presynaptic trace
    Apost::VFT = zero(W) # postsynaptic trace
    fireI::VBT # postsynaptic firing
    fireJ::VBT # presynaptic firing
    g::VFT # postsynaptic conductance
    records::Dict = Dict()
end

"""
[Spking Synapse](https://brian2.readthedocs.io/en/2.0b4/resources/tutorials/2-intro-to-brian-synapses.html)
"""
SpikingSynapse

function SpikingSynapse(pre, post, sym; σ = 0.0, p = 0.0, kwargs...)
    w = σ * sprand(post.N, pre.N, p)
    #spy(w)
    #savefig("random_wiring.png")
    rowptr, colptr, I, J, index, W = dsparse(w)
    fireI, fireJ = post.fire, pre.fire
    g = getfield(post, sym)
    SpikingSynapse(;@symdict(rowptr, colptr, I, J, index, W, fireI, fireJ, g)..., kwargs...)
end


function SpikingSynapse(Lee::SparseMatrixCSC{Float32, Int64},Lei::SparseMatrixCSC{Float32, Int64},Lii::SparseMatrixCSC{Float32, Int64},Lie::SparseMatrixCSC{Float32, Int64}, kwargs...)
    inhib = Lii+Lie
    cnt=0
    for i in eachrow(inhib)
        if sum(i[:]) != 0.0
            cnt+=1
        end
    end

    exc = Lee+Lei
    cnt=0
    for i in eachrow(exc)
        if sum(i[:]) != 0.0
            cnt+=1
        end
    end
    total = inhib+exc
    cnt = length(total)
    Ipop = SNN.IF16(;N = cnt, param = SNN.IFParameter16())
    Epop = SNN.IF16(;N = cnt, param = SNN.IFParameter16())

    @assert maximum(Lii) <= 0.0
    @assert maximum(Lei) >= 0.0
    @assert maximum(Lie) <= 0.0
    @assert maximum(Lee) >= 0.0
    Lee = 0.001*Lee

    rowptr, colptr, I, J, index, W = dsparse(Lee)
    fireI, fireJ = Epop.fire, Epop.fire
    g = getfield(Epop, :ge)
    LeeSyn = SpikingSynapse(;@symdict(rowptr, colptr, I, J, index, W, fireI, fireJ, g)..., kwargs...)

    rowptr, colptr, I, J, index, W = dsparse(Lei)
    fireI, fireJ = Epop.fire, Ipop.fire
    g = getfield(Epop, :ge)
    LeiSyn = SpikingSynapse(;@symdict(rowptr, colptr, I, J, index, W, fireI, fireJ, g)..., kwargs...)

    rowptr, colptr, I, J, index, W = dsparse(Lii)
    fireI, fireJ = Ipop.fire, Ipop.fire
    g = getfield(Ipop, :gi)
    LiiSyn = SpikingSynapse(;@symdict(rowptr, colptr, I, J, index, W, fireI, fireJ, g)..., kwargs...)

    rowptr, colptr, I, J, index, W = dsparse(Lie)
    fireI, fireJ = Ipop.fire, Epop.fire
    g = getfield(Ipop, :gi)
    LieSyn = SpikingSynapse(;@symdict(rowptr, colptr, I, J, index, W, fireI, fireJ, g)..., kwargs...)


    (LeeSyn,LeiSyn,LiiSyn,LieSyn,Epop,Ipop)
end


function forward!(c::SpikingSynapse, param::SpikingSynapseParameter)
    @unpack colptr, I, W, fireJ, g = c
    @inbounds for j in 1:(length(colptr) - 1)
        if fireJ[j]
            for s in colptr[j]:(colptr[j+1] - 1)
                g[I[s]] += W[s]
            end
        end
    end
end

function plasticity!(c::SpikingSynapse, param::SpikingSynapseParameter, dt::Float32, t::Float32)
    @unpack rowptr, colptr, I, J, index, W, tpre, tpost, Apre, Apost, fireI, fireJ, g = c
    @unpack τpre, τpost, Wmax, ΔApre, ΔApost = param
    @inbounds for j in 1:(length(colptr) - 1)
        if fireJ[j]
            for s in colptr[j]:(colptr[j+1] - 1)
                Apre[s] *= exp32(- (t - tpre[s]) / τpre)
                Apost[s] *= exp32(- (t - tpost[s]) / τpost)
                Apre[s] += ΔApre
                tpre[s] = t
                W[s] = clamp(W[s] + Apost[s], 0f0, Wmax)
            end
        end
    end
    @inbounds for i in 1:(length(rowptr) - 1)
        if fireI[i]
            for st in rowptr[i]:(rowptr[i+1] - 1)
                s = index[st]
                Apre[s] *= exp32(- (t - tpre[s]) / τpre)
                Apost[s] *= exp32(- (t - tpost[s]) / τpost)
                Apost[s] += ΔApost
                tpost[s] = t
                W[s] = clamp(W[s] + Apre[s], 0f0, Wmax)
            end
        end
    end
end
#plasticity16!(c::SpikingSynapse, param::SpikingSynapseParameter, dt::Float16, t::Float32) = plasticity!(c::SpikingSynapse, param::SpikingSynapseParameter, dt::Float32, t::Float32)