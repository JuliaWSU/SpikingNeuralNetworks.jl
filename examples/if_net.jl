
using SpikingNeuralNetworks
SNN.@load_units
using SparseArrays
import LinearAlgebra.normalize!
using OnlineStats
using Plots
using JLD2
using UnicodePlots
#unicodeplots()

function makeNetGetTimes()
    scale = 1/10.0#1.0/200.0
    if !isfile("potjans_full_scale.jld2")

        @time (Lee,Lie,Lei,Lii),Ne,Ni = SNN.potjans_layer(scale)
        @save "potjans_full_scale.jld2" Lee Lie Lei Lii Ne Ni

    else
        @load "potjans_full_scale.jld2" Lee Lie Lei Lii Ne Ni
    end
    @time (NoisyInputSynInh,NoisyInputSyn,LeeSyn,LeiSyn,LiiSyn,LieSyn,E,I,Noisy) = SNN.SpikingSynapse(Lee,Lei,Lii,Lie)
    print("wiring done")    
    P = [E,I,Noisy] # populations     
    C = [NoisyInputSynInh,NoisyInputSyn,LeeSyn,LeiSyn,LiiSyn,LieSyn] # connections
    SNN.monitor([E,I], [:fire])
    @time SNN.sim!(P, C; duration = 0.2second)
    print("simulation done !")
    (times,nodes) = SNN.get_trains([E,I])
    @time o1 = HeatMap(zip(minimum(times):maximum(times)/100.0:maximum(times),minimum(nodes):maximum(nodes/100.0):maximum(nodes)) )
    @time fit!(o1,zip(times,convert(Vector{Float64},nodes)))
    plot(o1, marginals=false, legend=true) #|>display 
    Plots.savefig("default_heatmap.png")
    return (P,C,times,nodes)
end

(P,C,times,nodes) = makeNetGetTimes()
