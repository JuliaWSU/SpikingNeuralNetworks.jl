using SparseArrays
using SpikingNeuralNetworks
using Test

w = sprand(10,10,0.25)
rowptr, colptr, I, _, _, WW = SNN.dsparse(w)
g =   Vector{Bool}([true,false,false,false,false,true,false,true,false,false])#,0,0,1};
@show(g)

@show(colptr)

#Test.@test At.nzval == A.nzval[index]
    
#g = Array(Boolean)
#sprand(10,10,0.5)


#=
@testset "Constructors" begin
    @testset "Type parameters" begin
        for Model in (SNN.HH, SNN.IF, SNN.IF2, SNN.IZ, SNN.NoisyIF, SNN.Poisson, SNN.Rate)
            test_typeparams(Model)
        end
        test_typeparams(SNN.RateSynapse; args=(SNN.Rate(),SNN.Rate()))
        test_typeparams(SNN.PINningSynapse; args=(SNN.Rate(),SNN.Rate()))
        test_typeparams(SNN.FLSynapse; args=(SNN.Rate(),SNN.Rate()))
        test_typeparams(SNN.SpikingSynapse; args=(SNN.IF(),SNN.IF(),:ge))
    end # Type parameters
    end # Constructors
=#