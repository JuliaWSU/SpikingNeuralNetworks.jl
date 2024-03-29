function index_assignment!(item::NTuple{4, Int64}, w0Weights::SparseMatrixCSC{Float64, Int64}, g_strengths::Vector{Float64}, lee::SparseMatrixCSC{Float32, Int64}, lie::SparseMatrixCSC{Float32, Int64}, lii::SparseMatrixCSC{Float32, Int64}, lei::SparseMatrixCSC{Float32, Int64})
    # excitatory weights.
    (jee,_,jei,_) = g_strengths 
    # Relative inhibitory synaptic weight
    wig = -20*4.5
    (src,tgt,syn0,syn1) = item
    if syn0==1
        if syn1==1            
            setindex!(w0Weights, jee, src,tgt)    
            setindex!(lee,jee, src,tgt)
            @assert lee[src,tgt]>=0.0

        elseif syn1==0# meaning if the same as a logic: Inhibitory post synapse  is true                   
            setindex!(w0Weights, jei, src,tgt)
            setindex!(lei, jei, src,tgt)
            @assert lei[src,tgt]>=0.0

        end
    elseif syn0==0# meaning if the same as a logic: Inhibitory post synapse  is true   
        if syn1==1
            setindex!(w0Weights, wig, src,tgt)
            setindex!(lie, wig, src,tgt)
            @assert w0Weights[src,tgt]<=0.0
            @assert lie[src,tgt]<=0.0

        elseif syn1==0# eaning meaning if the same as a logic: if occursin("I",k1)      is true               
            @assert syn1==0
            setindex!(w0Weights, wig, src,tgt)
            setindex!(lii,wig, src,tgt)
            @assert w0Weights[src,tgt]<=0.0
            @assert syn1==0
            @assert lii[src,tgt]<=0.0

        end
    end
end
"""
This function contains iteration logic seperated from synapse selection logic for readability only.
Hi code reuse with pattern above, need to some how re-write as one method.
"""
function index_assignment!(item::NTuple{4, Int64},w0Weights::SparseMatrixCSC{Float64, Int64}, lee::SparseMatrixCSC{Float32, Int64}, lie::SparseMatrixCSC{Float32, Int64}, lii::SparseMatrixCSC{Float32, Int64}, lei::SparseMatrixCSC{Float32, Int64})
    (src,tgt,syn0,syn1) = item

    if syn0==1
        if syn1==1
            setindex!(lee, w0Weights[src,tgt], src,tgt)
            @assert lee[src,tgt]>=0.0
        elseif
            syn1==0# meaning if the same as a logic: occursin("I",k1) is true                   
            setindex!(lei, w0Weights[src,tgt], src,tgt)
            @assert lei[src,tgt]>=0.0
        end
    elseif syn0==0
        # meaning meaning if the same as a logic: elseif occursin("I",k) is true  
        if syn1==1 
            setindex!(lie, w0Weights[src,tgt], src,tgt)
            @assert lie[src,tgt]<=0.0
       elseif
            syn1==0# eaning meaning if the same as a logic: if occursin("I",k1)      is true               
            setindex!(lii, w0Weights[src,tgt], src,tgt)
            @assert lii[src,tgt]<=0.0
        end
    end  
    #      
end

