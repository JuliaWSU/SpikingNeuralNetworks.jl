                    EE = SNN.IFNF(;N = length(cumvalues), param = SNN.IFParameter())
                    synEE = make_proj(Lee[v,v1],EE)
                elseif syn1==0# meaning if the same as a logic: Inhibitory post synapse  is true                   
                    EI = SNN.IFNF(;N = length(cumvalues), param = SNN.IFParameter())
                    synEI = make_proj(Lei[v,v1],EI)
                end
            elseif syn0==0         
                if syn1==1 
                    IE = SNN.IFNF(;N = length(cumvalues), param = SNN.IFParameter())
                    synIE = make_proj(Lie[v,v1],IE)
                elseif syn1==0
                    II = SNN.IFNF(;N = length(cumvalues), param = SNN.IFParameter())
                    synII = make_proj(Lii[v,v1],II)
                end
            end 
        end
                    
    end
    (EE,EI,IE,II,synII,synIE,synEI,synEE)
end

function build_matrix(Ncells::Int32,just_iterator,g_strengths::Vector{Float64})
    #index_assignment!(::NTuple{4, Int64}, ::Vector{Float64}, ::SparseMatrixCSC{Float32, Int64}, ::Vector{Vector{Tuple{Int64, Int64}}}, ::Vector{Vector{Tuple{Int64, Int64}}}, ::Vector{Vector{Tuple{Int64, Int64}}}, ::Vector{Vector{Tuple{Int64, Int64}}})
    #Lee = spzeros(Boolean, (Ncells, Ncells))

    Lee = Vector{Vector{Tuple{Int64, Int64}}}[]
    Lie = Vector{Vector{Tuple{Int64, Int64}}}[]
    Lii = Vector{Vector{Tuple{Int64, Int64}}}[]
    Lei = Vector{Vector{Tuple{Int64, Int64}}}[]
    @showprogress for i in just_iterator
        index_assignment!(i[:],g_strengths,Lxx)#,Lee,Lie,Lii,Lei)
    end
    Lee = Lxx[(i[1],i[2]) for i in Lee]
    Lie = Lxx[Lie[1,:],Lie[2,:]]
    Lii = Lxx[Lii[1,:],Lii[2,:]]
    Lei = Lxx[Lei[1,:],Lei[2,:]]

    #Iterators.map(f, iterators...)
   # map(index_assignment!(item,w0Weights,g_strengths,Lee,Lie,Lii,Lei,just_iterator)
    #output = map(x -> samplesmallGram(L), 1:1:10)
    #map(samplesmallGram(L), just_iterator)
    Lexc = Lee+Lei
    Linh = Lie+Lii

    #map!(index_assignment!, item for iter_item)
    @assert maximum(Lexc[:])>=0.0
    @assert maximum(Linh[:])<=0.0


    #@show(just_iterator)
    #ploop(f, itr,w0Weights,g_strengths,Lee,Lie,Lii,Lei)
    #ploop(index_assignment!,just_iterator,w0Weights,g_strengths,Lee,Lie,Lii,Lei)
    
    ## this works
    ##
    Lee = Lxx[(i[1],i[2]) for i in Lee]
    Lie = Lxx[Lie[1,:],Lie[2,:]]
    Lii = Lxx[Lii[1,:],Lii[2,:]]
    Lei = Lxx[Lei[1,:],Lei[2,:]]

    #Iterators.map(f, iterators...)
    #map(index_assignment!(item,w0Weights,g_strengths,Lee,Lie,Lii,Lei,just_iterator)
    #output = map(x -> samplesmallGram(L), 1:1:10)
    #map(samplesmallGram(L), just_iterator)
    Lexc = Lee+Lei
    Linh = Lie+Lii

    #map!(index_assignment!, item for iter_item)
    @assert maximum(Lexc[:])>=0.0
    @assert maximum(Linh[:])<=0.0

    #Lee_ = MArray{Tuple{Ncells,Ncells},Float32}(Lee)
    #Lee_ = {Tuple{Ncells,Ncells},Float32}(Lee) #,2,9}

    return Lee,Lie,Lei,Lii
end


"""
Build the matrix from the Potjans parameterpotjans_layers.
"""
function potjans_weights(args)
    Ncells, g_strengths, ccu, scale = args
    #(;Ncells,g_strengths,ccu,scale)
    (cumulative,ccu,layer_names,conn_probs,syn_pol) = potjans_params(ccu,scale)    
    cumvalues = values(cumulative)
    #cumvalues = convert(Vector{Vector{Float32}},cumvalues)

    #just_iterator = []
    #Lxx = spzeros(Float32, (Ncells, Ncells))
    Lee = spzeros(Float32, (Ncells, Ncells))
    Lie = spzeros(Float32, (Ncells, Ncells))
    Lei = spzeros(Float32, (Ncells, Ncells))
    Lii = spzeros(Float32, (Ncells, Ncells))

    #rv = spzeros(Float32, (Ncells, Ncells))

    build_matrix_prot!(Lee,Lie,Lei,Lii,cumvalues,conn_probs,Ncells,syn_pol,g_strengths)
    (EE,EI,IE,II,synII,synIE,synEI,synEE) = build_neurons_connections(Lee,Lei,Lie,Lii,cumvalues, Ncells,syn_pol)
    #Lee,Lie,Lei,Lii = build_matrix(Ncells,just_iterator,g_strengths)
    #Lee,Lie,Lei,Lii
end


function auxil_potjans_param(scale=1.0::Float64)
	ccu = Dict{String, Int32}("23E"=>20683,
		    "4E"=>21915, 
		    "5E"=>4850, 
		    "6E"=>14395, 
		    "6I"=>2948, 
		    "23I"=>5834,
		    "5I"=>1065,
		    "4I"=>5479)
	ccu = Dict{String, Int32}((k,ceil(Int64,v*scale)) for (k,v) in pairs(ccu))
	Ncells = Int32(sum([i for i in values(ccu)])+1)
	Ne = Int32(sum([ccu["23E"],ccu["4E"],ccu["5E"],ccu["6E"]]))
    Ni = Int32(Ncells - Ne)
    Ncells, Ne, Ni, ccu

end

function potjans_layer(scale)
    
    Ncells,Ne,Ni, ccu = auxil_potjans_param(scale)    
    pree = 0.1
    K = round(Int, Ne*pree)
    sqrtK = sqrt(K)
    g = 1.0
    tau_meme = 10   # (ms)
    je = 2.0 / sqrtK * tau_meme * g
    ji = 2.0 / sqrtK * tau_meme * g 
    jee = 0.15je 
    jei = je 
    jie = -0.75ji 
    jii = -ji
    g_strengths = Vector{Float64}([jee,jie,jei,jii])

    genStaticWeights_args = (;Ncells,g_strengths,ccu,scale)
    potjans_weights(genStaticWeights_args),Ne,Ni
end

#=
SGet in-degrees for each connection for the full-scale (1 mm^2) model
function get_indegrees()
    K = np.zeros([n_layers * n_pops_per_layer, n_layers * n_pops_per_layer])
    for target_layer in layers:
        for target_pop in pops:
            for source_layer in layers:
                for source_pop in pops:
                    target_index = structure[target_layer][target_pop]
                    source_index = structure[source_layer][source_pop]
                    n_target = N_full[target_layer][target_pop]
                    n_source = N_full[source_layer][source_pop]
                    K[target_index][source_index] = round(np.log(1. -
                        conn_probs[target_index][source_index]) / np.log(
                        (n_target * n_source - 1.) / (n_target * n_source))) / n_target
                end
            end
        end
    end
    return K
   
end
=#


#=
"""
Adjust synaptic weights and external drive to the in-degrees
to preserve mean and variance of inputs in the diffusion approximation
function adjust_w_and_ext_to_K(K_full, K_scaling, w, DC)
    K_ext_new = Dict()
    I_ext = Dict()
    for target_layer in layers:
        K_ext_new[target_layer] = {}
        I_ext[target_layer] = {}
        for target_pop in pops:
            target_index = structure[target_layer][target_pop]
            x1 = 0
            for source_layer in layers:
                for source_pop in pops:
                source_index = structure[source_layer][source_pop]
                x1 += w[target_index][source_index] * K_full[target_index][source_index] * \
                        full_mean_rates[source_layer][source_pop]
                end
            end
            if input_type == 'poisson'
                x1 += w_ext*K_ext[target_layer][target_pop]*bg_rate
                K_ext_new[target_layer][target_pop] = K_ext[target_layer][target_pop]*K_scaling
            end
            I_ext[target_layer][target_pop] = 0.001 * neuron_params['tau_syn_E'] * \
                (1. - np.sqrt(K_scaling)) * x1 + DC[target_layer][target_pop]
            w_new = w / np.sqrt(K_scaling)
            w_ext_new = w_ext / np.sqrt(K_scaling)
        end
    end
    return w_new, w_ext_new, K_ext_new, I_ext
end
"""
=#
#=
# Create cortical populations
self.pops = {}
layer_structures = {}
total_cells = 0 
x_dim_scaled = x_dimension * math.sqrt(N_scaling)
z_dim_scaled = z_dimension * math.sqrt(N_scaling)
default_cell_radius = 10 # for visualisation 
default_input_radius = 5 # for visualisation 
for layer in layers:
    self.pops[layer] = {}
    for pop in pops:
        
        y_offset = 0
        if layer == 'L6': y_offset = layer_thicknesses['L6']/2
        if layer == 'L5': y_offset = layer_thicknesses['L6']+layer_thicknesses['L5']/2
        if layer == 'L4': y_offset = layer_thicknesses['L6']+layer_thicknesses['L5']+layer_thicknesses['L4']/2
        if layer == 'L23': y_offset = layer_thicknesses['L6']+layer_thicknesses['L5']+layer_thicknesses['L4']+layer_thicknesses['L23']/2
        
        layer_volume = Cuboid(x_dim_scaled,layer_thicknesses[layer],z_dim_scaled)
        layer_structures[layer] = RandomStructure(layer_volume, origin=(0,y_offset,0))
https://github.com/OpenSourceBrain/PotjansDiesmann2014/blob/master/PyNN/scaling.py
=#
