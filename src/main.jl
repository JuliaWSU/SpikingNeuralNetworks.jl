

function sim!(P, C; dt = 0.1ms, duration = 10ms)
    for t = 0ms:dt:(duration - dt)
        sim!(P, C, dt)
    end
end

function sim16!(P, C, dt::Float16)
    for p in P
        integrate!(p, p.param, Float16(dt))
        record!(p)
    end
    for c in C
        forward!(c, c.param)
        record!(c)
    end
end

function sim!(P, C, dt::Float64)
    for p in P
        #integrate!(::SpikingNeuralNetworks.IF{Vector{Float32}, Vector{Bool}}, ::SpikingNeuralNetworks.IFParameter, ::Float16)

        integrate!(p, p.param, dt)
        record!(p)
    end
    for c in C
        forward!(c, c.param)
        record!(c)
    end
end

function sim16!(P, C; dt::Float64 = 0.1ms, duration::Float64 = 10ms)
    for t = 0ms:dt:(duration - dt)
        sim!(P, C, Float16(dt))
    end
end



function train16!(P, C, dt::Float16, t = 0.0)
    for p in P
        integrate!(p, p.param, Float16(dt))
        record!(p)
    end
    for c in C
        forward!(c, c.param)
        plasticity!(c, c.param, Float16(dt), Float32(t))
        record!(c)
    end
end

function train16!(P, C; dt = 0.1ms, duration = 10ms)
    for t = 0ms:dt:(duration - dt)
        train!(P, C, Float16(dt), t)
    end
end


function train!(P, C, dt, t = 0)
    for p in P
        integrate!(p, p.param, Float32(dt))
        record!(p)
    end
    for c in C
        forward!(c, c.param)
        plasticity!(c, c.param, Float32(dt), Float32(t))
        record!(c)
    end
end

function train!(P, C; dt = 0.1ms, duration = 10ms)
    for t = 0ms:dt:(duration - dt)
        train!(P, C, Float32(dt), Float32(t))
    end
end
