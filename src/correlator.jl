function current_current_correlator(model, dτ, betas; parallel=false, mthreads=Threads.nthreads(), kwargs...)
    mcs = Vector(undef, length(betas))
    if !parallel
        for i in eachindex(betas)
            mc = DQMC(model, beta=betas[i], delta_tau=dτ; kwargs...)
            mc[:ccc] = MonteCarlo.Measurement(
                mc, model, GreensAt((betas[i]/2)/dτ, 0), EachBondPairByBravaisDistance(mc), MonteCarlo.FlavorIterator(mc, 2), MonteCarlo.cc_kernel
                )
            run!(mc)
            mcs[i] = mc
        end
    else
        idx = Threads.Atomic{Int}(1)
        n = length(betas)
        Threads.@sync for _ in 1:mthreads
            Threads.@spawn while true
                i = Threads.atomic_add!(idx, 1) 
                i > n && break  
                mc = DQMC(model, beta=betas[i], delta_tau=dτ; kwargs...)
                mc[:ccc] = MonteCarlo.Measurement(
                    mc, model, GreensAt(round((betas[i]/2)/dτ), 0), EachBondPairByBravaisDistance(mc), MonteCarlo.FlavorIterator(mc, 2), MonteCarlo.cc_kernel
                    )
                run!(mc)
                mcs[i] = mc
            end
        end
    end
    return mcs
end