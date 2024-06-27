### PREPARE ###

using Pkg

Pkg.activate(joinpath(homedir(), ".julia", "environments", "ABC"))

# Uncomment for very first execution (to install packages; may take a while...):
#Pkg.add("Distributions")
#Pkg.add("DifferentialEquations")
#Pkg.add("Plots")
#Pkg.add("CairoMakie")
#Pkg.add("PairPlots")

# define kinetic model
include("KM_simple_Julia.jl")
output_model = (all_params) -> simple(all_params)

### EXPERIMENTAL DATA ###

using Distributions

# "correct" parameters - Db_OL, H_O3, Td_O3, k_SLR
kp_truth = [1.0E-9, 8.0E-4, 4.0E-6, 1.0E+2]

# uncertainty ranges for prior
kp_logrange = [[-11.0, -7.0],
    [-5.0, -3.0],
    [-7.0, -3.0],
    [1.0, 5.0]
]

# particle radius (cm), ozone concentration (cm-3)
env_sets_exp = [[0.00225, 6.36E+12], 
[0.001125, 6.36E+12], 
[0.00225, 5.0E+13], 
[0.008, 1E+13],
]

# simulate each experiment
exp_t = []
exp_y = []
for env in env_sets_exp

    # add systematic error
    env1 = rand(Normal(1.0, 0.08), 1)*env[1] 
    env2 = rand(Normal(1.0, 0.08), 1)*env[2]

    # run model
    inpt = [kp_truth; env1; env2]
    out_t_raw, out_y_raw = output_model(inpt)

    # save results
    out_y = []
    out_t = []
    for step in [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
        for i in size(out_t_raw)
            index_of_closest = findmin(abs.(out_y_raw .- step))[2]
            push!(out_y, out_y_raw[index_of_closest])
            push!(out_t, out_t_raw[index_of_closest])
        end
    end

    # add absolute random error
    out_y = out_y .+ 0.02 .* randn(length(out_y))
    append!(exp_y, [out_y])
    append!(exp_t, [out_t])

end

# distance function: mean-squared log error
function logabs_error_metric((p1, p2, p3, p4))
    N = []
    for (env_i, (env1, env2)) in enumerate(env_sets_exp)
        inpt = [p1, p2, p3, p4, env1, env2]
        t, y = output_model(inpt)
        Nsub = []
        for (i, xt) in enumerate(exp_t[env_i])
            xy = exp_y[env_i][i]
            y_step = y[findmin(abs.(t .- xt))[2]]
            push!(Nsub, (log10(xy) - log10(abs(y_step)))^2)
        end
        push!(N, mean(Nsub))
    end
    MSLE = mean(N)
    return MSLE
end

### SAMPLING ###

# priors
priors = [
    LogUniform(10^kp_logrange[1][1], 10^kp_logrange[1][2]),
    LogUniform(10^kp_logrange[2][1], 10^kp_logrange[2][2]),
    LogUniform(10^kp_logrange[3][1], 10^kp_logrange[3][2]),
    LogUniform(10^kp_logrange[4][1], 10^kp_logrange[4][2])
]

# sampling function
function sample_and_evaluate(target_successes, threshold)

    successful_params = []
    errors = []
    total_attempts = 0
    last_reported_progress = 0  # Track the last reported progress

    while length(successful_params) < target_successes
        total_attempts += 1
        sampled_params = [rand(p) for p in priors]
        error = logabs_error_metric(sampled_params)
        if error < threshold
            push!(successful_params, sampled_params)
            push!(errors, error)
        end
        
        # progress report
        current_progress = floor((length(successful_params) / target_successes) * 10) * 10
        if current_progress > last_reported_progress
            acceptance_rate = length(successful_params) / total_attempts * 100
            println("Progress: $(current_progress)% complete. Acceptance Rate: $(round(acceptance_rate, digits=2))%")
            last_reported_progress = current_progress  # update the last reported progress
        end
    end
    
    # final report
    acceptance_rate = length(successful_params) / total_attempts * 100
    println("Completed. Total Attempts: $total_attempts. Acceptance Rate: $(round(acceptance_rate, digits=2))%")
    
    # sort the parameter sets by their errors
    sort_indices = sortperm(errors)
    sorted_params = [successful_params[i] for i in sort_indices]
    
    return successful_params, errors
end

n_samples = 30
threshold = 0.010
successful_params, errors = sample_and_evaluate(n_samples, threshold)

### PLOTTING ###

using Plots

prior_samples = [rand(p, 100000) for p in priors]

for i in 1:length(priors)
    p = plot(layout = (1, 2), legend = false)
    
    # log-transform the samples
    log_prior_samples = log10.(prior_samples[i])
    log_posterior_samples = log10.([x[i] for x in successful_params])
    
    # determine x-axis limits based on the log-transformed prior distribution
    xlims_log = (log10(minimum(support(priors[i]))), log10(maximum(support(priors[i]))))

    # histogram of the log-transformed prior samples with specified x-axis limits
    histogram!(p[1], log_prior_samples, color=:gray, bins=range(kp_logrange[i][1], kp_logrange[i][2], 10), title="Prior", xlabel="log10(Parameter $(i))", ylabel="frequency", xlims=xlims_log)

    # histogram of the log-transformed posterior samples with the same x-axis limits
    histogram!(p[2], log_posterior_samples, bins=range(kp_logrange[i][1], kp_logrange[i][2], 10), title="Posterior", xlabel="log10(Parameter $(i))", ylabel="frequency", xlims=xlims_log)
    
    # vertical line for the true parameter value in the posterior plot
    vline!(p[2], [log10(kp_truth[i])], color=:red, lw = 4)
    
    display(p)
end

# plot ensemble spread
p2 = []
for j in 1:length(env_sets_exp)

    # synthetic data
    p = plot(exp_t[j], exp_y[j], 
            seriestype=:scatter, 
            markercolor=:black, 
            ylimits=(0,1), 
            xlimits=(0, 1.3*maximum(exp_t[j])),
            legend=false)

    # kinetic model solutions
    for i in 1:size(successful_params, 1)

        inpt = [successful_params[i]; env_sets_exp[j]]
        t, y = output_model(inpt)

        plot!(t, y, seriescolor=:indianred1, lw = 2.5)
    end

    # "truth"
    truth = [kp_truth; env_sets_exp[j]]
    t, y = output_model(truth)
    plot!(t, y, seriescolor=:black, linestyle=:dash, lw = 2.5)
    plot!(exp_t[j], exp_y[j], seriestype=:scatter, markercolor=:black)

    push!(p2, p) # Add the plot to the array

end
p3 = plot(p2..., layout = (3, 2), legend = false)
display(p3)

# pair plot
using CairoMakie
using PairPlots

df_posterior = DataFrame()

for i in 1:length(successful_params[1])
    df_posterior[!, Symbol("Param$(i)")] = [x[i] for x in successful_params]
end

transform!(df_posterior, names(df_posterior) .=> (x -> log10.(x)) .=> names(df_posterior))

pairplot(df_posterior, fullgrid=true, axis=(;
Param1=(; lims=(;low=kp_logrange[1][1], high=kp_logrange[1][2])),
Param2=(; lims=(;low=kp_logrange[2][1], high=kp_logrange[2][2])),
Param3=(; lims=(;low=kp_logrange[3][1], high=kp_logrange[3][2])),
Param4=(; lims=(;low=kp_logrange[4][1], high=kp_logrange[4][2])),
))