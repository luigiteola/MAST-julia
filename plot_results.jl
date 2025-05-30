# This code manages the results of the optimization with rolling horizon

using Plots
using Measures
gr()

# Retrieve optimal variables from results dictionary

network_detail = string()
if network_model in Cu_plate
    network_detail = "Copper plate"
elseif network_model in Regional
    network_detail = "Regional"
else
    network_detail = "Nodal"
end

# Initialize arrays for global timeline (1:T)
v_pwr_gen = Dict(g => zeros(Float64, T) for g in UGen)
v_pwr_curtailed = Dict(g => zeros(Float64, T) for g in GenT2)
v_batc = Dict(s => zeros(Float64, T) for s in UStorage)
v_batd = Dict(s => zeros(Float64, T) for s in UStorage)
v_unserved = Dict(n => zeros(Float64, T) for n in UBus)

# Map results to global timeline using abs_subhorizon
for sub in 1:N_subhorizons
    # Calculate abs_subhorizon for this subhorizon (same as in main code)
    local abs_step_hours = (H - overlap_days) * hours_per_day
    local abs_t_start = (sub - 1) * abs_step_hours + 1
    local abs_t_end = min(sub * abs_step_hours, T)
    local abs_subhorizon = abs_t_start:abs_t_end

    # Map the results to the global timeline
    for (i, t) in enumerate(abs_subhorizon)
        for g in UGen
            v_pwr_gen[g][t] = results["Pwr_Gen_var"][sub][g][i]
        end
        for g in GenT2
            v_pwr_curtailed[g][t] = results["Pwr_curtailed"][sub][g][i]
        end
        for s in UStorage
            v_batc[s][t] = -1 * results["P_chrg"][sub][s][i]
            v_batd[s][t] = results["P_dischrg"][sub][s][i]
        end
        for n in UBus
            v_unserved[n][t] = results["unserved_demand"][sub][n][i]
        end
    end
end

# Get demand and solar generation
csmDemand_tot = Dict(t => sum(csmDemand[(n, t)] for n in UBus) for t in 1:T)
sorted_times = sort(collect(keys(csmDemand_tot)))

p_csmDemand = [csmDemand_tot[t] for t in sorted_times]
p_sysDemand = (p_csmDemand)*(1 + Loss_factor)
p_sysDemand_max = Int(round(maximum(p_sysDemand), digits=0))


# Plot generation per technology
p_blc = zeros(Float64, T)
p_brc = zeros(Float64, T)
p_gas = zeros(Float64, T)
p_hyd = zeros(Float64, T)
p_sol = zeros(Float64, T)
p_wnd = zeros(Float64, T)
p_geo = zeros(Float64, T)
p_batc = zeros(Float64, T)
p_batd = zeros(Float64, T)
p_unserved = zeros(Float64, T)

for t in 1:T
    p_batc[t] = sum(v_batc[s][t] for s in UStorage if haskey(v_batc, s); init=0.0)
    p_batd[t] = sum(v_batd[s][t] for s in UStorage if haskey(v_batd, s); init=0.0)
    p_unserved[t] = sum(v_unserved[n][t] for n in UBus if haskey(v_unserved, n); init=0.0)
end


for (gen_name, tech) in Gen_Tech_links
    for t in Time
        if tech in BlCoal_Tech
            p_blc[t] += v_pwr_gen[gen_name][t]
        elseif tech in BrCoal_Tech
            p_brc[t] += v_pwr_gen[gen_name][t]
        elseif tech in Gas_Tech
            p_gas[t] += v_pwr_gen[gen_name][t]
        elseif tech in Hydro_Tech
            p_hyd[t] += v_pwr_gen[gen_name][t]
        elseif tech in Solar_Tech
            p_sol[t] += v_pwr_gen[gen_name][t]
        elseif tech in Wind_Tech
            p_wnd[t] += v_pwr_gen[gen_name][t]
        end
    end
end

tech_data = [p_blc p_brc p_hyd p_wnd p_sol p_batd p_gas p_unserved]
tech_labels = ["BlackCoal" "BrownCoal" "Hydro" "Wind" "UtilitySolar" "UtilityStorage" "Gas" "UnservedDemand"]
tech_colors = [:black :brown :lightblue :green :orange :purple :cyan :gray]

p = areaplot(1:T,
    tech_data,
    label=tech_labels, 
    stack=:stack,
    fillcolor=tech_colors,
    # title="$(splitext(ModelFile)[1]) $(selected_year) $(network_detail) model \nPeak Demand: $(p_sysDemand_max) MW | Optimal Cost: $(round(total_cost/1e6, digits=2)) M",
    titlefont=(12, "Helvetica", :black),
    xlabel="Hours", 
    ylabel="Power (MW)",
    xguidefont=font(10, "Helvetica", :black),
    yguidefont=font(10, "Helvetica", :black), 
    legendfont=font(8, "Helvetica", :black), 
    legend=:outerright,
    top_margin=10mm,
    right_margin=10mm,
    left_margin=10mm,
    bottom_margin=10mm,
    tickformat=:plain,
    linecolor=:transparent,
    linewidth=0,
    size=(1200, 600)
)

areaplot!(1:T, p_batc,
    label="", 
    color=:purple) 


plot!(1:T, p_sysDemand, 
    label="SystemDemand", 
    color=:blue, 
    linewidth=3)


# Save the plot with a unique numbered suffix
let
    base_plot_path = joinpath(RESULTS_DIR, "system_generation")
    plot_ext = ".png"
    plot_path = base_plot_path * plot_ext
    counter = 1

    # Increment counter until an unused filename is found
    while isfile(plot_path)
        plot_path = "$(base_plot_path)_$(counter)$(plot_ext)"
        counter += 1
    end

    savefig(p, plot_path)

    # Confirm save
    if isfile(plot_path)
        @info "Plot successfully saved to: $plot_path"
    else
        @warn "Failed to save plot to: $plot_path"
    end
end