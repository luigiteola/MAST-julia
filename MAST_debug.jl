#  Copyright Notice
#  Copyright ©2025. Luigi Teola, Mohsen Aldaadi, Muhammad Adnan, and Gregor Verbič. Based on previous work done in MATLAB and AMPL by Shariq Riaz, Archie Chapman, and Gregor Verbič (2017). All Rights Reserved.

#  Permission to use, copy, modify, and distribute this software and 
#  its documentation for educational, research, and not-for-profit purposes,
#  without fee and without a signed licensing agreement, is hereby granted,
#  provided that the above copyright notice, this paragraph and the following
#  two paragraphs appear in all copies, modifications, and distributions.
#  Contact: Gregor Verbic, School of Electrical and Computer Engineering
#  Centre for Future Energy Networks, J03 - Electrical Engineering Building,
#  The University of Sydney, Ph +61 2 9351 8136,
#  gregor.verbic@sydney.edu.au.
#  
#  IN NO EVENT SHALL AUTHORS BE LIABLE TO ANY PARTY FOR DIRECT,
#  INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
#  LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS 
#  DOCUMENTATION, EVEN IF AUTHORS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
#  DAMAGE.
#  
#  AUTHORS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT 
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
#  PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, 
#  PROVIDED HEREUNDER IS PROVIDED "AS IS". AUTHORS HAVE NO OBLIGATION TO 
#  PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.


##############################################################################
# MAST_debug.jl - for debugging purposes. Same code with MAST.jl but with no command-line prompt.
# Just run/execute this file in the current Julia environment
##############################################################################

# List of required packages
required_packages = ["CSV", "DataFrames", "DataFramesMeta", "JuMP", "Gurobi", "HiGHS", "Plots", "VegaLite", "XLSX", "JSON", "Logging", "Colors",  "Measures"]

# Function to check and install missing packages
using Pkg

function install_required_packages(packages)
    # Get the list of packages in the current environment
    installed_pkgs = [dep.name for dep in values(Pkg.dependencies()) if !dep.is_direct_dep]
    installed_pkgs_direct = [dep.name for dep in values(Pkg.dependencies()) if dep.is_direct_dep]
    
    for pkg in packages
        if pkg in installed_pkgs_direct
            println("$pkg is already installed as a direct dependency.")
        elseif pkg in installed_pkgs
            println("$pkg is already installed as a transitive dependency.")
            # Optionally upgrade to direct dependency if needed
            # Pkg.add(pkg)
        else
            println("Installing $pkg...")
            Pkg.add(pkg)
        end
    end
end

# Install packages if missing
install_required_packages(required_packages)

using CSV, DataFrames, DataFramesMeta, JuMP, Gurobi, Plots, VegaLite, XLSX, JSON, Logging, ArgParse, Measures, Plots


# Load config file with fallback
function load_config(config_path::String)
    if isfile(config_path)
        return JSON.parsefile(config_path)
    else
        println("Warning: Config file '$config_path' not found. Using default values.")
        return Dict(
            "data_directory" => "models",
            "model_file" => "Test.xlsx",
            "demand_directory" => "ESOO_2013_Load_Traces",
            "solar_directory" => "Solar_Trace",
            "wind_directory" => "0910_Wind_Traces",
            "trace" => Dict("year" => 2020, "month" => 7, "day" => 1),
            "planning" => Dict("horizon_days" => 3, "rolling_horizon_days" => 2, "overlap_days" => 1, "voll" => 18600, "curtailment_penalty" => 600),
            "network_model" => "Nodal",
            "loss_factor" => 0.1,
            "reserve_margin" => 0.1,
            "solver_name" => "Gurobi",
            "mipgap" => 0.01,
            "plot_horizon_days": 3
        )
    end
end

function store_subhorizon_variables(model::JuMP.Model, sub::Int, subhorizon, abs_subhorizon)
    df = DataFrame(Subhorizon=Int[], Time=Int[], VariableName=String[], Index=String[], Value=Float64[], VarIndex=String[])
    vars = all_variables(model)
    
    for var in vars
        var_name = name(var)
        if !isempty(var_name)
            base_name = split(var_name, "[")[1]  # e.g., "Pwr_Gen_var"
            index_str = occursin("[", var_name) ? split(var_name, "[")[2][1:end-1] : "scalar"  # e.g., "gen1,73"
            indices = split(index_str, ",")
            if length(indices) == 2
                id, t_local = indices  # "generator_name", "local subhorizon time index"
                t_local_global = parse(Int, t_local)  # Global time index
                t_local_idx = t_local_global - first(subhorizon) + 1  # Convert to local subhorizon time index
                if t_local_idx < 1 || t_local_idx > length(subhorizon)
                    println("Warning: t_local $t_local_global (local idx $t_local_idx) out of bounds for $subhorizon in $var_name")
                    continue  # Skip if out of bounds
                end
                t_global = subhorizon[t_local_idx]  # Map local index to global time
                # Only store if t_global is within the abstracted subhorizon
                if t_global in abs_subhorizon
                    val = value(var)
                    var_index = string(base_name, "_", id)  # e.g., "S_Up_var_G1"
                    push!(df, (Subhorizon=sub, Time=t_global, VariableName=base_name, Index=id, Value=val, VarIndex=var_index))
                end
            elseif length(indices) == 1
                val = value(var)
                # Only store scalar variables (e.g., Build_line) in the first subhorizon
                if sub == 1
                    var_index = string(base_name, "_", index_str)
                    push!(df, (Subhorizon=sub, Time=0, VariableName=base_name, Index=index_str, Value=val, VarIndex=var_index))
                end
            end
        else
            println("Unnamed variable found: ", var)
        end
    end

    obj_value = objective_value(model)
    push!(df, (Subhorizon=sub, Time=first(abs_subhorizon), VariableName="Objective", Index="Subhorizon_$sub", Value=obj_value, VarIndex="Objective_Subhorizon_$sub"))

    fixed_cost = sum(Generator_data_dic[g]["Fix_Cost"] * value(model[:Status_var][g, t]) for g in UGen for t in abs_subhorizon)
    startup_cost = sum(Generator_data_dic[g]["Start_up_Cost"] * value(model[:S_Up_var][g, t]) for g in UGen for t in abs_subhorizon)
    shutdown_cost = sum(Generator_data_dic[g]["Shut_down_Cost"] * value(model[:S_Down_var][g, t]) for g in UGen for t in abs_subhorizon)
    variable_cost = sum(Generator_data_dic[g]["Variable_Cost"] * value(model[:Pwr_Gen_var][g, t]) for g in UGen for t in abs_subhorizon)
    unserved_demand_cost = sum(value(model[:unserved_demand][n, t]) * voll for n in UBus for t in abs_subhorizon)

    # Store cost components
    push!(df, (Subhorizon=sub, Time=first(abs_subhorizon), VariableName="Fixed_Cost", Index="Subhorizon_$sub", Value=fixed_cost, VarIndex="Fixed_Cost_Subhorizon_$sub"))
    push!(df, (Subhorizon=sub, Time=first(abs_subhorizon), VariableName="Startup_Cost", Index="Subhorizon_$sub", Value=startup_cost, VarIndex="Startup_Cost_Subhorizon_$sub"))
    push!(df, (Subhorizon=sub, Time=first(abs_subhorizon), VariableName="Shutdown_Cost", Index="Subhorizon_$sub", Value=shutdown_cost, VarIndex="Shutdown_Cost_Subhorizon_$sub"))
    push!(df, (Subhorizon=sub, Time=first(abs_subhorizon), VariableName="Variable_Cost", Index="Subhorizon_$sub", Value=variable_cost, VarIndex="Variable_Cost_Subhorizon_$sub"))
    push!(df, (Subhorizon=sub, Time=first(abs_subhorizon), VariableName="Unserved_Demand_Cost", Index="Subhorizon_$sub", Value=unserved_demand_cost, VarIndex="Unserved_Demand_Cost_Subhorizon_$sub"))
    
    return df
end

#-------------------------------------------------------------------------------------------------------------------------------------------

 # Initialize

 # Load config file
 const config_path = joinpath(@__DIR__, "config.json")
 config = load_config(config_path)

 # Select model and input files
 const DATA_DIR = joinpath(@__DIR__, config["data_directory"])
 const Data_Demand_DIR = joinpath(@__DIR__, config["demand_directory"])
 const Data_Solar_DIR = joinpath(@__DIR__, config["solar_directory"])
 const Data_Wind_DIR = joinpath(@__DIR__, config["wind_directory"])
 const ModelFile = config["model_file"]
 const network_model = config["network_model"] # Network model detail - Cu_plate, Nodal or Regional

 # Identification of different generator types
 const BlCoal_Tech = Set(["BlCT", "Coal", "Black Coal", "Sub Critical"])
 const BrCoal_Tech = Set(["BrCT", "Brown Coal"])
 const Hydro_Tech = Set(["HYD", "HYDR", "Hydro", "Water", "Hydro RoR", "Hydropower", "Hydroelectric", "Run-of-River", "Large Hydro", "Small Hydro", "Pumped Hydro", "Dam-type Hydro"])
 const Gas_Tech = Set(["OCGT", "CCGT", "Gas", "Natural Gas", "NatGas", "Cogen", "NG", "CHP", "Gas Turbine", "Gas Engine", "GT"])
 const Solar_Tech = Set(["PV", "Solar", "SOLR", "Solar PV", "CST", "CSP", "Utility Solar", "Solar Thermal", "Photovoltaic"])
 const Wind_Tech = Set(["WND", "Wind", "Onshore Wind", "Offshore Wind", "Wind Turbine"])

 # Identification of network detail options
 const Cu_plate = Set(["Cu_plate", "Cu Plate" ,"Copper Plate", "Copperplate", "Copper plate", "Single-bus", "Single Bus"])
 const Nodal = Set(["Nodal", "Node", "Full Network", "Whole Network", "Full"])
 const Regional = Set(["Regional", "Region", "Zonal", "Zone"])

 # Trace date selection
 const selected_year = config["trace"]["year"]
 const selected_month = config["trace"]["month"]
 const selected_day = config["trace"]["day"]

 # Planning horizon parameters
 const D = config["planning"]["horizon_days"]
 const H = config["planning"]["rolling_horizon_days"]
 const overlap_days = config["planning"]["overlap_days"]
 const hours_per_day = 24
 const T = D * hours_per_day  # Calculated value
 const subhorizon_total = 1:T  # Calculated range
 const abs_step_days = H - overlap_days  # Effective step size in days after accounting for overlap
 const abs_step_hours = abs_step_days * hours_per_day  # Effective step size in hours
 const N_subhorizons = ceil(Int, D / (H - overlap_days))  # Adjusted number of subhorizons
 const voll = config["planning"]["voll"] # Value of lost load (VOLL) in $/MWh for unserved demand
 const curtailment_penalty = config["planning"]["curtailment_penalty"] # Cost of curtailment in $/MWh

 # Plotting parameters
 const plot_horizon_days = config["plot_horizon_days"] # Days for plotting
 const T_plot = plot_horizon_days * hours_per_day  # Total hours for plotting
 
 # Additional parameters
 const Loss_factor = config["loss_factor"]
 const Reserve_margin = config["reserve_margin"]

 # Solver parameters
 const solver_name = config["solver_name"]  # Solver selection (e.g., Gurobi, CPLEX, etc.)
 const mipgap = config["mipgap"]  # MIP gap tolerance

 # Performance Metrics
 global solver_time = 0.0       # Total Gurobi solve time
 global solver_nodes = 0        # Total MIP nodes explored
 global solver_iterations = 0   # Total simplex/barrier iterations
 global total_cost = 0.0        # Total cost of the optimization

  # Results metrics
  global total_cost = 0.0        # Total cost of the optimization
  global gen_capacity = 0.0
  global total_unserved_demand = 0.0 # Total unserved demand
  global total_unserved_demand_cost = 0.0
  global total_fixed_cost = 0.0
  global total_startup_cost = 0.0
  global total_shutdown_cost = 0.0
  global total_variable_cost = 0.0
  global total_curtailment = 0.0
  global total_curtailment_cost = 0.0


  # Results storage (for post-processing)
  results = Dict("Status_var" => [],
                "S_Up_var" => [],
                "S_Down_var" => [],
                "Pwr_Gen_var" => [], 
                "SOC" => [], 
                "P_chrg" => [], 
                "P_dischrg" => [], 
                "U_bin" => [], 
                "psm_batt_egy" => [], 
                "psm_grid_pwr" => [], 
                "psm_feedin_pwr" => [], 
                "psm_battc_pwr" => [], 
                "psm_battd_pwr" => [], 
                "psm_batt_bin" => [],
                "psm_pv_pwr" => [], 
                "psm_pv_spill" => [],  
                "psm_feedin_priceratio" => [],
                "Pwr_line_var" => [],
                "Angle_bus_var" => [],
                "unserved_demand" => [],
                "Fixed_Cost" => [],
                "Startup_Cost" => [],
                "Shutdown_Cost" => [],
                "Variable_Cost" => [],
                "Pwr_curtailed" => [],
                "Curtailment_Cost" => [],
)
 
const RESULTS_DIR = joinpath(@__DIR__, "results")

function main()

    # Data input handling
    include("slack_bus.jl") # Slack bus function
    include("data_handling.jl") # Data input files handling
    include("demand_traces.jl") # Demand traces handling
    include("solar_traces.jl") # Solar PV and prosumer solar PV traces handling
    include("wind_traces.jl") # Wind traces handling
    include("GenT2.jl") # Generator type 2 (Wind/Solar) handling
   

    # DataFrame to collect all variables across subhorizons
    all_vars_df = DataFrame(Time=Int[])
    pivoted_all_vars_df = DataFrame(Time=Int[])
    
    #-------------------------------------------------------------------------------------------------------------------------------------------

    # Rolling horizon optimization loop
        opt_time = @elapsed begin
            for sub in 1:N_subhorizons

                println("\nStarting Subhorizon $sub")

                # Define the time window for the current subhorizon in hours
                t_start = max(1, (sub - 1) * (H - overlap_days) * hours_per_day + 1)  # Start of subhorizon, including overlap
                t_end = min(t_start + (H * hours_per_day) - 1, T)  # End of subhorizon
                subhorizon = t_start:t_end
                T_sub = length(subhorizon)  # Length of current subhorizon in hours

                # Define the effective subhorizon (non-overlapping portion) for storing results
                abs_t_start = (sub - 1) * abs_step_hours + 1  # Start of non-overlapping portion
                abs_t_end = min(sub * abs_step_hours, T)  # End of non-overlapping portion
                abs_subhorizon = abs_t_start:abs_t_end  # Non-overlapping portion to store

                println("Subhorizon $sub: t_start = $t_start, t_end = $t_end, subhorizon = $subhorizon")
                println("Effective subhorizon = $abs_subhorizon")

                if solver_name == "HiGHS"
                    model = Model(HiGHS.Optimizer)
                    set_optimizer_attribute(model, "mip_rel_gap", mipgap)
                elseif solver_name == "Gurobi"
                    model = Model(Gurobi.Optimizer)
                    set_optimizer_attribute(model, "MIPGap", mipgap)
                elseif solver_name == "CPLEX"
                    model = Model(CPLEX.Optimizer)
                    set_optimizer_attribute(model, "mipgap", mipgap)
                else
                    error("Unsupported solver: $solver_name")
                end

                # set_silent(model)

                M0 = 1e6 # Large constant

                # Generator decision variables
                @variable(model, S_Up_var[UGen, subhorizon] >= 0, Int)
                @variable(model, 0 <= Status_var[g in UGen, t in subhorizon] <= Generator_data_dic[g]["Number_Units"], Int)
                @variable(model, S_Down_var[UGen, subhorizon] >= 0, Int)
                @variable(model, Pwr_Gen_var[UGen, subhorizon] >= 0)

                # Transmission decision variables
                @variable(model, Pwr_line_var[ULine, subhorizon])
                @variable(model, 3.1416 >= Angle_bus_var[UBus, subhorizon] >= -3.1416)

                # Utility storage decision variables
                @variable(model, P_chrg[s in UStorage, t in subhorizon] >= 0)
                @variable(model, P_dischrg[s in UStorage, t in subhorizon] >= 0)
                @variable(model, SOC[s in UStorage, t in subhorizon])
                @variable(model, U_bin[s in UStorage, t in subhorizon], Bin)

                # Other decision variables
                @variable(model, unserved_demand[UBus, subhorizon] >= 0) # Unserved demand for each bus in the subhorizon
                @variable(model, Pwr_curtailed[g in GenT2, t in subhorizon] >= 0) # Power curtailed for each Type 2 generator in the subhorizon
                @constraint(model, [n in UBus, t in subhorizon], unserved_demand[n, t] <= csmDemand[n, t]) # Unserved demand constraint
                @constraint(model, [s in UStorage, t in subhorizon], sum(unserved_demand[n, t] for n in UBus) <= M0 * (1 - U_bin[s, t])) # Unserved demand constraint for storage

                # Constraints

                # Active Power reserve constraints
                @constraint(model, u_Power_Reserve[r in URegion, t in subhorizon], 
                    sum(Status_var[g, t] * (Generator_data_dic[g]["Maximum_Real_Power"]) for (g, r) in GenT1_Region_links) -
                    sum(Pwr_Gen_var[g, t] for (g, r) in GenT1_Region_links) >= sum((Reserve_margin) * csmDemand[n, t] for (n, r) in Bus_Region_links))

                # Generator limits
                @constraint(model, u_Gen_max_pwr[g in GenT1, t in subhorizon], Pwr_Gen_var[g, t] <= (Generator_data_dic[g]["Maximum_Real_Power"]) * Status_var[g, t])
                @constraint(model, u_Gen_min_pwr[g in GenT1, t in subhorizon], (Generator_data_dic[g]["Minimum_Real_Power"]) * Status_var[g, t] <= Pwr_Gen_var[g, t])

                # On/Off constraints
                @constraint(model, u_On_Off_initial[g in GenT1], S_Up_var[g, t_start] - S_Down_var[g, t_start] == Status_var[g, t_start] - prev_Status_var[g])

                if T_sub > 1
                    @constraint(model, u_On_Off[g in GenT1, t in subhorizon[2:end]], 
                    S_Up_var[g, t] - S_Down_var[g, t] == Status_var[g, t] - Status_var[g, t-1])
                end

                @constraint(model, u_max_ONunits[g in UGen, t in subhorizon], Status_var[g, t] <= (Generator_data_dic[g]["Number_Units"]))

                # Ramping constraints
                @constraint(model, u_ramp_up_ini[g in GenT1; Generator_data_dic[g]["Ramp_Up_Rate"] < Generator_data_dic[g]["Maximum_Real_Power"]], 
                    Pwr_Gen_var[g, t_start] - prev_Pwr_Gen_var[g] <= 
                    Status_var[g, t_start] * Generator_data_dic[g]["Ramp_Up_Rate"] + 
                    (Generator_data_dic[g]["Minimum_Real_Power"] - Generator_data_dic[g]["Ramp_Up_Rate"]) * S_Up_var[g, t_start])
                @constraint(model, u_ramp_down_ini[g in GenT1; (Generator_data_dic[g]["Ramp_Down_Rate"]) < (Generator_data_dic[g]["Maximum_Real_Power"])], 
                    prev_Pwr_Gen_var[g] - Pwr_Gen_var[g, t_start] <= 
                    Status_var[g, t_start] * (Generator_data_dic[g]["Ramp_Down_Rate"]) + 
                    ((Generator_data_dic[g]["Minimum_Real_Power"]) - (Generator_data_dic[g]["Ramp_Down_Rate"])) * S_Down_var[g, t_start])

                if T_sub > 1
                    @constraint(model, u_ramp_up[g in GenT1, t in subhorizon[2:end]; Generator_data_dic[g]["Ramp_Up_Rate"] < Generator_data_dic[g]["Maximum_Real_Power"]], 
                        Pwr_Gen_var[g, t] - Pwr_Gen_var[g, t-1] <= 
                        Status_var[g, t] * Generator_data_dic[g]["Ramp_Up_Rate"] + 
                        (Generator_data_dic[g]["Minimum_Real_Power"] - Generator_data_dic[g]["Ramp_Up_Rate"]) * S_Up_var[g, t])

                    @constraint(model, u_ramp_down[g in GenT1, t in subhorizon[2:end]; (Generator_data_dic[g]["Ramp_Down_Rate"]) < (Generator_data_dic[g]["Maximum_Real_Power"])], 
                        Pwr_Gen_var[g, t-1] - Pwr_Gen_var[g, t] <= 
                        Status_var[g, t-1] * (Generator_data_dic[g]["Ramp_Down_Rate"]) +
                        ((Generator_data_dic[g]["Minimum_Real_Power"]) - (Generator_data_dic[g]["Ramp_Down_Rate"])) * S_Down_var[g, t-1])
                end

                # Generator Minimum Up/Down Time Constraints
                for g in GenT1
                    if Generator_data_dic[g]["MUT_h"] > 1
                        # Minimum Up Time
                        for t in subhorizon
                            if t >= t_start + Generator_data_dic[g]["MUT_h"] - 1
                                @constraint(model, Status_var[g, t] >= sum(S_Up_var[g, t - t1] for t1 in 0:Generator_data_dic[g]["MUT_h"]-1))
                            elseif t_start > 1 && t < t_start + Generator_data_dic[g]["MUT_h"] - 1
                                # Use prev_S_Up_var for times before t_start
                                past_hours = (t_start - 1):-1:max(t_start - Generator_data_dic[g]["MUT_h"], 1)
                                @constraint(model, Status_var[g, t] >= 
                                    sum(t - t1 >= t_start ? S_Up_var[g, t - t1] : prev_S_Up_var[g] for t1 in 0:(t - t_start)) + 
                                    sum(MUT_ini[(g, t1)] for t1 in past_hours if t1 <= t - t_start))
                            elseif t_start == 1 && t < Generator_data_dic[g]["MUT_h"]
                                # Initial subhorizon
                                @constraint(model, Status_var[g, t] >= sum(S_Up_var[g, t - t1] for t1 in 0:(t - t_start)) + 
                                    sum(MUT_ini[(g, t1)] for t1 in 1:(t - 1) if t1 <= t - t_start))
                            end
                        end
                    end

                    if Generator_data_dic[g]["MDT_h"] > 1
                        # Minimum Down Time
                        for t in subhorizon
                            if t >= t_start + Generator_data_dic[g]["MDT_h"] - 1
                                @constraint(model, Status_var[g, t] <= Generator_data_dic[g]["Number_Units"] - 
                                    sum(S_Down_var[g, t - t1] for t1 in 0:Generator_data_dic[g]["MDT_h"]-1))
                            elseif t_start > 1 && t < t_start + Generator_data_dic[g]["MDT_h"] - 1
                                # Use prev_S_Down_var for times before t_start
                                past_hours = (t_start - 1):-1:max(t_start - Generator_data_dic[g]["MDT_h"], 1)
                                @constraint(model, Status_var[g, t] <= Generator_data_dic[g]["Number_Units"] - 
                                    sum(t - t1 >= t_start ? S_Down_var[g, t - t1] : prev_S_Down_var[g] for t1 in 0:(t - t_start)) - 
                                    sum(MDT_ini[(g, t1)] for t1 in past_hours if t1 <= t - t_start))
                            elseif t_start == 1 && t < Generator_data_dic[g]["MDT_h"]
                                # Initial subhorizon
                                @constraint(model, Status_var[g, t] <= Generator_data_dic[g]["Number_Units"] - 
                                    sum(S_Down_var[g, t - t1] for t1 in 0:(t - t_start)) - 
                                    sum(MDT_ini[(g, t1)] for t1 in 1:(t - 1) if t1 <= t - t_start))
                            end
                        end
                    end
                end

                # Resource availability constraint for Type 2 generators (Wind and Solar)
                @constraint(model,u_Resource_availability_T2[g in GenT2, t in subhorizon], Pwr_Gen_var[g, t] + Pwr_curtailed[g, t] == Resource_trace_T2[(g, t)])
                @constraint(model,u_G_T2_min_pwr[g in GenT2, t in subhorizon], Generator_data_dic[g]["Minimum_Real_Power"] <= Pwr_Gen_var[g, t])
                @constraint(model, [g in GenT2, t in subhorizon], Pwr_curtailed[g,t] <= Resource_trace_T2[(g, t)])

                # Transmission constraints

                # DC Optimal Power Flow Constraints 
                @constraint(model, u_power_flow[l in ULine, t in subhorizon],
                    Pwr_line_var[l, t] == Base_power * (1 / Line_data_dic[l]["Reactance_pu"]) * 
                    (sum(Angle_bus_var[n1, t] for (ll, n1) in Line_end1_Bus_links if ll == l) -
                    sum(Angle_bus_var[n2, t] for (ll, n2) in Line_end2_Bus_links if ll == l)))

                # Thermal limits
                @constraint(model, u_thermal_limit_ub[l in ULine, t in subhorizon], Pwr_line_var[l, t] <= (Line_data_dic[l]["Thermal_Limit_MVA"]))
                @constraint(model, u_thermal_limit_lb[l in ULine, t in subhorizon], -(Line_data_dic[l]["Thermal_Limit_MVA"]) <= Pwr_line_var[l, t])

                # AC line angle stablility
                @constraint(model,
                    u_angle_limit_ub[l in ULine, t in subhorizon],
                    sum(Angle_bus_var[n1, t] for (ll, n1) in Line_end1_Bus_links if ll == l) -
                    sum(Angle_bus_var[n2, t] for (ll, n2) in Line_end2_Bus_links if ll == l) <=
                    Line_data_dic[l]["Maximum_Angle_Limit_degree"])

                @constraint(model,
                    u_angle_limit_lb[l in ULine, t in subhorizon],
                    -Line_data_dic[l]["Maximum_Angle_Limit_degree"] <=
                    sum(Angle_bus_var[n1, t] for (ll, n1) in Line_end1_Bus_links if ll == l) -
                    sum(Angle_bus_var[n2, t] for (ll, n2) in Line_end2_Bus_links if ll == l))


                # Reference angle for the slack bus
                @constraint(model, u_angle_zero[t in subhorizon], Angle_bus_var[Slack_bus, t] == 0)

                # Utility Scale Battery Constraints
                @constraint(model, [s in UStorage], SOC[s, t_start] == prev_SOC[s] + Charging_eff[s] * P_chrg[s, t_start] - P_dischrg[s, t_start] / Discharging_eff[s])
                if T_sub > 1
                    @constraint(model, [s in UStorage, t in subhorizon[2:end]], SOC[s, t] == SOC[s, t-1] + (Charging_eff[s] * P_chrg[s, t]) - (P_dischrg[s, t] / Discharging_eff[s]))
                end
                
                @constraint(model, [s in UStorage, t in subhorizon], P_chrg[s, t] <= Utility_storage_data_dic[s]["Maximum_Charge_Rate_MWh"] * U_bin[s, t])
                @constraint(model, [s in UStorage, t in subhorizon], P_dischrg[s, t] <= Utility_storage_data_dic[s]["Maximum_Discharge_Rate_MWh"] * (1 - U_bin[s, t]))
                @constraint(model, [s in UStorage, t in subhorizon], Utility_storage_data_dic[s]["Minimum_Storage_Capacity_MWh"] <= SOC[s, t] <= Utility_storage_data_dic[s]["Maximum_Storage_Capacity_MWh"])

                @constraint(model, u_Balance[n in UBus, t in subhorizon],
                    sum(Pwr_Gen_var[g, t] for (g, nb) in GenT1_Bus_links if nb == n) +
                    sum(Pwr_Gen_var[g, t] for (g, nb) in GenT2_Bus_links if nb == n) +
                    sum(Pwr_line_var[l2, t] for (l2, nb) in Line_end2_Bus_links if nb == n) +
                    sum(P_dischrg[s, t] for (s, nb) in Storage_Bus_links if nb == n) +
                    unserved_demand[n,t]
                    ==
                    csmDemand[n,t] * (1 + Loss_factor) +
                    sum(Pwr_line_var[l1, t] for (l1, nb) in Line_end1_Bus_links if nb == n) +
                    sum(P_chrg[s, t] for (s, nb) in Storage_Bus_links if nb == n)
                )

                #-------------------------------------------------------------------------------------------------------------------------------------------

                # Objecive Function
                @objective(model, Min,
                    sum(sum(
                            (Generator_data_dic[g]["Fix_Cost"]) * Status_var[g, t] +
                            (Generator_data_dic[g]["Start_up_Cost"]) * S_Up_var[g, t] +
                            (Generator_data_dic[g]["Shut_down_Cost"]) * S_Down_var[g, t] +
                            (Generator_data_dic[g]["Variable_Cost"]) * Pwr_Gen_var[g, t] +
                            (g in GenT2 ? curtailment_penalty * Pwr_curtailed[g, t] : 0)
                            for g in UGen) +
                            voll * sum(unserved_demand[n, t] for n in UBus) 
                        for t in subhorizon)
                )

                #-------------------------------------------------------------------------------------------------------------------------------------------

                optimize!(model)

                # Collect solver stats
                global solver_time += JuMP.solve_time(model)
                global solver_nodes += MOI.get(model, MOI.NodeCount())          
                global solver_iterations += MOI.get(model, MOI.SimplexIterations())

                # Check optimization status
                if termination_status(model) == MOI.OPTIMAL

                    # Get cost and scale to account only for non-overlapping hours
                    non_overlap_hours = length(abs_subhorizon)
                    total_hours = length(subhorizon)
                    cost_scaling_factor = non_overlap_hours / total_hours
                    subhorizon_cost = objective_value(model) * cost_scaling_factor
                    
                    unserved_demand_sub = sum((value(unserved_demand[n, t]) for n in UBus for t in abs_subhorizon), init=0.0)
                    curtailment_sub = sum((value(Pwr_curtailed[g, t]) for g in GenT2 for t in abs_subhorizon), init=0.0)

                    global total_cost += subhorizon_cost
                    global gen_capacity = sum(Generator_data_dic[g]["Maximum_Real_Power"] * Generator_data_dic[g]["Number_Units"] for g in UGen)
                    global total_unserved_demand += unserved_demand_sub
                    global total_curtailment += curtailment_sub

                    fixed_cost_sub = sum(Generator_data_dic[g]["Fix_Cost"] * value(Status_var[g, t]) for g in UGen for t in abs_subhorizon) * cost_scaling_factor
                    startup_cost_sub = sum(Generator_data_dic[g]["Start_up_Cost"] * value(S_Up_var[g, t]) for g in UGen for t in abs_subhorizon) * cost_scaling_factor
                    shutdown_cost_sub = sum(Generator_data_dic[g]["Shut_down_Cost"] * value(S_Down_var[g, t]) for g in UGen for t in abs_subhorizon) * cost_scaling_factor
                    variable_cost_sub = sum(Generator_data_dic[g]["Variable_Cost"] * value(Pwr_Gen_var[g, t]) for g in UGen for t in abs_subhorizon) * cost_scaling_factor
                    unserved_demand_cost_sub = sum((value(unserved_demand[n, t]) * voll for n in UBus for t in abs_subhorizon),init=0.0) * cost_scaling_factor
                    curtailment_cost_sub = sum((value(Pwr_curtailed[g, t]) * curtailment_penalty for g in GenT2 for t in abs_subhorizon),init=0.0) * cost_scaling_factor

                    global total_fixed_cost += fixed_cost_sub
                    global total_startup_cost += startup_cost_sub
                    global total_shutdown_cost += shutdown_cost_sub
                    global total_variable_cost += variable_cost_sub
                    global total_unserved_demand_cost += unserved_demand_cost_sub
                    global total_curtailment_cost +=  curtailment_cost_sub

                    println("Sub-horizon $sub (Hours $abs_t_start to $abs_t_end) optimized.")
                    println("   Optimal Cost (\$) =  $(round(subhorizon_cost, digits=0))")
                    println("   Installed Capacity (MW) = $(round(gen_capacity, digits=0 ))")
                    println("   Dispatched Generation (MWh) = $(round(sum(value(Pwr_Gen_var[g, t]) for g in UGen for t in abs_subhorizon), digits=0))")
                    println("   Unserved Demand (MWh) = $(round(unserved_demand_sub, digits=0))")  
                    println("   Curtailment (MWh) = $(round(curtailment_sub, digits=0))")

                    # Store results
                    push!(results["Status_var"], Dict(g => [value(Status_var[g, t]) for t in abs_subhorizon] for g in UGen))
                    push!(results["S_Up_var"], Dict(g => [value(S_Up_var[g, t]) for t in abs_subhorizon] for g in UGen))
                    push!(results["S_Down_var"], Dict(g => [value(S_Down_var[g, t]) for t in abs_subhorizon] for g in UGen))
                    push!(results["Pwr_Gen_var"], Dict(g => [value(Pwr_Gen_var[g, t]) for t in abs_subhorizon] for g in UGen))
                    push!(results["P_chrg"], Dict(s => [value(P_chrg[s, t]) for t in abs_subhorizon] for s in UStorage))
                    push!(results["P_dischrg"], Dict(s => [value(P_dischrg[s, t]) for t in abs_subhorizon] for s in UStorage))
                    push!(results["SOC"], Dict(s => [value(SOC[s, t]) for t in abs_subhorizon] for s in UStorage))
                    push!(results["U_bin"], Dict(s => [value(U_bin[s, t]) for t in abs_subhorizon] for s in UStorage))
                    push!(results["Pwr_curtailed"], Dict(g => [value(Pwr_curtailed[g, t]) for t in abs_subhorizon] for g in GenT2))
                    push!(results["Pwr_line_var"], Dict(l => [value(Pwr_line_var[l, t]) for t in abs_subhorizon] for l in ULine))
                    push!(results["Angle_bus_var"], Dict(n => [value(Angle_bus_var[n, t]) for t in abs_subhorizon] for n in UBus))
                    push!(results["unserved_demand"], Dict(n => [value(unserved_demand[n, t]) for t in abs_subhorizon] for n in UBus))
                    push!(results["Fixed_Cost"], fixed_cost_sub)
                    push!(results["Startup_Cost"], startup_cost_sub)
                    push!(results["Shutdown_Cost"], shutdown_cost_sub)
                    push!(results["Variable_Cost"], variable_cost_sub)
                    push!(results["Curtailment_Cost"], curtailment_cost_sub)
 
                    # # Store all decision variables for this subhorizon
                    sub_df = store_subhorizon_variables(model, sub, subhorizon, abs_subhorizon)

                    # Pivot the subhorizon DataFrame
                    sub_pivoted_df = unstack(sub_df, :Time, :VarIndex, :Value, combine=first, fill=0.0)
                    sort!(sub_pivoted_df, :Time)

                    # Ensure all time steps in abs_subhorizon are present
                    expected_times = DataFrame(Time=collect(abs_subhorizon))
                    sub_pivoted_df = leftjoin(expected_times, sub_pivoted_df, on=:Time)
                    for col in names(sub_pivoted_df)
                        if eltype(sub_pivoted_df[!, col]) <: Union{Missing, Float64}
                            sub_pivoted_df[!, col] = coalesce.(sub_pivoted_df[!, col], 0.0)
                        end
                    end

                    # Append to all_vars_df, ensuring consistent columns
                    if isempty(all_vars_df)
                        all_vars_df = sub_pivoted_df
                    else
                        # Ensure all columns from sub_pivoted_df are in all_vars_df
                        for col in names(sub_pivoted_df)
                            if col ∉ names(all_vars_df) && col != "Time"
                                all_vars_df[!, col] = zeros(Float64, nrow(all_vars_df))
                            end
                        end
                        # Ensure all columns from all_vars_df are in sub_pivoted_df
                        for col in names(all_vars_df)
                            if col ∉ names(sub_pivoted_df) && col != "Time"
                                sub_pivoted_df[!, col] = zeros(Float64, nrow(sub_pivoted_df))
                            end
                        end
                        # Reorder columns in sub_pivoted_df to match all_vars_df
                        if !isempty(names(all_vars_df))
                            sub_pivoted_df = select(sub_pivoted_df, names(all_vars_df))
                        end
                        # Append vertically
                        append!(all_vars_df, sub_pivoted_df)
                    end

                    # Update prev_ variables
                    for g in UGen
                        prev_Status_var[g] = value(Status_var[g, t_end])
                        prev_S_Up_var[g] = value(S_Up_var[g, t_end])  # Update previous startup
                        prev_S_Down_var[g] = value(S_Down_var[g, t_end])  # Update previous shutdown
                        if g in GenT1
                            prev_Pwr_Gen_var[g] = value(Pwr_Gen_var[g, t_end])
                        end
                    end
                    for s in UStorage
                        prev_SOC[s] = value(SOC[s, t_end])
                    end
                else
                    println("Sub-horizon $sub failed to solve optimally. Status: ", termination_status(model))
                    println("Primal Status: ", primal_status(model))
                    println("Dual Status: ", dual_status(model))
                    error("Optimization failed for Sub-horizon $sub")
                    break
                end
            end

            csmDemand_tot = Dict(t => sum(csmDemand[(n, t)] for n in UBus) for t in 1:T)
            sorted_times = sort(collect(keys(csmDemand_tot)))

            # Convert dictionary values to array (in sorted order)
    
            p_csmDemand = [csmDemand_tot[t] for t in sorted_times]
            p_sysDemand = (p_csmDemand * (1 + Loss_factor)) # total demand
            p_sysDemand_max = Int(round(maximum(p_sysDemand), digits=0))

            # System metrics
            println("\nSystem Summary:")
            println(" Model: $ModelFile")
            println(" Planning Horizon: $D days ($T hours) from $selected_year-$selected_month-$selected_day")
            println(" Rolling Horizon: $H days") 
            println(" Installed Generation Capacity (MW): $gen_capacity")
            println(" System Peak Demand (MW): $p_sysDemand_max")
            println(" Capacity Margin (%): $(round((gen_capacity - p_sysDemand_max) / gen_capacity * 100, digits=2))")
            println(" Total Curtailment (MWh) = $(round(total_curtailment, digits=0))")
            println(" Total Unserved Demand (MWh) = $(round(total_unserved_demand, digits=0))")

            println("\nCost Summary:")
            println(" Total Fixed Cost (\$): $(round(total_fixed_cost, digits=2))")
            println(" Total Startup Cost (\$): $(round(total_startup_cost, digits=2))")
            println(" Total Shutdown Cost (\$): $(round(total_shutdown_cost, digits=2))")
            println(" Total Variable Cost (\$): $(round(total_variable_cost, digits=2))")
            println(" Total Curtailment Cost (\$): $(round(total_curtailment_cost, digits=2))")
            println(" Total Unserved Demand Cost (\$): $(round(total_unserved_demand_cost, digits=2))")
            println(" Total Cost (\$): $(round(total_cost, digits=2))")
            

           # Storing of solution
            @info "Storing solution..."
            local base_output_filename = "solution"  # Hardcode the filename
            local base_output_path = joinpath(RESULTS_DIR, base_output_filename) # Base output path
            local output_ext = ".csv"
            local output_path = joinpath(RESULTS_DIR, base_output_filename * output_ext)  # Use joinpath for the full path
            local counter = 1

            # Increment counter until an unused filename is found
            while isfile(output_path)
                output_path = joinpath(RESULTS_DIR, "$(base_output_filename)_$(counter)$(output_ext)")
                counter += 1
            end

            # Sort all_vars_df by Time for consistency
            sort!(all_vars_df, :Time)

            # Save solution to CSV with confirmation
            CSV.write(output_path, all_vars_df)
            if isfile(output_path)
                @info "Solution saved successfully to: $output_path"
            else
                @warn "Failed to save decision variables to: $output_path (file does not exist after write)"
            end
            # Plotting of results
            @info "Plotting results..."
            include("plot_results.jl") # Results handling - system generation
           
        end
    end

main()