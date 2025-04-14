# This code manages the demand trace data for use in the optimization

# Define directory containing demand trace files (assumed to be defined earlier as Data_Demand_DIR)
N = length(UBus_orig)
Time = 1:T

Demand = Matrix{Float64}(undef, N, T)
psmDemand_matrix = Matrix{Float64}(undef, N, T)
csmDemand_matrix = Matrix{Float64}(undef, N, T)

for (index, bus_name) in enumerate(UBus_orig)
    demand_trace_name = Bus_data_dic[bus_name]["Demand_Trace_Name"]
    
    # Search for files in Data_Demand_DIR
    files = readdir(Data_Demand_DIR; join=true)  # Full paths
    
    # New fallback: Look for the first file starting with bus_name
    matching_files_bus_fallback = [f for f in files if startswith(lowercase(basename(f)), lowercase(bus_name)) && 
                                                      (endswith(f, ".csv") || endswith(f, ".xlsx"))]
    
    # If no explicit Demand_Trace_Name, search based on Bus_Region
    if ismissing(demand_trace_name) || demand_trace_name == "" || isnothing(demand_trace_name)
        bus_region = Bus_data_dic[bus_name]["Bus_Region"]
        matching_files = [f for f in files if occursin(lowercase(bus_region), lowercase(basename(f))) && 
                                            (endswith(f, ".csv") || endswith(f, ".xlsx"))]
        if !isempty(matching_files)
            file_path = matching_files[1]
            demand_trace_name = basename(file_path)[1:(end - (endswith(file_path, ".csv") ? 4 : 5))]  # Remove .csv or .xlsx
            println("\nProcessing demand for bus: ", bus_name, " (Demand Trace: ", demand_trace_name, " - fallback from Region)")
        elseif !isempty(matching_files_bus_fallback)
            # Use the first file starting with bus_name
            file_path = matching_files_bus_fallback[1]
            demand_trace_name = basename(file_path)[1:(end - (endswith(file_path, ".csv") ? 4 : 5))]  # Remove .csv or .xlsx
            println("\nNo region match, using first file starting with bus (", bus_name, ") for bus: ", bus_name, " (Demand Trace: ", demand_trace_name, ")")
        else
            println("\nNo matching demand trace file found for bus ", bus_name, " (Region: ", bus_region, ")")
            Demand[index, :] .= 0.0
            psmDemand_matrix[index, :] .= 0.0
            csmDemand_matrix[index, :] .= 0.0
            continue
        end
    else
        # Search for files containing demand_trace_name (e.g., region code) in the filename
        matching_files = [f for f in files if occursin(lowercase(demand_trace_name), lowercase(basename(f))) && 
                                            (endswith(f, ".csv") || endswith(f, ".xlsx"))]
        if !isempty(matching_files)
            file_path = matching_files[1]
            demand_trace_name = basename(file_path)[1:(end - (endswith(file_path, ".csv") ? 4 : 5))]  # Remove .csv or .xlsx
            # println("\nProcessing demand for bus: ", bus_name, " (Demand Trace: ", demand_trace_name, ")")
        elseif !isempty(matching_files_bus_fallback)
            # Use the first file starting with bus_name
            file_path = matching_files_bus_fallback[1]
            demand_trace_name = basename(file_path)[1:(end - (endswith(file_path, ".csv") ? 4 : 5))]  # Remove .csv or .xlsx
            println("\nNo trace match, using first file starting with bus (", bus_name, ") for bus: ", bus_name, " (Demand Trace: ", demand_trace_name, ")")
        else
            println("\nFile not found for bus: ", bus_name, " (Demand Trace: ", demand_trace_name, ")")
            Demand[index, :] .= 0.0
            psmDemand_matrix[index, :] .= 0.0
            csmDemand_matrix[index, :] .= 0.0
            continue
        end
    end

    # Check if the file exists before proceeding
    if !isempty(file_path) && isfile(file_path)
        # Determine file type and read accordingly
        if endswith(file_path, ".csv")
            demand_df = CSV.read(file_path, DataFrames.DataFrame)
        elseif endswith(file_path, ".xlsx")
            demand_xf = XLSX.readxlsx(file_path)
            sheet_name = XLSX.sheetnames(demand_xf)[1]  # Assume data is in the first sheet
            demand_df = DataFrame(XLSX.readtable(file_path, sheet_name)...)
        else
            println("Unsupported file format for bus: ", bus_name, " (Demand Trace: ", demand_trace_name, ")")
            continue
        end

        println("\n\nProcessing demand for bus: ", bus_name, " (Demand Trace: ", demand_trace_name, ")")

        # Filter the dataframe based on user-selected year, month, and day
        start_idx = findfirst((demand_df.Year .== selected_year) .&
                            (demand_df.Month .== selected_month) .&
                            (demand_df.Day .== selected_day))

        # Check if start_idx is not missing (ensure the starting date is found)
        if isnothing(start_idx)
            println("Start date not found for ", selected_year, "-", selected_month, "-", selected_day)
            continue
        end

        # Check if there are enough rows to cover all planning days
        end_idx = start_idx + D-1
        if end_idx > nrow(demand_df)
            println("Insufficient data for planning days starting from ", selected_year, "-", selected_month, "-", selected_day)
            continue
        end

        # Extract demand data for each day starting from the user-selected date
        for day in 0:D-1
            row = start_idx + day  # Row index starting from the selected date
            start_col = 4          # Start extracting from the 4th column (H1)
            end_col = 27           # End at the 27th column (H24)

            # Extract the demand data for the day
            global Demand_temp = Vector(demand_df[row, start_col:end_col]) * (Bus_data_dic[bus_name]["Demand_Trace_Weightage"])

            # Print demand for current bus and day
            println("DAY ", day + 1, "\nPeak demand ", maximum(Demand_temp))
            # println("Day ", day + 1, " Demand: ", Demand_temp)

            Prosumer_Demand_Weightage = (Bus_data_dic[bus_name]["Prosumer_Demand_p"]) / 100
            global psmDemand_temp = Prosumer_Demand_Weightage * Demand_temp
            global csmDemand_temp = (1 - Prosumer_Demand_Weightage) * Demand_temp

            println("Consumer demand\n", csmDemand_temp)
            println("Prosumer demand\n", psmDemand_temp)

            # Assign values to the correct slice of the Demand matrix
            Demand[index, day*24+1:(day+1)*24] = Demand_temp
            psmDemand_matrix[index, day*24+1:(day+1)*24] = psmDemand_temp
            csmDemand_matrix[index, day*24+1:(day+1)*24] = csmDemand_temp
        end
    else
        println("\nFile not found for bus: ", bus_name, " (Demand Trace: ", demand_trace_name, ")") # File not found, set to zero all data
        Demand[index, :] .= 0.0
        psmDemand_matrix[index, :] .= 0.0
        csmDemand_matrix[index, :] .= 0.0
    end
end

# Assign data to dict to use the same indexes of the optimization
global psmDemand = Dict()
global csmDemand = Dict()
for (index, bus_name) in enumerate(UBus_orig)
    for t in Time
        psmDemand[(bus_name, t)] = psmDemand_matrix[index, t]
        csmDemand[(bus_name, t)] = csmDemand_matrix[index, t]
    end
end

if network_model in Cu_plate
    csmDemand_cu_plate = Dict()
    psmDemand_cu_plate = Dict()
    println("Copper plate aggregating demand across all buses: ", UBus_orig)
    for t in subhorizon_total
        csm_demand_sum = sum(csmDemand[(b, t)] for b in UBus_orig if (b, t) in keys(csmDemand))
        csmDemand_cu_plate[("System", t)] = csm_demand_sum
        # println("Aggregated csmDemand for System at time $t: $csm_demand_sum")

        psm_demand_sum = sum(psmDemand[(b, t)] for b in UBus_orig if (b, t) in keys(psmDemand))
        psmDemand_cu_plate[("System", t)] = psm_demand_sum
        # println("Aggregated psmDemand for System at time $t: $psm_demand_sum")
    end
    global csmDemand = csmDemand_cu_plate
    global psmDemand = psmDemand_cu_plate
    # println("Redefined csmDemand for Cu_plate: ", keys(csmDemand))
    # println("Redefined psmDemand for Cu_plate: ", keys(psmDemand))

elseif network_model in Regional
    csmDemand_regional = Dict()
    psmDemand_regional = Dict()
    println("Regional model aggregating demand across buses per region:")
    for r in URegion
        region_buses = [b for (b, reg) in Bus_Region_links_orig if reg == r]
        println("Region $r includes buses: ", region_buses)
        for t in subhorizon_total
            csm_demand_sum = sum((csmDemand[(b, t)] for b in region_buses if (b, t) in keys(csmDemand)), init=0.0)
            csmDemand_regional[(r, t)] = csm_demand_sum
            # println("Aggregated csmDemand for region $r at time $t: $csm_demand_sum")

            psm_demand_sum = sum((psmDemand[(b, t)] for b in region_buses if (b, t) in keys(psmDemand)), init=0.0)
            psmDemand_regional[(r, t)] = psm_demand_sum
            # println("Aggregated psmDemand for region $r at time $t: $psm_demand_sum")
        end
    end
    global csmDemand = csmDemand_regional
    global psmDemand = psmDemand_regional
    # println("Redefined csmDemand for Regional: ", keys(csmDemand))
    # println("Redefined psmDemand for Regional: ", keys(psmDemand))

end

# Convert psmDemand dictionary to a DataFrame
psmDemand_df = DataFrame(
    Bus=[key[1] for key in keys(psmDemand)],   # Extract bus names from keys
    Time=[key[2] for key in keys(psmDemand)],  # Extract Time steps from keys
    Demand=collect(values(psmDemand))          # Convert demand values to a vector and assign
)

# Convert csmDemand dictionary to a DataFrame
csmDemand_df = DataFrame(
    Bus=[key[1] for key in keys(csmDemand)],   # Extract bus names from keys
    Time=[key[2] for key in keys(csmDemand)],  # Extract Time steps from keys
    Demand=collect(values(csmDemand))          # Convert demand values to a vector and assign
)

# Sort the DataFrames by Bus and Time using @orderby from DataFramesMeta
psmDemand_sorted = @orderby(psmDemand_df, :Bus, :Time)
csmDemand_sorted = @orderby(csmDemand_df, :Bus, :Time)

# Save sorted DataFrames as CSV files
CSV.write("psmDemand_sorted.csv", psmDemand_sorted)
CSV.write("csmDemand_sorted.csv", csmDemand_sorted)

println("\npsmDemand and csmDemand have been successfully sorted and saved as CSV files.")
