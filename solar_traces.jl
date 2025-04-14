# This code manages the solar trace data for use in the optimization

# Define directory containing solar trace files (assumed to be defined earlier as Data_Solar_DIR)

# Convert Bus_Region_links (Set of tuples) to a Dict for easier lookup
Bus_to_region = Dict(b => r for (b, r) in Bus_Region_links_orig)

N = length(UBus_orig)

# Filter only solar generators
solar_generators = filter(g -> Generator_data_dic[g]["Generation_Tech"] in Solar_Tech, UGen)

# Initialize a matrix for solar trace data
Solar_trace_DR_Matrix = Matrix{Float64}(undef, N, T)

# Iterate over buses and fetch corresponding Solar trace and generator real power for 7 days
for (index, bus_name) in enumerate(UBus_orig)
    PV_trace_name = Bus_data_dic[bus_name]["PV_Trace_Name"]  # Fetch the solar trace name for the bus
    region = Bus_to_region[bus_name]  # Get the region for this bus from Bus_to_region
    
    # Search for files in Data_Solar_DIR
    files = readdir(Data_Solar_DIR; join=true)  # Full paths
    
    # Step 1: Look for an exact match of the trace name (without extension)
    matching_files_exact = [f for f in files if lowercase(basename(f)[1:(end - (endswith(f, ".csv") ? 4 : 5))]) == lowercase(string(PV_trace_name)) && 
                                              (endswith(f, ".csv") || endswith(f, ".xlsx"))]
    
    # Step 2: If no exact match, look for files containing the region name and "sat"
    matching_files_region_sat = [f for f in files if occursin(lowercase(region), lowercase(string(basename(f)))) && 
                                                  occursin("sat", lowercase(string(basename(f)))) && 
                                                  (endswith(f, ".csv") || endswith(f, ".xlsx"))]
    
    # Step 3: New fallback - Look for files starting with bus_name and containing "Solar Real PV"
    matching_files_bus_fallback = [f for f in files if startswith(lowercase(basename(f)), lowercase(bus_name)) && 
                                                      occursin("solar real pv", lowercase(basename(f))) && 
                                                      (endswith(f, ".csv") || endswith(f, ".xlsx"))]
    
    # Choose the file to use
    if !isempty(matching_files_exact)
        # Use the exact match
        file_path = matching_files_exact[1]
        solar_trace_name = basename(file_path)[1:(end - (endswith(file_path, ".csv") ? 4 : 5))]  # Remove .csv or .xlsx
        println("\nProcessing solar trace for bus: ", bus_name, " (Solar Trace: ", solar_trace_name, ")")
    elseif !isempty(matching_files_region_sat)
        # Use the region + SAT file
        file_path = matching_files_region_sat[1]
        solar_trace_name = basename(file_path)[1:(end - (endswith(file_path, ".csv") ? 4 : 5))]  # Remove .csv or .xlsx
        println("\nNo exact match, using region (", region, ") + SAT solar trace for bus: ", bus_name, " (Solar Trace: ", solar_trace_name, ")")
    elseif !isempty(matching_files_bus_fallback)
        # Use the bus_name + "Solar Real PV" fallback
        file_path = matching_files_bus_fallback[1]
        solar_trace_name = basename(file_path)[1:(end - (endswith(file_path, ".csv") ? 4 : 5))]  # Remove .csv or .xlsx
        println("\nNo region match, using bus (", bus_name, ") + 'Solar Real PV' fallback for bus: ", bus_name, " (Solar Trace: ", solar_trace_name, ")")
    else
        Solar_trace_DR_Matrix[index, :] = zeros(T)  # No matching file found, set to zero for all planning days
        println("No matching solar trace file found for bus ", bus_name, " (PV_Trace_Name: ", PV_trace_name, ")")
        continue
    end

    # Check if the file exists before proceeding
    if isfile(file_path)
        # Determine file type and read accordingly
        if endswith(file_path, ".csv")
            solar_df = CSV.read(file_path, DataFrames.DataFrame)
        elseif endswith(file_path, ".xlsx")
            solar_xf = XLSX.readxlsx(file_path)
            sheet_name = XLSX.sheetnames(solar_xf)[1]  # Assume data is in the first sheet
            solar_df = DataFrame(XLSX.readtable(file_path, sheet_name)...)
        else
            println("Unsupported file format for bus: ", bus_name, " (Solar Trace: ", solar_trace_name, ")")
            continue
        end

        # Find the corresponding solar generator linked to this bus using Gen_Bus_links_orig
        gen_name = nothing
        for (g, b) in Gen_Bus_links_orig
            if b == bus_name && g in solar_generators
                gen_name = g  # Found the solar generator for this bus
                break
            end
        end

        if gen_name != nothing
            max_real_power = Generator_data_dic[gen_name]["Maximum_Real_Power"]

            # Filter the dataframe based on user-selected year, month, and day
            start_idx = findfirst((solar_df.Year .== selected_year) .&
            (solar_df.Month .== selected_month) .&
            (solar_df.Day .== selected_day))

            # Check if start_idx is not missing (ensure the starting date is found)
            if isnothing(start_idx)
                println("Start date not found for ", selected_year, "-", selected_month, "-", selected_day)
                continue
            end

            # Check if there are enough rows to cover all planning days
            end_idx = start_idx + D-1
            if end_idx > nrow(solar_df)
                println("Insufficient data for planning days starting from ", selected_year, "-", selected_month, "-", selected_day)
                continue
            end

            # Fetch the appropriate solar trace and scale by maximum real power for 7 days
            for day in 0:D-1
                row = start_idx + day
                start_col = 4
                end_col = 27

                PV_temp = Vector(solar_df[row, start_col:end_col])
                println("Solar trace for ", PV_trace_name, " before scaling: ", PV_temp)
                println("Maximum real power for generator ", gen_name, ": ", max_real_power)
                PV_temp = PV_temp * max_real_power
                println("Solar trace for ", PV_trace_name, " after scaling: ", PV_temp)

                # Store the solar trace data for this bus for the 24 hours of the corresponding day
                Solar_trace_DR_Matrix[index, day*24+1:(day+1)*24] = PV_temp
            end
        else
            Solar_trace_DR_Matrix[index, :] = zeros(T)  # No solar generator found, set to zero for all planning days
            println("No solar generator found for bus ", bus_name)
        end
    else
        Solar_trace_DR_Matrix[index, :] = zeros(T)  # File not found, set to zero for all planning days
        println("Solar trace file not found for bus ", bus_name, " (Solar Trace: ", solar_trace_name, ")")
        Solar_trace_DR[index, :] .= 0.0
    end
end

# Assign solar trace data to a dictionary for use in optimization
Solar_trace_DR = Dict()

for (index, value) in enumerate(UBus_orig)
    for t in 1:T
        Solar_trace_DR[(value, t)] = Solar_trace_DR_Matrix[index, t]  # Map solar data to dictionary
    end
end

#----------------------------------------------------------------------------------------------------------------

# Prosumer Solar PV parameters

# Initialize a matrix for solar trace data
Rooftop_solar_trace_DR_Matrix = Matrix{Float64}(undef, N, T)

# Iterate over buses and fetch corresponding Solar trace and generator real power for all planning days
for (index, bus_name) in enumerate(UBus_orig)
    Rooftop_PV_trace_name = Bus_data_dic[bus_name]["Rooftop_PV_Trace_Name"]  # Fetch the solar trace name for the bus
    region = Bus_to_region[bus_name]  # Get the region for this bus from Bus_to_region
    
    # Search for files in Data_Solar_DIR
    files = readdir(Data_Solar_DIR; join=true)  # Full paths
    
    # Step 1: Look for an exact match of the trace name (without extension)
    matching_files_exact = [f for f in files if lowercase(basename(f)[1:(end - (endswith(f, ".csv") ? 4 : 5))]) == lowercase(string(Rooftop_PV_trace_name)) && 
                                              (endswith(f, ".csv") || endswith(f, ".xlsx"))]
    
    # Step 2: If no exact match, look for files containing the region name and "sat"
    matching_files_region_sat = [f for f in files if occursin(lowercase(region), lowercase(string(basename(f)))) && 
                                                  occursin("sat", lowercase(string(basename(f)))) && 
                                                  (endswith(f, ".csv") || endswith(f, ".xlsx"))]
    
    # Step 3: New fallback - Look for files starting with bus_name and containing "Solar Real PV"
    matching_files_bus_fallback = [f for f in files if startswith(lowercase(basename(f)), lowercase(bus_name)) && 
                                                      occursin("solar real pv", lowercase(basename(f))) && 
                                                      (endswith(f, ".csv") || endswith(f, ".xlsx"))]
    
    # Choose the file to use
    if !isempty(matching_files_exact)
        # Use the exact match
        file_path = matching_files_exact[1]
        rooftop_solar_trace_name = basename(file_path)[1:(end - (endswith(file_path, ".csv") ? 4 : 5))]  # Remove .csv or .xlsx
        println("\nProcessing rooftop solar trace for bus: ", bus_name, " (Rooftop solar Trace: ", rooftop_solar_trace_name, ")")
    elseif !isempty(matching_files_region_sat)
        # Use the region + SAT file
        file_path = matching_files_region_sat[1]
        rooftop_solar_trace_name = basename(file_path)[1:(end - (endswith(file_path, ".csv") ? 4 : 5))]  # Remove .csv or .xlsx
        println("\nNo exact match, using region (", region, ") + SAT rooftop solar trace for bus: ", bus_name, " (Rooftop solar Trace: ", rooftop_solar_trace_name, ")")
    elseif !isempty(matching_files_bus_fallback)
        # Use the bus_name + "Solar Real PV" fallback
        file_path = matching_files_bus_fallback[1]
        rooftop_solar_trace_name = basename(file_path)[1:(end - (endswith(file_path, ".csv") ? 4 : 5))]  # Remove .csv or .xlsx
        println("\nNo region match, using bus (", bus_name, ") + 'Solar Real PV' fallback for bus: ", bus_name, " (Rooftop solar Trace: ", rooftop_solar_trace_name, ")")
    else
        # No matching file found
        Rooftop_solar_trace_DR_Matrix[index, :] = zeros(T)  # Set to zero for all planning days
        println("No matching rooftop solar trace file found for bus ", bus_name, " (Rooftop_PV_Trace_Name: ", Rooftop_PV_trace_name, ")")
        continue
    end

    # Check if the file exists before proceeding
    if isfile(file_path)
        # Determine file type and read accordingly
        if endswith(file_path, ".csv")
            rooftop_solar_df = CSV.read(file_path, DataFrames.DataFrame)
        elseif endswith(file_path, ".xlsx")
            rooftop_solar_xf = XLSX.readxlsx(file_path)
            sheet_name = XLSX.sheetnames(rooftop_solar_xf)[1]  # Assume data is in the first sheet
            rooftop_solar_df = DataFrame(XLSX.readtable(file_path, sheet_name)...)
        else
            println("Unsupported file format for bus: ", bus_name, " (Rooftop solar Trace: ", rooftop_solar_trace_name, ")")
            continue
        end

        rooftop_pv_max_real_power = Bus_data_dic[bus_name]["Rooftop_PV_Capacity_MW"]

        # Filter the dataframe based on user-selected year, month, and day
        start_idx = findfirst((rooftop_solar_df.Year .== selected_year) .&
        (rooftop_solar_df.Month .== selected_month) .&
        (rooftop_solar_df.Day .== selected_day))

        # Check if start_idx is not missing (ensure the starting date is found)
        if isnothing(start_idx)
            println("Start date not found for ", selected_year, "-", selected_month, "-", selected_day)
            continue
        end

        # Check if there are enough rows to cover all planning days
        end_idx = start_idx + D-1
        if end_idx > nrow(rooftop_solar_df)
            println("Insufficient data for planning days starting from ", selected_year, "-", selected_month, "-", selected_day)
            continue
        end

        # Fetch the appropriate solar trace and scale by maximum real power for all planning days
        for day in 0:D-1
            row = start_idx + day
            start_col = 4
            end_col = 27

            rooftop_PV_temp = Vector(rooftop_solar_df[row, start_col:end_col])
            println("Rooftop solar trace for ", Rooftop_PV_trace_name, " before scaling: ", rooftop_PV_temp)
            println("Maximum real power for rooftop solar at bus ", bus_name, ": ", rooftop_pv_max_real_power)
            rooftop_PV_temp = rooftop_PV_temp * rooftop_pv_max_real_power
            println("Rooftop solar trace for ", Rooftop_PV_trace_name, " after scaling: ", rooftop_PV_temp)

            # Store the solar trace data for this bus for the 24 hours of the corresponding day
            Rooftop_solar_trace_DR_Matrix[index, day*24+1:(day+1)*24] = rooftop_PV_temp
        end
    else
        Rooftop_solar_trace_DR_Matrix[index, :] = zeros(T)  # File not found, set to zero for all planning days
        println("Rooftop solar trace file not found for bus ", bus_name, " (Rooftop solar trace: ", rooftop_solar_trace_name, ")")
    end
end

# Assign solar trace data to a dictionary for use in optimization
Rooftop_solar_trace_DR = Dict()

for (index, value) in enumerate(UBus_orig)
    for t in Time
        Rooftop_solar_trace_DR[(value, t)] = Rooftop_solar_trace_DR_Matrix[index, t]  # Map solar data to dictionary
    end
end

if network_model in Cu_plate
    Rooftop_solar_trace_cu_plate = Dict()
    println("Copper plate aggregating rooftop solar across all buses: ", UBus_orig)
    for t in subhorizon_total
        rooftop_solar_sum = sum(Rooftop_solar_trace_DR[(b, t)] for b in UBus_orig if (b, t) in keys(Rooftop_solar_trace_DR); init=0.0)
        Rooftop_solar_trace_cu_plate[("System", t)] = rooftop_solar_sum
        # println("Aggregated Rooftop_solar_trace for System at time $t: $rooftop_solar_sum")
    end
    global Rooftop_solar_trace_DR = Rooftop_solar_trace_cu_plate
    # println("Redefined Rooftop_solar_trace_DR for Cu_plate: ", keys(Rooftop_solar_trace_DR))

elseif network_model in Regional
    Rooftop_solar_trace_regional = Dict()
    println("Regional model aggregating rooftop solar across buses per region:")
    for r in URegion
        region_buses = [b for (b, reg) in Bus_Region_links_orig if reg == r]
        println("Region $r includes buses: ", region_buses)
        for t in subhorizon_total
            rooftop_solar_sum = sum(Rooftop_solar_trace_DR[(b, t)] for b in region_buses if (b, t) in keys(Rooftop_solar_trace_DR); init=0.0)
            Rooftop_solar_trace_regional[(r, t)] = rooftop_solar_sum
            # println("Aggregated Rooftop_solar_trace for region $r at time $t: $rooftop_solar_sum")
        end
    end
    global Rooftop_solar_trace_DR = Rooftop_solar_trace_regional
    # println("Redefined Rooftop_solar_trace_DR for Regional: ", keys(Rooftop_solar_trace_DR))

else  # Nodal
    # No aggregation needed; keep original bus-level traces
    println("Nodal model: No aggregation for rooftop solar traces")
    # println("Rooftop_solar_trace_DR keys (unchanged): ", keys(Rooftop_solar_trace_DR))
end