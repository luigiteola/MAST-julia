# This code manages the wind trace data for use in the optimization

# Define directory containing wind trace files (assumed to be defined earlier as Data_Wind_DIR)

N = length(UBus_orig)
Bus_to_region = Dict(b => r for (b, r) in Bus_Region_links_orig)

wind_generators = filter(g -> Generator_data_dic[g]["Generation_Tech"] in Wind_Tech, UGen)
# Initialize a matrix for wind trace data
Wind_trace_DR_Matrix = Matrix{Float64}(undef, N, T)

# Iterate over buses and fetch corresponding Wind trace and generator real power for 7 days
for (index, bus_name) in enumerate(UBus_orig)
    Wind_trace_name = Bus_data_dic[bus_name]["Wind_Trace_Name"]  # Fetch the wind trace name for the bus
    region = Bus_to_region[bus_name]  # Get the region for this bus from Bus_to_region
    
    # Search for files in Data_Wind_DIR
    files = readdir(Data_Wind_DIR; join=true)  # Full paths
    
    # Step 1: Look for an exact match of the trace name (without extension)
    matching_files_exact = [f for f in files if lowercase(basename(f)[1:(end - (endswith(f, ".csv") ? 4 : 5))]) == lowercase(string(Wind_trace_name)) && 
                                              (endswith(f, ".csv") || endswith(f, ".xlsx"))]
    
    # Step 2: If no exact match, look for files containing the region name
    matching_files_region = [f for f in files if occursin(lowercase(region), lowercase(string(basename(f)))) && 
                                               (endswith(f, ".csv") || endswith(f, ".xlsx"))]
    
    # Step 3: New fallback - Look for the first file starting with bus_name
    matching_files_bus_fallback = [f for f in files if startswith(lowercase(basename(f)), lowercase(bus_name)) && 
                                                      (endswith(f, ".csv") || endswith(f, ".xlsx"))]
    
    # Choose the file to use
    if !isempty(matching_files_exact)
        # Use the exact match
        file_path = matching_files_exact[1]
        Wind_trace_name = basename(file_path)[1:(end - (endswith(file_path, ".csv") ? 4 : 5))]  # Remove .csv or .xlsx
        println("\nProcessing wind trace for bus: ", bus_name, " (Wind Trace: ", Wind_trace_name, ")")
    elseif !isempty(matching_files_region)
        # Use the region-based file
        file_path = matching_files_region[1]
        Wind_trace_name = basename(file_path)[1:(end - (endswith(file_path, ".csv") ? 4 : 5))]  # Remove .csv or .xlsx
        println("\nNo exact match, using region (", region, ") wind trace for bus: ", bus_name, " (Wind Trace: ", Wind_trace_name, ")")
    elseif !isempty(matching_files_bus_fallback)
        # Use the first file starting with bus_name
        file_path = matching_files_bus_fallback[1]
        Wind_trace_name = basename(file_path)[1:(end - (endswith(file_path, ".csv") ? 4 : 5))]  # Remove .csv or .xlsx
        println("\nNo region match, using first file starting with bus (", bus_name, ") for bus: ", bus_name, " (Wind Trace: ", Wind_trace_name, ")")
    else
        # No matching file found
        Wind_trace_DR_Matrix[index, :] = zeros(T)  # Set to zero for all planning days
        println("No matching wind trace file found for bus ", bus_name, " (Wind_Trace_Name: ", Wind_trace_name, ")")
        continue
    end

    # Check if the file exists before proceeding
    if isfile(file_path)
        # Determine file type and read accordingly
        if endswith(file_path, ".csv")
            wind_df = CSV.read(file_path, DataFrames.DataFrame)
        elseif endswith(file_path, ".xlsx")
            wind_xf = XLSX.readxlsx(file_path)
            sheet_name = XLSX.sheetnames(wind_xf)[1]  # Assume data is in the first sheet
            wind_df = DataFrame(XLSX.readtable(file_path, sheet_name)...)
        else
            println("Unsupported file format for bus: ", bus_name, " (Wind Trace: ", Wind_trace_name, ")")
            continue
        end

        # Find the corresponding wind generator linked to this bus using Gen_Bus_links_orig
        gen_name = nothing
        for (g, b) in Gen_Bus_links_orig
            if b == bus_name && g in wind_generators
                gen_name = g  # Found the wind generator for this bus
                break
            end
        end

        if gen_name != nothing
            max_real_power = Generator_data_dic[gen_name]["Maximum_Real_Power"]

            # Filter the dataframe based on user-selected year, month, and day
            start_idx = findfirst((wind_df.Year .== selected_year) .&
            (wind_df.Month .== selected_month) .&
            (wind_df.Day .== selected_day))

            # Check if start_idx is not missing (ensure the starting date is found)
            if isnothing(start_idx)
                println("Start date not found for ", selected_year, "-", selected_month, "-", selected_day)
                continue
            end

            # Check if there are enough rows to cover all planning days
            end_idx = start_idx + D-1
            if end_idx > nrow(wind_df)
                println("Insufficient data for planning days starting from ", selected_year, "-", selected_month, "-", selected_day)
                continue
            end

            # Fetch the appropriate wind trace and scale by maximum real power for all planning days
            for day in 0:D-1
                row = start_idx + day
                start_col = 4
                end_col = 27

                Wind_temp = Vector(wind_df[row, start_col:end_col])
                println("Wind trace for ", Wind_trace_name, " before scaling: ", Wind_temp)
                println("Maximum real power for generator ", gen_name, ": ", max_real_power)
                Wind_temp = Wind_temp * max_real_power
                println("Wind trace for ", Wind_trace_name, " after scaling: ", Wind_temp)

                # Store the wind trace data for this bus for the corresponding day
                Wind_trace_DR_Matrix[index, day*24+1:(day+1)*24] = Wind_temp
            end
        else
            Wind_trace_DR_Matrix[index, :] = zeros(T)  # No wind generator found, set to zero for all planning days
            println("No wind generator found for bus ", bus_name)
        end
    else
        Wind_trace_DR_Matrix[index, :] = zeros(T)  # File not found, set to zero for all planning days
        println("Wind trace file not found for bus ", bus_name, " (Wind Trace: ", Wind_trace_name, ")")
    end
end

# Assign wind trace data to a dictionary for use in optimization
Wind_trace_DR = Dict()

for (index, value) in enumerate(UBus_orig)
    for t in 1:T
        Wind_trace_DR[(value, t)] = Wind_trace_DR_Matrix[index, t]  # Map wind data to dictionary
    end
end
