#  Copyright Notice
#  Copyright ©2025. Luigi S. Teola, Mohsen S. Aldaadi, Muhammad Adnan, and Gregor Verbič. Based on previous work done in MATLAB and AMPL by Shariq Riaz, Archie C. Chapman, and Gregor Verbič (2017) . All Rights Reserved.
#  This code manages the input data for use in the optimization

worksheet = joinpath(DATA_DIR, ModelFile)

# Filter until end of dataset
function filter_data_end(df)
    for (i, row) in enumerate(eachrow(df))
        # Check if "END OF DATA" is in any column AND there’s at least one missing value
        row_values = collect(row)  # Convert row to array for easier checking
        has_end_of_data = any(x -> x isa String && x == "END OF DATA", row_values)
        has_missing = any(ismissing, row_values)
        if has_end_of_data && has_missing
            # Return rows up to (but not including) this row
            return df[1:i-1, :]
        end
    end
    # If "END OF DATA" is not found, return the whole DataFrame
    return df
end


# Create data frames
Generator_data_df = filter_data_end(DataFrame(XLSX.readtable(worksheet, "Generator Data")))
Bus_data_df = filter_data_end(DataFrame(XLSX.readtable(worksheet, "Bus Data")))
Line_data_df = filter_data_end(DataFrame(XLSX.readtable(worksheet, "Branch Data")))
Utility_data_df = filter_data_end(DataFrame(XLSX.readtable(worksheet, "Utility Storage Data")))


# Create dictionaries per data frame to assign data per column. Column names should striclty be followed

Generator_data_dic = Dict(
    Generator_data_df[!, "Generator Name"][i] => Dict(
        "Location_Bus" => Generator_data_df[!, "Location Bus"][i],
        "Number_Units" => Generator_data_df[!, "Number of Units"][i],
        "Connected_Grid" => Generator_data_df[!, "Connected to Grid"][i],
        "Apparent_Power_Rating" => Generator_data_df[!, "Apparent Power Rating (MVA)"][i],
        "Fix_Cost" => Generator_data_df[!, "Fix Cost (\$)"][i],
        "Start_up_Cost" => Generator_data_df[!, "Start up Cost (\$)"][i],
        "Shut_down_Cost" => Generator_data_df[!, "Shut down Cost (\$)"][i],
        "Variable_Cost" => Generator_data_df[!, "Variable Cost (\$/MW)"][i],
        "Maximum_Real_Power" => Generator_data_df[!, "Maximum Real Power (MW)"][i],
        "Minimum_Real_Power" => Generator_data_df[!, "Minimum Real Power (MW)"][i],
        "Maximum_Reactive_Power" => Generator_data_df[!, "Maximum Reactive Power (MVar)"][i],
        "Minimum_Reactive_Power" => Generator_data_df[!, "Minimum Reactive Power (MVar)"][i],
        "Ramp_Up_Rate" => Generator_data_df[!, "Ramp Up Rate (MW/h)"][i],
        "Ramp_Down_Rate" => Generator_data_df[!, "Ramp Down Rate (MW/h)"][i],
        "MUT_h" => Generator_data_df[!, "MUT (hour)"][i],
        "MDT_h" => Generator_data_df[!, "MDT (hour)"][i],
        "Generation_Type" => Generator_data_df[!, "Generation Type"][i],
        "Plot_Color" => Generator_data_df[!, "Plot Color"][i],
        "Generation_Tech" => Generator_data_df[!, "Generation Tech"][i]
    )
    for i in 1:nrow(Generator_data_df)
)

# Create Bus_data_dic from Bus_data_df
Bus_data_dic = Dict()
for i in 1:nrow(Bus_data_df)
    # Base dictionary with all mandatory fields
    bus_dict = Dict(
        "Bus_Region" => Bus_data_df[!, "Bus Region"][i],
        "Bus_Type" => Bus_data_df[!, "Bus Type"][i],
        "Demand_Trace_Weightage" => Bus_data_df[!, "Demand Trace Weightage"][i],
        "Power_Factor" => Bus_data_df[!, "Power Factor"][i],
        "Prosumer_Demand_p" => Bus_data_df[!, "Prosumer Demand (%)"][i],
        "Rooftop_PV_Capacity_MW" => Bus_data_df[!, "Rooftop PV Capacity (MW)"][i],
        "Feedin_Price_Ratio" => Bus_data_df[!, "Feedin Price Ratio"][i],
        "Maximum_Battery_Capacity_MWh" => Bus_data_df[!, "Maximum Battery Capacity (MWh)"][i],
        "Minimum_Battery_Capacity_MWh" => Bus_data_df[!, "Minimum Battery Capacity (MWh)"][i],
        "Maximum_Charge_Rate" => Bus_data_df[!, "Maximum Charge Rate (MW/h)"][i],
        "Maximum_Discharge_Rate" => Bus_data_df[!, "Maximum Discharge Rate (MW/h)"][i],
        "Battery_Efficiency" => Bus_data_df[!, "Battery Efficiency (%)"][i],
        "Minimum_Voltage" => Bus_data_df[!, "Minimum Voltage (pu)"][i],
        "Maximum_Voltage" => Bus_data_df[!, "Maximum Voltage (pu)"][i],
        "Base_kV" => Bus_data_df[!, "Base_kV"][i],
        "Demand_Trace_Name" => Bus_data_df[!, "Demand Trace Name"][i],
        "Wind_Trace_Name" => Bus_data_df[!, "Wind Trace Name"][i],
        "PV_Trace_Name" => Bus_data_df[!, "PV Trace Name"][i],
        "CST_Trace_Name" => Bus_data_df[!, "CST Trace Name"][i],
        "Rooftop_PV_Trace_Name" => Bus_data_df[!, "Rooftop PV Trace Name"][i]
    )

    # Conditionally add Bus_Subregion if the "Sub-region" column exists
    if "Sub-region" in names(Bus_data_df)
        bus_dict["Bus_Subregion"] = Bus_data_df[!, "Sub-region"][i]
    end

    # Conditionally add Bus_REZ if the "REZ" column exists
    if "REZ" in names(Bus_data_df)
        bus_dict["Bus_REZ"] = Bus_data_df[!, "REZ"][i]
    end

    # Assign the dictionary to the bus name
    Bus_data_dic[Bus_data_df[!, "Bus Name"][i]] = bus_dict
end

Line_data_dic = Dict(
    Line_data_df[!, "Line Name"][i] => Dict(
        "End_1_Bus" => Line_data_df[!, "End 1 Bus Name"][i],
        "End_2_Bus" => Line_data_df[!, "End 2 Bus Name"][i],
        "In_Service" => Line_data_df[!, "In Service"][i],
        "Resistance_pu" => Line_data_df[!, "Resistance (pu)"][i],
        "Reactance_pu" => Line_data_df[!, "Reactance (pu)"][i],
        "Susceptance_pu" => Line_data_df[!, "Susceptance (pu)"][i],
        "Thermal_Limit_MVA" => Line_data_df[!, "Thermal Limit (MVA)"][i],
        "Nature_of_Line" => Line_data_df[!, "Nature of Line"][i],
        "Maximum_Angle_Limit_degree" => Line_data_df[!, "Maximum Angle Limit (degree)"][i],
    )
    for i in 1:nrow(Line_data_df)
)

# Create Utility_storage_data_dic from Utility_data_df
Utility_storage_data_dic = Dict()
for i in 1:nrow(Utility_data_df)
    # Base dictionary with all mandatory fields
    utility_dict = Dict(
        "Location_Bus" => Utility_data_df[!, "Location Bus"][i],
        "Connected_to_Grid" => Utility_data_df[!, "Connected to Grid"][i],
        "Maximum_Storage_Capacity_MWh" => Utility_data_df[!, "Maximum Storage Capacity (MWh)"][i],
        "Minimum_Storage_Capacity_MWh" => Utility_data_df[!, "Minimum Storage Capacity (MWh)"][i],
        "Maximum_Charge_Rate_MWh" => Utility_data_df[!, "Maximum Charge Rate (MW/h)"][i],
        "Maximum_Discharge_Rate_MWh" => Utility_data_df[!, "Maximum Discharge Rate (MW/h)"][i],
        "Plot_Color" => Utility_data_df[!, "Plot Color"][i]
    )

    # Check for efficiency columns
    if "Storage Efficiency (%)" in names(Utility_data_df)
        # Use Storage Efficiency for both charging and discharging
        utility_dict["Charging_Efficiency"] = Utility_data_df[!, "Storage Efficiency (%)"][i]
        utility_dict["Discharging_Efficiency"] = Utility_data_df[!, "Storage Efficiency (%)"][i]
    elseif "Charging Efficiency (%)" in names(Utility_data_df) && "Discharging Efficiency (%)" in names(Utility_data_df)
        # Use separate Charging and Discharging Efficiency
        utility_dict["Charging_Efficiency"] = Utility_data_df[!, "Charging Efficiency (%)"][i]
        utility_dict["Discharging_Efficiency"] = Utility_data_df[!, "Discharging Efficiency (%)"][i]
    else
        # Neither set of columns is present; throw an error or set a default
        error("No efficiency columns found in Utility_data_df for utility storage $(Utility_data_df[!, "Utility Storage Name"][i])")
        # Default values if not found:
        # utility_dict["Charging_Efficiency"] = 0.9 (90%)
        # utility_dict["Discharging_Efficiency"] = 0.9 (90%)
    end

    # Assign the dictionary to the utility storage name
    Utility_storage_data_dic[Utility_data_df[!, "Utility Storage Name"][i]] = utility_dict
end


Generator_data_keys = keys(Generator_data_dic)
Bus_data_keys = keys(Bus_data_dic)
Line_data_keys = keys(Line_data_dic)
Utility_data_keys = keys(Utility_storage_data_dic)

# Filter generators not connected to grid (Connected_Grid != 1)
UGen = Set(key for key in Generator_data_keys if Generator_data_dic[key]["Connected_Grid"] == 1)

# Filter lines that are not in service (In_Service != 1)
# ULine = Set(key for key in Line_data_keys if Line_data_dic[key]["In_Service"] == 1)
ULine = Set(Line_data_keys) # For testing purposes, all lines are in service

# Filter utility storage not connected to grid (Connected_to_Grid != 1)
UStorage = Set(key for key in Utility_data_keys if Utility_storage_data_dic[key]["Connected_to_Grid"] == 1)


UBus = Set(Bus_data_keys)
UBus_orig = UBus
N = length(UBus)


# CROSS SETS LINKS

# Generator to bus links
Gen_Bus_links = Set()
for value in UGen
    Location_bus = Generator_data_dic[value]["Location_Bus"]
    println("Generator ", value, " connected to Bus ", Location_bus)
    push!(Gen_Bus_links, (value, Location_bus))
end
Gen_Bus_links_orig = Gen_Bus_links # Save original Bus_Region_links for use in different network model details

# Generator to region links
Gen_Region_links = Set()
for value in UGen
    Location_bus = Generator_data_dic[value]["Location_Bus"]
    Bus_region = Bus_data_dic[Location_bus]["Bus_Region"]
    println("Generator ", value, " in Region ", Bus_region)
    push!(Gen_Region_links, (value, Bus_region))
end

# Storage to bus links
Storage_Bus_links = Set()
for value in UStorage
    Location_bus = Utility_storage_data_dic[value]["Location_Bus"]
    println("Utility Storage ", value, " connected to Bus ", Location_bus)
    push!(Storage_Bus_links, (value, Location_bus))
end

# Generator type 1 set
GenT1 = Set()
for value in UGen
    Generation_type = Generator_data_dic[value]["Generation_Type"]
    if Generation_type == 1
        println("Generator ", value, " is Type ", Generation_type)
        push!(GenT1, (value))
    end
end

# Define GenT2 for wind and solar generators
GenT2 = Set()
for value in UGen
    Generation_type = Generator_data_dic[value]["Generation_Type"]
    if Generation_type == 2  # Wind and Solar
        println("Generator ", value, " is Type ", Generation_type)
        push!(GenT2, (value))
    end
end

# Generator type 1 to region links
GenT1_Bus_links = Set()
GenT1_Region_links = Set()
for value in UGen
    Generation_type = Generator_data_dic[value]["Generation_Type"]
    if Generation_type == 1
        Location_bus = Generator_data_dic[value]["Location_Bus"]
        Bus_region = Bus_data_dic[Location_bus]["Bus_Region"]
        # println("Type 1 Generator ", value, " is in Region ", Bus_region)
        push!(GenT1_Bus_links, (value, Location_bus))
        push!(GenT1_Region_links, (value, Bus_region))
    end
end

# Generator type 2 to region links
GenT2_Bus_links = Set()
GenT2_Region_links = Set()
for value in UGen
    Generation_type = Generator_data_dic[value]["Generation_Type"]
    if Generation_type == 2
        Location_bus = Generator_data_dic[value]["Location_Bus"]
        Bus_region = Bus_data_dic[Location_bus]["Bus_Region"]
        # println("Type 2 Generator ", value, " is in Region ", Bus_region)
        push!(GenT2_Bus_links, (value, Location_bus))
        push!(GenT2_Region_links, (value, Bus_region))
    end
end

# Generator type 1 to region links
GenT1_Region_links = Set()
for value in UGen
    Generation_type = Generator_data_dic[value]["Generation_Type"]
    if Generation_type == 1
        Location_bus = Generator_data_dic[value]["Location_Bus"]
        Bus_region = Bus_data_dic[Location_bus]["Bus_Region"]
        println("Type 1 Generator ", value, " is in Region ", Bus_region)
        push!(GenT1_Region_links, (value, Bus_region))
    end
end


# Generator type 2 to region links
GenT2_Region_links = Set()
for value in UGen
    Generation_type = Generator_data_dic[value]["Generation_Type"]
    if Generation_type == 2
        Location_bus = Generator_data_dic[value]["Location_Bus"]
        Bus_region = Bus_data_dic[Location_bus]["Bus_Region"]
        println("Type 2 Generator ", value, " is in Region ", Bus_region)
        push!(GenT2_Region_links, (value, Bus_region))
    end
end

# Region links
URegion = Set()
for value in UBus
    Bus_region = Bus_data_dic[value]["Bus_Region"]
    push!(URegion, Bus_region)
end

# Bus to region links
Bus_Region_links = Set()
for value in UBus
    Bus_region = Bus_data_dic[value]["Bus_Region"]
    println("Bus ", value, " is in Region ", Bus_region)
    push!(Bus_Region_links, (value, Bus_region))
end
Bus_Region_links_orig = Bus_Region_links # Save original Bus_Region_links for use in different network model details

# Line end bus links
Line_end1_Bus_links = Dict{String,String}()
for value in ULine
    End_1_Bus = Line_data_dic[value]["End_1_Bus"]
    Line_end1_Bus_links[value] = End_1_Bus
    println("Line ", value, " is connected to Bus ", End_1_Bus)
end

Line_end2_Bus_links = Dict{String,String}()
for value in ULine
    End_2_Bus = Line_data_dic[value]["End_2_Bus"]
    Line_end2_Bus_links[value] = End_2_Bus
    println("Line ", value, " is connected to Bus ", End_2_Bus)
end


# Bus to prosumer links
Bus_Prosumer_links = Set()
for value in UBus
    Bus_prosumer = Bus_data_dic[value]["Prosumer_Demand_p"]
    println("Bus ", value, " has prosumer demand ", Bus_prosumer)
    if Bus_prosumer != 0 # Gets only buses with prosumer demand
        push!(Bus_Prosumer_links, value) 
    end
end

# Generator to plot color links
Gen_Tech_links = Set()
for value in UGen
    Gen_tech = Generator_data_dic[value]["Generation_Tech"]
    Plot_color = Generator_data_dic[value]["Plot_Color"]
    # println("Generator: ", value, " Technology: ", Gen_tech)
    push!(Gen_Tech_links, (value, Gen_tech)) 
end

# Generator to plot color links
Gen_Color_links = Set()
for value in UGen
    Plot_color = Generator_data_dic[value]["Plot_Color"]
    push!(Gen_Color_links, (value, Plot_color))
end

# Storage to plot color links
Storage_Color_links = Set()
for value in UStorage
    Plot_color = Utility_storage_data_dic[value]["Plot_Color"]
    push!(Storage_Color_links, (value, Plot_color))
end

# Technology to plot color links
Tech_Color_links = Set()
for value in UGen
    Gen_tech = Generator_data_dic[value]["Generation_Tech"]
    Plot_color = Generator_data_dic[value]["Plot_Color"]
    # println("Generator: ", value, " Technology: ", Gen_tech, " PlotColor: ", Plot_color)
    push!(Tech_Color_links, (value, Gen_tech, Plot_color)) 
end

# Initial values

# Global technical parameters
SOC_psm_ini = 0.0 # Initial state of charge of batteries (Prosumer)
SOC_utl_ini = 0.0 # Initial state of charge of batteries (Utility-scale)

# Generation intial values
Slack_bus = determine_slack_bus(Generator_data_dic, Bus_data_dic, UGen, UBus) # Calls function to determine slack bus
Base_power = Bus_data_dic[Slack_bus]["Base_kV"] # Determine base power from slack bus
Pwr_Gen_ini_v = zeros(length(UGen), 1)
Status_ini_v = zeros(length(UGen), 1)
S_Down_ini_v = zeros(length(UGen), 1)
MUT_ini_v = zeros(length(UGen), T)
MDT_ini_v = zeros(length(UGen), T)
Status_ini = Dict()
Pwr_Gen_ini = Dict()
S_Down_ini = Dict()
MUT_ini = Dict()
MDT_ini = Dict()

for (index, value) in enumerate(UGen)
    Status_ini[value] = Status_ini_v[index]
    Pwr_Gen_ini[value] = Pwr_Gen_ini_v[index]
    S_Down_ini[value] = S_Down_ini_v[index]
end
for (index, value) in enumerate(UGen)
    for t in subhorizon_total
        MUT_ini[(value, t)] = MUT_ini_v[index, t]
        MDT_ini[(value, t)] = MDT_ini_v[index, t]
    end
end

# Utility storage initial values
Chrg_rate_strg = Dict(value => Utility_storage_data_dic[value]["Maximum_Charge_Rate_MWh"] for value in UStorage)
Dchrg_rate_strg = Dict(value => Utility_storage_data_dic[value]["Maximum_Discharge_Rate_MWh"] for value in UStorage)
Min_SOC_strg = Dict(value => Utility_storage_data_dic[value]["Minimum_Storage_Capacity_MWh"] for value in UStorage)
Max_SOC_strg = Dict(value => Utility_storage_data_dic[value]["Maximum_Storage_Capacity_MWh"] for value in UStorage)
Charging_eff = Dict(value => Utility_storage_data_dic[value]["Charging_Efficiency"] / 100 for value in UStorage)
Discharging_eff = Dict(value => Utility_storage_data_dic[value]["Discharging_Efficiency"] / 100 for value in UStorage)
Enrg_Strg_ini = Dict(value => SOC_utl_ini for value in UStorage)

# Prosumer storage initial values
psm_batt_eff = Dict(value => Bus_data_dic[value]["Battery_Efficiency"]/100 for value in UBus_orig)
psm_batt_maxcap = Dict(value => Bus_data_dic[value]["Maximum_Battery_Capacity_MWh"] for value in UBus_orig)
psm_batt_mincap = Dict(value => Bus_data_dic[value]["Minimum_Battery_Capacity_MWh"] for value in UBus_orig)
psm_batt_maxdischarge = Dict(value => Bus_data_dic[value]["Maximum_Discharge_Rate"] for value in UBus_orig)
psm_batt_maxcharge = Dict(value => Bus_data_dic[value]["Maximum_Charge_Rate"] for value in UBus_orig)
psm_feedin_priceratio = Dict(value => Bus_data_dic[value]["Feedin_Price_Ratio"] for value in UBus_orig)
psm_batt_init = Dict(value => SOC_psm_ini for value in UBus_orig)


# Initialize dictionaries for passing states
prev_Status_var = Dict(g => Status_ini[g] for g in UGen)
prev_SOC = Dict(s => Enrg_Strg_ini[s] for s in UStorage)
prev_Pwr_Gen_var = Dict(g => Pwr_Gen_ini[g] for g in GenT1)
prev_psm_batt_egy = Dict(n => SOC_psm_ini for n in UBus_orig)
prev_S_Up_var = Dict(g => 0 for g in UGen)  # Track previous startups
prev_S_Down_var = Dict(g => 0 for g in UGen)  # Track previous shutdowns

# Expansion
existing_lines = [l for l in ULine if Line_data_dic[l]["In_Service"] == 1]
println("Existing lines: ", existing_lines)


# Define the network model

if network_model in Cu_plate
    # Redefine UBus as a single "system" bus
    global UBus = Set(["System"])
    global N = length(UBus)
    global Slack_bus = "System"
    println("Copper plate aggregating all buses into system bus " , UBus)
    println("Number of buses: ", N)

    # Map all generators to the single "System" bus
    global GenT1_Bus_links = Set()
    for value in GenT1
        push!(GenT1_Bus_links, (value, "System"))
        println("Mapped generator $value to system bus")
    end
    # println("GenT1_Bus_links for Cu_plate: ", GenT1_Bus_links)

    global GenT2_Bus_links = Set()
    for value in GenT2
        push!(GenT2_Bus_links, (value, "System"))
        println("Mapped generator $value to system bus")
    end
    # println("GenT2_Bus_links for Cu_plate: ", GenT2_Bus_links)

    # Map all storage to the single "System" bus
    global Storage_Bus_links = Set()
    for value in UStorage
        push!(Storage_Bus_links, (value, "System"))
        println("Mapped storage $value to system bus")
    end
    # println("Storage_Bus_links for copper plate: ", Storage_Bus_links)

    # No transmission lines in copper plate
    global ULine = Set{String}()  # Empty set
    global URegion = copy(UBus)
    global Line_end1_Bus_links = Dict{String,String}()
    global Line_end2_Bus_links = Dict{String,String}()
    global Bus_Region_links = Set((b, "System") for b in UBus)
    global GenT1_Region_links = Set((g, "System") for g in GenT1 for gb in GenT1_Bus_links if gb[1] == g)
    global GenT2_Region_links = Set((g, "System") for g in GenT2 for gb in GenT2_Bus_links if gb[1] == g)
    println("Redefined ULine as empty for copper plate: ", ULine)
    println("Line_end1_Bus_links: ", Line_end1_Bus_links)
    println("Line_end2_Bus_links: ", Line_end2_Bus_links)
    println("Bus_Region_links: ", Bus_Region_links)
    println("GenT1_Region_links: ", GenT1_Region_links)
    println("GenT2_Region_links: ", GenT2_Region_links)
    println("Existing lines: ", existing_lines)

    global psm_batt_eff["System"] = maximum(psm_batt_eff[b] for b in UBus_orig if haskey(psm_batt_eff, b); init=0.0)
    global psm_batt_maxcap["System"] = sum(psm_batt_maxcap[b] for b in UBus_orig if haskey(psm_batt_maxcap, b); init=0.0)
    global psm_batt_mincap["System"] = sum(psm_batt_mincap[b] for b in UBus_orig if haskey(psm_batt_mincap, b); init=0.0)
    global psm_batt_maxdischarge["System"] = sum(psm_batt_maxdischarge[b] for b in UBus_orig if haskey(psm_batt_maxdischarge, b); init=0.0)
    global psm_batt_maxcharge["System"] = sum(psm_batt_maxcharge[b] for b in UBus_orig if haskey(psm_batt_maxcharge, b); init=0.0)
    global psm_feedin_priceratio["System"] = maximum(psm_feedin_priceratio[b] for b in UBus_orig if haskey(psm_feedin_priceratio, b); init=0.0)
    global psm_batt_init["System"] = maximum(psm_batt_init[b] for b in UBus_orig if haskey(psm_batt_init, b); init=0.0)
    global prev_psm_batt_egy = Dict(n => SOC_psm_ini for n in UBus)

elseif network_model in Regional
    # Redefine UBus as URegion
    global UBus = copy(URegion)
    global N = length(UBus)
    global Slack_bus = Bus_data_dic[Slack_bus]["Bus_Region"]
    println("Regional model aggregating buses per region ", UBus)
    println("Number of regions: ", N)

    # Map generators to regions
    global GenT1_Bus_links = Set()
    for value in GenT1
        bus = Generator_data_dic[value]["Location_Bus"]
        region = Bus_data_dic[bus]["Bus_Region"]
        push!(GenT1_Bus_links, (value, region))
        println("Mapped generator $value from bus $bus to region $region")
    end
    # println("GenT1_Bus_links for Regional: ", GenT1_Bus_links)

    global GenT2_Bus_links = Set()
    for value in GenT2
        bus = Generator_data_dic[value]["Location_Bus"]
        region = Bus_data_dic[bus]["Bus_Region"]
        push!(GenT2_Bus_links, (value, region))
        println("Mapped generator $value from bus $bus to region $region")
    end
    # println("GenT2_Bus_links for Regional: ", GenT2_Bus_links)


    # Map storage to regions
    global Storage_Bus_links = Set()
    for value in UStorage
        bus = Utility_storage_data_dic[value]["Location_Bus"]
        region = Bus_data_dic[bus]["Bus_Region"]
        push!(Storage_Bus_links, (value, region))
        println("Mapped storage $value from bus $bus to region $region")
    end
    # println("Storage_Bus_links for Regional: ", Storage_Bus_links)
    
    global Bus_Region_links = Set{String}() # empty set

    # Define inter-regional lines (ULine_region)
    global ULine_region = Set(l for l in ULine if Bus_data_dic[Line_data_dic[l]["End_1_Bus"]]["Bus_Region"] != Bus_data_dic[Line_data_dic[l]["End_2_Bus"]]["Bus_Region"])
    global ULine = ULine_region  # Replace ULine with regional subset
    global Line_end1_Bus_links = Dict(l => Bus_data_dic[Line_data_dic[l]["End_1_Bus"]]["Bus_Region"] for l in ULine)
    global Line_end2_Bus_links = Dict(l => Bus_data_dic[Line_data_dic[l]["End_2_Bus"]]["Bus_Region"] for l in ULine)
    global existing_lines = [l for l in ULine if Line_data_dic[l]["In_Service"] == 1]

    println("Redefined ULine as ULine_region (inter-regional lines): ", ULine)
    println("Line_end1_Bus_links: ", Line_end1_Bus_links)
    println("Line_end2_Bus_links: ", Line_end2_Bus_links)
    println("Existing inter-regional lines: ", existing_lines)

    for r in URegion
        region_buses = [b for (b, reg) in Bus_Region_links_orig if reg == r]

        # Aggregate prosumer variables for the region
        global psm_batt_eff[r] = maximum(psm_batt_eff[b] for b in region_buses if haskey(psm_batt_eff, b); init=0.0)
        global psm_batt_maxcap[r] = sum(psm_batt_maxcap[b] for b in region_buses if haskey(psm_batt_maxcap, b); init=0.0)
        global psm_batt_mincap[r] = sum(psm_batt_mincap[b] for b in region_buses if haskey(psm_batt_mincap, b); init=0.0)
        global psm_batt_maxdischarge[r] = sum(psm_batt_maxdischarge[b] for b in region_buses if haskey(psm_batt_maxdischarge, b); init=0.0)
        global psm_batt_maxcharge[r] = sum(psm_batt_maxcharge[b] for b in region_buses if haskey(psm_batt_maxcharge, b); init=0.0)
        global psm_feedin_priceratio[r] = maximum(psm_feedin_priceratio[b] for b in region_buses if haskey(psm_feedin_priceratio, b); init=0.0)
        global psm_batt_init[r] = maximum(psm_batt_init[b] for b in region_buses if haskey(psm_batt_init, b); init=0.0)

        println("Aggregated prosumer variables for region $r: eff=$(psm_batt_eff[r]), maxcap=$(psm_batt_maxcap[r])")
    end
end

