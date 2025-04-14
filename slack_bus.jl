function determine_slack_bus(Generator_data_dic::Dict, Bus_data_dic::Dict, UGen::Set, UBus::Set)
    # Check for Bus_Type == 3 among buses in UBus
    for bus in UBus
        data = Bus_data_dic[bus]
        bus_type = get(data, "Bus_Type", nothing)
        if !isnothing(bus_type) && bus_type == 3
            println("Slack bus selected: $bus (Bus_Type = 3 detected)")
            println("Base_kV of slack bus: ", Bus_data_dic[bus]["Base_kV"])
            return bus
        end
    end

    max_gen = nothing
    max_power = -Inf
    max_bus = nothing

    # Only consider generators in UGen (connected to grid)
    for gen in UGen
        data = Generator_data_dic[gen]
        power = coalesce(data["Maximum_Real_Power"], 0.0)
        if power > max_power
            max_power = power
            max_gen = gen
            max_bus = data["Location_Bus"]
        end
    end

    if max_bus !== nothing && max_bus in UBus
        slack_bus = max_bus
        println("\n\nNo Bus_Type = 3 found. Slack bus selected: $slack_bus (largest generator)")
        println("Connected generator: $max_gen with Maximum_Real_Power: $max_power MW")
        println("Base_kV of slack bus: ", Bus_data_dic[slack_bus]["Base_kV"])
        return slack_bus
    else
        slack_bus = first(UBus)
        println("\n\nNo Bus_Type = 3 or valid generators found. Slack bus arbitrarily set to: $slack_bus")
        return slack_bus
    end
end