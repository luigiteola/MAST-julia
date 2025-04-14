# This code manages the Type 2 Generator (Solar PV/Wind) data for use in the optimization

# Initialize Resource_trace_T2 for both solar and wind generators
Resource_trace_T2 = Dict()

# Loop over Type 2 generators (wind and solar)
for gen_name in GenT2
    # Find the bus the generator is located at
    bus_name = Generator_data_dic[gen_name]["Location_Bus"]

    # Check if the generator is solar (PV) or wind (WND) based on its technology type
    gen_type = Generator_data_dic[gen_name]["Generation_Tech"]

    if gen_type in Solar_Tech
        # Solar generator: Get the solar trace and scale by maximum power rating
        for t in Time
            trace_value = Solar_trace_DR[(bus_name, t)]
            Resource_trace_T2[(gen_name, t)] = trace_value
        end
    elseif gen_type in Wind_Tech
        # Wind generator: Get the wind trace and scale by maximum power rating
        for t in Time
            trace_value = Wind_trace_DR[(bus_name, t)]
            Resource_trace_T2[(gen_name, t)] = trace_value
        end
    else
        # Unsupported generator type
        println("Unsupported technology for generator: ", gen_name)
    end
end

# Print Resource_trace_T2 for verification
# println("Resource_trace_T2:")
# for (key, value) in Resource_trace_T2
#     println("Generator: ", key[1], ", Time: ", key[2], " => Resource Trace: ", value)
# end



