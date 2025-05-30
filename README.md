This is a Market Simulation Tool (MAST) coded in the Julia Programming Language.
See MAST_Julia_guide_v1.0.pdf for installation and initialization instructions.


Added functionality for v1.4
1. MAST_debug.jl - new file for debugging purposes. Same as MAST.jl but without the command line prompts. Open the file and run/execute using your IDE.
2. Improved handling of output CSV file - instead of a single column, now outputs a table of decision variables at each time step.
3. New Demand, Solar, and Wind traces added - 'Demand/Solar/Wind Traces Hourly Resolution' Folders
4. New Model added - 'NEM_14gen_59bus.xlsx'
5. Unserved Demand added as slack variable - for test cases when supply does not meet demand to avoid infeasibility. It is multiplied by 'voll' set in config.json as a penalty.
6. Curtailment added as slack variable for VRE - for cases when there is oversupply of VRE, optimization can curtail excess. It is multiplied by 'curtailment_penalty' set in config.json as a penalty.
7. Minor improvements in system and cost summary outputs via terminal
8. Bug fixes in storage initialization


Added functionality for v1.3
Network model detail
Determine network model detail in config.json
1. Copper plate - single bus model
2. Regional - buses grouped per region
3. Nodal (default) - all buses and branches considered

Allowed values for "network_model"
 1. Copper plate = "Cu_plate", "Cu Plate" ,"Copper Plate", "Copperplate", "Copper plate", "Single-bus", "Single Bus"
 2. Regional = "Regional", "Region", "Zonal", "Zone"
 3. Nodal = "Nodal", "Node", "Full Network", "Whole Network", "Full"


