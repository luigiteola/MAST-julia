This is a Market Simulation Tool (MAST) coded in the Julia Programming Language.
See MAST_Julia_guide_v1.0.pdf for installation and initialization instructions.


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
