
Monte Carlo eXtreme (MCXLAB) - CUDA - with Brillouin scattering
=========================
This is a Master's thesis project at the University of Maryland Optics Biotech Lab to add Brillouin scattering behavior to the existing CUDA MCXLAB code. For information about MCX refer to README_MCX.txt. The thesis documenting this work is titled "Monte Carlo Simulations of Brillouin Scattering in Turbid Media" by Stephanie Lashley. 

-   Author: Stephanie Lashley (stephanielashley13@gmail.com)

- For example files refer to guassian_source_basic.m and gaussian_source_DOT.m
- Usage: [flux,detpt,vol,seed,trajectory,~,all_brillouin_data]=mcxlab(cfg); 

- The detected Brillouin data is stored in detpt.data(end-1,:), this is the cosine of the Brillouin scattering angle
- The Brillouin shift of every simulated photon (not just detected) is stored in all_brillouin_data
- When multiple materials are simulated, the material IDs are stored in detpt.data(end,:) to calculated the Brillouin shift for each material type
- Do not edit cfg.savedetflags
- The ~ output is a legacy variable from development that was included as an example of how to add an additional output to this code. It is called test_var in mcx_core.cu. 

