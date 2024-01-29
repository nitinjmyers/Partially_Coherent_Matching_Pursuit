# Partially_Coherent_Matching_Pursuit
This algorithm called Partially Coherent Matching Pursuit (PC-MP) is developed to perform sparse recovery, i.e., compressive sensing (CS), under phase errors. These phase errors are constant for each batch of CS measurements and vary across the batches.

Description of files and steps to execute our algorithm:
- "PC_MP.m" contains an implementation of the PC-MP algorithm discussed in our paper [1].
- "comparison.m" is the main file that calls our PC-MP function and also the benchmarks.
   Set "simulation_id=1" in line 3 to get performance with SNR, "simulation_id=2" for performance with CS measurements, "simulation_id=3" for performance with sparsity 
   level.
   The variable "aveNumber" denotes number of realizations for which the performance is averaged. For the results in [1], we used  "aveNumber=400".
  
[1] Weijia Yi, Nitin Jonathan Myers, Geethu Joseph, "Sparse Millimeter Wave Channel Estimation From Partially Coherent Measurements," to appear in Proc. of the IEEE Intl. Conf. on Communications 2024. 

