# Hi-BDiSCO

#### System Requirement ####

The binary executables (under bin/ directory) are compiled to run on the popular Linux CentOS 7, RedHat/Roky 7 and 8, Ubuntu 20 and 22 platforms. They work with both Intel and AMD CPUs. 

#### Dependencies ####

1) MC simulation:

The binary bin/chrom_ap1.x for MC simulation requires OpenMPI --version >= 4

2) BD simulation:

The binary bin/code for BD simulation requires 

3) Initial struction generation:

Python scripts under src/ requires following libs:

pyBigWig 0.3.18  
cooler 0.8.11
numpy 1.18.5  

#### Data Input Preparation ####

1) Experimental data:

The experimental data are under data/ directory (example data are provided).

For the example data:

-> nucleosome_position.bed is a bed file generated by DANPOS using MNase-seq data that contains the nucleosome positioning information.

-> tail_acetylation.narrowPeak is a narrowPeak file generated by MACS from Chip-seq data that contains the tail acetylation information.

-> LH.bw is a bigwig file of Chip-seq data that contains the linker histone information.

-> Hi-C-sample.mcool is a mcool file contains Micro-C contact data.

The users can replace these data with their own data of interest.

2) Inputs:

The inputs can be modified via input.txt:

For the default example:

chr = chr18                                                           -> chromosome number
start = 58110000                                                      -> chromosome start
end = 58150000                                                        -> chromosome end
nucleosome_position = data/nucleosome_position.bed                    -> nucleosome position file location (optional, if not provided, life-like fiber will be generated)
tail_acetylation = data/tail_acetylation.narrowPeak                   -> tail acetylation file location (optional, if not provided, wildtype tail will be used)
LH = data/LH.bw                                                       -> linker histone file location (optional, if not provided, random occupancy of LH will be assigned based on the LH density)
LH_ratio = 0.5                                                        -> linker histone density
Hi-C = data/Hi-C-sample.mcool                                         -> Hi-C-like map location
Hi-C-path = /resolutions/50                                           -> path name of the Hi-C-like map file
N_rep = 1000                                                          -> number of replicas to distribute restraints
N_sim = 100                                                           -> number of simulations to be performed

#### Run simulation ####

After modify the data inputs, the users can run the simulation by execute the following shell scripts:

1) ./file_preparation.sh                                              -> Prepare the input files for MC and BD code based on the provided experimental files

2) ./run_part1_initial_structure_generation.sh                        -> Perform MC simulations to generate N_sim numbers of random initial structures (submitting HPC jobs, output/org_sys/mc/sub_batch.s should be modified accordingly)

3) ./run_part2_BD_reconstruction.sh                                   -> Perform BD simulations to fold the target region of chromosome based on the Hi-C restraints (submitting HPC jobs, output/org_sys/bd/run-code.sbatch should be modified accordingly)

4) ./run_part3_MC_subsequent_simulation.sh                            -> Perform MC simulations to taking account the affects of tail acetylation, linker histion, etc., and solve spatial problems (submitting HPC jobs, output/org_sys/mc/sub_batch.s should be modified accordingly)



