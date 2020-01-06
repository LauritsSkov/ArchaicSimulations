#!/usr/bin/python
import msprime as msp
import numpy as np
import random



'''
Files you will need:

    recombination map (hg19)
    bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip




'''


# Generation time, mutation rate and recomination rate
gen_time = 29.0 
rec_rate = 1.2e-8 
mu = 1.2e-8
chrom = 19
Overall_scenario = 'test'


# Populations sizes (number of individuals)
Ne_Africa = 27122
Ne_Europe = 5000
Modernhumans = 7000
Denisovasize = 5000
Neanderthal_recentsize = 1000
Neanderthal_latersize = 2000
Neanderthal_Denisovasize = 5000

# Split times (years)
OutOfafrica = 62000/gen_time
Denisova_split = 350000 
IntrogressingVindijaSplit = 90000
DenisovaNeanderthal = 420000

# Admix parameters (years and admixture proportions (in percent))
Denisova_admix_time = 45000 
DenisovaProportion = 0.08
admix_into = 1 # (1 for humans, 5 for neanderthals)


# Number of samples
n_ingroup = 2
African_samples = 1000
samples = [msp.Sample(0, 0)]*African_samples + [msp.Sample(1, 0)]*n_ingroup + [msp.Sample(3, 80000/gen_time)]*n_ingroup + [msp.Sample(4, 120000/gen_time)]*n_ingroup + [msp.Sample(6, 60000/gen_time)]*n_ingroup

population_configurations = [
    msp.PopulationConfiguration(initial_size = Ne_Africa), #0
    msp.PopulationConfiguration(initial_size = Ne_Europe), #1
    msp.PopulationConfiguration(initial_size = Denisovasize), #2
    msp.PopulationConfiguration(initial_size = Denisovasize), #3
    msp.PopulationConfiguration(initial_size = Neanderthal_recentsize), #4
    msp.PopulationConfiguration(initial_size = Neanderthal_recentsize), #5
    msp.PopulationConfiguration(initial_size = Neanderthal_recentsize), #6
]


demographic_events_dict = {


    # admixture times
    1551.06896552 + random.uniform(0,1)/100000.0: msp.MassMigration(time = 1551.06896552, source = 1, destination = 5,proportion = 0.02),
    Denisova_admix_time/gen_time + random.uniform(0,1)/100000.0: msp.MassMigration(time = Denisova_admix_time/gen_time, source = admix_into, destination = 2 ,proportion = DenisovaProportion),

    # Human parameters
    OutOfafrica - 300 + random.uniform(0,1)/100000.0: msp.PopulationParametersChange(time = OutOfafrica - 300, initial_size = 1305, growth_rate = 0, population_id = 1),
    OutOfafrica - 200 + random.uniform(0,1)/100000.0: msp.PopulationParametersChange(time = OutOfafrica - 200, initial_size = 5000, growth_rate = 0, population_id = 1),
    OutOfafrica - 100 + random.uniform(0,1)/100000.0: msp.PopulationParametersChange(time = OutOfafrica - 100, initial_size = 250, growth_rate = 0, population_id = 1),
    OutOfafrica + random.uniform(0,1)/100000.0: msp.MassMigration(time = OutOfafrica, source = 1, destination = 0, proportion = 1.0),
    OutOfafrica + 0.0001 + random.uniform(0,1)/100000.0: msp.PopulationParametersChange(time = OutOfafrica  + 0.0001, initial_size = Modernhumans, growth_rate = 0, population_id = 0),

    # archaic population merges parameters
    IntrogressingVindijaSplit/gen_time + random.uniform(0,1)/100000.0: msp.MassMigration(time = IntrogressingVindijaSplit/gen_time, source = 6, destination = 5, proportion = 1.0),
    130000/gen_time - 0.0001 + random.uniform(0,1)/100000.0: msp.MassMigration(time = 130000/gen_time - 0.0001, source = 5, destination = 4, proportion = 1.0), 
    Denisova_split/gen_time - 0.0001 + random.uniform(0,1)/100000.0: msp.MassMigration(time = Denisova_split/gen_time, source = 3, destination = 2, proportion = 1.0),
    DenisovaNeanderthal/gen_time - 0.0001 + random.uniform(0,1)/100000.0: msp.MassMigration(time = DenisovaNeanderthal/gen_time, source = 4, destination = 2, proportion = 1.0),


    # psms pop sizes
    130000/gen_time + random.uniform(0,1)/100000.0: msp.PopulationParametersChange(time = 130000/gen_time + 0.0001, initial_size = Neanderthal_latersize, growth_rate = 0, population_id = 4),


    # Neanderthal and denisova merges and have a bigger pop size
    DenisovaNeanderthal/gen_time + random.uniform(0,1)/100000.0: msp.PopulationParametersChange(time = DenisovaNeanderthal/gen_time, initial_size = Neanderthal_Denisovasize, growth_rate = 0, population_id = 2),


    # Humans and archaics are the same population
    580000/gen_time + random.uniform(0,1)/100000.0: msp.MassMigration(time = 580000/gen_time, source = 2, destination = 0,proportion = 1.0),
    580000/gen_time + 0.0001 + random.uniform(0,1)/100000.0: msp.PopulationParametersChange(time = 580000/gen_time + 0.0001, initial_size = 7000, growth_rate = 0, population_id = 0),


    }


demographic_events = []

for i, key in enumerate(sorted(demographic_events_dict)):
    demographic_events.append(demographic_events_dict[key])



# --------------------------------------------------------------------------------------
# Simulate
# --------------------------------------------------------------------------------------

# Run simulations
map_file = '/faststorage/home/laurits2/Denisova_introgression/helperfiles/Recombination_map/genetic_map_GRCh37_chr{}.txt'.format(chrom) #"path to recombination map"
recomb_map = msp.RecombinationMap.read_hapmap(map_file)


ts = msp.simulate(
    samples = samples,
    population_configurations = population_configurations,
    demographic_events = demographic_events,
    recombination_map = recomb_map,  
    #length = 5000000, recombination_rate = rec_rate,      
    mutation_rate = mu
)


#--------------------------------------------------------------------------------------
# Get archaic segments
#--------------------------------------------------------------------------------------
print ('removing outgroup')

Testpopulation = ts.get_samples(1)
Outgroup = ts.get_samples(0)
Archaics = np.concatenate((ts.get_samples(3), ts.get_samples(4), ts.get_samples(6)))



# Write genotype output (only keep sites which are not varying in Africa)
print('\t'.join(['chrom','pos','human haplotype','Denisova','Altai','Vindija']))
for variant in ts.variants():

    pos = str(variant.site.position)

    if np.sum(variant.genotypes[Outgroup]) == 0:

        if variant.genotypes[Testpopulation][0] == 1:
            print('\t'.join([str(chrom), pos, '{0}|{1}'.format(variant.genotypes[Testpopulation][0],variant.genotypes[Testpopulation][1]), '{0}|{1}'.format(variant.genotypes[Archaics][0],variant.genotypes[Archaics][1]), '{0}|{1}'.format(variant.genotypes[Archaics][2],variant.genotypes[Archaics][3]), '{0}|{1}'.format(variant.genotypes[Archaics][4],variant.genotypes[Archaics][5])]))
