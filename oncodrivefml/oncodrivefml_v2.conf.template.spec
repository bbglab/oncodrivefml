[genome]
build = option('hg18', 'hg19', 'hg38', 'mm10', 'c3h', default='hg19')

[signature]
method = option('none', 'full', 'complement', 'bysample', 'file', default='full')

classifier = option('CANCER_TYPE', 'SAMPLE', 'SIGNATURE', default='CANCER_TYPE')
include_mnp = boolean(default=True)
normalize_by_sites = option('whole_genome', 'wgs', 'whole_exome', 'wxs', 'wes', default=None)
only_mapped_mutations = boolean(default=False)

path = string(default=None)
column_ref = string(default=None)
column_alt = string(default=None)
column_probability = string(default=None)

[score]
file = string
format = option('tabix', 'pack')
chr = integer
chr_prefix = string
pos = integer
ref = integer(default=None)
alt = integer(default=None)
score = integer
element = integer(default=None)
extra = integer(default=None)

minimum_number_of_stops = integer(default=3)
mean_to_stop_function = string(default=None)


[statistic]
method = option('amean', 'gmean', default='amean')
discard_mnp = boolean(default=False)

sampling = integer(default=100000)
sampling_max = integer(default=1000000)
sampling_chunk = integer(default=100)
sampling_min_obs = integer(default=10)

per_sample_analysis = option('amean', 'gmean', 'max', default=None)

    [[indels]]
        include = boolean(default=True)
        method = option('stop', 'max')
        max_consecutive = integer(default=0)
        simulate_with_signature =  boolean(default=True)

        gene_exomic_frameshift_ratio = boolean(default=False)
        stops_function = option('mean', 'median', 'random', 'random_choice', default='mean')



[settings]
cores = integer(default=None)

[logging]
version = integer(default=1)