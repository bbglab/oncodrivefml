import json
import gzip
import glob
import os

MAIN_FOLDER = "/shared/projects/fmdrivers"
mutations_files = glob.glob("{}/output/*/*/*/mutations.tsv.gz".format(MAIN_FOLDER))
for mutations_file in mutations_files:
    project, cancer_type, feature, file = mutations_file.replace(MAIN_FOLDER + "/output/", '').split(os.sep)


    mutations = json.load(gzip.open(mutations_file, 'rt'))
    indels_file = os.path.join(MAIN_FOLDER, 'output', project, cancer_type, feature, "indels.txt")
    if os.path.exists(indels_file):
        print("{} - {} - {} [skip]".format(project, cancer_type, feature))
        continue
    else:
        print("{} - {} - {} [process]".format(project, cancer_type, feature))

    with open(indels_file, 'w') as fd:
        for element, muts in mutations.items():
            indels = list(filter(lambda v: v['TYPE'] == 'indel', muts))
            if len(indels) > 0:
                for i in indels:
                    print("{}\t{}\t{}\t{}\t{}\t{}".format(i['CHROMOSOME'], i['POSITION'], i['REF'], i['ALT'], i['SAMPLE'], element), file=fd)