from os.path import join

# Constants
# TODO: remove and pass as tool arguments

DEFAULT_SAMPLING_SIZE = 10000
TABIX = "/soft/bio/sequence/tabix-0.2.3/tabix"

SIGNATURES = {
    'pilot': '/projects_bg/bg/shared/projects/fmdrivers/input/signature_probabilities/signature_probabilities_pilotproject.tsv',
    'tcga': '/projects_bg/bg/shared/projects/fmdrivers/input/signature_probabilities/signature_probabilities_tcga.tsv',
    'cll': '/projects_bg/bg/shared/projects/fmdrivers/input/signature_probabilities/signature_probabilities_cll.tsv',
    'pwcag': '/projects_bg/bg/shared/projects/fmdrivers/input/signature_probabilities/signature_probabilities_pwcag.tsv'
}

HG19_DIR = "/projects_bg/bg/soft/intogen_home/gencluster/software/mutsigCV/reffiles/chr_files_hg19"

SAMPLING_INPUT_FOLDER = "/projects_bg/bg/shared/projects/fmdrivers/input"

SCORES = {
    'cadd': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'cadd', 'whole_genome_SNVs.tsv.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 5, 'element': None
    },
    'cadd_11': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'cadd_11', 'whole_genome_SNVs.tsv.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 5, 'element': None
    },
    'cadd_12': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'cadd_12', 'whole_genome_SNVs.tsv.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 5, 'element': None
    },
    'funseq2': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'funseq2', 'hg19_wg_score.tsv.gz'),
        'chr': 0, 'chr_prefix': 'chr', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 4, 'element': None
    },
    'rnasnp': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'rnasnp', 'hg19_rnasnp_scores.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 5, 'element': 6
    },
    'tfbs_creation': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'tfbs_creation', 'tfbs_creation.tsv.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 4, 'element': None
    },
    'tfbs_disruption': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'tfbs_disruption', 'tfbs_disruption.tsv.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 2, 'alt': 3, 'score': 4, 'element': None
    },
    'mirnats_bladder': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 6, 'element': 5
    },
    'mirnats_blood': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 7, 'element': 5
    },
    'mirnats_brain': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 8, 'element': 5
    },
    'mirnats_breast': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 9, 'element': 5
    },
    'mirnats_colon': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 10, 'element': 5
    },
    'mirnats_heart': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 11, 'element': 5
    },
    'mirnats_kidney': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 12, 'element': 5
    },
    'mirnats_liver': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 13, 'element': 5
    },
    'mirnats_lung': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 14, 'element': 5
    },
    'mirnats_ovary': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 15, 'element': 5
    },
    'mirnats_pancreas': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 16, 'element': 5
    },
    'mirnats_prostate': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 17, 'element': 5
    },
    'mirnats_spleen': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 18, 'element': 5
    },
    'mirnats_stomach': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 19, 'element': 5
    },
    'mirnats_testicle': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 20, 'element': 5
    },
    'mirnats_uterus': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 21, 'element': 5
    },
    'mirnats_all': {
        'file': join(SAMPLING_INPUT_FOLDER, 'scores', 'mirnats', 'hg19_mirnaTS_score.txt.gz'),
        'chr': 0, 'chr_prefix': '', 'pos': 1, 'ref': 3, 'alt': 4, 'score': 22, 'element': 5
    }
}

class Configuration():

    MUTS_FILE_NAME = "mutations.json.gz"
    INDELS_FILE_NAME = "indels.txt"
    SAMPLING_RESULTS_FILE_NAME = "sampling_results_{}.json.gz"
    PVALUES_FILE_NAME = 'pvalues_{}_{}.tsv.gz'

    # Define indel score strategies
    INDEL_STRATEGY_PERCENTILES = 'percentiles'
    INDEL_STRATEGY_MAX = 'max'
    INDEL_STRATEGY_USER_INPUT = 'user-input'
    INDEL_STRATEGIES = [INDEL_STRATEGY_PERCENTILES, INDEL_STRATEGY_USER_INPUT, INDEL_STRATEGY_MAX]


    def __init__(self, project, tissue, score, feature, output_dir, cache_dir, input_dir_or_file, indel_strategy=None,
                 cores=4, testing=False, expression_file=None, indel_file=None, signature="",
                 ensembl_file="/projects_bg/bg/shared/datasets/non-coding/ensemble_genes_75.txt"):

        self.expression_file = expression_file
        self.testing = testing
        self.input = input_dir_or_file
        self.project = project
        self.tissue = tissue
        self.score = score
        self.feature = feature
        self.output_dir = output_dir
        self.cache_dir = cache_dir
        self.ensembl_file = ensembl_file
        self.cores = cores
        self.indel_strategy = indel_strategy
        self.indel_file = indel_file
        self.signature = self.create_signature(signature)

    def create_signature(self, signature):
        if signature == "":
            return "{}_{}".format(self.project, self.tissue) if self.tissue != 'pan' else 'none'
        else:
            return signature
