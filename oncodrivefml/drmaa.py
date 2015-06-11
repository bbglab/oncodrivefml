import gzip
import logging
import os
import pickle


def drmaa_run(variants_dict, signature_dict, task, size, figures=True):

    try:
        import drmaa
    except RuntimeError:
        raise RuntimeError("It's not possible to import 'drmaa' python package. May be it's not installed or you"
                           "don't have the DRMAA_LIBRARY_PATH environment variable defined.")

    # Split variants file into several chunks
    variants_list = list(variants_dict.items())
    variants_list_split = [variants_list[i:i+size] for i in range(0, len(variants_list), size)]
    variants_dict_split = [{k: v for k, v in i} for i in variants_list_split]
    logging.info("Splitting the input in {} jobs".format(len(variants_dict_split)))
    for i, split in enumerate(variants_dict_split):
        with gzip.open(os.path.join(task.output_folder, "split_{}.pickle.gz".format(i)), 'wb') as fd:
            pickle.dump(split, fd)

    # Save signature dict
    logging.info("Store signature dictionary")
    with gzip.open(os.path.join(task.output_folder, "signature.pickle.gz"), 'wb') as fd:
        pickle.dump(signature_dict, fd)

    # QMap chuncks
    logging.info("Submit {} jobs to the cluster".format(len(variants_dict_split)))
    # TODO

    # Join results
    logging.info("Joining jobs output")
    # TODO

    # Plot figures
    if figures:
        logging.info("Plotting figures")
        # TODO

    logging.info("Done")
    return 1
