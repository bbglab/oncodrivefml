import gzip
import logging
import os
import pickle
import shutil
from bgqmap import QMapExecutor
from oncodrivefml.compute import multiple_test_correction
from oncodrivefml.qqplot import qqplot_png, qqplot_html, add_symbol


def drmaa_run(variants_dict, signature_dict, task, size, figures=True):

    try:
        import drmaa
    except RuntimeError:
        raise RuntimeError("It's not possible to import 'drmaa' python package. May be it's not installed or you"
                           "don't have the DRMAA_LIBRARY_PATH environment variable defined.")

    # Save signature dict
    logging.info("Store signature dictionary")
    signature_file = os.path.join(task.output_folder, "signature.pickle.gz")
    with gzip.open(signature_file, 'wb') as fd:
        pickle.dump(signature_dict, fd)

    # Split variants file into several chunks
    variants_list = list(variants_dict.items())
    variants_list_split = [variants_list[i:i+size] for i in range(0, len(variants_list), size)]
    variants_dict_split = [{k: v for k, v in i} for i in variants_list_split]
    logging.info("Splitting the input in {} jobs".format(len(variants_dict_split)))
    arguments = []
    partial_results = []
    partial_inputs = []
    for i, split in enumerate(variants_dict_split):
        split_file = os.path.join(task.output_folder, "split_in_{}.pickle.gz".format(i))
        with gzip.open(split_file, 'wb') as fd:
            pickle.dump(split, fd)
            arguments.append("-s {} -i {} -t none:{} -r {} -n split_out_{} -o {}".format(task.score_file, split_file, signature_file, task.regions_file, i, task.output_folder))
            partial_results.append("split_out_{}.pickle.gz".format(i))
            partial_inputs.append(split_file)

    # QMap chuncks
    logging.info("Submit {} jobs to the cluster".format(len(variants_dict_split)))

    retry = 1
    jobs_fail = 0
    logs_dir = os.path.join(task.output_folder, 'logs')
    while retry <= 5:

        # Create the qmap executor
        executor = QMapExecutor(
            ['normal', 'long', 'short-high', 'short-low', 'bigmem'],
            task.max_jobs,
            task.cores,
            output_folder=logs_dir,
            interactive=False,
            adaptative=False,
            commands_per_job=1
        )

        # Run all
        jobs_done, jobs_fail, jobs_skip = executor.run(
            "oncodrivefml",
            arguments,
            len(arguments),
            job_name="ofml"
        )

        # Close the executor
        executor.exit()

        if jobs_fail == 0:
            if not task.debug:
                shutil.rmtree(logs_dir)
                os.remove(signature_file)
                for partial_input in partial_inputs:
                    os.remove(partial_input)
            break
        else:
            retry += 1
            logging.info("Some jobs fail, retry {} of maximum 5".format(retry))

    if jobs_fail > 0:
        logging.error("%d jobs fail. Check the logs at '%s'.", jobs_fail, logs_dir)
        return -1

    # Join results
    logging.info("Joining jobs output")
    results = {}
    for partial_result_file in partial_results:
        partial_result_path = os.path.join(task.output_folder, partial_result_file)
        with gzip.open(partial_result_path, 'rb') as fd:
            partial_result = pickle.load(fd)
            for k, v in partial_result.items():
                results[k] = v
        if not task.debug:
            os.remove(partial_result_path)

    # Run multiple test correction
    logging.info("Computing multiple test correction")
    results_concat = multiple_test_correction(results, num_significant_samples=2)

    # Sort and store results
    results_concat.sort('pvalue', 0, inplace=True)
    fields = ['muts', 'muts_recurrence', 'samples_mut', 'pvalue', 'qvalue']
    df = results_concat[fields].copy()
    df.reset_index(inplace=True)
    df = add_symbol(df)
    with open(task.results_file, 'wt') as fd:
        df.to_csv(fd, sep="\t", header=True, index=False)

    if figures:
        logging.info("Creating figures")
        qqplot_png(task.results_file, task.qqplot_file + ".png")
        qqplot_html(task.results_file, task.qqplot_file + ".html")

    logging.info("Done")
    return 0
