import os
import re
import argparse
import pandas as pd
from tqdm import tqdm
import scipy.stats as stat
import matplotlib.pyplot as plt


# Script to bootstrap samples based on the output of generate_sensitivity_data.py in order to determine an optimal
# parameter set for plasmidcounts.py.
#
# Note that the true plasmid counts file must have the following columns:
#   - names_unique: the unique identifier of the genome that must correspond to the filename
#   - shortread_genome_group: the list of plasmid present on the genome. This should be formatted as follows:
#       + each line contains ALL plasmids detected on the given genome
#       + plasmids on the same replicon should be separated by a pipe (|)
#       + plasmids on different replicons should be separated by a semicolon (;)
#   - shortread_genome_num: the count of plasmids on the genome
#   - contains_true_value: binary flag indicating whether the genome should be processed (1) or not (0)
#
# All plasmidcount_results.csv files should be formatted as output by plasmidcounts.py and/or generate_sensitivity.py


# function to count valid matches of a plasmidcount result to a truth file
def count_plasmid_matches(truth_df, check_df):
    iteration_contigs = list(check_df["contig_name"].drop_duplicates())

    match_count = 0
    for contig in iteration_contigs:
        current_df = check_df.loc[check_df["contig_name"] == contig]
        current_truth = truth_df['shortread_genome_group'].loc[contig]
        current_truth = re.sub(r'\([^)]*\)', '', current_truth)
        split_truth = current_truth.split(";")

        for index, row in current_df.iterrows():
            curr_row = row.plasmid_name
            curr_row = re.sub(r'\([^)]*\)', '', curr_row)
            all_plasmids = curr_row.split(";")

            for plasmid in all_plasmids:
                if "|" in plasmid:
                    has_mult = True
                    plasmid = plasmid.split("|")
                else:
                    has_mult = False

                if not has_mult:
                    is_match = (plasmid in current_truth)
                else:
                    is_match = False
                    already_matched = False
                    matched_against = None
                    previous_match = None

                    for p in plasmid:
                        if p in current_truth:
                            is_match = True
                            matched_against = [x for x in split_truth if p in x][0]

                        if is_match & (not already_matched):
                            already_matched = True
                            previous_match = matched_against
                        elif is_match & already_matched & (matched_against == previous_match):
                            continue
                        elif is_match & already_matched & (matched_against != previous_match):
                            is_match = False
                            break

                if is_match:
                    match_count += 1

    return match_count / truth_df['shortread_genome_num'].sum()


def run_bootstrap(root_dir, true_count_file, extension, match_cutoff, num_bootstrap=100):

    # script setup
    df = pd.read_csv(os.path.join(root_dir, true_count_file)).set_index("names_unique")
    true_df = df[df['contains_true_value'] == 1]
    all_strains = true_df.index.tolist()
    calibrations = range(3, len(all_strains)-5)  # the number of strains selected in each loop
    full_path = os.path.join(root_dir, true_count_file)
    all_dirs = os.listdir(root_dir)
    strain_plasmid_df = pd.DataFrame(columns=None, index=all_strains)
    plasmid_matches = {}

    # build the strain plasmid count dataframe
    for d in all_dirs:
        if not os.path.isdir(os.path.join(root_dir, d)):
            continue

        # container vars
        genome_names = []
        plasmid_names = []

        # check that results file exists
        if not os.path.isfile(os.path.join(root_dir, d, "plasmidcount_results.csv")):
            print(f"Warning: no plasmid count found for {d}. Continuing...")
        else:
            with open(os.path.join(root_dir, d, "plasmidcount_results.csv"), 'r') as f:
                lines = f.readlines()[1:]
                for line in lines:
                    stripped_line = line.rstrip()
                    genome_names.append(stripped_line.split(",")[1].split("/")[-1].split(extension)[0])
                    plasmid_names.append(stripped_line.split(",")[5])

            # build out the strain-to-plasmid dataframe
            condition_name = d
            plasmid_dict = dict(zip(genome_names, plasmid_names))
            new_df = pd.DataFrame(data={condition_name: plasmid_names}, index=genome_names)
            new_df.index.name = "genome_names"
            new_df = new_df.groupby(new_df.index)[condition_name].apply(lambda x: ';'.join(x))
            strain_plasmid_df = pd.concat([strain_plasmid_df, new_df], axis=1)

            # get the plasmid replicon match count for the current iteration
            iteration_df = pd.read_csv(os.path.join(root_dir, d, "plasmidcount_results.csv"))
            iteration_df["contig_name"] = iteration_df["filename"].str.split("/").str[-1].str.split(extension).str[0]
            plasmid_matches[d] = count_plasmid_matches(true_df, iteration_df)

    # loop through calibrations and bootstraps
    all_conditions = list(strain_plasmid_df.columns)
    plasmid_match_per_calibration = []

    for n in calibrations:
        print(f"\nCalibrating with {n} strains...")

        # container vars
        plasmid_match_collect = []
        plasmid_match_condition = []
        plasmid_match_thresh = []

        for b in tqdm(range(num_bootstrap)):
            # sample the strain/plasmid df
            sampled_df = strain_plasmid_df.sample(n).sort_index()

            # get the list of genomes
            sampled_genomes = list(sampled_df.index)

            # get the list of holdout genomes
            holdout_genomes = list(set(strain_plasmid_df.index) - set(sampled_genomes))

            # get the same genomes from the true df
            sampled_true_df = true_df.loc[sampled_genomes]

            # calculate pct match for sampled genomes
            pct_match = [count_plasmid_matches(
                sampled_true_df, sampled_df[x].reset_index().rename(columns={"index": "contig_name", x: "plasmid_name"}))
                for x in list(sampled_df.columns)]

            # overindexed conditions should be ignored
            pct_match = [float("-inf") if i > 1 else i for i in pct_match]

            # get the condition for the max calculated pct match
            max_match_idx = pct_match.index(max(pct_match))
            max_match_condition = all_conditions[max_match_idx]

            # calculate pct match for holdout genomes for the condition chosen above and save it
            holdout_df = (strain_plasmid_df[max_match_condition]
                          .loc[holdout_genomes]
                          .reset_index()
                          .rename(columns={"index": "contig_name", max_match_condition: "plasmid_name"}))

            # perform pct match calculation
            pct_match_holdout = count_plasmid_matches(true_df.loc[holdout_genomes], holdout_df)

            # only add pct_match to the tracking list if it is valid (ie <1)
            if pct_match_holdout <= 1:
                plasmid_match_collect.append(pct_match_holdout)
                plasmid_match_condition.append(max_match_condition)

            # get all conditions with match pct > cutoff
            thresh_idx = [x for x, val in enumerate(pct_match) if val > match_cutoff]
            thresh_condition = [all_conditions[x] for x in thresh_idx]

            # calculate the match pct for ALL strains for the conditions found above
            thresh_vals = [
                count_plasmid_matches(true_df.loc[all_strains],
                                      strain_plasmid_df[x]
                                      .reset_index()
                                      .rename(columns={"index": "contig_name", x: "plasmid_name"}))
                for x in thresh_condition]

            # again, remove conditions where pct >1
            thresh_vals = [x for x in thresh_vals if x <= 1]
            plasmid_match_thresh.append(stat.tmean(thresh_vals))

        # append to the tracking dict
        plasmid_match_per_calibration.append({'calibration': n,
                                              'raw': plasmid_match_collect,
                                              'conditions': plasmid_match_condition,
                                              'thresh_mean': stat.tmean(plasmid_match_thresh),
                                              'thresh_sem': stat.sem(plasmid_match_thresh),
                                              'thresh_std': stat.tstd(plasmid_match_thresh),
                                              'mean': stat.tmean(plasmid_match_collect),
                                              'std': stat.tstd(plasmid_match_collect),
                                              'sem': stat.sem(plasmid_match_collect)})

    # plot count match threshold
    match_mean = [x['thresh_mean'] for x in plasmid_match_per_calibration]
    match_err = [x['thresh_sem'] for x in plasmid_match_per_calibration]
    plt.figure(figsize=(8, 5))
    plt.errorbar(calibrations, match_mean, yerr=match_err, fmt='-o', capsize=5, color='blue', label='Mean +- SEM')
    plt.xlabel("Calibration Size")
    plt.ylabel("Thresh Match Pct")
    plt.title("Bootstrap Calibration Performance - Threshold Count Match Percentage")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("bootstrap_performance_plot.pdf")

    # get condition with best match
    all_conditions = [i for sl in [x['conditions'] for x in plasmid_match_per_calibration] for i in sl]
    all_pct = [i for sl in [x['raw'] for x in plasmid_match_per_calibration] for i in sl]
    pct_mapping_df = pd.DataFrame({"condition": all_conditions, "pct": all_pct})
    grouped_pct_df = pct_mapping_df.groupby("condition", as_index=False)["pct"].mean()
    grouped_pct_df['pct'] = grouped_pct_df['pct'].apply(lambda x: x if x < 1 else float('-inf'))

    # print best match and related info
    best_match_condition = grouped_pct_df.loc[grouped_pct_df['pct'] == max(grouped_pct_df['pct'])]["condition"].iloc[0]
    plasmid_match_pct = plasmid_matches[best_match_condition]
    print(f"Best match condition: {best_match_condition}; matched plasmids: {plasmid_match_pct * 100}%")


def get_parser():
    parser = argparse.ArgumentParser(description="Bootstrap to determine optimal parameter set for plasmidcounts.")

    parser.add_argument("--root_dir", required=True,
                        help="The base directory; this should contain output from generate_sensitivity_data.py")
    parser.add_argument("--true_count_file", required=True,
                        help="The file containing the actual count of plasmids for each genome.")
    parser.add_argument("--extension", required=True,
                        help="The extension of fasta files that follows the unique identifier of each genome.")
    parser.add_argument("--match_cutoff", required=True, help="The match pct threshold.")
    parser.add_argument("--num_bootstrap", required=True, default=100, help="The number of bootstrap iterations.")

    return parser


def main():
    args = get_parser().parse_args()
    run_bootstrap(args.root_dir,
                  args.true_count_file,
                  args.extension,
                  args.match_cutoff,
                  args.num_bootstrap)


if __name__ == "__main__":
    main()