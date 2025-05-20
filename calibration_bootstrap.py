import os
import scipy.stats as stat
import pandas as pd
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt

def count_matches(df1, df2):
    try:
        match_count = (df1 == df2).value_counts()[True]
        return match_count
    except KeyError:
        return 0


# script params
root_dir = "/path/to/data/root"
true_count_file = "true_counts.csv"
true_df = pd.read_csv(os.path.join(root_dir, true_count_file)).set_index("names_unique")
full_path = os.path.join(root_dir, true_count_file)
all_dirs = os.listdir(root_dir)
all_strains = ['AL_100', 'AL_12', 'AL_16', 'AL_25', 'AL_37', 'AL_49', 'AL_61', 'AL_99', 'KH2_17', 'KH2_22',
              'KH2_37', 'KH2_7', 'KH_10', 'KH_11', 'KH_12', 'KH_13', 'KH_24', 'KH_35', 'KH_3', 'KH_40',
              'KH_45', 'KH_48', 'KH_4', 'KH_5', 'KH_6', 'KH_7', 'KH_8', 'KH_26', 'KH_39', 'KH_41', 'KH_9',
              'AL_98', 'KH_31']

calibrations = range(3,31)  # the number of strains selected in each loop
num_bootstrap = 100 # the number of bootstraps run for each strain set

strain_plasmid_df = pd.DataFrame(columns=None, index=all_strains)

# build the strain plasmid count dataframe
for d in all_dirs:
    if not os.path.isdir(os.path.join(root_dir, d)):
        continue
    elif d == "base_params":
        continue

    genome_names = []
    plasmid_count = []

    if not os.path.isfile(os.path.join(root_dir, d, "plasmidcount_aggregated.csv")):
        print(f"Warning: no aggregated plasmid count found for {d}. Continuing...")
    else:
        with open(os.path.join(root_dir, d, "plasmidcount_aggregated.csv"), 'r') as f:
            lines = f.readlines()[1:]
            for line in lines:
                stripped_line = line.rstrip()
                genome_names.append(stripped_line.split(",")[-2].split("/")[-1].split("_genome.fasta")[0])
                plasmid_count.append(stripped_line.split(",")[-1])

        plasmid_count = [int(x) for x in plasmid_count]
        condition_name = d
        plasmid_dict = dict(zip(genome_names, plasmid_count))
        new_df = pd.DataFrame(data={condition_name: plasmid_count}, index=genome_names)
        strain_plasmid_df = pd.concat([strain_plasmid_df, new_df], axis=1)

# loop through calibrations and bootstraps
all_conditions = list(strain_plasmid_df.columns)
r2_per_calibration = []
match_per_calibration = []

for n in calibrations:
    print(f"Calibrating with {n} strains...")

    r2_collect = []
    r2_condition = []
    r2_thresh = []
    match_collect = []
    match_condition = []
    match_thresh = []

    for b in range(num_bootstrap):
        # sample the strain/plasmid df
        sampled_df = strain_plasmid_df.sample(n).sort_index()

        # get the list of genomes
        sampled_genomes = list(sampled_df.index)

        # get the list of holdout genomes
        holdout_genomes = list(set(strain_plasmid_df.index) - set(sampled_genomes))

        # get the same genomes from the true df
        sampled_true_df = true_df.loc[sampled_genomes]

        # calculate r2 for the sampled genomes
        r2s = [r2_score(sampled_true_df['shortread_genome_num'], sampled_df[x]) for x in list(sampled_df.columns)]

        # get the condition for the max calculated r2
        max_r2_idx = r2s.index(max(r2s))
        max_r2_condition = all_conditions[max_r2_idx]

        # calculate r2 for the holdout genomes for the condition chosen above and save it
        r2_collect.append(r2_score(strain_plasmid_df[max_r2_condition].loc[holdout_genomes],
                                   true_df['shortread_genome_num'].loc[holdout_genomes]))
        r2_condition.append(max_r2_condition)

        # get all conditions with r2 > 0.9
        thresh_idx = [x for x, val in enumerate(r2s) if val > 0.9]
        thresh_condition = [all_conditions[x] for x in thresh_idx]

        # calculate r2 for ALL strains for the conditions found above
        thresh_vals = [r2_score(strain_plasmid_df[x], true_df['shortread_genome_num'].loc[all_strains]) for x in thresh_condition]
        r2_thresh.append(stat.tmean(thresh_vals))

        # calculate pct match for sampled genomes
        pct_match = [count_matches(sampled_true_df['shortread_genome_num'], sampled_df[x])/n for x in list(sampled_df.columns)]

        # get the condition for the max calculated pct match
        max_match_idx = pct_match.index(max(pct_match))
        max_match_condition = all_conditions[max_match_idx]

        # calculate pct match for holdout genomes for the condition chosen above and save it
        match_collect.append(count_matches(strain_plasmid_df[max_match_condition].loc[holdout_genomes],
                                           true_df['shortread_genome_num'].loc[holdout_genomes]) / len(holdout_genomes))
        match_condition.append(max_match_condition)

        # get all conditions with match pct > 0.9
        thresh_idx = [x for x, val in enumerate(pct_match) if val > 0.9]
        thresh_condition = [all_conditions[x] for x in thresh_idx]

        # calculate the match pct for ALL strains for the conditions found above
        thresh_vals = [count_matches(true_df['shortread_genome_num'].loc[all_strains], strain_plasmid_df[x])/len(all_strains) for x in thresh_condition]
        match_thresh.append(stat.tmean(thresh_vals))

    r2_per_calibration.append({'calibration':n,
                               'raw':r2_collect,
                               'conditions': r2_condition,
                               'thresh_mean': stat.tmean(r2_thresh),
                               'thresh_sem': stat.sem(r2_thresh),
                               'mean': stat.tmean(r2_collect),
                               'sem': stat.sem(r2_collect)})

    match_per_calibration.append({'calibration': n,
                                  'raw': match_collect,
                                  'conditions': match_condition,
                                  'thresh_mean': stat.tmean(match_thresh),
                                  'thresh_sem': stat.sem(match_thresh),
                                  'mean': stat.tmean(match_collect),
                                  'sem': stat.sem(match_collect)})

# plot R2
r2_mean = [x['mean'] for x in r2_per_calibration]
r2_sem = [x['sem'] for x in r2_per_calibration]
plt.figure(figsize=(8, 5))
plt.errorbar(calibrations, r2_mean, yerr=r2_sem, fmt='-o', capsize=5, color='blue', label='Mean +- SEM')
plt.xlabel("Calibration Size")
plt.ylabel("R2")
plt.title("Bootstrap Calibration Performance - R2")
plt.grid(True)
plt.tight_layout()

# plot pct match
match_mean = [x['mean'] for x in match_per_calibration]
match_sem = [x['sem'] for x in match_per_calibration]
plt.figure(figsize=(8, 5))
plt.errorbar(calibrations, match_mean, yerr=match_sem, fmt='-o', capsize=5, color='blue', label='Mean +- SEM')
plt.xlabel("Calibration Size")
plt.ylabel("Match Pct")
plt.title("Bootstrap Calibration Performance - Match Percentage")
plt.grid(True)
plt.tight_layout()

# plot r2 threshold
r2_mean = [x['thresh_mean'] for x in r2_per_calibration]
r2_sem = [x['thresh_sem'] for x in r2_per_calibration]
plt.figure(figsize=(8, 5))
plt.errorbar(calibrations, r2_mean, yerr=r2_sem, fmt='-o', capsize=5, color='blue', label='Mean +- SEM')
plt.xlabel("Calibration Size")
plt.ylabel("Thresh R2")
plt.title("Bootstrap Calibration Performance - Threshold R2")
plt.grid(True)
plt.tight_layout()

# plot match threshold
match_mean = [x['thresh_mean'] for x in match_per_calibration]
match_sem = [x['thresh_sem'] for x in match_per_calibration]
plt.figure(figsize=(8, 5))
plt.errorbar(calibrations, match_mean, yerr=match_sem, fmt='-o', capsize=5, color='blue', label='Mean +- SEM')
plt.xlabel("Calibration Size")
plt.ylabel("Thresh Match Pct")
plt.title("Bootstrap Calibration Performance - Threshold Match Percentage")
plt.grid(True)
plt.tight_layout()
plt.show()