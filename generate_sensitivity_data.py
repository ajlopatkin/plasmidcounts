import random
import numpy as np
from plasmidcounts import plasmidcounts

# script vars
sens_range = 0.25
num_slices = 50
num_sets = 100
base_dir = "/path/to/data/root"
file_dir = "input_data"
read_type = "short"
db_path = "/path/to/plasmid/db"
output_base = "/path/to/output/base"

#base params
mash_cutoff = 0.05
plasfinder_cov = 80
plasfinder_pid = 90
blast_id_cutoff = 90
blast_len_thresh = 27000
overlap_cutoff_short = 0.01
overlap_cutoff_long = 0.15

# generate sensitivty ranges
mash_cutoff_range = np.linspace(mash_cutoff*(1-sens_range), mash_cutoff*(1+sens_range), num_slices).tolist()
plasfinder_cov_range = np.linspace(plasfinder_cov*(1-sens_range), plasfinder_cov*(1+sens_range), num_slices).tolist()
plasfinder_pid_range = np.linspace(plasfinder_pid*(1-sens_range), plasfinder_pid*(1+sens_range), num_slices).tolist()
blast_id_cutoff_range = np.linspace(blast_id_cutoff*(1-sens_range), blast_id_cutoff*(1+sens_range), num_slices).tolist()
blast_len_thresh_range = np.linspace(blast_len_thresh*(1-sens_range), blast_len_thresh*(1+sens_range), num_slices).tolist()
overlap_cutoff_short_range = np.linspace(overlap_cutoff_short*(1-sens_range), overlap_cutoff_short*(1+sens_range), num_slices).tolist()
overlap_cutoff_long_range = np.linspace(overlap_cutoff_long*(1-sens_range), overlap_cutoff_long*(1+sens_range), num_slices).tolist()

param_sets = []
for i in range(72, 76):
    # randomly choose a parameter value
    current_params = {"mash_cutoff":random.choice(mash_cutoff_range),
                      "plasfinder_cov":random.choice(plasfinder_cov_range),
                      "plasfinder_pid":random.choice(plasfinder_pid_range),
                      "blast_id_cutoff":random.choice(blast_id_cutoff_range),
                      "blast_len_thresh":random.choice(blast_len_thresh_range),
                      "overlap_cutoff_short":random.choice(overlap_cutoff_short_range),
                      "overlap_cutoff_long":random.choice(overlap_cutoff_long_range)}

    # save the current param set
    param_sets.append(current_params)
    current_output = f"{output_base}/iteration_{i}"

    # run plasmidcounts with the current params
    plasmidcounts(base_dir, [file_dir], [read_type], current_output, db_path,
                     mash_cutoff=current_params['mash_cutoff'],
                     plasfinder_cov=current_params['plasfinder_cov'],
                     plasfinder_pid=current_params['plasfinder_pid'],
                     blast_id_cutoff=current_params['blast_id_cutoff'],
                     blast_len_thresh=current_params['blast_len_thresh'],
                     overlap_cutoff_short=current_params['overlap_cutoff_short'],
                     overlap_cutoff_long=current_params['overlap_cutoff_long'])

    # write the current params to a .txt
    with open(f"{current_output}/params.txt", "w") as f:
        f.write(f"mash_cutoff:{current_params['mash_cutoff']}\n")
        f.write(f"plasfinder_cov:{current_params['plasfinder_cov']}\n")
        f.write(f"plasfinder_pid:{current_params['plasfinder_pid']}\n")
        f.write(f"blast_id_cutoff:{current_params['blast_id_cutoff']}\n")
        f.write(f"blast_len_thresh:{current_params['blast_len_thresh']}\n")
        f.write(f"overlap_cutoff_short:{current_params['overlap_cutoff_short']}\n")
        f.write(f"overlap_cutoff_long:{current_params['overlap_cutoff_long']}\n")

# run the base params last
current_output = f"{output_base}/base_params"
plasmidcounts(base_dir, [file_dir], [read_type], current_output, db_path)