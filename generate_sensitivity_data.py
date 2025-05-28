import random
import numpy as np
from plasmidcounts import plasmidcounts

# # script vars
# sens_range = 0.75
# num_slices = 50
# num_sets = 100
# base_dir = "./"
# file_dir = ["input_data_simulatedFastas"]
# read_type = ["short"]
# db_path = "./db"
# output_base = "./processing_052021_simulatedFastas_rnd2"
# python_path = "/Users/gregory.wood/miniconda3/envs/plasmidcounts/bin/python"
# plasmidfinder_path = "/Users/gregory.wood/miniconda3/envs/plasmidcounts/bin/plasmidfinder.py"
def run_generate_data(python_path, plasmidfinder_path, base_dir, file_dir, read_type, output_base, db_path):
    # base params (round 0)
    # mash_cutoff = 0.05
    # plasfinder_cov = 80
    # plasfinder_pid = 90
    # blast_id_cutoff = 90
    # blast_len_thresh = 27000
    # overlap_cutoff_short = 0.01
    # overlap_cutoff_long = 0.15

    # base params (round 1)
    # script vars
    sens_range = 0.75
    num_slices = 50
    num_sets = 100
    mash_cutoff = 0.026
    plasfinder_cov = 60
    plasfinder_pid = 75
    blast_id_cutoff = 90
    blast_len_thresh = 10000
    overlap_cutoff_short = 0.01
    overlap_cutoff_long = 0.24


    def bounded_linspace(center, srange, slices, min_val, max_val):
        lower = max(center * (1 - srange), min_val)
        upper = min(center * (1 + srange), max_val)
        return np.linspace(lower, upper, slices).tolist()


    # update each parameter with limits
    mash_cutoff_range = bounded_linspace(mash_cutoff, sens_range, num_slices, 0, 1)
    plasfinder_cov_range = bounded_linspace(plasfinder_cov, sens_range, num_slices, 0, 100)
    plasfinder_pid_range = bounded_linspace(plasfinder_pid, sens_range, num_slices, 65, 100)
    blast_id_cutoff_range = bounded_linspace(blast_id_cutoff, sens_range, num_slices, 70, 100)
    blast_len_thresh_range = bounded_linspace(blast_len_thresh, sens_range, num_slices, 0, 100000)
    overlap_cutoff_short_range = bounded_linspace(overlap_cutoff_short, sens_range, num_slices, 0, 1)
    overlap_cutoff_long_range = bounded_linspace(overlap_cutoff_long, sens_range, num_slices, 0, 1)

    param_sets = []
    for i in range(num_sets):
        print(i)
        # randomly choose a parameter value
        current_params = {"mash_cutoff": random.choice(mash_cutoff_range),
                        "plasfinder_cov": random.choice(plasfinder_cov_range),
                        "plasfinder_pid": random.choice(plasfinder_pid_range),
                        "blast_id_cutoff": random.choice(blast_id_cutoff_range),
                        "blast_len_thresh": random.choice(blast_len_thresh_range),
                        "overlap_cutoff_short": random.choice(overlap_cutoff_short_range),
                        "overlap_cutoff_long": random.choice(overlap_cutoff_long_range)}

        # save the current param set
        param_sets.append(current_params)
        current_output = f"{output_base}/iteration_{i}"
        print(current_output)
        
        # run plasmidcounts with the current params
        plasmidcounts(python_path, plasmidfinder_path, output_base, base_dir, file_dir, read_type, current_output, db_path,
                    mash_cutoff=current_params['mash_cutoff'],
                    plasfinder_cov=current_params['plasfinder_cov'],
                    plasfinder_pid=current_params['plasfinder_pid'],
                    blast_id_cutoff=current_params['blast_id_cutoff'],
                    blast_len_thresh=current_params['blast_len_thresh'],
                    overlap_cutoff_short=current_params['overlap_cutoff_short'],
                    overlap_cutoff_long=current_params['overlap_cutoff_long'],
                    skip_interactive=True)

        # write the current params to a .txt
        # os.makedirs(current_output, exist_ok=True)
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
    plasmidcounts(python_path, plasmidfinder_path, output_base, base_dir, file_dir,
                read_type, current_output, db_path, skip_interactive=True)
