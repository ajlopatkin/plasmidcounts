import random
import numpy as np
from plasmidcounts import plasmidcounts
import os

def run_generate_data(python_path, plasmidfinder_path, base_dir, file_dir, read_type, output_base, db_path, output_dir=None):
    # script vars
    sens_range = 0.75
    num_slices = 50
    num_sets = 100

    # base params (round 1)
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

    for i in range(num_sets):
        print(i)
        current_params = {
            "mash_cutoff": random.choice(mash_cutoff_range),
            "plasfinder_cov": random.choice(plasfinder_cov_range),
            "plasfinder_pid": random.choice(plasfinder_pid_range),
            "blast_id_cutoff": random.choice(blast_id_cutoff_range),
            "blast_len_thresh": random.choice(blast_len_thresh_range),
            "overlap_cutoff_short": random.choice(overlap_cutoff_short_range),
            "overlap_cutoff_long": random.choice(overlap_cutoff_long_range),
        }

        full_output_base = os.path.join(output_base, output_dir) if output_dir else output_base
        os.makedirs(full_output_base, exist_ok=True)
        
        current_output = os.path.join(full_output_base, f"iteration_{i}")
        os.makedirs(current_output, exist_ok=True)

        plasmidcounts(
            python_path, plasmidfinder_path, output_base, base_dir,
            file_dirs=[file_dir], read_type=[read_type],
            output_dir=current_output, db_path=db_path,
            mash_cutoff=current_params['mash_cutoff'],
            plasfinder_cov=current_params['plasfinder_cov'],
            plasfinder_pid=current_params['plasfinder_pid'],
            blast_id_cutoff=current_params['blast_id_cutoff'],
            blast_len_thresh=current_params['blast_len_thresh'],
            overlap_cutoff_short=current_params['overlap_cutoff_short'],
            overlap_cutoff_long=current_params['overlap_cutoff_long'],
            skip_interactive=True
        )

        with open(os.path.join(current_output, "params.txt"), "w") as f:
            for k, v in current_params.items():
                f.write(f"{k}:{v}\n")

    base_output = os.path.join(full_output_base, "base_params")
    
    os.makedirs(base_output, exist_ok=True)
    plasmidcounts(
        python_path, plasmidfinder_path, output_base, base_dir,
        file_dirs=[file_dir], read_type=[read_type],
        output_dir=base_output, db_path=db_path,
        mash_cutoff=mash_cutoff,
        plasfinder_cov=plasfinder_cov,
        plasfinder_pid=plasfinder_pid,
        blast_id_cutoff=blast_id_cutoff,
        blast_len_thresh=blast_len_thresh,
        overlap_cutoff_short=overlap_cutoff_short,
        overlap_cutoff_long=overlap_cutoff_long,
        skip_interactive=True
    )

if __name__ == "__main__":
    pass
