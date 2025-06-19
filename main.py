import os
import yaml
from generate_sensitivity_data import run_generate_data
from plasmidcounts import plasmidcounts

def load_config():
    config_path = os.path.join(os.path.dirname(__file__), "configuration_example.yml")
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def main():
    config = load_config()

    # Extract values from config
    mode = config.get("mode", "plasmidcounts")
    base_dir = config["base_dir"]
    file_dirs = config["file_dirs"].split(",") if isinstance(config["file_dirs"], str) else config["file_dirs"]
    read_types = config["read_types"].split(",") if isinstance(config["read_types"], str) else config["read_types"]
    output_base = config["output_base"]
    output_dir = config.get("output_dir", "output")
    db_path = config["db_path"]
    python_path = config["python_path"]
    plasmidfinder_path = config["plasmidfinder_path"]

    print(f"Running in mode: {mode}")

    if mode == "plasmidcounts":
        for file_dir, read_type in zip(file_dirs, read_types):
            full_output_dir = os.path.join(output_base, output_dir)
            os.makedirs(full_output_dir, exist_ok=True)

            plasmidcounts(
                python_path=python_path,
                plasmidfinder_path=plasmidfinder_path,
                output_base=output_base,
                base_dir=base_dir,
                file_dirs=[file_dir],         
                read_type=[read_type],   
                output_dir=output_dir,
                db_path=db_path,
                skip_interactive=True
            )

    elif mode == "sensitivity":
        for file_dir, read_type in zip(file_dirs, read_types):
            run_generate_data(
                base_dir=base_dir,
                file_dir=file_dir,
                read_type=read_type,
                output_base=output_base,
                output_dir=output_dir,   # ← 추가
                db_path=db_path,
                python_path=python_path,
                plasmidfinder_path=plasmidfinder_path
            )
    else:
        raise ValueError(f"Invalid mode '{mode}' in configuration file. Must be 'plasmidcounts' or 'sensitivity'.")

if __name__ == "__main__":
    main()
