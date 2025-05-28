import os
from plasmidcounts import plasmidcounts
from generate_sensitivity_data import run_generate_data 

def run_plasmidcounts():
    print("\n--- Run plasmidcounts ---")
    base_dir = input("Enter base_dir path: ")
    file_dir = input("Enter file_dir name(s) (comma-separated if multiple): ").split(",")
    read_type = input("Enter read_type(s) (comma-separated): ").split(",")
    output_base = input("Enter output_base path: ")
    output_dir = input("Enter output_dir path: ")
    db_path = input("Enter db_path: ")
    tmp_path = None
    python_path = input("Enter full path to python executable: ")
    plasmidfinder_path = input("Enter full path to plasmidfinder.py: ")

    # Run
    plasmidcounts(
        python_path, plasmidfinder_path,
        output_base, base_dir, file_dir, read_type,
        output_dir, db_path,
        tmp_path=tmp_path,
        skip_interactive=True
    )


def run_sensitivity_analysis():
    print("\n--- Run sensitivity analysis ---")
    base_dir = input("The base directory to run plasmidcounts: ")
    file_dir = input("The directories where files are stored; can be specified multiple times.(comma-separated): ").split(",")
    read_type = input("The read type for each directory in file_dirs. (comma-separated): ").split(",")
    output_base = input("The base path to write output summary files like plasmidfinder_results.csv: ")
    db_path = input("The path to the plasmid database: ")
    python_path = input("Full path to the Python interpreter for running plasmidfinder: ")
    plasmidfinder_path = input("Full path to the plasmidfinder.py script: ")

    # Run
    run_generate_data(
        python_path=python_path,
        plasmidfinder_path=plasmidfinder_path,
        base_dir=base_dir,
        file_dir=file_dir,
        read_type=read_type,
        output_base=output_base,
        db_path=db_path
    )


def main():
    print("Choose which analysis to run:")
    print("1: Run plasmidcounts")
    print("2: Run sensitivity analysis")

    choice = input("Enter 1 or 2: ").strip()
    if choice == "1":
        run_plasmidcounts()
    elif choice == "2":
        run_sensitivity_analysis()
    else:
        print("Invalid choice. Exiting.")


if __name__ == "__main__":
    main()
