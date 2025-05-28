# env: plasmidCount_env
import os
import re
import json
import time
import shutil
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import OrderedDict
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline


def get_reads(desc):
    m = re.search('reads=(.+?) ', desc)
    try:
        return int(m.group(1))
    except:
        return 1000000


def get_coverage(desc):
    m = re.search('cov_(.+?)_', desc)
    try:
        return float(m.group(1))
    except:
        return 1000000


def get_contig_len(desc):
    m1 = re.search('len=(.+?) ', desc)
    m2 = re.search('_length_(.+?)_', desc)
    try:
        return float(m1.group(1))
    except:
        try:
            return float(m2.group(1))
        except:
            return 1000000


def build_overlap_matrix(blast_results):
    mat_size = len(blast_results)
    if mat_size < 2:
        return np.zeros((0, 0))
    else:
        overlap_matrix = np.zeros((mat_size, mat_size))

    for x in range(mat_size):
        for y in range(mat_size):
            if x == y:
                continue
            set1 = set(blast_results[x])
            set2 = set(blast_results[y])
            overlap_matrix[x, y] = len(set1.intersection(set2)) / len(set1)

    return overlap_matrix


# main function for PlasmidCounts pipeline.
def plasmidcounts(python_path, plasmidfinder_path, output_base, base_dir, file_dirs, read_type, output_dir, db_path,
                  tmp_path=None,
                  contig_len_cutoff=500,
                  mash_cutoff=0.05,
                  plasfinder_cov=80,
                  plasfinder_pid=90,
                  read_threshold=5,
                  cover_threshold=1,
                  blast_id_cutoff=90,
                  blast_len_thresh=27000,
                  overlap_cutoff_short=0.01,
                  overlap_cutoff_long=0.15,
                  blast_len_cutoff=1000,
                  filter_col=True,
                  filter_std=True,
                  num_contig_thresh=3,
                  skip_interactive=False):
    if not isinstance(file_dirs, list):
        print("Error: file_dirs must be a list.")
        raise ValueError
    elif not isinstance(read_type, list):
        print("Error: read_type must be a list.")
        raise ValueError

    # initialize containers
    existing_unique_files = []
    skip_header = 0

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # stage 1: run plasmidfinder on all reads
    # this produces two files: pfinder_raw.json and plasmidfinder_results.csv
    # pfinder_raw is a collection of all output from PlasmidFinder, which by default outputs JSON objects
    # plasmidfinder_results is the this data transformed into a CSV for easier downstream processing in Pandas
    # if an existing file is found, the script will just add net new strains to the file instead of recreating everything
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # if tmp_path not provided, put it in output_dir
    if not tmp_path:
        tmp_path = os.path.join(output_dir, "tmp_dir")
        os.mkdir(tmp_path)

    pfinder_result_path = output_base + "/plasmidfinder_results.csv"
    # "os.path.join(output_base, "/plasmidfinder_results.csv")
    print(pfinder_result_path)

    if os.path.isfile(pfinder_result_path):
        if not skip_interactive:
            print("Warning: plasmidfinder_results.csv already exists. Run plasmidfinder anyways?")
            ans = input("Enter [y/n]:")
            if ans.lower() not in ["y", "n"]:
                print("Input not recognized. Not running...")
                do_plasfinder = 0
            elif ans.lower() == "y":
                print("Appending new files to existing results...")
                existing_pfinder_results = pd.read_csv(pfinder_result_path)
                existing_unique_files = [s.split("/")[-1] for s in list(existing_pfinder_results["filename"].unique())]
                skip_header = 1
                do_plasfinder = 1
            else:
                print("Not running...")
                do_plasfinder = 0
        else:
            print("Found existing plasmidfinder_results.csv. Skipping plasmidfinder...")
            do_plasfinder = 0
    else:
        do_plasfinder = 1

    if do_plasfinder:
        # obj_container will collect all plasmidfinder json objects
        obj_container = []

        # loop through file directories
        print(tmp_path)
        json_path = os.path.join(tmp_path, "data.json")
        pf_tmp = os.path.join(tmp_path, "tmp")

        for cdir, filetype in zip(file_dirs, read_type):
            full_path = os.path.join(base_dir, cdir)
            # loop through files in each directory
            for (dirpath, dirnames, filenames) in os.walk(full_path):
                for cfile in filenames:
                    if cfile.split(".")[-1] != "fasta":
                        continue
                    elif cfile[0] == ".":
                        continue

                    print(f"Processing sample {cfile}...")

                    if cfile in existing_unique_files:
                        print(f"File {cfile} already exists in plasmidfinder results. Skipping...")
                        continue

                    if os.path.isdir(pf_tmp):
                        try:
                            shutil.rmtree(pf_tmp)
                        except FileNotFoundError:
                            print(f"Warning: could not delete {pf_tmp}; skipping...")

                    if os.path.isfile(json_path):
                        os.remove(json_path)

                    time.sleep(0.5)

                    curr_fname = os.path.join(dirpath, cfile).replace("//", "/")
                    ret = os.system(f"{python_path} {plasmidfinder_path} -i {curr_fname} -o {tmp_path} >/dev/null 2>&1")
                    # print(f"Return code: {ret}")

                    # get the results in json object and then append to the container
                    with open(json_path) as f:
                        obj = json.load(f)

                    obj["filetype"] = filetype
                    obj_container.append(obj)

        # dump the raw json results for potential reuse
        raw_json_path = os.path.join(tmp_path, "pfinder_raw.json")
        with open(raw_json_path, "w") as f:
            json.dump(obj_container, f)

        # write header to plasmidfinder results
        if not skip_header:
            with open(pfinder_result_path, "w") as f:
                f.write("filename,file_type,plasmid_name,plasmid_accession,contig,pctId,coverage,hit_len,template_length,subj_start,subj_end,query_start,query_end,read_len,contig_cover,contig_len")
                f.write("\n")

        # build the CSV file with plasmidfinder results
        for i in obj_container:
            filename = i["plasmidfinder"]["user_input"]["filename(s)"][0]
            results = i["plasmidfinder"]["results"]
            filetype = i["filetype"]
            result_keys = results.keys()
            file_seqs = list(SeqIO.parse(filename, "fasta"))
            seqs_to_reads = {x.name: get_reads(x.description) for x in file_seqs}
            seqs_to_cover = {x.name: get_coverage(x.description) for x in file_seqs}
            seqs_to_len = {x.name: get_contig_len(x.description) for x in file_seqs}
            for k in result_keys:
                org = results[k]
                org_keys = org.keys()
                for kk in org_keys:
                    org_l2 = org[kk]
                    if isinstance(org_l2, str):
                        if k != "Gram Positive":
                            with open(pfinder_result_path, "a") as f:
                                f.write(
                                    f"{filename},{filetype},no plasmid found,n/a,n/a,1000,1000,0,0,0,0,0,0,1000000,1000000,1000000")
                                f.write("\n")
                        continue
                    hits = org_l2.keys()
                    for contig in hits:
                        curr_hit = org_l2[contig]
                        curr_plasmid = curr_hit["plasmid"]
                        curr_contig = contig.split(" ")[0]
                        curr_accession = curr_hit["accession"]
                        curr_start_subj = curr_hit["positions_in_contig"].split("..")[0]
                        curr_end_subj = curr_hit["positions_in_contig"].split("..")[1]
                        curr_start_query = curr_hit["position_in_ref"].split("..")[0]
                        curr_end_query = curr_hit["position_in_ref"].split("..")[1]
                        curr_pctid = curr_hit["identity"]
                        curr_coverage = curr_hit["coverage"]
                        curr_len = curr_hit["HSP_length"]
                        curr_templen = curr_hit["template_length"]
                        curr_readlen = seqs_to_reads[curr_contig.split(":")[0]]
                        curr_cover = seqs_to_cover[curr_contig.split(":")[0]]
                        contig_len = seqs_to_len[curr_contig.split(":")[0]]

                        with open(pfinder_result_path, "a") as f:
                            f.write(f"{filename},{filetype},{curr_plasmid},{curr_accession},{curr_contig},{curr_pctid},\
                                      {curr_coverage},{curr_len},{curr_templen},{curr_start_subj},{curr_end_subj},\
                                      {curr_start_query},{curr_end_query},{curr_readlen},{curr_cover},{contig_len}")
                            f.write("\n")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # stage 2: filter and collapse using pandas
    # we now take all results from plasmidfinder, and first drop any duplicate rows (i.e., same file, contig, and plasmid)
    # we then group by filename (i.e., strain) and contig, which gives us one line per contig
    # each contig potentially has multiple plasmid replicons, but we de-duplicate the plasmid accession and name
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # get the plasmidfinder results and make sure plasmid acc is a string
    df_raw = pd.read_csv(pfinder_result_path)
    df_raw['plasmid_accession'] = df_raw['plasmid_accession'].astype(str)

    # drop entries that did not have a valid plasmid
    df_raw = df_raw.loc[df_raw['pctId'] < 1000]

    # filter col plasmids
    if filter_col:
        df_raw = df_raw[~df_raw['plasmid_name'].str.contains("Col", na=False)]

    # drop entries that do not meet the read length criteria
    df_raw = df_raw.loc[df_raw['read_len'] > read_threshold]

    # drop entries that do not meet the contig coverage threshold
    df_raw = df_raw.loc[df_raw['contig_cover'] > cover_threshold]

    # drop entries that do not meet the coverage threshold
    df_raw = df_raw.loc[df_raw['coverage'] >= plasfinder_cov]

    # drop any exact duplicates where rows have the same filename, contig, plasmid acc and name
    df_raw.drop_duplicates(subset=['filename', 'contig', 'plasmid_accession', 'plasmid_name'], inplace=True)

    # if enabled, filter high-variance contigs in strains
    if filter_std:
        strains_to_filter = df_raw.assign(counts=1)[['filename', 'counts']].groupby(['filename']).count()
        strains_to_filter = list(strains_to_filter.loc[strains_to_filter["counts"] > num_contig_thresh].index)
        all_files = df_raw.groupby('filename')['contig_len']
        contig_means = all_files.transform('mean')
        contig_stds = all_files.transform('std')
        m = df_raw['contig_len'].lt(contig_means.sub(contig_stds.mul(1)))
        contigs_to_drop = df_raw[(df_raw['contig_len'] < blast_len_thresh)
                                 & (df_raw["filename"].isin(strains_to_filter))
                                 & m].index
        df_raw.drop(contigs_to_drop, inplace=True)

    # group the resulting dataframe by filename and contig, collapsing contigs with >1 replicon to a single line
    grouped_df = (df_raw[['filename', 'file_type', 'contig', 'plasmid_name', 'plasmid_accession', 'contig_len']]
                  .assign(contig=lambda x: x.contig.astype(str))
                  .assign(contig=lambda x: x.contig.str.split(":").str[0])
                  .groupby(['filename', 'file_type', 'contig', 'contig_len'])
                  .agg({'plasmid_name': '|'.join, 'plasmid_accession': '|'.join})
                  .reset_index())

    # drop contigs that are < contig_cutoff length
    ind_to_drop = []
    file_list = list(set(list(grouped_df.filename)))
    for file in file_list:
        short_contigs = [x.id for x in list(SeqIO.parse(file, "fasta")) if len(x.seq) < contig_len_cutoff]
        ind_to_drop = list(grouped_df[(grouped_df["filename"] == file)
                                      & (grouped_df.contig.isin(short_contigs))].index)
        grouped_df.drop(index=ind_to_drop, inplace=True)

    # remove duplicate accessions from a single line
    grouped_df["plasmid_accession"] = (grouped_df["plasmid_accession"]
                                       .apply(lambda x: "|".join(list(OrderedDict.fromkeys(x.split("|"))))))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # stage 3: remove repeat replicons that exist on the same strain
    # in this stage, we compare all contigs of a strain that have the same replicon using the mash distance
    # if any two contigs fall within the distance defined by mash_cutoff, they are considered duplicate, and one is removed
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # contig count df gets any files that have >1 contigs with duplicate plasmid accessions
    # these will be the contigs that we compare via mash below
    contig_count_df = (grouped_df[['filename', 'plasmid_accession', 'contig']]
                       .assign(counts=1)
                       .groupby(['filename', 'plasmid_accession'])
                       .agg({'contig': '|'.join, 'counts': 'sum'})
                       .reset_index())
    contig_count_df = contig_count_df[contig_count_df.counts > 1]

    # remove the tmp mash directory if it already exists
    msh_tmp = os.path.join(output_dir, "mash_folder")
    if os.path.exists(msh_tmp):
        try:
            shutil.rmtree(msh_tmp)
        except FileNotFoundError:
            print("Warning: could not delete mash tmp folder. Skipping...")

    similarity_df = pd.DataFrame([], columns=['filename', 'plasmid_acc', 'contig1', 'contig2'])
    for index, row in contig_count_df.iterrows():
        print(f"Checking contigs for {row['filename'].split('/')[-1]} plasmid {row['plasmid_accession']}...")

        # set up the directories and loop variables for mash
        try:
            os.mkdir(msh_tmp)
        except FileExistsError:
            print("Warning: mash directory already exists. Reusing...")

        contig_list = row['contig'].split('|')
        filename = row['filename']
        records = [x for x in list(SeqIO.parse(filename, "fasta")) if x.id in contig_list]

        if not records:
            try:
                shutil.rmtree(msh_tmp)
            except FileNotFoundError:
                print("Warning: could not delete mash tmp folder. Skipping...")
            continue

        # write individual fasta files into the mash folder
        for record in records:
            mash_target = os.path.join(msh_tmp, record.id + ".fasta")
            SeqIO.write(record, mash_target, "fasta")

        # run the mash sketch and dist commands
        os.system(f"mash sketch -o {msh_tmp}/reference {msh_tmp}/*.fasta >/dev/null 2>&1")
        result_path = os.path.join(msh_tmp, 'mash_results.txt')
        ref_path = os.path.join(msh_tmp, 'reference.msh')
        for file in os.listdir(msh_tmp):
            if ".fasta" not in file:
                continue
            else:
                file_path = os.path.join(msh_tmp, file)
                os.system(f"mash dist {ref_path} {file_path} >> {result_path}")

        # read mash results into a dataframe
        mash_df = pd.read_csv(result_path, delimiter="\t", header=None)
        mash_df.rename(columns={0: 'seq1', 1: 'seq2', 2: 'dist', 3: 'pVal', 4: 'hashNum'}, inplace=True)

        # loop through the mash result dataframe
        for index2, row2 in mash_df.iterrows():
            if row2['seq1'] == row2['seq2']:
                continue
            # if the distance between the two contigs is less than the cutoff, add a row to similarity df
            elif row2['dist'] < mash_cutoff:
                similarity_df = similarity_df.append({'filename': filename,
                                                      'plasmid_acc': row['plasmid_accession'],
                                                      'contig1': row2['seq1'].split("/")[1],
                                                      'contig2': row2['seq2'].split("/")[1]}, ignore_index=True)
        # clean up tmp folder
        try:
            shutil.rmtree(msh_tmp)
        except FileNotFoundError:
            print("Warning: could not delete Mash tmp folder. Skipping...")

    if not similarity_df.empty:
        # check string eliminates duplicate rows with reversed order (i.e., p1,p2 and p2,p1)
        similarity_df['check_string'] = similarity_df.apply(lambda x: ''.join(sorted([x['contig1'], x['contig2']])),
                                                            axis=1)
        similarity_df.drop_duplicates('check_string', inplace=True)

        # loop through the similarity dataframe, where each row is a duplicate to be dropped
        for index, row in similarity_df.iterrows():
            # get the filename, plasmid accession and contig for each row
            file = row["filename"]
            plasmid_acc = row["plasmid_acc"]
            contig = row["contig1"].split(".fasta")[0]

            # drop the rows that are determined to be duplicate
            grouped_df.drop(grouped_df[(grouped_df["filename"] == file)
                                       & (grouped_df["plasmid_accession"] == plasmid_acc)
                                       & (grouped_df["contig"].str.contains(contig))].index,
                            inplace=True)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # stage 4: blast each identified plasmid contig against PLSdb
    # now, we blast the contig in each row of grouped_df against the PLSdb database and record the top hits for each
    # hits that are kept are those that meet the %ID and length thresholds defined by the cutoffs at the top of the script
    # for a given strain, if the top hits of two contigs overlap by >90%, those two contigs are considered duplicate
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # make the blast DB if it doesn't exist
    plsdb_fna = os.path.join(db_path, "plsdb.fna")
    if not os.path.isfile(os.path.join(db_path, "plsdb.fna.nhr")):
        print("Creating PLSdb database...")
        cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=os.path.join(base_dir, plsdb_fna))
        cline()

    # for each row in grouped df, get all the records in the corresponding fasta
    grouped_df["blast_plasmids"] = ""
    blast_tmp_dir = os.path.join(tmp_path, "blast_tmp")
    if not os.path.isdir(blast_tmp_dir):
        os.mkdir(blast_tmp_dir)

    for index, row in grouped_df.iterrows():
        for seq_record in SeqIO.parse(row.filename, "fasta"):
            if seq_record.id in row.contig:
                print(f"Blasting {row.filename.split('/')[-1]} against PLSdb...")

                # write the fasta file and then execute the blast
                tmp_fna = os.path.join(blast_tmp_dir, "tmp.fasta")
                tmp_txt = os.path.join(blast_tmp_dir, "blast_tmp.txt")
                SeqIO.write(seq_record, tmp_fna, "fasta")
                cline = NcbiblastnCommandline(query=tmp_fna,
                                              db=plsdb_fna,
                                              out=tmp_txt,
                                              outfmt="6 sseqid pident length evalue")
                cline()

                # remove the beginning lines, since they are not structured, then read the blast results
                os.system(f"sed -i '' '/^#/d' {tmp_txt}")
                blast_df = pd.read_csv(tmp_txt, sep='\t', header=None)
                blast_df.rename(columns={0: 'seqID', 1: 'pctId', 2: 'length', 3: 'eval'}, inplace=True)

                # drop rows that do not have enough coverage or pctID
                blast_df.drop(blast_df[(blast_df["pctId"] < blast_id_cutoff)
                                       | (blast_df["length"] < blast_len_cutoff)].index,
                              inplace=True)

                # get unique hits and insert them into the dataframe
                blast_df.drop_duplicates(subset=['seqID'], inplace=True)
                unique_hits = list(blast_df["seqID"])
                if not blast_df.empty:
                    top_hit = blast_df.sort_values(by=["pctId", "length"], ascending=False).reset_index()["seqID"][0]
                    grouped_df.at[index, "top_hit"] = top_hit
                grouped_df.at[index, "blast_plasmids"] = unique_hits

                # remove temp blast files
                if os.path.isfile(tmp_fna):
                    os.remove(tmp_fna)
                if os.path.isfile(tmp_txt):
                    os.remove(tmp_txt)

    # pickle the final DF to be used later (faster than having to re-run blast, etc.)
    pkl_path = os.path.join(output_dir, "plasmid_results.pkl")
    grouped_df.to_pickle(pkl_path)

    # loop through all files in grouped_df
    for filename in list(grouped_df["filename"].unique()):

        # get the lines corresponding to the current filename that have at least 2 blast results
        tmp_df = grouped_df[(grouped_df["filename"] == filename)
                            & (grouped_df["blast_plasmids"].apply(lambda x: len(x)) > 1)]

        # all plasmids contains a list of lists with blast results for each contig
        all_plasmids = [x for x in tmp_df["blast_plasmids"]]
        if len(all_plasmids) < 2:
            continue

        # build a matrix of overlap percentages
        overlap_mat = build_overlap_matrix(all_plasmids)

        # build list of cutoffs; for plasmids with <75 blast hits, don't combine (since percentages could be skewed)
        overlap_list = [overlap_cutoff_long if x >= blast_len_thresh else overlap_cutoff_short for x in
                        tmp_df["contig_len"]]
        overlap_list = [overlap_list[x] if len(all_plasmids[x]) > 75 else 1 for x in range(len(overlap_list))]

        # get initial combination list based on overlap
        combine_list = {}
        for x in range(len(all_plasmids)):
            for y in range(len(all_plasmids)):
                if overlap_mat[x, y] > overlap_list[y]:
                    combine_list[x] = y

        # reduce redundant entries
        for key, val in combine_list.items():
            d_to_change = [x for x in combine_list if combine_list[x] == key]
            for ind in d_to_change:
                combine_list[ind] = val

        if combine_list:
            combine_list = {key: val for key, val in combine_list.items() if val != key}

        # remove collapsed entries from grouped_df
        for key, val in combine_list.items():
            key_idx = tmp_df.index[key]
            val_idx = tmp_df.index[val]

            # collapse the rows that are determined to be duplicates
            grouped_df.at[val_idx, "plasmid_name"] = f"{grouped_df['plasmid_name'].loc[val_idx]}|" \
                                                     f"{grouped_df['plasmid_name'].loc[key_idx]}"
            grouped_df.at[val_idx, "plasmid_accession"] = f"{grouped_df['plasmid_accession'].loc[val_idx]}|" \
                                                          f"{grouped_df['plasmid_accession'].loc[key_idx]}"
            grouped_df.at[val_idx, "contig"] = f"{grouped_df['contig'].loc[val_idx]}|" \
                                               f"{grouped_df['contig'].loc[key_idx]}"
            grouped_df.drop(index=key_idx, inplace=True)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # stage 5: generate multi-fastas (not currently used)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # print("Writing combined contig multi-fastas...")
    # unique_files = list(grouped_df["filename"].unique())
    # combined_fasta_dir = os.path.join(output_dir, "combined_fastas")
    # if not os.path.isdir(combined_fasta_dir):
    #     os.mkdir(combined_fasta_dir)

    # for file in unique_files:
    #     base_file = file.split("/")[-1].split(".fasta")[0]
    #     target_file = base_file + "_plasmidcontigs.fasta"
    #     out = os.path.join(combined_fasta_dir, target_file)

    #     filtered_df = grouped_df[grouped_df["filename"] == file]
    #     raw_fasta = list(SeqIO.parse(file, "fasta"))
    #     full_seq = []

    #     for index, row in filtered_df.iterrows():
    #         contigs = row.contig.split("|")
    #         seq_to_add = Seq("")
    #         header = f"sample={base_file}; contigs={row.contig}; replicons={row.plasmid_name}; contig_len="
    #         for contig in contigs:
    #             if contig == 'nan':
    #                 continue
    #             curr_fasta = [x for x in raw_fasta if x.name == contig][0]
    #             seq_to_add += curr_fasta.seq
    #             header += f"{len(curr_fasta.seq)}|"

    #         full_seq.append(SeqRecord(seq_to_add,
    #                                   id=header.rstrip("|"),
    #                                   name=header.rstrip("|"),
    #                                   description=""))
    #     SeqIO.write(full_seq, out, "fasta")

    # print("Writing combined PLSdb multi-fastas...")
    # combined_plsdb_dir = os.path.join(output_dir, "combined_plsdb")
    # if not os.path.isdir(combined_plsdb_dir):
    #     os.mkdir(combined_plsdb_dir)

    # for file in unique_files:
    #     full_plsdb = SeqIO.parse(plsdb_fna, "fasta")
    #     base_file = file.split("/")[-1].split(".fasta")[0]
    #     target_file = base_file + "_allplsdb.fasta"
    #     out = os.path.join(combined_plsdb_dir, target_file)

    #     filtered_df = grouped_df[grouped_df["filename"] == file]
    #     all_plasmids = [x for x in filtered_df["top_hit"]]
    #     plas_seqs = [f for f in full_plsdb if f.id in all_plasmids]
    #     SeqIO.write(plas_seqs, out, "fasta")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # cleanup & output
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # final_count_df contains all files and the count of distinct plasmids on each
    final_count_df = (grouped_df[['filename']]
                      .assign(counts=1)
                      .groupby(['filename'])
                      .agg({'counts': 'sum'})
                      .reset_index())

    # add lines for any samples that did not have plasmids
    all_dirs = [os.path.join(base_dir, x) for x in file_dirs]
    files_in_df = [x for x in list(grouped_df["filename"].unique())]
    for d in all_dirs:
        all_files = [os.path.join(d, x) for x in os.listdir(d) if (x.split(".")[-1] == "fasta")]
        missing_files = [x for x in all_files if (x not in files_in_df) and ("plasmid.fasta" not in x)]
        for f in missing_files:
            new_row = {'filename': f,
                       'file_type': 'n/a',
                       'contig': 'n/a',
                       'contig_len': 1000000,
                       'plasmid_name': 'no plasmid found',
                       'plasmid_accession': 'n/a',
                       'blast_plasmids': []}
            grouped_df = pd.concat([grouped_df, pd.DataFrame([new_row])])

            new_row = {'filename': f,
                       'counts': 0}
            final_count_df = pd.concat([final_count_df, pd.DataFrame([new_row])])

    # write the outputs
    grouped_df.drop(labels="blast_plasmids", axis=1).to_csv(os.path.join(output_dir, "plasmidcount_results.csv"))
    print(output_dir)
    final_count_df.to_csv(os.path.join(output_dir, "plasmidcount_aggregated.csv"))

    return grouped_df, final_count_df


def get_parser():
    parser = argparse.ArgumentParser(description="Count plasmids in long/short read sequences.")
    parser.add_argument("--python_path", required=True, help="Full path to the Python interpreter for running plasmidfinder.")
    parser.add_argument("--plasmidfinder_path", required=True, help="Full path to the plasmidfinder.py script.")

    parser.add_argument("--base_dir", required=True, help="The base directory to run plasmidcounts.")
    parser.add_argument("--file_dir", action="append", required=True,
                        help="The directories where files are stored; can be specified multiple times.")
    parser.add_argument("--read_type", action="append", required=True,
                        help="The read type for each directory in file_dirs; specify once for each file_dir.")
    parser.add_argument("--output_base", required=True, help="The base path to write output summary files like plasmidfinder_results.csv")
    parser.add_argument("--output_dir", required=True, help="The directory to write outputs.")
    parser.add_argument("--db_path", required=True, help="The path to the plasmid database.")
    parser.add_argument("--tmp_path", action="store", default=None,
                        help="the path to store temporary files; defaults to <output_dir>/tmp_dir/")
    parser.add_argument("--contig_len_cutoff", action="store", default=500,
                        help="the minimum length of contigs to keep")
    parser.add_argument("--mash_cutoff", action="store", default=0.05,
                        help="mash distance to consider two contigs duplicate")
    parser.add_argument("--plasfinder_cov", action="store", default=80,
                        help="coverage threshold for plasmidfinder results")
    parser.add_argument("--plasfinder_pid", action="store", default=90,
                        help="pct ID threshold for plasmidfinder results")
    parser.add_argument("--read_threshold", action="store", default=5,
                        help="threshold for number of reads in a contig (must be contained in fasta header)")
    parser.add_argument("--cover_threshold", action="store", default=1, help="threshold for coverage of a contig")
    parser.add_argument("--blast_id_cutoff", action="store", default=90,
                        help="cutoff for plasmid identity in PLSdb blast")
    parser.add_argument("--blast_len_thresh", action="store", default=27000,
                        help="threshold for contig length to blast")
    parser.add_argument("--overlap_cutoff_short", action="store", default=0.01,
                        help="cutoff for overlap of blast results for short contigs")
    parser.add_argument("--overlap_cutoff_long", action="store", default=0.15,
                        help=" cutoff for overlap of blast results for long contigs")
    parser.add_argument("--blast_len_cutoff", action="store", default=1000,
                        help="cutoff for plasmid length in PLSdb blast")
    parser.add_argument("--filter_col", action="store_true", default=False, help="flag to filter Col plasmids")
    parser.add_argument("--filter_std", action="store_true", default=False, help="flag to filter high-variance strains")
    parser.add_argument("--num_contig_thresh", action="store", default=3,
                        help="only filter high-variance strains with at least this many contigs")
    parser.add_argument("--plasfinder_loc", action="store", help="location of plasmidfinder executable, if not on path")
    parser.add_argument("--mash_loc", action="store", help="location of mash executable, if not on path")
    parser.add_argument("--skip_interactive", action="store_true", default=False,
                        help="skip interactive prompts while running")

    return parser


def main():
    args = get_parser().parse_args()
    plasmidcounts(args.python_path,
                  args.plasmidfinder_path,
                  args.output_base,
                  args.base_dir,
                  args.file_dir,
                  args.read_type,
                  args.output_dir,
                  args.db_path,
                  tmp_path=args.tmp_path,
                  contig_len_cutoff=args.contig_len_cutoff,
                  mash_cutoff=args.mash_cutoff,
                  plasfinder_cov=args.plasfinder_cov,
                  plasfinder_pid=args.plasfinder_pid,
                  read_threshold=args.read_threshold,
                  cover_threshold=args.cover_threshold,
                  blast_id_cutoff=args.blast_id_cutoff,
                  blast_len_thresh=args.blast_len_thresh,
                  overlap_cutoff_short=args.overlap_cutoff_short,
                  overlap_cutoff_long=args.overlap_cutoff_long,
                  blast_len_cutoff=args.blast_len_cutoff,
                  filter_col=args.filter_col,
                  filter_std=args.filter_std,
                  num_contig_thresh=args.num_contig_thresh,
                  skip_interactive=args.skip_interactive)

# python_path, plasmidfinder_path, output_base, base_dir, file_dirs, read_type, output_dir, db_path,
if __name__ == "__main__":
    main()
