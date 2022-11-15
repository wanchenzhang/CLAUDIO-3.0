import pandas as pd
import click
import os


@click.command()
@click.option("-i", "--input-filepath", default="data/out/dist_reeval/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv.sqcs.csv")
def main(input_filepath):
    data = pd.read_csv(input_filepath, index_col=0)
    num_of_nans = sum(pd.isna(data["eucl_dist"]))
    print("Number of Nans in dist column:", num_of_nans)
    print("Number of Nans in topolink dist column:", sum(pd.isna(data["eucl_dist_tplk"]) & ~pd.isna(data["eucl_dist"])))
    print("Number of datapoints with dist = 0:", sum(data["eucl_dist"] == 0))
    print("Number of computed distances in dist column:", len(data.index) - num_of_nans)
    print(data.columns)
    print(data.head)
    print(data.eucl_dist.values)
    print(f"number computed eucl_dist: {len(data[~pd.isna(data.eucl_dist)].index)} of {len(data[data.path != '-'].index)}")
    print(f"number fulfilled criteria: {len(data[data.lysin_criteria_fulfilled].index)} of {len(data[data.path != '-'].index)}")
    print(f"number of unique proteins: {len(data[data.lysin_criteria_fulfilled]['unip_id_a'].unique())}")
    print(f"number of new inter interactions: {len(data[data.new_XL_type == 'inter'].index)}")
    project_path = '/'.join(os.path.abspath(__file__).split('/')[:-3])
    print(f"number of used methods:\n\talphafold:\n\t\ta:{len(data[data.method_a == 'alphafold'].index)}\n\t\t\t"
          f"dist_true:{len(data[(data.method_a == 'alphafold') & (~pd.isna(data.eucl_dist))].index)}\n\t\t\t\t"
          f"lys_c_true:{len(data[(data.method_a == 'alphafold') & data.lysin_criteria_fulfilled].index)}\n\t\t"
          f"b:{len(data[data.method_b == 'alphafold'].index)}\n\t\t\t"
          f"dist_true:{len(data[(data.method_b == 'alphafold') & (~pd.isna(data.eucl_dist))].index)}\n\t\t\t\t"
          f"lys_c_true:{len(data[(data.method_b == 'alphafold') & data.lysin_criteria_fulfilled].index)}\n\t"
          f"realigning:\n\t\ta:{len(data[data.method_a == 'realigning'].index)}\n\t\t\t"
          f"dist_true:{len(data[(data.method_a == 'realigning') & (~pd.isna(data.eucl_dist))].index)}\n\t\t\t\t"
          f"lys_c_true:{len(data[(data.method_a == 'realigning') & data.lysin_criteria_fulfilled].index)}\n\t\t"
          f"b:{len(data[data.method_b == 'realigning'].index)}\n\t\t\t"
          f"dist_true:{len(data[(data.method_b == 'realigning') & (~pd.isna(data.eucl_dist))].index)}\n\t\t\t\t"
          f"lys_c_true:{len(data[(data.method_b == 'realigning') & data.lysin_criteria_fulfilled].index)}\n\t"
          f"pdb_chain_uniprot:\n\t\ta:{len(data[data.method_a == 'pdb_chain_uniprot'].index)}\n\t\t\t"
          f"dist_true:{len(data[(data.method_a == 'pdb_chain_uniprot') & (~pd.isna(data.eucl_dist))].index)}\n\t\t\t\t"
          f"lys_c_true:{len(data[(data.method_a == 'pdb_chain_uniprot') & data.lysin_criteria_fulfilled].index)}\n\t\t"
          f"b:{len(data[data.method_a == 'pdb_chain_uniprot'].index)}\n\t\t\t"
          f"dist_true:{len(data[(data.method_a == 'pdb_chain_uniprot') & (~pd.isna(data.eucl_dist))].index)}\n\t\t\t\t"
          f"lys_c_true:{len(data[(data.method_b == 'pdb_chain_uniprot') & data.lysin_criteria_fulfilled].index)}")
    print(f"\n\nSuspicious datapoints:\n\talphafold + lysin_crit_false (total number = {len(data[(data.method_a == 'alphafold') & ~data.lysin_criteria_fulfilled].index)}):\n\t\t"
          f"{data[(data.method_a == 'alphafold') & ~data.lysin_criteria_fulfilled][['unip_id_a', 'pdb_id', 'chain', 'method_a', 'pub', 'lysin_criteria_fulfilled', 'eucl_dist', 'XL_type', 'new_XL_type']]}\n\t(publication + lysin_criteria_true)/publication:\n\t\t"
          f"Schweppe et al. 2017: {round((len(data[(data.pub == 'Schweppe et al. 2017') & data.lysin_criteria_fulfilled].index) / (len(data[data.pub == 'Schweppe et al. 2017'].index))) * 100, 2)}%"
          f"\n\t\tLiu et al. 2018: {round((len(data[(data.pub == 'Liu et al. 2018') & data.lysin_criteria_fulfilled].index) / (len(data[data.pub == 'Liu et al. 2018'].index))) * 100, 2)}%"
          f"\n\tLys_criteria unfullfilled (total number = {len(data[~data.lysin_criteria_fulfilled].index)}): {data[~data.lysin_criteria_fulfilled][['unip_id_a', 'pos_a', 'pos_b', 'lys_crit_a', 'lys_crit_b']]}"
          f"\n\tNo distance computed (total number = {len(data[pd.isna(data.eucl_dist) & (data.path != '-')].index)}): {data[pd.isna(data.eucl_dist) & (data.path != '-')][['unip_id_a', 'pos_a', 'pos_b', 'eucl_dist']]}"
          f"\n\tNot found in databases (total number = {len(data[(data.path == '-')].index)}): {data[(data.path == '-')][['unip_id_a']]}"
          f"\n\tVery high distances (total number = {len(data[data.eucl_dist >= 100].index)}): {data[data.eucl_dist >= 100][['unip_id_a', 'eucl_dist']]}"
          f"\n\tNo distance computed but lysin_crit_true (total number = {len(data[pd.isna(data.eucl_dist) & data.lysin_criteria_fulfilled].index)}): {data[pd.isna(data.eucl_dist) & data.lysin_criteria_fulfilled][['pdb_id', 'pos_a', 'pos_b', 'eucl_dist', 'lysin_criteria_fulfilled']]}"
          f"\n\tVery low distances (total number = {len(data[(0 < data.eucl_dist) & (data.eucl_dist < 5)].index)}): {data[(0 < data.eucl_dist) & (data.eucl_dist < 5)][['unip_id_a', 'eucl_dist']]}"
          f"\n\tTopolink is Nan but self computed is not zero or Nan (total number = {sum((~(data.eucl_dist == 0) & ~(pd.isna(data.eucl_dist))) & (pd.isna(data.eucl_dist_tplk)))}): {data[(~(data.eucl_dist == 0) & ~(pd.isna(data.eucl_dist))) & (pd.isna(data.eucl_dist_tplk))][['unip_id_a', 'pdb_pos_a', 'pdb_pos_b', 'pdb_id', 'eucl_dist', 'eucl_dist_tplk']]}")
    final_data = pd.read_csv(f'{project_path}/data/out/new_inter/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant_final.csv', index_col=0)
    print(f"\n\nFinal results:"
          f"\n\ttotal number of interactions: {len(final_data.index)}"
          f"\n\ttotal number of proteins: {len(final_data.unip_id_a.unique())}"
          f"\n\tnumber of new inters: {sum(final_data.final_XL_type == 'inter')}"
          f"\n\tnumber of discovered homomers: {len(os.listdir(f'{project_path}/data/out/new_inter/homomers'))}")

# nohup sh -c "python3 start01b_liu_uniprot_search.py -s True && python3 claudio_mod02_02.py -i data/out/uniprot_search/supp_RA117.000470_133922_0_supp_23973_3zf3c5_table1.csv.sqcs -s True && python3 claudio_mod02_03.py -i2 data/out/uniprot_search/supp_RA117.000470_133922_0_supp_23973_3zf3c5_table1.csv.sqcs && python3 read_out.py" > logs/log111121_0344.log 2>&1 &

# nohup sh -c "python3 claudio_mod02_01.py && python3 claudio_mod02_02.py && python3 claudio_mod02_03.py && python3 read_out.py" > logs/log220323_1400.log 2>&1 &

# nohup sh -c "python3 claudio_mod02_02.py && python3 claudio_mod02_03.py && python3 read_out.py" > logs/log220204_1300.log 2>&1 &

# nohup sh -c "python3 claudio_mod02_03.py && python3 read_out.py" > logs/log220207_2000.log 2>&1 &

# nohup sh -c "python3 claudio_lists.py" > nohup.log 2>&1 &

# nohup sh -c "python3 claudio_ops.py" > nohup.log 2>&1 &

# nohup sh -c "python3 claudio_structdi.py && python3 read_out.py -i 'data/out/module02/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv.sqcs.csv'" > logs/log220413_1330.log 2>&1 &

# nohup sh -c "python3 claudio_lists.py -o 'data/out/full' && python3 claudio_structdi.py -o 'data/out/full' && python3 claudio_ops.py -o 'data/out/full' && python3 claudio_xl.py -i 'data/out/full/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv.sqcs.csv' -i2 'data/out/full/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv_homosig.csv' -o 'data/out/full'" > logs/log220505_1230.log 2>&1 &

# nohup sh -c "python3 claudio.py -rt True" > logs/log221109_1830.log 2>&1 &
