import argparse
from multiprocessing.pool import MaybeEncodingError
import requests

import numpy as np
import pandas as pd
from sklearn.metrics import pairwise
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

from opencadd.databases.klifs import setup_remote
from kissim.api import encode

def parse_args():
    parser = argparse.ArgumentParser(description="Rank kinase off-targets based on similarity to known targets.")
    
    parser.add_argument("-j", type=int, default=1, help="Number of parallel jobs to run.")
    parser.add_argument('--use_save', type=str, default=None, help="Use saved fingerprints if available.")
    parser.add_argument("--target", type=str, required=True, help="Output file path for ranked off-targets.")
    parser.add_argument("--output", type=str, default="ranked_offtargets.csv", help="Output file path for ranked off-targets.")
    
    return parser.parse_args()

def download_kinase_data(args):
    r = requests.get("https://klifs.net/api/kinase_names?species=HUMAN")
    r.raise_for_status()
    data = r.json()
    all_ids = [record["kinase_ID"] for record in data]
    
    r = requests.get(f"https://klifs.net/api/kinase_information?kinase_ID={','.join(map(str, all_ids))}")
    r.raise_for_status()
    data = r.json()
    
    items = []
    for record in data:
        items.append({
            "KLIFS Name": record["name"],
            "HGNC Name": record["HGNC"],
            "Family": record["family"],
            "Group": record["group"],
            "Class": record["kinase_class"],
            "Kinase Name": record["full_name"],
            "UniprotID": record["uniprot"],
            "IUPHAR": record["iuphar"],
            "Pocket": record["pocket"],
            "KLIFS ID": record["kinase_ID"]
        })
    klifs = pd.DataFrame.from_dict(items).sort_values(by="KLIFS Name")
    
    cols = ['HGNC Name', 'KLIFS Name', 'UniprotID', 'Group', 'Kinase Name', 'Pocket']
    kinase_selection_df = klifs[cols]
    kinase_selection_df.columns = ['kinase', 'kinase_klifs', 'uniprot_id', 'group', 'full_kinase_name', 'pocket']
    
    return kinase_selection_df

def preprocess_kinase_data(kinase_selection_df):
    klifs_session = setup_remote()
    
    # Get list of kinase names
    kinase_names = kinase_selection_df["kinase_klifs"].to_list()

    # Get all available structures for these kinases
    structures_df = klifs_session.structures.by_kinase_name(kinase_names=kinase_names)
    print(f"Number of structures: {len(structures_df)}")
    print("Kinases:", *structures_df["kinase.klifs_name"].unique())

    structures_df = structures_df[
        (structures_df["species.klifs"] == "Human")
        & (structures_df["structure.dfg"] == "in")
        & (structures_df["structure.resolution"] <= 3)
        & (structures_df["structure.qualityscore"] >= 6)
    ]
    print(f"Number of structures: {len(structures_df)}")
    print("Kinases:", *structures_df["kinase.klifs_name"].unique())

    structure_klifs_ids = structures_df["structure.klifs_id"].to_list()
    print(f"Number of structures: {len(structure_klifs_ids)}")
    
    return structure_klifs_ids

def cal_KiSSim_fps(structure_klifs_ids, N_CORES=1):
    kissim_fingerprints_df_list = []
    error_ids = []

    for i in tqdm(range(len(structure_klifs_ids)//N_CORES+1)):
        try:
            kissim_fingerprints = encode(structure_klifs_ids[N_CORES*i:N_CORES*(i+1)], n_cores=N_CORES)
        except ValueError as e:
            print(i, N_CORES*i, N_CORES*(i+1), e)
            error_ids.extend(structure_klifs_ids[N_CORES*i:N_CORES*(i+1)])
            continue
        except MaybeEncodingError as e:
            print(i, N_CORES*i, N_CORES*(i+1), e)
            error_ids.extend(structure_klifs_ids[N_CORES*i:N_CORES*(i+1)])
            continue
        
        # Save fingerprints in csv file
        _structure_klifs_ids = list(kissim_fingerprints.data.keys())
        kissim_fingerprints_array = [
            fingerprint.values_array().tolist()
            for _structure_klifs_id, fingerprint in kissim_fingerprints.data.items()
        ]
        
        kissim_fingerprints_array = np.array(kissim_fingerprints_array)
        kissim_fingerprints_df = pd.DataFrame(kissim_fingerprints_array, index=_structure_klifs_ids)
        kissim_fingerprints_df_list.append(kissim_fingerprints_df)

    print('Try again for error ids')
    for i in tqdm(error_ids):
        try:
            kissim_fingerprints = encode([i], n_cores=1)
        except ValueError as e:
            print(i, e)
            continue
        except MaybeEncodingError as e:
            print(i, e)
            continue
        
        # Save fingerprints in csv file
        _structure_klifs_ids = list(kissim_fingerprints.data.keys())
        kissim_fingerprints_array = [
            fingerprint.values_array().tolist()
            for _structure_klifs_id, fingerprint in kissim_fingerprints.data.items()
        ]
        
        kissim_fingerprints_array = np.array(kissim_fingerprints_array)
        kissim_fingerprints_df = pd.DataFrame(kissim_fingerprints_array, index=_structure_klifs_ids)
        kissim_fingerprints_df_list.append(kissim_fingerprints_df)
    
    kissim_fingerprints_df = pd.concat(kissim_fingerprints_df_list, axis=0)
    print(f"Matrix shape: {kissim_fingerprints_df.shape}")
    print(f"Number of fingerprints: {kissim_fingerprints_df.shape[0]}")
    print(f"Number of fingerprint bits: {kissim_fingerprints_df.shape[1]}")
    
    return kissim_fingerprints_df

def standarize_fingerprints(kissim_fingerprints_df):
    scaler = StandardScaler()
    kissim_fingerprints_df_z = pd.DataFrame(scaler.fit_transform(kissim_fingerprints_df), index=kissim_fingerprints_df.index, columns=kissim_fingerprints_df.columns)
    return kissim_fingerprints_df_z

def compare_structures(kissim_fingerprints_df_z):
    structure_distance_matrix_array = pairwise.nan_euclidean_distances(kissim_fingerprints_df_z.values)
    
    # Create DataFrame with structure KLIFS IDs as index/columns
    structure_klifs_ids = kissim_fingerprints_df.index.to_list()
    structure_distance_matrix_df = pd.DataFrame(
        structure_distance_matrix_array, index=structure_klifs_ids, columns=structure_klifs_ids
    )
    print(f"Structure distance matrix size: {structure_distance_matrix_df.shape}")
    
    return structure_distance_matrix_df

def map_structure(structure_distance_matrix_df):
    # Copy distance matrix to kinase matrix
    kinase_distance_matrix_df = structure_distance_matrix_df.copy()
    # Replace structure KLIFS IDs with the structures' kinase names
    kinase_names = structures_df.set_index("structure.klifs_id").loc[
        structure_klifs_ids, "kinase.klifs_name"
    ]
    kinase_distance_matrix_df.index = kinase_names
    kinase_distance_matrix_df.columns = kinase_names
    
    # We unstack the matrix (each pairwise comparison in a single row)
    # We group by kinase names (level=[0, 1] ensures that the order of the kinases is ignored
    # We take the minimum value in each kinase pair group
    # We unstack the remaining data points
    kinase_distance_matrix_df = (
        kinase_distance_matrix_df.unstack().groupby(level=[0, 1]).min().unstack(level=1)
    )
    # Cosmetics: Remove the index and column names
    kinase_distance_matrix_df.index.name = None
    kinase_distance_matrix_df.columns.name = None
    
    print(
        f"Structure matrix of shape {structure_distance_matrix_df.shape} "
        f"reduced to kinase matrix of shape {kinase_distance_matrix_df.shape}."
    )
    
    return kinase_distance_matrix_df

def rank_kinase_offtargets(kinase_distance_matrix_df, target, output):
    df = kinase_distance_matrix_df.loc[target].sort_values()
    df.name = 'distance'
    df.index.name = 'kinase'
    df.to_csv(output)
    
def main(args):
    if args.use_save is None:
        kinase_selection_df = download_kinase_data(args)
        structure_klifs_ids = preprocess_kinase_data(kinase_selection_df)
        kissim_fingerprints_df = cal_KiSSim_fps(structure_klifs_ids, N_CORES=args.j)
        kissim_fingerprints_df_z = standarize_fingerprints(kissim_fingerprints_df)
        structure_distance_matrix_df = compare_structures(kissim_fingerprints_df_z)
        kinase_distance_matrix_df = map_structure(structure_distance_matrix_df)
    else:
        kinase_distance_matrix_df = pd.read_csv(args.use_save, index_col=0)
        print(f"Loaded kinase distance matrix of shape {kinase_distance_matrix_df.shape} from {args.use_save}.")
        
    rank_kinase_offtargets(kinase_distance_matrix_df, args.target, args.output)
    
if __name__ == "__main__":
    args = parse_args()
    main(args)