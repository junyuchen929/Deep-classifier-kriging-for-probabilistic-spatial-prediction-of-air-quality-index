import numpy as np
import pandas as pd
from tqdm import tqdm
import os

regions = ["CA", "NE", "SE"]

def minmax_normalize(coords):
    min_vals = coords.min(axis=0)
    max_vals = coords.max(axis=0)
    return (coords - min_vals) / (max_vals - min_vals), min_vals, max_vals

def compute_phi_matrix(coords, num_basis_list):
    N = coords.shape[0]
    phi = np.zeros((N, sum(num_basis_list)))
    knots_1d = [np.linspace(0, 1, int(np.sqrt(i))) for i in num_basis_list]
    K = 0

    for res in range(len(num_basis_list)):
        theta = 1 / np.sqrt(num_basis_list[res]) * 2.5
        knots_s1, knots_s2 = np.meshgrid(knots_1d[res], knots_1d[res])
        knots = np.column_stack((knots_s1.flatten(), knots_s2.flatten()))

        for i in range(num_basis_list[res]):
            d = np.linalg.norm(coords - knots[i, :], axis=1) / theta
            mask = (d >= 0) & (d <= 1)
            phi[mask, i + K] = ((1 - d[mask]) ** 6) * (35 * d[mask] ** 2 + 18 * d[mask] + 3) / 3

        K += num_basis_list[res]
    return phi

def generate_embeddings_for_region(reg):
    print(f"\n===============================")
    print(f"Generating embeddings for {reg}")
    print("===============================\n")

    df_train = pd.read_csv(f"projection_matrix_{reg}_train.csv")
    df_test = pd.read_csv(f"projection_matrix_{reg}_test.csv")

    coords_train = df_train[["uj_lat", "uj_lon"]].values
    coords_test = df_test[["uj_lat", "uj_lon"]].values
    all_coords = np.vstack([coords_train, coords_test])

    # Normalize
    all_coords_norm, coord_min, coord_max = minmax_normalize(all_coords)
    coords_train_norm = all_coords_norm[:len(df_train)]
    coords_test_norm = all_coords_norm[len(df_train):]

    # Same basis resolution for all regions
    num_basis = [5**2, 9**2, 15**2, 21**2]  

    print(f"Computing φ_train for {reg}...")
    phi_train = compute_phi_matrix(coords_train_norm, num_basis)

    print(f"Computing φ_test for {reg}...")
    phi_test = compute_phi_matrix(coords_test_norm, num_basis)

    # Remove all-zero columns
    idx_zero = np.where(np.all(phi_train == 0, axis=0))[0]
    phi_train = np.delete(phi_train, idx_zero, axis=1)
    phi_test = np.delete(phi_test, idx_zero, axis=1)

    print(f"{reg} final φ shape: {phi_train.shape}")

    os.makedirs(f"embeddings_{reg}", exist_ok=True)
    np.save(f"embeddings_{reg}/phi_train_{reg}.npy", phi_train.astype(np.float16))
    np.save(f"embeddings_{reg}/phi_test_{reg}.npy", phi_test.astype(np.float16))
    np.save(f"embeddings_{reg}/idx_zero_{reg}.npy", idx_zero)

    print(f"✔ Saved φ embeddings for region: {reg}")

if __name__ == '__main__':
    for reg in regions:
        generate_embeddings_for_region(reg)

    print("\n🎯 All region embeddings successfully generated!")
