import numpy as np
import pandas as pd
from tqdm import tqdm
import os

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

def main():
    # Load datasets
    df_train = pd.read_csv("projection_matrix_train.csv")
    df_test = pd.read_csv("projection_matrix_test.csv")

    coords_train = df_train[["uj_lat", "uj_lon"]].values
    coords_test = df_test[["uj_lat", "uj_lon"]].values
    all_coords = np.vstack([coords_train, coords_test])

    # Normalize coordinates
    all_coords_norm, coord_min, coord_max = minmax_normalize(all_coords)
    coords_train_norm = all_coords_norm[:len(df_train)]
    coords_test_norm = all_coords_norm[len(df_train):]

    # Set basis function resolutions
    num_basis = [5**2, 7**2, 11**2]

    print("Generating phi_train...")
    phi_train = compute_phi_matrix(coords_train_norm, num_basis)

    print("Generating phi_test...")
    phi_test = compute_phi_matrix(coords_test_norm, num_basis)

    # Remove all-zero columns (same for train and test)
    print("Removing all-zero columns...")
    idx_zero = np.where(np.all(phi_train == 0, axis=0))[0]
    phi_train = np.delete(phi_train, idx_zero, axis=1)
    phi_test = np.delete(phi_test, idx_zero, axis=1)
    print(f"Final phi shape: {phi_train.shape}")

    # Save all outputs
    os.makedirs("embeddings", exist_ok=True)
    np.save("embeddings/phi_train.npy", phi_train.astype(np.float16))
    np.save("embeddings/phi_test.npy", phi_test.astype(np.float16))

    print("Embedding matrices saved to 'embeddings/' folder.")

if __name__ == '__main__':
    main()
