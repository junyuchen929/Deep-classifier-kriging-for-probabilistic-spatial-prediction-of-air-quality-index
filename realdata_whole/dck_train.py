import numpy as np
import pandas as pd
import time
#from keras.utils import np_utils
from keras.utils import to_categorical

from keras.models import load_model
from model_function import model_function  # 假设你把模型结构封装到了一个文件里

# === Load data ===
df_train = pd.read_csv("projection_matrix_train.csv")
df_test = pd.read_csv("projection_matrix_test.csv")

df_train["class"] = df_train["class"] - 1
#dummy_y = np_utils.to_categorical(df_train["class"])
dummy_y = to_categorical(df_train["class"])
n = dummy_y.shape[1]
print('Total number of classes: %d' % n)

# === Load precomputed embedding matrices ===
phi_train = np.load("embeddings/phi_train.npy")  # shape: (N_train, D)
phi_test = np.load("embeddings/phi_test.npy")    # shape: (N_test, D)

print("Training embedding shape:", phi_train.shape)
print("Test embedding shape:", phi_test.shape)

# === Train the model ===
time_records = []
train_start = time.time()
model = model_function(df_train, phi_train, dummy_y, n)
train_end = time.time()

# === Predict ===
pred = model.predict(phi_test)
pred_df = pd.DataFrame(pred)
df_test_preds = pd.concat([df_test.reset_index(drop=True), pred_df], axis=1)

# === Save predictions ===
df_test_preds.to_csv("predictions_results.csv", index=False)

# === Save training time ===
df_times = pd.DataFrame([{"train_time_sec": train_end - train_start}])
df_times.to_csv("DNN_time_records.csv", index=False)
print("Average training time (sec):", df_times['train_time_sec'].mean())
