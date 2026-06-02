import os
import json
import pandas as pd
import matplotlib.pyplot as plt
from math import ceil, floor
import math
import numpy as np
import re

def b_to_B(x):
    return x // 8

def b_to_KB(x):
    return x // 8192

def b_to_MB(x):
    return x // (1024 * 1024)

def b_to_GB(x):
    return x // (1024 * 1024 * 1024)

def B_to_KB(x):
    return np.round(x / 1024, 2)

def B_to_MB(x):
    return x // (1024 * 1024)

def B_to_GB(x):
    return x // (1024 * 1024 * 1024)

def MB_to_GB(x):
    return x / 1024

def read_and_flatten(results_dir, name=None):
    data = []
    # Iterate through all files in the "results" directory
    for filename in os.listdir(results_dir):
        if filename.endswith(".json"):

            # Read the content of the file
            with open(os.path.join(results_dir, filename), 'r') as file:
                log_data = json.load(file)
                
                # Add prefixes to the keys of "offline" and "online" data
                offline_data = {f"offline-{key}": value for key, value in log_data["offline"].items()}
                online_data = {f"online-{key}": value for key, value in log_data["online"].items()}
                
                # Flatten the nested JSON structure for easier analysis
                if "specs" in log_data:
                    flattened_data = {
                        **log_data['specs'],
                        **offline_data,
                        **online_data,
                    }
                else:
                    flattened_data = {
                        **offline_data,
                        **online_data,
                    }
                # Append the flattened data to the list
                data.append(flattened_data)
    # Create a DataFrame from the collected data
    df = pd.DataFrame(data)
    
    #create column called 'name' and set it all to the arguement given in the inputs
    if name:
        df['name'] = name
        
    return preprocess(df)

def preprocess(df):
    """
    Applies common and specific preprocessing steps to a DataFrame,
    checking for column existence before applying transformations.
    """

    if 'offline-serverTimeMs' in df.columns and 'offline-encodeTimeMs' in df.columns:
        df['offline-totalTimeMs'] = df['offline-serverTimeMs'] + df['offline-encodeTimeMs']
    elif 'offline-serverTimeMs' in df.columns:
        df['offline-totalTimeMs'] = df['offline-serverTimeMs']
    else:
        df['offline-totalTimeMs'] = 0 

    if 'offline-totalTimeMs' in df.columns:
        df['offline-totalTimeS'] = (df['offline-totalTimeMs']/1000).apply(lambda x: np.round(x))

    if 'inputDatabaseSizeMb' in df.columns and 'online-serverTimeMs' in df.columns:
        df['online-throughputGBs'] = (1000 * df['inputDatabaseSizeMb'] / (1024 * df['online-serverTimeMs'])).apply(lambda x: np.round(x))
        df['online-throughputMBs'] = (1000 * df['inputDatabaseSizeMb'] / df['online-serverTimeMs']).apply(lambda x: np.round(x, -1))

    if 'online-rgswTimeUs' in df.columns:
        df["online-rgswTimeMs"] = df["online-rgswTimeUs"] // 1000

    if 'online-firstPassTimeUs' in df.columns:
        df["online-firstPassTimeMs"] = (df["online-firstPassTimeUs"] / 1000).apply(lambda x: round(x, -1) if x < 1000 else round(x, -2))

    if 'online-firstPackTimeUs' in df.columns:
        df["online-firstPackTimeMs"] = (df["online-firstPackTimeUs"] / 1000).apply(lambda x: round(x, -1))

    # These are primarily for inspire-pareto, but common preprocess should handle them if they exist
    if 'online-secondPassTimeUs' in df.columns:
        df["online-secondPassTimeMs"] = (df["online-secondPassTimeUs"] / 1000).apply(lambda x: round(x, 1))
    if 'online-packingKeyRotationsTimeUs' in df.columns:
        df["online-packingKeyRotationsTimeMs"] = (df["online-packingKeyRotationsTimeUs"] / 1000).apply(round)
    if 'online-packingMaskOnlineTimeUs' in df.columns:
        df["online-packingMaskOnlineTimeMs"] = (df["online-packingMaskOnlineTimeUs"] / 1000).apply(round)
    if 'online-packingBodyOnlineTimeUs' in df.columns:
        df["online-packingBodyOnlineTimeMs"] = (df["online-packingBodyOnlineTimeUs"] / 1000).apply(round)

    # These are common to both dataframes after the user's latest preprocess structure
    if 'online-serverTimeMs' in df.columns:
        df["online-serverTimeMs"] = (df["online-serverTimeMs"]).apply(lambda x: round(x, -1))
    if 'online-totalBytesKB' in df.columns:
        df["online-totalBytesKB"] = (df["online-totalBytesKB"]).apply(lambda x: round(x))

    for col in [
        'online-uploadKeys',
        'online-uploadQuery',
        'online-uploadBytes',
        'online-downloadBytes',
        'online-totalBytes',
    ]:
        if col in df.columns:
            df[f'{col}KB'] = np.round(B_to_KB(df[col]))
    if 'offline-serverTimeMs' in df.columns: 
        df['offline-serverTimeS'] = df['offline-serverTimeMs'].apply(lambda x: np.round(x/1000, 1))
    # df['online-throughputGBs'] = (MB_to_GB(df['inputDatabaseSizeMb']) / (df['online-serverTimeMs']/1000)).apply(lambda x: np.round(x, 4)) 
    # df['online-throughputMBs'] = (df['inputDatabaseSizeMb'] / (df['online-serverTimeMs']/1000)).apply(lambda x: np.round(x)) 


    return df

def pareto_frontier(df, x_col, y_col, tolerance=0.00):
  """Removes rows that are not on or within tolerance of the Pareto frontier.
  
  This version correctly handles ties by keeping the point with the minimum x_col
  for any given y_col on the frontier.

  Args:
      df (pd.DataFrame): The DataFrame to process.
      x_col (str): The name of the column to minimize.
      y_col (str): The name of the column to minimize.
      tolerance (float): The tolerance percentage for including points near the
        Pareto frontier. Defaults to 0.0.

  Returns:
      pd.DataFrame: A new DataFrame with Pareto frontier rows and rows within
      tolerance.
  """
  # Create a copy to avoid modifying the original DataFrame
  df_copy = df.copy()
  
  # --- Step 1: Sort values and drop duplicates ---
  # Drop any fully identical rows first
  df_copy = df_copy.drop_duplicates(subset=[x_col, y_col])
  # Sort by the column you want to minimize first
  df_copy = df_copy.sort_values(by=x_col)

  # --- Step 2: Find the initial frontier candidates (your original logic) ---
  df_copy["pareto_front"] = df_copy[y_col].cummin()
  df_copy["tolerance_threshold"] = df_copy["pareto_front"] * (1 + tolerance)
  pareto_candidates = df_copy[df_copy[y_col] <= df_copy["tolerance_threshold"]]

  # --- Step 3: Refine the frontier to handle ties in the y-dimension ---
  # For any points with the same y_col value, we only want to keep the one 
  # with the minimum x_col value.
  # We use groupby() on the y_col and find the index of the minimum x_col.
  idx_to_keep = pareto_candidates.groupby(y_col)[x_col].idxmin()
  
  final_pareto_df = pareto_candidates.loc[idx_to_keep]

  return final_pareto_df.sort_values(by=x_col).drop(columns=["pareto_front", "tolerance_threshold"])


def transpose_and_round(df_to_transpose, int_columns):
    df_final_transposed = pd.DataFrame(index=df_to_transpose.columns, columns=df_to_transpose.index)
    # int_columns = [
    #     'gamma0', 'gamma1', 'gamma2', 'payloadB', 'resizedDbFirstDim',
    #     'online-firstPassTimeMs', 'online-firstPackTimeMs',
    #     'online-packingKeyRotationsTimeMs', 'online-packingMaskOnlineTimeMs',
    #     'online-packingBodyOnlineTimeMs', 'online-serverTimeMs'
    # ]

    # Iterate through each column of the original DataFrame.
    # In the transposed DataFrame, these columns will become rows.
    for col_name in df_to_transpose.columns:
        # Get the series (column) from the original DataFrame
        series_to_transpose = df_to_transpose[col_name]

        if col_name in int_columns:
            if all(val == int(val) for val in series_to_transpose if pd.notna(val)):
                df_final_transposed.loc[col_name] = series_to_transpose.astype(int)
            else:
                df_final_transposed.loc[col_name] = series_to_transpose
        else:
            df_final_transposed.loc[col_name] = series_to_transpose

    # Print the final transposed DataFrame with the corrected data types
    return df_final_transposed
