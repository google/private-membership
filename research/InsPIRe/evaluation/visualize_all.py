from visualize_base import *

# Directory containing log files
results_dir_ypir = f"results/ypir"
results_dir_simpleypir = f"results/simpleypir"
results_dir_inspire_0 = f"results/inspire_0"
results_dir_rgswpir = f"results/inspire"
results_dir_inspire_squared = f"results/inspire-squared"
results_dir_kspir = f"results/kspir"
results_dir_hintlesspir = f"results/hintlesspir"

figures_dir_all = f"figures/all"
os.makedirs(figures_dir_all, exist_ok=True)

df_inspire_0_raw = read_and_flatten(results_dir=results_dir_inspire_0, name='inspire')
df_inspire_0_raw['name'] = 'inspire_0'
df_inspire_0_raw = df_inspire_0_raw[df_inspire_0_raw["gamma0"] == 2048]

# Read and preprocess inspire-pareto data
df_inspire_squared_raw = read_and_flatten(results_dir=results_dir_inspire_squared, name='inspire^2')
# add column names if they don't exists
for col in df_inspire_0_raw.columns:
    if col not in df_inspire_squared_raw.columns:
        df_inspire_squared_raw[col] = 0
df_inspire_squared_raw = df_inspire_squared_raw[df_inspire_squared_raw["gamma1"] == 1024]
df_inspire_squared_raw = df_inspire_squared_raw[df_inspire_squared_raw["gamma0"] == df_inspire_squared_raw["gamma2"]]
df_inspire_squared_raw['payloadB'] = df_inspire_squared_raw['gamma0'] * 16 // 8
df_inspire_squared_raw['resizedItemSizeBits'] = df_inspire_squared_raw['gamma0'] * 16


# Read and preprocess rgswpir data
df_rgswpir_raw = read_and_flatten(results_dir=results_dir_rgswpir, name='inspire')
df_rgswpir_raw = df_rgswpir_raw[df_rgswpir_raw['polyLen'] == 2048]
df_rgswpir_raw['payloadBits'] = df_rgswpir_raw['resizedItemSizeBits'] / df_rgswpir_raw['interpolateDegree']
df_rgswpir_raw['payloadKB'] = df_rgswpir_raw['resizedItemSizeBits'] / df_rgswpir_raw['interpolateDegree'] / 8192

# kspir data
df_kspir_raw = read_and_flatten(results_dir=results_dir_kspir)
df_kspir_raw['name'] = 'kspir'
df_kspir_raw['resizedItemSizeBits'] = 32*1024*8

# doublepir data
df_ypir_raw = read_and_flatten(results_dir=results_dir_ypir, name='ypir')
df_ypir_raw['name'] = 'ypir'
df_ypir_raw['resizedItemSizeBits'] = 8

df_simpleypir_raw = read_and_flatten(results_dir=results_dir_simpleypir, name='simpleypir')
df_simpleypir_raw['name'] = 'simpleypir'
df_simpleypir_raw['resizedItemSizeBits'] = 32*1024*8

df_hintlesspir_raw = read_and_flatten(results_dir=results_dir_hintlesspir, name='hintlesspir')
df_hintlesspir_raw['name'] = 'hintlesspir'
df_hintlesspir_raw['resizedItemSizeBits'] = 32*1024*8

# Get all unique dbSizes from both dataframes

all_db_sizes_mb = [1024.0]
# # Uncomment this line to perform full experiments
# all_db_sizes_mb = [1024.0, 8192.0, 32768.0]

all_entry_sizes_bits = [1, 512, 32*1024*8]

for input_item_size_bits in all_entry_sizes_bits:

    this_name = f'final-entry={input_item_size_bits}-bits'
    with open(os.path.join(figures_dir_all, f"{this_name}.txt"), 'w') as f:
        f.write("")

    df_inspire_0_entry_selected = df_inspire_0_raw[df_inspire_0_raw["resizedItemSizeBits"] / input_item_size_bits == df_inspire_0_raw["resizedItemSizeBits"] // input_item_size_bits]
    df_inspire_squared_entry_selected = df_inspire_squared_raw[df_inspire_squared_raw["resizedItemSizeBits"] / input_item_size_bits == df_inspire_squared_raw["resizedItemSizeBits"] // input_item_size_bits]
    df_rgswpir_entry_selected = df_rgswpir_raw[df_rgswpir_raw["payloadBits"] / input_item_size_bits == df_rgswpir_raw["payloadBits"] // input_item_size_bits]
    df_kspir_entry_selected = df_kspir_raw[df_kspir_raw["resizedItemSizeBits"] / input_item_size_bits == df_kspir_raw["resizedItemSizeBits"] // input_item_size_bits]
    df_ypir_entry_selected = df_ypir_raw[df_ypir_raw["resizedItemSizeBits"] / input_item_size_bits == df_ypir_raw["resizedItemSizeBits"] // input_item_size_bits]
    df_simpleypir_entry_selected = df_simpleypir_raw[df_simpleypir_raw["resizedItemSizeBits"] / input_item_size_bits == df_simpleypir_raw["resizedItemSizeBits"] // input_item_size_bits]
    df_hintlesspir_entry_selected = df_hintlesspir_raw[df_hintlesspir_raw["resizedItemSizeBits"] / input_item_size_bits == df_hintlesspir_raw["resizedItemSizeBits"] // input_item_size_bits]

    for dbSize in all_db_sizes_mb:

        df_kspir_selected = df_kspir_entry_selected[df_kspir_entry_selected["inputDatabaseSizeMb"] == dbSize]
        df_ypir_selected = df_ypir_entry_selected[df_ypir_entry_selected["inputDatabaseSizeMb"] == dbSize]
        df_simpleypir_selected = df_simpleypir_entry_selected[df_simpleypir_entry_selected["inputDatabaseSizeMb"] == dbSize]
        df_hintlesspir_selected = df_hintlesspir_entry_selected[df_hintlesspir_entry_selected["inputDatabaseSizeMb"] == dbSize]

        df_inspire_0_selected = df_inspire_0_entry_selected[df_inspire_0_entry_selected["inputDatabaseSizeMb"] == dbSize]

        df_rgswpir_selected = df_rgswpir_entry_selected[df_rgswpir_entry_selected["inputDatabaseSizeMb"] == dbSize].sort_values(
            by=[
                "online-totalBytesKB",
                "online-serverTimeMs",
            ]
        )

        # remove points that are not on the pareto optimal of (online-serverTimeMs, online-totalBytesKB)
        df_rgswpir_selected = pareto_frontier(df_rgswpir_selected, "online-totalBytesKB", "online-serverTimeMs")

        # part1 = transpose_and_round(df_rgswpir_selected[[
        #     'name',
        #     'interpolateDegree',
        #     'offline-totalTimeS',
        #     "online-uploadKeysKB",
        #     "online-uploadQueryKB",
        #     "online-downloadBytesKB",
        #     "online-totalBytesKB",
        #     "online-serverTimeMs",
        #     'online-throughputMBs',
        # ]],
        # int_columns=[
        #     'interpolateDegree',
        # ])
        # print(part1.to_string(header=False))

        # print(
        #     f"\n-------------------------------------------------------- INSPIRE-PARETO --------------------------------------------------------"
        # )

        df_inspire_squared_selected = df_inspire_squared_entry_selected[df_inspire_squared_entry_selected["inputDatabaseSizeMb"] == dbSize].sort_values(
            by=[
                "gamma0",
                "online-totalBytesKB",
                "online-serverTimeMs",
                "gamma1",
                "gamma2",
            ]
        )

        df_inspire_squared_selected = pareto_frontier(df_inspire_squared_selected, "online-totalBytesKB", "online-serverTimeMs")

        our_work_combined = pd.concat([df_inspire_0_selected, df_inspire_squared_selected, df_rgswpir_selected], ignore_index=True, axis=0).sort_values(
            by=[
                "online-totalBytesKB",
                "online-serverTimeMs",
            ]
        )
        our_work_combined_pareto = pareto_frontier(our_work_combined, "online-totalBytesKB", "online-serverTimeMs")
        
        # seperate 'our_work_combined_pareto' based on name
        df_inspire_0_selected = our_work_combined_pareto[our_work_combined_pareto["name"] == "inspire_0"]
        df_inspire_squared_selected = our_work_combined_pareto[our_work_combined_pareto["name"] == "inspire^2"]
        df_rgswpir_selected = our_work_combined_pareto[our_work_combined_pareto["name"] == "inspire"]

        # df_inspire_squared_selected = df_inspire_squared_selected[[
        #     'name',
        #     "gamma0",
        #     "gamma1",
        #     "gamma2",
        #     'payloadB',
        #     "resizedDbFirstDim",
        #     "offline-totalTimeS",
        #     "online-uploadKeysKB",
        #     "online-uploadQueryKB",
        #     "online-downloadBytesKB",
        #     "online-totalBytesKB",
        #     "online-firstPassTimeMs",
        #     "online-secondPassTimeMs",
        #     "online-firstPackTimeMs",
        #     "online-packingKeyRotationsTimeMs",
        #     "online-packingMaskOnlineTimeMs",
        #     "online-packingBodyOnlineTimeMs",
        #     "online-serverTimeMs",
        #     "online-throughputMBs"
        # ]]

        # int_columns = [
        #     'gamma0', 'gamma1', 'gamma2', 'payloadB', 'resizedDbFirstDim',
        #     "offline-totalTimeS",
        #     'online-firstPassTimeMs', 'online-firstPackTimeMs',
        #     'online-packingKeyRotationsTimeMs', 'online-packingMaskOnlineTimeMs',
        #     'online-packingBodyOnlineTimeMs', 'online-serverTimeMs',
        #     'online-throughputMBs'
        # ]

        # part2 = transpose_and_round(df_inspire_squared_selected, int_columns)
        # print(part2.to_string(header=False))
        
        # append the columns of part1 and part2
        combined = pd.concat([df_ypir_selected, df_simpleypir_selected, df_kspir_selected, df_hintlesspir_selected, df_inspire_0_selected, df_inspire_squared_selected, df_rgswpir_selected], ignore_index=True, axis=0)
        combined['name'] = combined['name'].replace({
            'inspire_0' : "$\\ours_0$",
            'inspire^2' : "$\\oursdouble$",
            'inspire' : "$\\ours$",
        })

        ready = transpose_and_round(combined[[
            'name',
            "offline-totalTimeS",
            "online-uploadKeysKB",
            "online-uploadQueryKB",
            "online-downloadBytesKB",
            "online-totalBytesKB",
            "online-serverTimeMs",
            "online-throughputMBs",
        ]],[
            "offline-totalTimeS",
            "online-uploadKeysKB",
            "online-uploadQueryKB",
            "online-downloadBytesKB",
            "online-totalBytesKB",
            "online-serverTimeMs",
            "online-throughputMBs",
        ])

        ready = ready.rename(index={
            # 'name' : '',
            'offline-totalTimeS': 'Offline Time',
            'online-uploadKeysKB': 'Upload (Keys)',
            'online-uploadQueryKB': 'Upload (Query)',
            'online-downloadBytesKB': 'Download',
            'online-totalBytesKB': 'Total Comm.',
            'online-serverTimeMs': 'Server Time',
            'online-throughputMBs': 'Throughput'
        })

        this_name = f'final-entry={input_item_size_bits}-bits'
        # ready.to_csv(os.path.join(figures_dir_all, f"{this_name}.csv"))
        # print(ready.to_string(header=False))
        with open(os.path.join(figures_dir_all, f"{this_name}.txt"), 'a') as f:
            f.write(f"------------------------------ DB Size = {dbSize} MB ------------------------------\n")
            f.write(ready.to_string(header=False))
            f.write("\n")

print(f"Tables written to {figures_dir_all}")