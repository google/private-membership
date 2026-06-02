from visualize_base import *

# Directory containing log files
results_dir = f"results/inspire-squared"
figures_dir = f"figures/inspire-squared"
os.makedirs(figures_dir, exist_ok=True)

df = read_and_flatten(results_dir=results_dir)
df = df[df["gamma1"] == 1024]
df = df[df["gamma0"] == df["gamma2"]]

df['payloadB'] = df['gamma0'] * 16 // 8

for dbSize in sorted(df["inputDatabaseSizeMb"].unique()):
    df_db = df[df["inputDatabaseSizeMb"] == dbSize].sort_values(
        by=[
            "gamma0",
            "online-totalBytesKB",
            "online-serverTimeMs",
            "gamma1",
            "gamma2",
        ]
    )

    print(
        f"-------------------------------------------------------- {dbSize} MB --------------------------------------------------------"
    )
    df_db = df_db[[
        "gamma0",
        "gamma1",
        "gamma2",
        'payloadB',
        "resizedDbFirstDim",
        "online-uploadKeysKB",
        "online-uploadQueryKB",
        "online-downloadBytesKB",
        "online-totalBytesKB",
        "online-firstPassTimeMs",
        "online-secondPassTimeMs",
        "online-firstPackTimeMs",
        # "online-packingKeyRotationsTimeMs",
        "online-packingMaskOnlineTimeMs",
        "online-packingBodyOnlineTimeMs",
        "online-serverTimeMs",
    ]]        

    int_columns = [
        'gamma0', 'gamma1', 'gamma2', 'payloadB', 'resizedDbFirstDim',
        "online-uploadKeysKB",
        "online-uploadQueryKB",
        "online-downloadBytesKB",
        "online-totalBytesKB",
        'online-firstPassTimeMs', 'online-firstPackTimeMs',
        'online-packingKeyRotationsTimeMs', 'online-packingMaskOnlineTimeMs',
        'online-packingBodyOnlineTimeMs', 'online-serverTimeMs'
    ]

    # Print the final transposed DataFrame with the corrected data types
    df_final_transposed = transpose_and_round(df_db, int_columns).rename(index={
            # 'name' : '',
            'gamma0': '$\\gamma_0$',
            'gamma1': '$\\gamma_1$',
            'gamma2': '$\\gamma_2$',
            'payloadB': 'Entry Size',
            'resizedDbFirstDim': '$t$',
            'offline-totalTimeS': 'Offline Time',
            'online-uploadKeysKB': 'Upload (Keys)',
            'online-uploadQueryKB': 'Upload (Query)',
            'online-downloadBytesKB': 'Download',
            'online-totalBytesKB': 'Total Comm.',
            "online-firstPassTimeMs": 'First Mat Mul.',
            "online-firstPackTimeMs": 'First Pack',
            "online-secondPassTimeMs": 'Second Mat Mul.',
            "online-packingMaskOnlineTimeMs": 'Second Pack',
            "online-packingBodyOnlineTimeMs": 'Third Pack',
            'online-serverTimeMs': 'Server Time',
            'online-throughputMBs': 'Throughput'
        })


    print(df_final_transposed.to_string(header=False))