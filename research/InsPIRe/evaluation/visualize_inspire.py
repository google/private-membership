from visualize_base import *
import math

# Directory containing log files
results_dir = f"results/inspire"
figures_dir = f"figures/inspire"
os.makedirs(figures_dir, exist_ok=True)
# os.makedirs(os.path.join(figures_dir, "poly-degree-ablation"), exist_ok=True)

df = read_and_flatten(results_dir=results_dir)
df = df[df['polyLen'] == 2048]

df['printable1'] = "(" + df['resizedDbFirstDim'].astype(str) + "," + df['online-uploadKeysKB'].astype(int).astype(str) + ")"
df['printable2'] = "(" + df['resizedDbFirstDim'].astype(str) + "," + df['online-uploadQueryKB'].astype(int).astype(str) + ")"
df['printable3'] = "(" + df['resizedDbFirstDim'].astype(str) + "," + df['online-downloadBytesKB'].astype(int).astype(str) + ")"

df['printable4'] = "(" + df['resizedDbFirstDim'].astype(str) + "," + df['online-firstPassTimeMs'].astype(int).astype(str) + ")"
df['printable5'] = "(" + df['resizedDbFirstDim'].astype(str) + "," + df['online-firstPackTimeMs'].astype(int).astype(str) + ")"
df['printable6'] = "(" + df['resizedDbFirstDim'].astype(str) + "," + df['online-rgswTimeMs'].astype(int).astype(str) + ")"

for dbSize in [1024.0, 8192.0, 32768.0]:
    df_raw = df[df["inputDatabaseSizeMb"] == dbSize]
    print(
        f"-------------------------------------------------------- {dbSize} MB --------------------------------------------------------\n"
    )

    for input_item_size_bits in sorted(df_raw["inputItemSizeBits"].unique()):
        df_db = df_raw[df_raw["inputItemSizeBits"] == input_item_size_bits]
        print(f"----------------------------------------------------- {input_item_size_bits} bits ---------------------------------------------\n")

        df_db = df_db.sort_values(
            by=[
                "online-totalBytesKB",
                "online-serverTimeMs",
            ]
        )

        # remove points that are not on the pareto optimal of (online-serverTimeMs, online-totalBytesKB)
        pareto=pareto_frontier(df_db, "online-totalBytesKB", "online-serverTimeMs")

        # Create a scattar plot (with no lines) with "online-serverTimeMs" on the x axis and "online-totalBytesKB" on the y-axis
        plt.plot(
            pareto["online-serverTimeMs"],
            pareto["online-totalBytesKB"],
            '.-',
            label=f"DB={dbSize/1024:.1f} GB"
        )

        print(transpose_and_round(pareto[[
            'resizedDbFirstDim',
            'interpolateDegree',
            # 'payloadKB',
            # 'inputItemSizeBits',
            'offline-totalTimeS',
            'online-clientQueryGenTimeMs',
            "online-uploadKeysKB",
            "online-uploadQueryKB",
            "online-downloadBytesKB",
            "online-totalBytesKB",
            "online-firstPassTimeMs",
            "online-firstPackTimeMs",
            "online-rgswTimeMs",
            "online-serverTimeMs",
            # 'online-throughputMBs',
            # 'printable1',
            # 'printable2',
            # 'printable3',
            # 'printable4',
            # 'printable5',
            # 'printable6',
        ]], [
            'interpolateDegree',
            'offline-totalTimeS',
            'online-clientQueryGenTimeMs',
            "online-uploadKeysKB",
            "online-uploadQueryKB",
            "online-downloadBytesKB",
            "online-totalBytesKB",
            "online-firstPassTimeMs",
            "online-firstPackTimeMs",
            "online-serverTimeMs",
            'online-throughputMBs',
        ]).rename(index={
            'resizedDbFirstDim': '$N/t$',
            'offline-totalTimeS': 'Offline Time',
            'online-clientQueryGenTimeMs': 'Query Gen',
            'online-uploadKeysKB': 'Upload Keys',
            'online-uploadQueryKB': 'Upload Query',
            'online-downloadBytesKB': 'Download',
            'online-totalBytesKB': 'Total Comm.',
            "online-firstPassTimeMs": 'Mat Mul.',
            "online-firstPackTimeMs": 'Packing',
            "online-rgswTimeMs": 'Poly Eval',
            'online-serverTimeMs': 'Total Time',
            'online-throughputMBs': 'Throughput'
        }).to_string(header=False))