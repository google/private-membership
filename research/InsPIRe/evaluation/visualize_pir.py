from visualize_base import *

# Directory containing log files
figures_dir = f"figures/pir-experiments"
os.makedirs(figures_dir, exist_ok=True)

with open(os.path.join(figures_dir, "table.txt"), 'w') as table_file:
    for experiment_type in ['small', 'large']:
    # for experiment_type in ['small']:
        results_dir = f"results/pir-experiments/{experiment_type}"

        df = read_and_flatten(results_dir).sort_values(by=['inputDatabaseSizeMb'])
        # ###################################
        ## For filtering out, if necessary
        df = df.set_index(['protocolType', 'secondLevelPackingMask', 'secondLevelPackingBody', 'gamma0',]).sort_index()
        to_keep = [
            ('SimplePIR', 'CDKS', 'NoPacking', 2048),
            ('SimplePIR', 'InspiRING', 'NoPacking', 1024),
            ('SimplePIR', 'InspiRING', 'NoPacking', 2048),
            ('DoublePIR', 'CDKS', 'CDKS', 2048),
            ('DoublePIR', 'InspiRING', 'NoPacking', 1024),
            ('DoublePIR', 'InspiRING', 'NoPacking', 2048),
        ]
        existing_entries = []

        for entry in to_keep:
            if entry in df.index:
                existing_entries.append(entry)
        
        df = df.loc[existing_entries]
        df = df.reset_index()
        # ###################################

        df['name'] = df.apply(lambda row: f"{row['protocolType']}-{row['secondLevelPackingMask']}-{row['secondLevelPackingBody']}-{row['gamma0']}", axis=1)
        columns_to_print = [
            "inputDatabaseSizeMb",
            'name',
            'offline-serverTimeS',
            "online-uploadKeysKB",
            "online-uploadQueryKB",
            "online-downloadBytesKB",
            "online-totalBytesKB",
            "online-serverTimeMs",
            "online-throughputGBs"
        ]
        # print(df.sort_values(by=['inputDatabaseSizeMb'])[columns_to_print])

        # Group the data into categories
        grouped = df.groupby([
            'protocolType',
            'secondLevelPackingMask',
            'secondLevelPackingBody',
            'gamma0'
        ])

        # Plot serverTimeMs as a function of dbSizeMB
        plt.figure()
        for name, group in grouped:
            plt.plot(group['inputDatabaseSizeMb'], group['online-serverTimeMs'], label=f"{name}")
        plt.xlabel('DB Size (MB)')
        plt.xscale('log')
        # plt.yscale('log')
        plt.ylabel('online-serverTimeMs')
        plt.title('serverTimeMs as a function of inputDatabaseSizeMb')
        plt.legend()
        plt.savefig(os.path.join(figures_dir, 'serverTimeMs_vs_inputDatabaseSizeMb.png'))

        # Plot totalCommKB as a function of dbSizeMB
        plt.figure()
        for name, group in grouped:
            plt.plot(group['inputDatabaseSizeMb'], group['online-totalBytesKB'], label=f"{name}")
        plt.xlabel('DB Size (MB)')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel('Total Comm. (KB)')
        # plt.ylim(0, 1000)
        plt.title('Total Comm. vs. DB Size')
        plt.legend()
        plt.savefig(os.path.join(figures_dir, 'totalCommKB_vs_inputDatabaseSizeMb.png'))


        table_file.write(f"==================================================== {experiment_type} ================================================================\n")
        for db_size in df['inputDatabaseSizeMb'].unique():
            to_show = df[df['inputDatabaseSizeMb'] == db_size][columns_to_print].drop(["inputDatabaseSizeMb"], axis=1).reset_index(drop=True).transpose()
            # to_show['print'] = to_show[0].apply(lambda x: str(x)) + " & " + to_show[1].apply(lambda x: str(x)) + " & " + to_show[2].apply(lambda x: str(x)) + " \\\\"
            # print(to_show['print'])
            table_file.write(f"------------------------------------------------ DB Size = {db_size} MB -------------------------------------------------\n")
            table_file.write(to_show.to_string())
            table_file.write("\n")
            
            # for name, group in grouped:
            #     print(name)
            #     print(
            #         group[group['inputDatabaseSizeMb']==db_size]
            #         [columns_to_print]
            #     )
                