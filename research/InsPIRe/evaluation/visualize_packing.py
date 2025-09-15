from visualize_base import *

version = 'vtemp'
results_dir = f"results/packing"
figures_dir = f"figures/packing"
data_dir = f"figures/packing/data"
os.makedirs(figures_dir, exist_ok=True)
os.makedirs(data_dir, exist_ok=True)

def legend_name(group_name):
    (method, prerotate, gamma) = group_name
    if method.lower() == "cdks":
        return "CDKS"
    else:
        if bool(prerotate) == False:
            if gamma == 2048:
                return f"InspiRINGFull_I"
            elif gamma == 1024:
                return f"InspiRING_I"
            else:
                return f"InspiRING_I(gamma={gamma})"
        else:
            if gamma == 2048:
                return f"InspiRINGFull"
            elif gamma == 1024:
                return f"InspiRING"
            else:
                return f"InspiRING(gamma={gamma})"

# List to store the data from each log file
data = []

# Iterate through all files in the "results" directory
for filename in os.listdir(results_dir):
    if filename.endswith(".json"):

        # Read the content of the file
        with open(os.path.join(results_dir, filename), 'r') as file:
            log_data = json.load(file)
            
            # Flatten the nested JSON structure for easier analysis
            # Append the flattened data to the list
            data.append(log_data)

# Create a DataFrame from the collected data
df = pd.DataFrame(data).sort_values(by=["totalNumToPack", "packingType", "gamma"])
df['sizeAfterPackingKB'] = b_to_KB(
    (df['packingType'] != 'CDKS') * (df['totalNumToPack'] / df['gamma']) * ( 20 * 2048 + 28 * df['gamma']) +
    (df['packingType'] == 'CDKS') * (df['totalNumToPack'] / df['gamma']) * ( (20 + 28) * 2048)
)
# df['sizeAfterPackingKB'] = B_to_KB(df['outputSizeBytes'])
df['timeOnlineMs'] = df['timeOnlineUs'] / 1000
df['timeOnlineStdDevMs'] = df['timeOnlineStdDev'] / 1000
df['timeOfflineS'] = df['timeOfflineUs'] / 1000000
df['packingRate'] = (df['totalNumToPack'] / df['gamma']) * df['logP'] * df['gamma'] / (8 * 1024 * df['sizeAfterPackingKB'])
df['outputSizeKB'] = B_to_KB(df['outputSizeBytes'])
df['legendName'] = df.apply(lambda row: legend_name((row['packingType'], row['preRotate'], row['gamma'])), axis=1)

df['keysSizeBytes'] = (df['gamma'] == 1024) * (df['keysSizeBytes']/2) + (df['gamma'] != 1024) * (df['keysSizeBytes']) 
df['keysSizeKB'] = B_to_KB(df['keysSizeBytes'])

# # First plot: for all group, filter out the rows where "totalNumToPack" is
# # equal to 2048 and plot the totalOnlineTime and totalOfflineTime
# # as a function of "gamma"

# # Group the data into categories
# grouped = df.groupby([
#     'packingType',
#     'preRotate',
# ])

# plt.figure(figsize=(10, 5))
# color=['red', 'green', 'blue']
# for i, (group_name, group_data) in enumerate(grouped):
#     c = color[i]
#     group_name = (group_name[0], bool(group_name[1]))
#     print(group_name)
#     selected = group_data[group_data['totalNumToPack'] == 2048]
#     plt.plot(selected['gamma'], selected['timeOnlineMs'],
#              'o-', c=c,
#              label=f'{group_name}'
#             )
#     plt.plot(selected['gamma'], selected['timeOfflineMs'],
#              's-', c=c,
#             #  label=f'{group_name}-offline'
#              )

# plt.xlabel('gamma')
# # plt.xscale('log', base=2)
# plt.ylabel('time')
# # plt.yscale('log')
# plt.legend()
# plt.title("Packing time of 2048 LWE ciphertexts (Protocol, preRotate)")
# plt.savefig(f"figures/packing/{version}/variable-gamma-2048-lwe-packing-time.png")

# plt.figure(figsize=(10, 5))
# for i, (group_name, group_data) in enumerate(grouped):
#     c = color[i]
#     print(group_name)
#     group_name = (group_name[0], bool(group_name[1]))
#     selected = group_data[group_data['totalNumToPack'] == 2048]
#     plt.plot(selected['gamma'], selected['sizeAfterPackingKB'],
#              '^-', c=c,
#              label=f'{group_name}'
#              )

# plt.xlabel('gamma')
# # plt.xscale('log', base=2)
# plt.ylabel('sizeAfterPackingKB')
# # plt.yscale('log')
# plt.legend()
# plt.title("Packed size of 2048 LWE ciphertexts (Protocol, preRotate)")
# plt.savefig(f"figures/packing/{version}/variable-gamma-2048-lwe-packing-size.png")


# Second plot: keeo those with gamma>=256. group by packingType and preRotate
# and gamma and plot the timeOnlineUs and timeOfflineUs as a function of
# totalNumToPack

table_data = df[df['totalNumToPack'] == 2**12]
cols_to_show_in_table=[
    'legendName',
    'keysSizeKB',
    'outputSizeKB',
    'sizeAfterPackingKB',
    'packingRate',
    'timeOfflineS',
    'timeOnlineMs',
]
print(table_data[table_data['gamma'] >= 1024][cols_to_show_in_table].transpose())

# Group the data into categories
grouped = df[df['gamma'] >= 1024].groupby([
    'packingType',
    'preRotate',
    'gamma',
])

for name, col in [('Online', 'timeOnlineMs'), ('Offline', 'timeOfflineS')]:
    plt.figure(figsize=(10, 5))
    color=['red', 'green', 'blue', 'orange', 'yellow', 'cyan', 'purple', 'brown', 'magenta', 'black', ]
    for i, (group_name, group_data) in enumerate(grouped):
        c = color[i]
        # print(group_name)
        plt.plot(group_data['totalNumToPack'], group_data[col],
                '.-', c=c,
                label=f'{legend_name(group_name)}'
                )

    plt.xlabel('Total Number of LWEs')
    plt.ylabel('Time')
    plt.legend()
    plt.title(f"{name} packing time of LWE ciphertexts")
    plt.savefig(f"{figures_dir}/variable-totalNumToPack-lwe-packing-time-{name.lower()}.pdf")
    plt.savefig(f"{figures_dir}/variable-totalNumToPack-lwe-packing-time-{name.lower()}.png")

for i, (group_name, group_data) in enumerate(grouped):
    file_name = legend_name(group_name)
    group_data.to_csv(f"{data_dir}/{file_name}.csv", index=False)
    




