import json
import os

results_dir = 'results/ypir'
for filename in os.listdir(results_dir):
    if filename.endswith(".txt"):
        with open(os.path.join(results_dir, filename), 'r') as f:
            data = json.load(f)
            data['online']['uploadKeys'] = 11 * 3 * 2048 * 56 // 8
            data['online']['uploadQuery'] = data['online']['uploadBytes'] - data['online']['uploadKeys']                 
            data['online']['totalBytes'] = data['online']['uploadBytes'] + data['online']['downloadBytes']
            with open(f'{results_dir}/{filename}.json', 'w') as f:
                json.dump(data, f, indent=4)


results_dir = 'results/simpleypir'
for filename in os.listdir(results_dir):
    if filename.endswith(".txt"):
        with open(os.path.join(results_dir, filename), 'r') as f:
            data = json.load(f)
            data['online']['uploadKeys'] = 11 * 3 * 2048 * 56 // 8
            data['online']['uploadQuery'] = data['online']['uploadBytes'] - data['online']['uploadKeys']                 
            data['online']['totalBytes'] = data['online']['uploadBytes'] + data['online']['downloadBytes']
            with open(f'{results_dir}/{filename}.json', 'w') as f:
                json.dump(data, f, indent=4)