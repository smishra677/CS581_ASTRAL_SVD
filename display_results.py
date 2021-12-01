import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

'''
fig, axes = plt.subplots(1,5, figsize=(10,3))
fig.tight_layout()
fig.subplots_adjust(top=0.80)

groupby_results = {}
df = pd.read_csv('svd_output/FinalData.csv')
j = 0
for i in [50, 100, 250, 500, 1000]:
    #df1 = df.loc[df['nl'] == 26]
    df1 = df.loc[df['Method'] == 'usa']
    df1 = df1.loc[df1['dataset'] == i]

    for index, row in df.iterrows():
        method_type = row['Method']
        dataset_type = row['dataset']
        if method_type in groupby_results.keys():
            if dataset_type in groupby_results[method_type].keys():
                groupby_results[method_type][dataset_type].append(row['rf'])
            else:
                groupby_results[method_type][dataset_type] = []
        else:
            groupby_results[method_type] = {}

    print(groupby_results)

    sns.boxplot(x="Threshold", y="rf", data=df1, palette = "Blues", ax = axes[j]).set_title(str(i) +  " sites")
    j = j + 1
fig.suptitle('Comparing USA Across Thresholds (t = 0.8 vs t = 0.4)')
plt.savefig("AcrossThreshold", bbox_inches='tight')
'''

time_data = {'Methods':["usa_0.4", "usa_0.8", "Astral", "svd"], 'Time (s)':[419.17734932899476, 859.55630517005916, 247.43555879592896, 184.36803460121155]}
ax = sns.barplot(x='Methods', y='Time (s)', data=time_data, palette= "Blues").set_title("Runtime Comparison Across Methods")
plt.xlabel("Method")
plt.ylabel("Time (s)")

plt.savefig("Time")
