import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns





groupby_results = {}
df = pd.read_csv('stats_final_1_t.csv')
#print(df)

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

#print(groupby_results)

ax = sns.boxplot(x="Method", y="rf", hue = "dataset", data=df)
fig = ax.get_figure()
fig.savefig('preliminary_results')