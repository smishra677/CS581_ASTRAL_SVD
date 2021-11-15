
# Import libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

data = pd.read_csv('stats_final_1_t.csv',names=['Method','dataset','nl','ei1','ei2','fp','fn','rf'])

Astral	=data.query('Method == "Astral"')
svd =data.query('Method == "svd"')



fig, axes = plt.subplots(1,2)
fig.set_size_inches(17,17)



df1 = pd.DataFrame(data=svd)
print(df1.describe())

box1= sns.boxplot(x=df1['dataset'], y=df1['rf'], orient='v', ax=axes[0])
box1.set_xlabel("svd", fontsize=10)
box1.xaxis.set_label_position('top') 
box1.set_ylabel("Average Error from all the methods", fontsize=10)

'''
df = pd.DataFrame(data=Svdquartets, columns=["rf"])
print(df.describe())
box1= sns.boxplot(x="variable", y="value", data=pd.melt(df), orient='v', ax=axes[1])
box1.set_xlabel("Svdquartets", fontsize=10)
box1.xaxis.set_label_position('top') 
box1.set_ylabel("", fontsize=10)

'''


plt.savefig('box_ave_method.jpg', dpi=300)  



plt.show()


