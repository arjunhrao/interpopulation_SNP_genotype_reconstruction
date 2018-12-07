import numpy as np
import json
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

#Importing Sergei's data
with open("17q25_0.01_reconstruction.json") as f:
    data = json.load(f)

#Importing population list
with open("order_of_populations.txt") as f:
    populationorder=f.readlines()

#Strip \n
for i in range(0,len(populationorder)):
    populationorder[i]=populationorder[i].strip("\n")

#Making dictionary
dicdata={}
for i in range(0,26):
    dicdata[data[i][1]]={}

for k in range(0,26):
    for j in range(0,26):
        dicdata[data[k*26+j][0]][data[k*26+j][1]]=data[k*26+j][2]

# making reference matrix
references=populationorder
populations=populationorder[::-1]
harvest=[]

for i in range(0,26):
    harvest.append([])
    for j in range(0,26):
        harvest[i].append(float(dicdata[references[i]][populations[j]]))

popforlabel=[]
for i in populations:
    popforlabel.append(i[0])

harvest = np.array(harvest)

#Generating serborn graph
pallet = ["#ffffff", "#E2B1AC", "#C97E81", "#CD5E3D", "#D35220"]
pallet=pallet[::-1]
ax=sns.heatmap(harvest,linewidths=1.5,linecolor="#eeeeee",cmap=pallet,yticklabels=references,xticklabels=popforlabel,vmin=0,vmax=0.180,square=True,cbar_kws={"orientation": "horizontal","shrink": 0.3})
ax.yaxis.tick_right()
plt.yticks(rotation=0)
plt.xticks(rotation=0)
plt.title("17q25")
plt.show()

