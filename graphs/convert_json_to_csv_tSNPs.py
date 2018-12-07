import glob
import csv
import re
files=glob.glob('*.json')

#converting json files in the folder into csv file ["#of tSNPs","population"]
def readfile(file):
    with open(file) as f:
        lines=f.readlines()
    values=[len(lines[0].replace("[","").replace("]","").split(","))]
    population = re.search(r'_(.*)_0.1_' , file)
    values.append(population.group(1))
    return values

#making list file
snplist=[]
for i in files:
    snplist.append(readfile(i))

#output values
with open('extracted_values.csv', 'w') as f:
    for k in snplist:
        print(snplist)
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow(k)
        
