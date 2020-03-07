import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import hdbscan
import sys

# path of file
file = sys.argv[1]
output_path = sys.argv[2]
mcs = sys.argv[3]
cthres = sys.argv[4]

df_file = pd.read_csv(file)
plt.figure(figsize=(15,15))
sns_plot = sns.scatterplot(data=df_file,x='x',y='y')
filename = file.split('.txt')[0].split('data\\')[-1]
sns_plot.get_figure().savefig(r"{0}\{1}.png".format(output_path,filename),encoding='utf8',format='png')
plt.close()

clusterer = hdbscan.HDBSCAN(min_cluster_size=int(mcs))
clusterer.fit(df_file[['x','y']])

df_file['clusters'] = clusterer.labels_
df_file['prob'] = clusterer.probabilities_
df_file['clusters_thres'] = [int(x) if x > float(cthres) else -1 for x in clusterer.probabilities_]

plt.figure(figsize=(15,15))
sns_plot_clustering = sns.scatterplot(data=df_file,x='x',y='y',hue='clusters')
sns_plot_clustering.get_figure().savefig(r"{0}\{1}_clustering.png".format(output_path,filename),encoding='utf8',format='png')
plt.close()


plt.figure(figsize=(15,15))
sns_plot_clusters = sns.scatterplot(data=df_file[df_file['clusters']!=-1],x='x',y='y',hue='clusters')
sns_plot_clusters.get_figure().savefig(r"{0}\{1}_clusters.png".format(output_path,filename),encoding='utf8',format='png')
plt.close()

plt.figure(figsize=(15,15))
sns_plot_rest = sns.scatterplot(data=df_file[df_file['clusters']==-1],x='x',y='y',hue='clusters')
sns_plot_rest.get_figure().savefig(r"{0}\{1}_rest.png".format(output_path,filename),encoding='utf8',format='png')
plt.close()

plt.figure(figsize=(15,15))
sns_plot_thres = sns.scatterplot(data=df_file,x='x',y='y',hue='clusters_thres')
sns_plot_thres.get_figure().savefig(r"{0}\{1}_thres.png".format(output_path,filename),encoding='utf8',format='png')
plt.close()