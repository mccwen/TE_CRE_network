# visualizing network
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

def network_polt(filename, min_high_corr):
    db=pd.read_csv(filename)
    # select only highly correlated enhancers for visualization
    sub_db=db[db["coaccess"] > min_high_corr]
    adj=pd.crosstab(sub_db.Peak1, sub_db.Peak2)
    idx=adj.columns.union(adj.index)
    adj=adj.reindex(index=idx, columns=idx, fill_value=0)
    plt.rc('font', size=3) 
    plt.figure(figsize=(16,13), dpi=90)
    nx.draw_networkx(nx.from_pandas_adjacency(m))
    plt.savefig('highly_correlated_enhancer_network_coaccess.png')
    