# To add a new cell, type '#%%'
# To add a new markdown cell, type '#%% [markdown]'
#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataScience.changeDirOnImportExport setting
# ms-python.python added
import os
try:
	os.chdir(os.path.join(os.getcwd(), 'build'))
	print(os.getcwd())
except:
	pass

#%%
import pandas as pd


#%%
data = pd.read_csv("../results/28augf1S1/overall.csv")


#%%



#%%
data.head()


#%%
data['gap'] = (data['global_upper_bound'] - data['global_lower_bound']
               )/(data['global_lower_bound'] + 0.00001)
data['opt'] = data['global_lower_bound'] == data['global_upper_bound']
data['reduction'] = (data['first_size_graph'] - data['size_after_reduced_cost'])/(data['first_size_graph'] + 0.000001)
lst = []
import re
for l in data['NameInstance']:
    test = re.search(r'.*\_(\d+)', l).group(1)
    lst.append(int(test))

se = pd.Series(lst)
data['Id'] = se.values


#%%
grouped = data.groupby(['pricing_solver','n','m'])


#%%
import numpy as np
aggregation = {"tot_lb": {np.max, np.mean},
               "gap": {np.max, np.mean},
               "first_size_graph": {np.max, np.mean},
               "opt": np.sum,
               "reduction": {np.max,np.mean},
              "tot_cputime" : {np.max, np.mean}}
            #    "rel_error": {np.max,np.mean},
#%%
grouped["NameInstance"].count()
grouped2 = data.groupby("NameInstance")
grouped2["Id"].count().to_csv("allinstances.csv")

#%%
result = grouped.agg(aggregation)
result


#%%
to_write = result.pivot_table(index=['n', 'm'], values=[
                              'tot_lb', 'gap', 'first_size_graph','reduction','opt'], columns=['pricing_solver'])
to_write

#%%
pivot = data.pivot_table(values=['tot_lb', 'gap', 'first_size_graph','reduction','opt','rel_error','nb_generated_col','global_lower_bound','global_upper_bound','tot_cputime'], index=['n','m','Id'], columns=['pricing_solver'])


#%%
to_write.columns.set_levels(['AFFS','AFFC','AFBS','AFBC','AFZFS','AFZFC','AFZBS','AFZBC','TI','ATI'],level=2,inplace=True)
to_write

#%%
to_write.head()


#%%
to_write.columns = ["_".join(x) for x in to_write.columns.ravel()]
to_write

#%%
to_write.to_csv('results_tue2019aug_summary.csv')


#%%
pivot


#%%
pivot.columns.set_levels(['AF','TI','ATI'],level=1,inplace=True)
pivot


#%%
pivot.columns = ["_".join(x) for x in pivot.columns.ravel()]


#%%
pivot


#%%
pivot.to_csv('results_tue2019_latex_pivot.csv')


#%%



