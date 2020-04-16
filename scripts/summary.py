# To add a new cell, type '#%%'
# To add a new markdown cell, type '#%% [markdown]'
#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataScience.changeDirOnImportExport setting
# ms-python.python added
import os
try:
	os.chdir("/home/daniel/")
	print(os.getcwd())
except:
	pass

#%%
import numpy as np
import pandas as pd
import re


#%%
file_name = "CG_overall_20190829.csv"
data = pd.read_csv("/home/daniel/Dropbox/PapersOR/PM/implementation/results/CG_results_20190829/CG_overall_20190829.csv")
match = re.search(r'.*\_(\d{4})(\d{2})(\d{2})\.csv',file_name)
year = match.group(1)
month = match.group(2)
day = match.group(3)


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
summary_grouped = data.groupby(['pricing_solver','n','m'])


#%%
aggregation = {"tot_lb": {np.max, np.mean},
               "gap": {np.max, np.mean},
               "first_size_graph": {np.max, np.mean},
               "opt": np.sum,
               "reduction": {np.max,np.mean},
              "tot_cputime" : {np.max, np.mean}}

#%%
summary_write = summary_grouped.agg(aggregation).pivot_table(index=['n', 'm'], values=[
                              'tot_lb', 'gap', 'first_size_graph','reduction','opt'], columns=['pricing_solver'])
summary_write.columns.set_levels(['AFFS','AFFC','AFBS','AFBC','AFZFS','AFZFC','AFZBS','AFZBC','TI','ATI'],level=2,inplace=True)
summary_write.columns = ["_".join(x) for x in summary_write.columns.ravel()]
summary_write.to_csv("CG_summary_"+year+month+day+".csv")

#%%
all_instances = data.pivot_table(values=['tot_lb', 'gap', 'first_size_graph','reduction','opt','rel_error','nb_generated_col','global_lower_bound','global_upper_bound','tot_cputime'], index=['n','m','Id'], columns=['pricing_solver'])
all_instances.columns.set_levels(['AFFS','AFFC','AFBS','AFBC','AFZFS','AFZFC','AFZBS','AFZBC','TI','ATI'],level=1,inplace=True)
all_instances.columns = ["_".join(x) for x in all_instances.columns.ravel()]
all_instances.to_csv("CG_allinstances_"+year+month+day+".csv")


#%%



