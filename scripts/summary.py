#%%
import os
try:
	os.chdir("/home/danielk/")
	print(os.getcwd())
except:
	pass

#%%
import numpy as np
import pandas as pd
import re


#%%
file_name = "CG_overall_20200506.csv"
data = pd.read_csv("/home/danielk/Dropbox/PapersOR/PM/implementation/" + file_name)
match = re.search(r'.*\_(\d{4})(\d{2})(\d{2})\.csv',file_name)
year = match.group(1)
month = match.group(2)
day = match.group(3)

data['gap'] = (data['global_upper_bound'] - data['global_lower_bound']
               )/(data['global_lower_bound'] + 0.00001)
data['opt'] = data['global_lower_bound'] == data['global_upper_bound']
data['reduction'] = (data['first_size_graph'] - data['size_after_reduced_cost'])/(data['first_size_graph'] + 0.000001)
data['Inst'] = data.NameInstance.apply(
	lambda x:  int(re.search(r'.*\_(\d+)', x).group(1)))


#%%
summary_grouped = data.groupby(['pricing_solver','n','m'])
aggregation = {"tot_lb": {np.max, np.mean},
               "gap": {np.max, np.mean},
               "first_size_graph": {np.max, np.mean},
               "opt": np.sum,
               "reduction": {np.max,np.mean},
              "tot_cputime" : {np.max, np.mean}}
summary_write = summary_grouped.agg(aggregation).pivot_table(index=['n', 'm'], values=[
                              'tot_lb', 'gap', 'first_size_graph','reduction','opt'], columns=['pricing_solver'])
summary_write.columns.set_levels(['AFBC','TI','ATI','TI_BDD'],level=2,inplace=True)
summary_write.columns = ["_".join(x) for x in summary_write.columns.ravel()]
summary_write.to_csv("CG_summary_"+year+month+day+".csv")

#%%
all_instances = data.pivot_table(values=['tot_lb', 'gap', 'first_size_graph','reduction','opt','rel_error','nb_generated_col','global_lower_bound','global_upper_bound','tot_cputime'], index=['n','m','Inst'], columns=['pricing_solver'])
all_instances.columns.set_levels(['AFBC','TI','ATI','TI_BDD'],level=1,inplace=True)
all_instances.columns = ["_".join(x) for x in all_instances.columns.ravel()]
all_instances.to_csv("CG_allinstances_"+year+month+day+".csv")

#%%
df_pessoa = pd.read_csv("/home/danielk/Dropbox/PapersOR/PM/implementation/results/all_pessoa.csv")
df_pessoa.Opt= df_pessoa.Opt.apply(str)
df_pessoa['best'] = df_pessoa.apply(lambda x: re.search(r'[^0-9]?(\d+)', x['Opt']).group(1),axis=1)
df_pessoa.best = df_pessoa.best.apply(pd.to_numeric)

# %%
df_all = pd.merge(data, df_pessoa, on=['Inst', 'n','m'])

# %%
df_all.loc[df_all['global_lower_bound'] > df_all['best'],['NameInstance','m','best','global_lower_bound']]

#%%
df_all.info()

#%%
df_all['NameInstance'].value_counts()
