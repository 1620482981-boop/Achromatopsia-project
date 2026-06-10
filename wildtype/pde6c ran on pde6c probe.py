#!/usr/bin/env python
# coding: utf-8

# In[53]:


import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from itertools import combinations
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.formula.api import ols


Path("wildtype/results").mkdir(parents=True, exist_ok=True)


# In[54]:


class qPCRdf: # define a class that includes qPCR dataframe obtained from qPCR machine
    """
    class qPCRdf: qPCR dataframe exported from qPCR machine
    Finding summary statistics for the ct_values
    Plotting graphes such as delta_ct difference and computing statistical analysis
    """

    def __init__(self,file_name):
        self.file_name = file_name
        self.df = pd.read_csv(self.file_name)
        self.df = self.df.drop(['Unnamed: 0', 'Fluor', 'SQ','Well'], axis=1, errors='ignore')
        #if self.df.loc[self.df['Target'] == 'Control', 'Cq'].notna().any():
            #raise ValueError('Ct value for controls should be NaN,but some are not')
        self.df['TechRep']=self.df.groupby(['Target','Content', 'Sample']).cumcount()
        self.df = self.df[~self.df['Content'].isin(['Uninjected', 'No RNA', 'No enzyme', 'No water'])]
        self.df=self.df.dropna(subset=['Cq']).reset_index(drop=True) # drop samples with any missing values

    def remove_outliers(self):
        def clean_dataframe(group):
            grouped=group['Cq'].to_numpy()
            median=np.median(grouped) #find out the median of technical variable
            range_=np.ptp(grouped) #find out the range of technical variables
            if range_>1:
                distance=np.abs(grouped-median)
                outlier_index=np.argmax(distance)
                clean_group= group.drop(group.index[outlier_index])
            else:
                clean_group=group.copy()

            target, content, sample = group.name
            clean_group['Target'] = target
            clean_group['Content'] = content
            clean_group['Sample'] = sample

            return clean_group

        self.df = self.df.groupby(['Target', 'Content', 'Sample'], group_keys=False).apply(clean_dataframe).reset_index(drop=True)
        self.df = self.df[['Target', 'Content', 'Sample', 'Cq','TechRep']]
        return self



    def Target(self):
        df=self.df
        return df[df['Target']!='SDHA'].reset_index().drop(columns=['index'])

    def Housekeeping(self):
        df=self.df
        return df[df['Target']=='SDHA'].reset_index().drop(columns=['index'])

    def Delta_Ct(self):
        Ct_Housekeeping=self.Housekeeping()
        Ct_Target=self.Target()
        Ct_Housekeeping_mean = Ct_Housekeeping.groupby(["Content", "Sample"], as_index=False)["Cq"].mean().rename(columns={"Cq": "mean_cq_housekeeping"})
        Ct_Target_mean =Ct_Target.groupby(["Content", "Sample"], as_index=False)["Cq"].mean().rename(columns={"Cq": "mean_cq_target"})
        merged = Ct_Target_mean.merge(Ct_Housekeeping_mean,on=["Content", "Sample"],how="inner")

        merged["avg_delta_ct"]=merged['mean_cq_target']-merged['mean_cq_housekeeping']

        return merged 

    def Delta_Delta_Ct_from_plates(self, qpcr_objects, base_line):
        all_dfs = [self.Delta_Ct().copy()]
        for qpcr in qpcr_objects:
            all_dfs.append(qpcr.Delta_Ct().copy())
        combined = pd.concat(all_dfs, ignore_index=True)
        baseline_mean_delta_ct = combined.loc[combined["Content"] == base_line,"avg_delta_ct"].mean()
        combined["delta_delta_ct"] = combined["avg_delta_ct"] - baseline_mean_delta_ct
        combined["fold_change"] = 2 ** (-combined["delta_delta_ct"])
        return combined

    def normality_check(self,qpcr_objects,base_line,para):
        df = self.Delta_Delta_Ct_from_plates(qpcr_objects,base_line)

        model = ols(f"{para} ~ C(Content)", data=df).fit()

        residuals = model.resid
        anova_table = sm.stats.anova_lm(model, typ=1)

        # levene test for constant variance check
        groups=[group["delta_delta_ct"] for name,group in df.groupby('Content')]
        stat, p = stats.levene(*groups)
        print(f'The p value for levene test is {p}'), print(anova_table)
        fig=sm.qqplot(residuals, line="45")
        plt.title("Q-Q plot of residuals")
        plt.show()
        return fig

    def repoted_data(self,qpcr_objects,base_line):
        df = self.Delta_Delta_Ct_from_plates(qpcr_objects,base_line)
        summary = (
        df.groupby("Content").agg(mean_delta_delta_ct=("delta_delta_ct", "mean"),sd_delta_delta_ct=("delta_delta_ct", "std"),n=("Sample", "count"))
        .reset_index())

        summary["reported_fold_change"] = 2 ** (-summary["mean_delta_delta_ct"])
        summary["log2_fold_change"] = -summary["mean_delta_delta_ct"]
        return summary


    def qc_boxplot(self,qpcr_objects,base_line):
        fig, axes = plt.subplots(1, 2, figsize=(13, 8))
        ddct_df = self.Delta_Delta_Ct_from_plates(qpcr_objects,base_line)

        sns.boxplot(x= 'Content',y= 'delta_delta_ct',data=ddct_df,ax=axes[0],showmeans=True)
        axes[0].set_title('Delta_Delta_Ct comparasion')

        sns.barplot(x="Content",y="log2_fold_change",data=self.repoted_data(qpcr_objects,base_line),errorbar=None, ax=axes[1])
        axes[1].axhline(0, linestyle="--")
        axes[1].set_title('Fold_Change comparasion')

        plt.tight_layout()
        plt.show()
        return fig

    def plot_difference(self,qpcr_objects, base_line ,unique_group, value,statistic_test):
        '''
        Run and produce a pair-wise comparasion plot between subgroups using the given statistical tests 
        '''
        ct_mean=self.Delta_Delta_Ct_from_plates(qpcr_objects,base_line)
        fig, ax = plt.subplots(figsize=(6, 5))
        sns.boxplot(x=unique_group, y= value, data=ct_mean, ax=ax)
        sns.stripplot(x=unique_group, y= value, data=ct_mean, color='black', alpha=0.6, ax=ax)
        pairs = list(combinations(ct_mean[unique_group].unique(), 2))
        annotator = Annotator(ax, pairs, data=ct_mean, x=unique_group, y= value)
        annotator.configure(test=statistic_test,        # or 't-test_ind', 'Mann-Whitney'
                            text_format='star',    # shows *, **, ***
                            loc='inside')           # or 'outside'
        annotator.apply_and_annotate()
        plt.show()
        return fig


# In[55]:


my_qPCR=qPCRdf('Pde6c wild type rerun.csv')
my_qPCR.remove_outliers()
my_qPCR.df.head(5)


# In[56]:


my_qPCR.Delta_Ct()


# In[57]:


my_qPCR.Delta_Delta_Ct_from_plates([],'P12')


# In[58]:


my_qPCR.repoted_data([],'P12').to_csv("wildtype/results/pde6c_ran_on_pde6c_probe_summary_table.csv", index=False)


# In[59]:


fig=my_qPCR.normality_check([],'P12','delta_delta_ct')
fig.savefig('wildtype/results/pde6c_ran_on_pde6c_probe_linearity_proof.png', dpi=300, bbox_inches="tight")


# In[60]:


fig=my_qPCR.qc_boxplot([],'P12')
fig.savefig('wildtype/results/pde6c_ran_on_pde6c_probe_fold_change.png', dpi=300, bbox_inches="tight")


# In[61]:


fig=my_qPCR.plot_difference([],'P12','Content','delta_delta_ct', 't-test_ind')
fig.savefig('wildtype/results/pde6c_ran_on_pde6c_probe_statistical_test.png', dpi=300, bbox_inches="tight")

