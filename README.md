# Codes for RT-qPCR analysis

Python codes instead of Excels are used to analyse and plot the summary plot since they produce better reproducibility
than simple excel-mediated workflows. This README file describes two versions of the code suited for routine qPCR analysis, plotting summary statistics(`delta_delta_ct`, `fold_change`) with respect to the control groups and last but not least,
conducting statistical tests based on summary statistic values. 

## Overview



Import neccessary packages

The code starts by defining a `qPCRdf` class.

```  
 def __init__(self,file_name):
        self.file_name = file_name
        self.df = pd.read_csv(self.file_name)
        self.df = self.df.drop(['Unnamed: 0', 'Fluor', 'SQ','Well'], axis=1, errors='ignore')
        if self.df.loc[self.df['Target'] == 'Control', 'Cq'].notna().any():
            raise ValueError('Ct value for controls should be NaN,but some are not')
        self.df['TechRep']=self.df.groupby(['Target','Content', 'Sample']).cumcount()
        self.df = self.df[~self.df['Content'].isin(['Uninjected', 'No RNA', 'No enzyme', 'No water'])]
        self.df=self.df.dropna(subset=['Cq']).reset_index(drop=True)
 ```
Adding 
