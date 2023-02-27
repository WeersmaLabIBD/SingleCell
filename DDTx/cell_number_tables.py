### original script by R. Oelen

import tensorflow as tf
import pandas as pd
import csv
import pickle as pkl
import matplotlib.pyplot as plt
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
from sccoda.util import comp_ana as mod
import os

# create class we will use
class ProportionAnalyzer:
  """
    Object to load the count data, set parameters and fetch results
    """
def __init__(self, count_data_loc, treatment_column='donor_recipient', exclusion_columns='ACR_grade', covariates_columns=[], fdr=0.05, sep='\t', reference=None):
  """constructor
        
        Parameters
        ----------
        count_data_loc : str
            The location of the count data to load
        treatment_column : str, optional
            The column that contains the data upon which to compare conditions or treatments
        exclusion_columns : list, optional
            The columns that are removed from the count data (for example covariates that are not used, or cell types to exclude)
        covariates_columns : list, optional
            The columns that contain the covariates to use
        fdr : float, optional
            The FDR cutoff to use
        sep : str, optional
            The separator for reading the count data
        
        """
# set up class variables from constructor
self.__count_data_loc = count_data_loc
self.__treatment_column = treatment_column
self.__exclusion_columns = exclusion_columns
self.__covariates_columns = covariates_columns
self.__fdr = fdr
self.__sep = sep
# other variables that will be populated later
self.__raw_counts = None
self.__filtered_counts = None
self.__sccoda = None
self.__reference = reference
# set up the count data
self.__setup_count_data()
# set up the reference
if self.__reference is None:
  self.__setup_reference()
# of course we will at some point have a result
self.__has_run = False
self.__formula = None
self.__model = None
self.__result = None
# set up a seed so that the results are always the same
tf.random.uniform([1], seed=7777)
tf.random.set_seed(7777)

def __setup_reference(self):
  """set up the reference for the model
        
        """
# get the unique entries possible
levels = self.__filtered_counts[self.__treatment_column].unique()
# the first level is the reference
self.__reference = levels[0]


def __make_columns_safe(self):
  """rename columns so that they don't contain quotes
        
        """
# get the current column names
current_colnames = self.__raw_counts.columns
# create a new list of columns
new_colnames = []
# replace each occurence
for old_column in current_colnames:
  # replace
  new_column = old_column.replace('"', '')
# add to list
new_colnames.append(new_column)

# zip to create a mapping from the old to the new column names
columnname_mapping = dict(zip(current_colnames, new_colnames))
# and do the actual renaming
self.__raw_counts.rename(columnname_mapping, axis = 'columns', inplace = True)


def __setup_count_data(self):
  """setup the count data to be used for the analysis
        
        """
# read the count data
self.__raw_counts = pd.read_csv(self.__count_data_loc, sep = self.__sep, quoting = csv.QUOTE_NONE)
# make the columns safe
self.__make_columns_safe()
# make a copy from which we will remove some columns
self.__filtered_counts = self.__raw_counts.copy(deep = True)
# remove the columns that we don't need
self.__filtered_counts.drop(self.__exclusion_columns, axis = 'columns', inplace = True)
# set up the sccoda object
self.__sccoda = dat.from_pandas(self.__filtered_counts, covariate_columns = self.__covariates_columns + [self.__treatment_column])


def create_boxplot(self):
  """plot the counts
        
        """
# create the plot
viz.boxplots(self.__sccoda, feature_name = self.__treatment_column)
plt.tight_layout()


def run_model(self):
  """run the model
        
        """
# set up the covariates
covariates_formula = ''
for covariate in self.__covariates_columns:
  covariates_formula = ''.join([covariates_formula, covariate, ' + '])
# set up treatment
treatment_formula = "C(" + self.__treatment_column + ", Treatment('" + self.__reference + "'))"
# create formula
self.__formula = covariates_formula + treatment_formula
# create the model
self.__model = mod.CompositionalAnalysis(self.__sccoda, formula = self.__formula)
# run the analysis
self.__result = self.__model.sample_hmc()
# we ran the analysis, so let's say so
self.__has_run = True


def get_result(self):
  """return the result of the model
        
        Returns
        -------
        result
           The result of the compositional analysis
        """
# check if we have run the model
if self.__has_run:
  return self.__result
else:
  print('model has not been run, returning None')
return None


def get_summary(self, fdr=None):
  """get a summary of the results for a specific FDR
        
        Parameters
        ----------
        fdr : float, optional
            The FDR cutoff to use
        
        Returns
        -------
        result
           summary of the compositional analysis
        """
# see if an fdr was supplied
fdr_to_use=None
if fdr is None:
  fdr_to_use = self.__fdr
else:
  fdr_to_use = fdr

# check if we have run the model
if self.__has_run:
  # if so, return a summary
  summary = self.__result.summary_prepare([fdr_to_use])
return(summary)
else:
  print('model has not been run, returning None')
return None


def get_raw_counts(self):
  """get the raw counts
        
        Returns
        -------
        pandas.DataFrame
           The raw counts before filtering columns
        
        """
return self.__raw_counts


def get_counts(self):
  """get the cell counts
        
        Returns
        -------
        pandas.Dataframe
            The filtered cell counts
        
        """
return self.__filtered_counts


def get_sccoda(self):
  """get the sccoda input object
        
        Returns
        -------
        sccoda
            The filtered cell counts
        
        """
return self.__sccoda


def get_fdr(self):
  """get the fdr set
        
        Returns
        -------
        float
            The fdr used
        
        """
return self.__fdr




# the locations of the cell proportion files
proportions_table_loc = '/xx/ddtx/ongoing/differential_proportion/cell_number_tables/donor_recipient/'
# locations of the specific cell type tables
proportions_table_compartment_loc = ''.join([proportions_table_loc, 'cell_numbers_compartment_donor_recipient.tsv'])
proportions_table_adult_elmentaite_martin_immune_loc = ''.join([proportions_table_loc, 'cell_numbers_adult_elmentaite_martin_immune.tsv'])
# the column of the sample
sample_column = 'timepoint'
# the column of the inflammation status
inflammation_column = 'donor_recipient'
# the FDR to use
fdr = 0.05



# perform the analysis for the compartment
compartment_analysis = ProportionAnalyzer(count_data_loc = proportions_table_compartment_loc, exclusion_columns = ['sample_status_timepoint', 'ACR_grade'], covariates_columns = [], reference = '"d"')
# show the plot
compartment_analysis.create_boxplot()
# do the analysis
compartment_analysis.run_model()
# get the result
compartment_summary_005 = compartment_analysis.get_summary()

compartment_summary_005[0]

compartment_summary_005[1]

# perform the analysis for the compartment
celltype_analysis = ProportionAnalyzer(count_data_loc = proportions_table_adult_elmentaite_martin_immune_loc, exclusion_columns = ['sample_status_timepoint', 'ACR_grade'], covariates_columns = [], reference = '"d"')
# show the plot
celltype_analysis.create_boxplot()
# do the analysis
celltype_analysis.run_model()
# get the result
celltype_summary_005 = celltype_analysis.get_summary()

celltype_summary_005[0]


celltype_summary_005[1]

# save the results

result_output_loc = '/xx/ddtx/ongoing/differential_proportion/sccoda/'
# specific files
result_compartment_summary_005_general_summary_loc = ''.join([result_output_loc, 'ddtx_sccoda_compartment_005', '_general_summary.tsv'])
result_compartment_summary_005_covariate_summary_loc = ''.join([result_output_loc, 'ddtx_sccoda_compartment_005', '_covariate_summary.tsv'])
result_celltype_summary_005_general_summary_loc = ''.join([result_output_loc, 'ddtx_sccoda_adult_elmentaite_martin_immune_005', '_general_summary.tsv'])
result_celltype_summary_005_covariate_summary_loc = ''.join([result_output_loc, 'ddtx_sccoda_adult_elmentaite_martin_immune_005', '_covariate_summary.tsv'])
# write the results
compartment_summary_005[0].to_csv(result_compartment_summary_005_general_summary_loc, sep = '\t')
compartment_summary_005[1].to_csv(result_compartment_summary_005_covariate_summary_loc, sep = '\t')
celltype_summary_005[0].to_csv(result_celltype_summary_005_general_summary_loc, sep = '\t')
celltype_summary_005[1].to_csv(result_celltype_summary_005_covariate_summary_loc, sep = '\t')


plt.rcParams['font.size'] = '5'
celltype_analysis.create_boxplot()
plt.xticks(rotation=45, ha='right')
plt.savefig('ddtx_celltypes_donor_recipient.png',dpi=1000)

