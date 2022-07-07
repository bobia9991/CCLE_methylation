from sklearnex import patch_sklearn
patch_sklearn()

import collections
import os
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.model_selection import KFold
from taigapy import TaigaClient
from tqdm import *

import warnings
warnings.filterwarnings("ignore")

# Function to run K-fold cross validation given a feature matrix and the corresponding labels
K_FOLD_NUMBER = 10
def evaluate(features, labels):
    # Stores predicted expression values
    predicted_expression = dict()
    
    # Intialize k-fold cv
    kfold = KFold(K_FOLD_NUMBER, shuffle=True, random_state=0)
    for train, test in kfold.split(features):
        # Split dataset into train-test splits and impute missing values
        train_dataset = features.iloc[train].dropna(axis=1, how='all')
        train_labels = labels.iloc[train]
        
        imputer = SimpleImputer(strategy='mean')
        train_dataset = pd.DataFrame(imputer.fit_transform(train_dataset), columns=train_dataset.columns, index=train_dataset.index)
        
        test_dataset = features.iloc[test][train_dataset.columns.tolist()]
        test_dataset = pd.DataFrame(imputer.transform(test_dataset), columns=test_dataset.columns, index=test_dataset.index)
        
        # Initialize and train model
        model = RandomForestRegressor(n_estimators=500, max_depth=5, n_jobs=-1, random_state=0)
        model.fit(train_dataset, train_labels)

        # Predict on test set and store test set predictions
        prediction = model.predict(test_dataset)
        cell_lines = test_dataset.index
        for i in range(len(prediction)):
            cell_line = cell_lines[i]
            predicted_expression[cell_line] = prediction[i]
            
    # Calculate performance (Pearson's r)
    predicted_series = pd.Series(predicted_expression)
    cv_correlation = predicted_series.corr(labels)

    return cv_correlation, predicted_series


def main():
    # Get all genes in the feature matrix folder
    all_genes = [filename.split('.')[0] for filename in os.listdir('feature_matrices') if filename.endswith('.csv') and "coverage" not in filename]
    
    # Get Expression data
    tc = TaigaClient()
    CCLE_expression = tc.get(name='depmap-a0ab', version=116, file='CCLE_expression')
    CCLE_expression.columns = [element.split(' ')[0] for element in CCLE_expression.columns]

    # Run K-fold CV function on all genes found and record the predictions/performance
    gene_to_correlation = collections.defaultdict(dict)
    threshold_to_use = 0.1
    series_list = list()
    for gene in tqdm(all_genes):
        # Read in feature matrix for this gene
        data_matrix_file = f'feature_matrices/{gene}.csv'
        df = pd.read_csv(data_matrix_file, header=0, index_col=0)
        
        # Get expression data for this gene
        expression_df = CCLE_expression[gene]
        
        # Subset expression and feature matrix to cell lines they have in common
        common_cell_lines = sorted(list(set(df.index) & set(expression_df.index)))
        df = df.loc[common_cell_lines]
        expression_df = expression_df.loc[common_cell_lines]

        # Remove and average features that are very sparse into a single new feature
        min_feature_value = df.shape[1] * threshold_to_use
        non_na_entries_per_feature = df.count()
        to_keep_features = list()
        to_combine_features = list()
        for index_ in non_na_entries_per_feature.index:
            if non_na_entries_per_feature[index_] > min_feature_value:
                to_keep_features.append(index_)
            else:
                to_combine_features.append(index_)
                
        if len(to_keep_features) == 0:
            continue
                
        # Create final feature matrix
        to_combine_features_df = df[to_combine_features]
        new_feature = to_combine_features_df.mean(axis=1)
        new_feature = new_feature.rename('combined_feature')
        df = df[to_keep_features]
        df = pd.concat([df, new_feature], axis=1)

        # Run ML
        if df.shape[1] >= 1 and expression_df.shape[0] >= 10:
            cv_correlation, predicted_series = evaluate(df, expression_df)
            predicted_series = predicted_series.rename(gene)
            series_list.append(predicted_series)
            gene_to_correlation[gene] = cv_correlation

    # Save results to disk
    df = pd.concat(series_list, axis=1)
    df.to_csv(f"ML_predictions.csv")
    
    performance_series = pd.Series(gene_to_correlation)
    performance_series.to_csv(f"ML_performance.csv")


if __name__ == "__main__":
    main()
