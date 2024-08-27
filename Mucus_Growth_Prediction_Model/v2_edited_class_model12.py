import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from imblearn.over_sampling import SMOTE
from sklearn.metrics import precision_score
import seaborn as sns
from sklearn.metrics import (accuracy_score, precision_score, recall_score, f1_score, 
                             matthews_corrcoef, roc_auc_score, confusion_matrix, 
                             precision_recall_curve, roc_curve, auc)
from sklearn.calibration import calibration_curve
from sklearn.model_selection import cross_val_score, learning_curve
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
from sklearn.feature_selection import SelectFromModel, SelectKBest, mutual_info_classif

# Input files and classifier from command line arguments
input_file = sys.argv[1]

# Read the training DataFrame from a text file
df = pd.read_csv(input_file, sep='\t', header=0, index_col=0)

# Separate features (X) and target (y)
X = df.iloc[:, :-1]  # Abundances of microbial species (all columns except the last)
y = df.iloc[:, -1]   # Target variable (BMI_class - last column)

# Feature selection using RandomForestClassifier
# Ensure X is a DataFrame
X = pd.DataFrame(X)

# Feature selection using RandomForestClassifier
rf_clf = RandomForestClassifier(n_estimators=100, random_state=42)
rf_clf.fit(X, y)
rf_selector = SelectFromModel(rf_clf, prefit=True)
X_rf_selected = rf_selector.transform(X)

# Feature selection using SupportVectorMachine
svm_clf = SVC(kernel='linear', C=1.0, random_state=42)
svm_clf.fit(X, y)
svm_selector = SelectFromModel(svm_clf, prefit=True)
X_svm_selected = svm_selector.transform(X)

# Feature selection using NaiveBayes with SelectKBest
nb_clf = GaussianNB()
nb_clf.fit(X, y)
nb_selector = SelectKBest(score_func=mutual_info_classif, k='all')
nb_selector.fit(X, y)
X_nb_selected = nb_selector.transform(X)

# Identify common features selected by all three classifiers
selected_features_rf = rf_selector.get_support()
selected_features_svm = svm_selector.get_support()
selected_features_nb = nb_selector.get_support()

common_features = selected_features_rf & selected_features_svm & selected_features_nb
X_selected = X.loc[:, common_features]

# Print the shape of the selected features
print(f"Original shape of X: {X}")
print(f"Shape of X after feature selection: {X_selected.shape}")
#print(f"X after feature selection: {X_selected}")
# Print the column names of X_selected
print("Selected feature names:")
print(X_selected.columns)

# Function to ensure training set contains at least 3 samples from each class
def custom_train_test_split(X, y, test_size=0.2, random_state=None):
    """
    Custom train-test split function to ensure that both training and testing sets
    contain at least 2 samples from each class.
    
    Parameters:
    - X: Features array or DataFrame
    - y: Labels array or Series
    - test_size: Proportion of the dataset to include in the test split
    - random_state: Seed used by the random number generator
    
    Returns:
    - X_train: Training features
    - X_test: Testing features
    - y_train: Training labels
    - y_test: Testing labels
    """
    np.random.seed(random_state)
    unique_classes = np.unique(y)
    
    X_train_list = []
    X_test_list = []
    y_train_list = []
    y_test_list = []
    
    for cls in unique_classes:
        cls_indices = np.where(y == cls)[0]
        if len(cls_indices) < 4:
            raise ValueError(f"Class {cls} does not have enough samples to ensure at least 2 samples in both training and testing sets.")
        
        np.random.shuffle(cls_indices)
        split_idx = max(2, int(len(cls_indices) * (1 - test_size)))
        
        X_train_list.append(X.iloc[cls_indices[:split_idx]])
        X_test_list.append(X.iloc[cls_indices[split_idx:]])
        y_train_list.append(y.iloc[cls_indices[:split_idx]])
        y_test_list.append(y.iloc[cls_indices[split_idx:]])
    
    X_train = pd.concat(X_train_list)
    X_test = pd.concat(X_test_list)
    y_train = pd.concat(y_train_list)
    y_test = pd.concat(y_test_list)
    
    return X_train, X_test, y_train, y_test

# Example usage:
# X_train, X_test, y_train, y_test = custom_train_test_split(X, y, test_size=0.2, random_state=42)# Generate 10 random numbers between 10 and 100

def calculate_roc_auc(clf, X_test, y_test):
    """
    Calculate the ROC AUC score based on the number of unique classes in y_test.

    Parameters:
    - clf: Trained classifier with a predict_proba method
    - X_test: Test features
    - y_test: Test labels

    Returns:
    - roc_auc: ROC AUC score if applicable, otherwise a message
    """
    unique_classes = np.unique(y_test)
    num_classes = len(unique_classes)
    
    if num_classes == 2:
        # Binary classification
        roc_auc = roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1])
    elif num_classes > 2:
        # Multiclass classification
        roc_auc = roc_auc_score(y_test, clf.predict_proba(X_test), multi_class='ovr', average='weighted')
    else:
        # Only one class present
        return "The number of classes is only one"
    
    return roc_auc

# Example usage:
# roc_auc = calculate_roc_auc(clf, X_test, y_test)

def calculate_tnr(conf_matrix, y_test):
    """
    Calculate the True Negative Rate (Specificity) for binary and multiclass classification.

    Parameters:
    - conf_matrix: Confusion matrix
    - y_test: Test labels

    Returns:
    - tnr: True Negative Rate (Specificity)
    """
    unique_classes = np.unique(y_test)
    num_classes = len(unique_classes)
    
    if num_classes == 2:
        # Binary classification
        tn, fp, fn, tp = conf_matrix.ravel()
        tnr = tn / (tn + fp)
    else:
        # Multiclass classification
        tnr_list = []
        for i in range(conf_matrix.shape[0]):
            tn = np.sum(np.delete(np.delete(conf_matrix, i, axis=0), i, axis=1))
            fp = np.sum(np.delete(conf_matrix[i, :], i))
            tnr_list.append(tn / (tn + fp))
        tnr = np.mean(tnr_list)
    
    return tnr

# Example usage:
# tnr = calculate_tnr(conf_matrix, y_test)
results_list = []

#random_states = np.random.randint(10, 101, size=10)
random_states = np.random.choice(range(5, 101), size=10, replace=False)
# List of classifiers to evaluate
classifiers = {
    'RandomForest': RandomForestClassifier(
        n_estimators=100, 
        max_depth=None, 
        class_weight='balanced', 
        random_state=42
    ),
    'LogisticRegression': LogisticRegression(
        solver='lbfgs', 
        max_iter=1000, 
        random_state=42
    ),
    'SupportVectorMachine': SVC(
        kernel='rbf', 
        C=1.0, 
        gamma='scale', 
        probability=True, 
        class_weight='balanced', 
        random_state=42
    ),
    'KNeighbors': KNeighborsClassifier(
        n_neighbors=5, 
        weights='uniform', 
        algorithm='auto'
    ),
    'NaiveBayes': GaussianNB(
        var_smoothing=1e-9
    )
}

# Initialize dictionaries to store results and cross-validation scores
# Initialize a list to store results
results_list = []
for rand_state in random_states:
    print(f"\nProcessing random_state={rand_state}\n")

    # Create a directory for the current random state
    output_dir = f'random_state_{rand_state}'
    os.makedirs(output_dir, exist_ok=True)

    # Split the dataset into 80% training and 20% testing
    #X_train, X_test, y_train, y_test = custom_train_test_split(X, y, test_size=0.2, random_state=rand_state)
    X_train, X_test, y_train, y_test = custom_train_test_split(X_selected, y, test_size=0.2, random_state=rand_state)
    # Determine the minimum number of samples in any class
    min_samples = min(y_train.value_counts())

    # Set k_neighbors to be less than or equal to min_samples - 1
    k_neighbors = min(5, min_samples - 1)

    # Initialize SMOTE with the adjusted k_neighbors
    smote = SMOTE(k_neighbors=k_neighbors, random_state=rand_state)

    # Apply SMOTE to the training data
    X_train_res, y_train_res = smote.fit_resample(X_train, y_train)
    # Initialize dictionaries to store results and cross-validation scores
    cv_scores_all = {name: [] for name in classifiers}
    accuracy_scores_all = {name: [] for name in classifiers}
    results = {}
    for name, clf in classifiers.items():
        for _ in range(10):  # 100 permutations
            cv_scores = cross_val_score(clf, X_train_res, y_train_res, cv=4)
            mean_cv_scores = np.mean(cv_scores)
            accuracy = np.mean(cv_scores)  # Accuracy is the mean of cross-validation scores
            cv_scores_all[name].append(mean_cv_scores)
            accuracy_scores_all[name].append(accuracy)

        # Train the model on the entire training set
        clf.fit(X_train_res, y_train_res)
        
        # Predict on the test set
        y_pred = clf.predict(X_test)
        
        # Evaluate the model
        accuracy = accuracy_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred, average='weighted', zero_division=0)
        recall = recall_score(y_test, y_pred, average='weighted')
        f1 = f1_score(y_test, y_pred, average='weighted')
        mcc = matthews_corrcoef(y_test, y_pred)
        roc_auc = calculate_roc_auc(clf, X_test, y_test)

        # Confusion Matrix
        conf_matrix = confusion_matrix(y_test, y_pred)
        
        # True Positive Rate (Sensitivity, Recall)
        tpr = recall
        
        # True Negative Rate (Specificity) for multiclass
        tnr = calculate_tnr(conf_matrix, y_test)
        # Store the results
        results[name] = {
            'Cross-Validation Score': np.mean(cv_scores_all[name]),
            'Accuracy': accuracy,
            'Precision': precision,
            'Recall': recall,
            'F1 Score': f1,
            'MCC': mcc,
            'ROC AUC': roc_auc,
            'Confusion Matrix': conf_matrix,
            'TPR': tpr,
            'TNR': tnr
        }
    # Convert the cv_scores_all and accuracy_scores_all dictionaries to DataFrames for plotting
    cv_scores_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in cv_scores_all.items()]))
    accuracy_scores_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in accuracy_scores_all.items()]))

    # Melt the DataFrames for seaborn boxplot
    cv_scores_melted = cv_scores_df.melt(var_name='Classifier', value_name='Mean CV Score')
    accuracy_scores_melted = accuracy_scores_df.melt(var_name='Classifier', value_name='Accuracy')

    # Create the boxplot for mean cross-validation scores with jittered data points
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Classifier', y='Mean CV Score', data=cv_scores_melted)
    sns.stripplot(x='Classifier', y='Mean CV Score', data=cv_scores_melted, color='black', jitter=True, alpha=0.5)
    plt.title('Mean Cross-Validation Scores for Each Classifier')
    plt.xlabel('Classifier')
    plt.ylabel('Mean CV Score')
    plt.xticks(rotation=45)
    plt.tight_layout()

    # Save the plot for mean cross-validation scores
    plt.savefig(os.path.join(output_dir, 'classifier_mean_cv_scores_boxplot.pdf'))
    plt.close()

    # Create the boxplot for accuracy scores with jittered data points
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Classifier', y='Accuracy', data=accuracy_scores_melted)
    sns.stripplot(x='Classifier', y='Accuracy', data=accuracy_scores_melted, color='black', jitter=True, alpha=0.5)
    plt.title('Accuracy Scores for Each Classifier')
    plt.xlabel('Classifier')
    plt.ylabel('Accuracy')
    plt.xticks(rotation=45)
    plt.tight_layout()

    # Save the plot for accuracy scores
    plt.savefig(os.path.join(output_dir, 'classifier_accuracy_scores_boxplot.pdf'))
    plt.close()

    print("Boxplot of mean cross-validation scores has been saved to 'classifier_mean_cv_scores_boxplot.pdf'.")
    print("Boxplot of accuracy scores has been saved to 'classifier_accuracy_scores_boxplot.pdf'.")
    
    # Create a DataFrame from the results list
    # Print and save results for all classifiers
    results_file_path = os.path.join(output_dir, 'results.txt')
    with open(results_file_path, 'w') as f:
        for name, result in results.items():
            f.write(f"Results for {name}:\n")
            f.write("Cross-Validation Score (Training Set): {:.4f}\n".format(result['Cross-Validation Score']))
            f.write("Accuracy (%): {:.2f}\n".format(result['Accuracy'] * 100))
            f.write("Precision: {:.4f}\n".format(result['Precision']))
            f.write("F1 Score: {:.4f}\n".format(result['F1 Score']))
            f.write("MCC: {:.4f}\n".format(result['MCC']))
            f.write("ROC AUC: {:.4f}\n".format(result['ROC AUC']))
            f.write("TPR (Sensitivity): {:.4f}\n".format(result['TPR']))
            f.write("TNR (Specificity): {:.4f}\n".format(result['TNR']))
            f.write("\n" + "="*60 + "\n")

    # Append results to the list
    for name, result in results.items():
        results_list.append({
            'Random State': rand_state,
            'Classifier': name,
            'Cross-Validation Score': result['Cross-Validation Score'],
            'Accuracy (%)': result['Accuracy'] * 100,
            'Precision': result['Precision'],
            'F1 Score': result['F1 Score'],
            'MCC': result['MCC'],
            'ROC AUC': result['ROC AUC'],
            'TPR (Sensitivity)': result['TPR'],
            'TNR (Specificity)': result['TNR']
        })

results_df = pd.DataFrame(results_list)

# Save the DataFrame to a CSV file
results_df.to_csv('evaluation_results.csv', index=False)

print("Evaluation results have been saved to 'evaluation_results.csv'.")
