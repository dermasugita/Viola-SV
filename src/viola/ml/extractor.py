import random
import collections
import numpy as np
import pandas as pd
import scipy
from sklearn.decomposition import NMF, non_negative_factorization
from typing import List
from sklearn.metrics import silhouette_score
from sklearn.metrics import pairwise_distances
from viola.ml._constrained_kmeans import cop_kmeans

def SV_signature_extractor(X, n_iter=10, name='test', **sklearn_nmf_parameters):
    """
    Parameters
    -----------
    X: ndarray or pd.DataFrame
        (n_sample, n_features) shaped ndarray.
    n_iter: int
        Number of iteration for resampling X.
    name: str or int
        Name of this run. Do not use '_' in this argument.
    **sklearn_nmf_parameters: dict
        kwargs to pass the sklearn.decomposition.NMF (https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.NMF.html).
        Specify number of signatures to the 'n_components' argument.
    
    Returns
    --------
    result_silhouette: float
        Average silhouette score of the extracted SV signatures
    result_metrics: float
        beta_loss of the extracted matrice.
    exposure_matrix: (n_samples, n_signatures)
        Exposure matrix.
    signature_matrix: (n_signatures, n_features)
        SV signature matrix extracted.
    """
    if not isinstance(sklearn_nmf_parameters.get('n_components', None), int):
        raise TypeError('n_components argument of sklearn.decomposition.NMF should be int type but {} was passed.'.format(type(sklearn_nmf_parameters.get('n_components'))))

    if isinstance(X, pd.DataFrame):
        X = X.values
    ls_X = []
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            ls_X += [(i, j) for _ in range(X[(i, j)])] 
    
    
    # NMF unit
    def _unit(ls_X, unit_name, **sklearn_nmf_parameters) -> List[pd.Series]:
        # monte carlo bootstrap sampling
        ls_sample = random.choices(ls_X, k=len(ls_X))
        dict_count = collections.Counter(ls_sample)
        new_X = np.zeros_like(X) # need to fix later
        for key, value in dict_count.items():
            new_X[key] = value

        # NMF model
        model = NMF(**sklearn_nmf_parameters)
        W = model.fit_transform(new_X)
        H = model.components_
        #print(model.reconstruction_err_)

        # vector generation
        ls_return = []
        for i in range(H.shape[0]):
            ser_name = str(unit_name) + '_' + str(i)
            ls_return.append(pd.Series(H[i, :], name=ser_name))
        return ls_return
    # /NMF unit
    
    def _get_distance(x: pd.Series, y: pd.Series) -> float:
        x_unit_name = x.name.split('_')[0]
        y_unit_name = y.name.split('_')[0]
        arr_x = x.values
        arr_y = y.values
        if x_unit_name == y_unit_name:
            return 10**16 # for the purpose of constraints
        cos_sim = (x*y).sum()/(np.linalg.norm(arr_x)*np.linalg.norm(arr_y))
        return 1 - cos_sim
    
    # get signature vectors
    ls_ser = []
    for i_iter in range(n_iter):
        unit_name = str(name) + str(i_iter)
        ls_ser += _unit(ls_X, unit_name=unit_name, **sklearn_nmf_parameters)
    
    print(str(name) + ': finished NMF')
    # /get signature vectors
    
    # run constrained k-means clustering
    n_cluster = sklearn_nmf_parameters.get('n_components')
    arr_evaluate = np.array(ls_ser)
    cannot_link = []
    for i in range(len(ls_ser)):
        for j in range(len(ls_ser)):
            if ls_ser[i].name.split('_')[0] == ls_ser[j].name.split('_')[0]:
                cannot_link += [(i, j)] 
    ls_clusters, ls_centers = cop_kmeans(dataset=arr_evaluate, k=n_cluster, cl=cannot_link)
    print(str(name) + ': finished kmeans clustering')
    # /run constrained k-means clustering

    # calculate silhouette score
    result_silhouette = silhouette_score(arr_evaluate, ls_clusters, metric="cosine")
    # /calculate silhouette score
    
    # create average sianature matrix
    mean_H = np.array(ls_centers)
    # /create average sianature matrix
    
    # calculate metrics
    metrics = sklearn_nmf_parameters.get('beta_loss', 'frobenius')
    result_W, result_H, result_iter = non_negative_factorization(X=X, H=mean_H, update_H=False, **sklearn_nmf_parameters)
    result_X = np.dot(result_W, result_H)
    if metrics == "frobenius":
        result_metrics = np.linalg.norm(X - result_X, ord='fro')
    elif metrics == "kullback-leibler":
        result_metrics = scipy.special.kl_div(X, result_X).sum()

    unscaled_W = result_W
    unscaled_H = result_H
    # /calculate metrics

    # create signature matrix and exposure matrix as the final result
    signature_matrix = result_H / result_H.sum(axis=1)[:, np.newaxis]
    exposure_matrix, signature_matrix, result_iter = non_negative_factorization(X=X, H=signature_matrix, update_H=False, **sklearn_nmf_parameters)
    # /create signature matrix and exposure matrix as the final result

    print(str(name) + ': finished all steps')
    print('Silhouette Score: {0}, {1}: {2}'.format(result_silhouette, metrics, result_metrics))
    print('\n==================\n')

    return (result_silhouette, result_metrics, exposure_matrix, signature_matrix)
