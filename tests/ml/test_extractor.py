import viola
import numpy as np
import pandas as pd
import os
import pytest
HERE = os.path.abspath(os.path.dirname(__file__))

def test_SV_signature_extractor():
    df_feature = pd.read_csv(os.path.join(HERE, 'data/feature_matrix.csv'), index_col = 0)
    df_feature.drop('others', axis=1, inplace=True)
    sil, fro, ex, sig = viola.SV_signature_extractor(df_feature, 10, 'test', n_components=2, init='nndsvda', beta_loss='frobenius', max_iter=10000)
    sil, fro, ex, sig = viola.SV_signature_extractor(df_feature, 10, 'test', n_components=2, init='nndsvda', max_iter=10000)
    with pytest.raises(TypeError):
        sil, fro, ex, sig = viola.SV_signature_extractor(df_feature, 10, 'test', init='nndsvda', beta_loss='frobenius', max_iter=10000)
