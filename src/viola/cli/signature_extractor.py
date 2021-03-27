import click
import viola
import numpy as np
import pandas as pd
from io import StringIO, TextIOWrapper
import os
default_tol = 10**(-4)

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version='1.0.0')
@click.option('--n-iter', default=10, help='Number of iteration of NMF.')
@click.option('--name', default='trial', help='The name of this NMF trial.')
@click.option('--n-signatures', default=2, help='Number of SV signature.')
@click.option('--init', default=None, help='Method used to initialize the procedure. Following options including this are corresponding to sklearn.decomposition.NMF.  (See https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.NMF.html)')
@click.option('--solver', default='cd', help='Numerical solver to use: ‘cd’ is a Coordinate Descent solver. ‘mu’ is a Multiplicative Update solver.')
@click.option('--beta-loss', default='frobenius', help='Beta divergence to be minimized, measuring the distance between X and the dot product WH. Note that values different from ‘frobenius’ (or 2) and ‘kullback-leibler’ (or 1) lead to significantly slower fits. Note that for beta_loss <= 0 (or ‘itakura-saito’), the input matrix X cannot contain zeros. Used only in ‘mu’ solver.')
@click.option('--tol', default=default_tol, help='Tolerance of the stopping condition.')
@click.option('--max-iter', default=10000, help='Maximum number of iterations before timing out.')
@click.option('--random-state', default=None, help='Used for initialisation (when init == ‘nndsvdar’ or ‘random’), and in Coordinate Descent. Pass an int for reproducible results across multiple function calls.')
@click.option('--alpha', default=0, help='Constant that multiplies the regularization terms. Set it to zero to have no regularization.')
@click.option('--l1-ratio', default=0, help='The regularization mixing parameter, with 0 <= l1_ratio <= 1. For l1_ratio = 0 the penalty is an elementwise L2 penalty (aka Frobenius Norm). For l1_ratio = 1 it is an elementwise L1 penalty. For 0 < l1_ratio < 1, the penalty is a combination of L1 and L2.')
@click.option('--verbose', default=0, help='Whether to be verbose.')
@click.option('--shuffle', is_flag=True, help='If specified, randomize the order of coordinates in the CD solver.')
@click.option('--regularization', default='both', help='Select whether the regularization affects the components (H), the transformation (W), both or none of them.')
@click.argument('input', type=click.File('r'))
@click.argument('output', type=click.File('w'))

def extract_signature(n_iter, name, n_signatures, init, solver, beta_loss, tol, max_iter,
    random_state, alpha, l1_ratio, verbose, shuffle, regularization, input, output):
    infile = pd.read_csv(input, index_col=0)
    result_sil, result_met, mat_exposure, mat_signature = viola.SV_signature_extractor(
        infile, n_iter=n_iter, name=name, n_components=n_signatures, init=init, solver=solver,
        beta_loss=beta_loss, tol=tol, max_iter=max_iter, random_state=random_state,
        alpha=alpha, l1_ratio=l1_ratio, verbose=verbose, shuffle=shuffle,
        regularization=regularization
    )
    np.savetxt(output, mat_signature, delimiter='\t')
