from astropy.io import ascii, fits
from astropy.table import Table
import pandas as pd
import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--save', action='store_true')
parser.add_argument('--print', action='store_true')
args = parser.parse_args()

dat = ascii.read('../arxiv/hizqa_catalog.csv')

# name = dat['name']
# # keep only the first 4 numbers after J and +/-
# for i, nm in enumerate(name):
#     if '+' in nm:
#         coord = nm.split('+')
#         name[i] = coord[0][:5] + '+' + coord[1][:4]
#     elif '-' in nm:
#         coord = nm.split('-')
#         name[i] = coord[0][:5] + '-' + coord[1][:4]
# dat['name'] = name
# dat['ra'] = np.round(dat['ra'], 4)
# dat['dec'] = np.round(dat['dec'], 4)
# dat['JAB'] = np.round(dat['JAB'], 2)

# if_unknown = (dat['label'] == 'unknown')
# if_problematic = (dat['label'] == 'problematic')
# if_qso = (dat['label'] == 'qso')
# if_inconclusive = (dat['label'] == 'inconclusive')
# if_felobal = (dat['label'] == 'FeLoBAL')
# if_mdwarf = [l.startswith('M') for l in dat['label']]
# if_ldwarf = [l.startswith('L') for l in dat['label']]
# if_tdwarf = [l.startswith('T') for l in dat['label']]
# if_mlt = np.array(if_mdwarf) | np.array(if_ldwarf) | np.array(if_tdwarf)

# dat['type'] = np.empty_like(dat['label'])
# dat['spectral type'] = np.empty_like(dat['label'])
# dat['redshift'] = np.empty_like(dat['label'])

# dat['type'][if_unknown] = np.full_like(dat[if_unknown]['type'], 'uNQ')
# dat['type'][if_problematic] = np.full_like(dat[if_problematic]['type'], 'problematic')
# dat['type'][if_qso] = np.full_like(dat[if_qso]['type'], 'QSO')
# dat['type'][if_inconclusive] = np.full_like(dat[if_inconclusive]['type'], 'inconclusive')
# dat['type'][if_felobal] = np.full_like(dat[if_felobal]['type'], 'FeLoBAL')
# dat['type'][if_mlt] = np.full_like(dat[if_mlt]['type'], 'star')
# dat['spectral type'][if_mlt] = dat['label'][if_mlt]
# dat['redshift'][if_qso] = [float(n[2:]) for n in dat['note'][if_qso]]
# dat['label'] = dat['type']

# TODO: JAB from wildhunt

# dat = dat['name', 'ra', 'dec', 'JAB', 'label', 'spectral type', 'redshift', 'source', 'note', 'run']
# # select targets with multiple entries
# count = pd.value_counts(dat['name'])
# multiname = count[count > 1].index.values
# # select targets in dat whose name is in multiname
# multitarget = dat[np.in1d(dat['name'], multiname)]

if args.save:
    dat.write('hizqa_catalog.csv')
elif args.print:
    latex = dat['name', 'ra', 'dec', 'JAB', 'label', 'spectral type', 'redshift', 'source', 'run'].to_pandas().to_latex(index=False)
    print(latex)