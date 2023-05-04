import pandas as pd
import numpy as np
cnts = pd.read_csv('read-counts.txt', sep='\t', comment='#', index_col=0)
cnts['clip_enrichment'] = cnts['CLIP-35L33G.bam'] / cnts['RNA-control.bam']
cnts['rden_change'] = (cnts['RPF-siLin28a.bam'] / cnts['RNA-siLin28a.bam']) / (cnts['RPF-siLuc.bam'] / cnts['RNA-siLuc.bam'])

from matplotlib import pyplot as plt

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
nonzero_clip_idx = ~np.isnan(cnts['clip_enrichment'].values) & (cnts['rden_change'].values > 0)
nonzero_rden_idx = ~np.isnan(cnts['rden_change'].values) & (cnts['clip_enrichment'].values > 0)
nonzero_idx = nonzero_clip_idx & nonzero_rden_idx
# only get the values that nonzero_idx[i] is True
ax.scatter(np.log2(cnts['clip_enrichment'].values[nonzero_idx]),np.log2(cnts['rden_change'].values[nonzero_idx]), s=1, alpha=0.5)

ax.set_xlabel('LIN28A CLIP enrichment (log2)')
ax.set_ylabel('Ribosome density change\nupon Lin28a knockdown (log2)')
"""
# Tick labels
ax.set_xticks([-6, -4, -2, 0, 2, 4])
ax.set_yticks([-2, -1, 0, 1, 2])
# limit the range of the x and y axis
ax.set_xlim(-6, 4)
ax.set_ylim(-2, 2)
"""
fig.savefig("test.png", dpi=300)

cnts = pd.read_csv('read-counts-location.txt', sep='\t', comment='#', index_col=0)
cnts['clip_enrichment'] = cnts['CLIP-35L33G.bam'] / cnts['RNA-control.bam']
cnts['rden_change'] = (cnts['RPF-siLin28a.bam'] / cnts['RNA-siLin28a.bam']) / (cnts['RPF-siLuc.bam'] / cnts['RNA-siLuc.bam'])
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
nonzero_clip_idx = ~np.isnan(cnts['clip_enrichment'].values) & (cnts['rden_change'].values > 0)
nonzero_rden_idx = ~np.isnan(cnts['rden_change'].values) & (cnts['clip_enrichment'].values > 0)
nonzero_idx = nonzero_clip_idx & nonzero_rden_idx
# Color the dots by location
location_to_color = {'nucleus': 'blue', 'integral membrane': 'red', 'cytoplasm': 'green'}
cnts['color'] = cnts['location'].apply(lambda x: location_to_color[x])
ax.scatter(np.log2(cnts['clip_enrichment'].values[nonzero_idx]),np.log2(cnts['rden_change'].values[nonzero_idx]), s=1, alpha=0.5, c=cnts['color'].values[nonzero_idx])

ax.set_xlabel('LIN28A CLIP enrichment (log2)')
ax.set_ylabel('Ribosome density change\nupon Lin28a knockdown (log2)')
fig.savefig('clip-vs-rden-filtered_colored.png', dpi=300)
