import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
import numpy as np
import sys

fp = FontProperties(family="Arial", weight="bold") 
globscale = 1.35
LETTERS = { "T" : TextPath((-0.305, 0), "T", size=1, prop=fp),
            "G" : TextPath((-0.384, 0), "G", size=1, prop=fp),
            "A" : TextPath((-0.35, 0), "A", size=1, prop=fp),
            "C" : TextPath((-0.366, 0), "C", size=1, prop=fp) }
COLOR_SCHEME = {'G': 'orange', 
                'A': 'red', 
                'C': 'blue', 
                'T': 'darkgreen'}

def letterAt(letter, x, y, yscale=1, ax=None):
    text = LETTERS[letter]

    t = mpl.transforms.Affine2D().scale(1*globscale, yscale*globscale) + \
        mpl.transforms.Affine2D().translate(x,y) + ax.transData
    p = PathPatch(text, lw=0, fc=COLOR_SCHEME[letter],  transform=t)
    if ax != None:
        ax.add_artist(p)
    return p

# Read in the maf file and save seq
seqs = []
i = 0
with open(sys.argv[1]) as f:
    line = f.readline()
    while line:
        if line[0] == 's':
            line = line.strip().split()[-1].upper()
            if i != 0:
                # Add line into ith seq
                seqs[i] += line
                i += 1
            else:
                seqs.append(line)
            line = f.readline()
        elif line[0] == 'a':
            line = f.readline()
            line = line.strip().split()[-1].upper()
            seqs[0] += line
            i += 1
            line = f.readline()
        else:
            line = f.readline()
# Make all the seqs the same length
# Pad with '-'
max_len = max([len(seq) for seq in seqs])
for i in range(len(seqs)):
    seqs[i] += '-' * (max_len - len(seqs[i]))

# Only get first 16 from each seq
seqs = [seq[:16] for seq in seqs]
max_len = 16

# Get the ratio of each base at each position
# Save the ratio into a list
# Save the bit score into another list
cnt = np.zeros((4, max_len))
for seq in seqs:
    for i in range(len(seq)):
        if seq[i] == 'A':
            cnt[0][i] += 1
        elif seq[i] == 'C':
            cnt[1][i] += 1
        elif seq[i] == 'G':
            cnt[2][i] += 1
        elif seq[i] == 'T':
            cnt[3][i] += 1

# Calculate the ratio
cnt_sum = np.sum(cnt, axis=0)
cnt_ratio = cnt / cnt_sum

# Plot the cnt_sum by the position
# make the color of the line fancy
plt.plot(range(-8, 8), cnt_sum, color='black')
# Make the title bigger
plt.title("Mirlet7d MSA depth", fontsize=20)
plt.xlabel("Distance from cross-linked nucleotide")
plt.ylabel("Depth")
plt.ylim(0, 70)
plt.xticks(range(-8, 8))
# Vertical line at 0
plt.axvline(x=0, color='black', linestyle='--')
plt.savefig('mirlet7d_MSAdepth.png')
plt.close()

# Calculate the bit score
# first, get 2 - shannon entropy
# second, get the bit score
shannon_entropy = 2 + (cnt_ratio * np.log2(cnt_ratio + 1e-30)).sum(axis=0)
preds_bit = cnt_ratio * shannon_entropy

# Plot the sequence logo with preds_tot
all_scores = []
score_each = []
tmp_dict = {0 : 'A', 1 : 'C', 2 : 'G', 3 : 'T'}
for j in range(len(cnt_ratio[0])):
    score_each = []
    for i in range(len(cnt_ratio)):
        score_each.append((tmp_dict[i], cnt_ratio[i][j]))
    all_scores.append(score_each)

fig, ax = plt.subplots(figsize=(10,3))

x = -8
maxi = 0
for scores in all_scores:
    y = 0
    for base, score in scores:
        letterAt(base, x,y, score, ax)
        y += score
    x += 1
    maxi = max(maxi, y)

plt.xticks(range(-8,x))
plt.xlim((-9, x)) 
plt.ylim((0, maxi))
plt.tight_layout()
# X-axis label
plt.xlabel("Distance from cross-linked nucleotide")
plt.ylabel("Frequency")
# Save the figure
plt.savefig('mirlet7d_logo.png', bbox_inches='tight')
plt.close()

# One more plot for bit score
all_scores = []
score_each = []
for j in range(len(preds_bit[0])):
    score_each = []
    for i in range(len(preds_bit)):
        score_each.append((tmp_dict[i], preds_bit[i][j]))
    all_scores.append(score_each)

fig, ax = plt.subplots(figsize=(10,3))

x = -8
maxi = 0
for scores in all_scores:
    y = 0
    for base, score in scores:
        letterAt(base, x,y, score, ax)
        y += score
    x += 1
    maxi = max(maxi, y)

plt.xticks(range(-8,x))
plt.xlim((-9, x))
plt.ylim((0, 2))
#plt.tight_layout()
plt.xlabel("Distance from cross-linked nucleotide")
plt.ylabel("Bits")
# Save the figure
plt.savefig('mirlet7d_logo_bit.png', bbox_inches='tight')