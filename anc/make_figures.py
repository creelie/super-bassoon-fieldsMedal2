"""Figures of the paper (3D PNG, 300 dpi).  Usage: python3 -I make_figures.py <anc-dir> <out-dir>"""
import math, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpmath import mp, digamma, zeta, diff, log, pi

anc, out = sys.argv[1], sys.argv[2]
plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm', 'font.size': 10,
                     'axes.edgecolor': '#52514e', 'axes.labelcolor': '#0b0b0b',
                     'xtick.color': '#52514e', 'ytick.color': '#52514e'})
BLUE, ORANGE = '#2a78d6', '#eb6834'
SEQ = LinearSegmentedColormap.from_list('seq', ['#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b'])


def pane(ax):
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.set_facecolor('#fcfcfb'); a.pane.set_edgecolor('#d9d8d3')
        a._axinfo['grid']['color'] = '#e4e3de'


# --- Figure 1: margin of Proposition 4.2 at B = 7/8, log q0 = 22, c = 2
mp.dps = 20
def h0(s0):
    return float(-log(pi) / 2 + digamma((s0 + 1) / 2) / 2 - diff(zeta, s0) / zeta(s0))
def margin(B, L, s0, lq0=22.0, c=2.0):
    k = 1 / (1 - B); logy = k * math.log(c * lq0) - L; y = math.exp(logy); Y = math.exp(L + logy)
    Psi = ((1 - 1 / logy) * (math.exp(L) - 1) - 1.03883 * (1 - math.exp(-L))) * y / L
    R = lq0 + 2 * L * lq0 / math.log(2)
    E = 1.03883 * (math.sqrt(Y) + Y ** (1 / 3) * math.log(Y) / math.log(2))
    T = math.exp(L) / y * math.pi ** 2 / 8 / L
    H = lq0 / 2 + h0(s0); M = max(4 * (s0 - 0.5), 1 / (s0 - B))
    Z = 4 / L * H * (math.exp(B * logy) * math.cosh(B * L / 2) ** 2 * M
                     + math.sqrt(y) * math.cosh(L / 4) ** 2 * s0 / (1 - B) ** 2)
    return (Psi - R - 2 * E - T - Z) / y

Ls = np.linspace(5, 11, 41); ss = np.linspace(1.15, 1.35, 41)
H0 = {s: h0(s) for s in ss}
Zm = np.array([[margin(0.875, L, s) for L in Ls] for s in ss])
LL, SS = np.meshgrid(Ls, ss)
fig = plt.figure(figsize=(6.4, 4.8)); ax = fig.add_subplot(projection='3d'); pane(ax)
Zc = np.arcsinh(Zm)
ax.plot_surface(LL, SS, Zc, cmap=SEQ, linewidth=0.2, edgecolor='#ffffff', antialiased=True, alpha=0.95)
ax.contour(LL, SS, Zm, levels=[0], colors=[ORANGE], linewidths=2, offset=-8)
ax.contour(LL, SS, Zc, levels=[0], colors=[ORANGE], linewidths=2)
m0 = margin(0.875, 8, 1.225)
ax.scatter([8], [1.225], [math.asinh(m0)], s=40, color=ORANGE, edgecolor='#fcfcfb', linewidth=1.5, depthshade=False, zorder=10)
ax.text(7.2, 1.225, math.asinh(m0) + 1.2, r'$L=8,\ \sigma_0=1.225$', color='#0b0b0b', fontsize=9)
ax.set_xlabel(r'$L$'); ax.set_ylabel(r'$\sigma_0$'); ax.set_zlabel(r'$\operatorname{arsinh}(\mathrm{margin}/y)$')
ax.set_zlim(-8, 5); ax.view_init(elev=28, azim=-60)
fig.tight_layout(); fig.savefig(os.path.join(out, 'fig_margin.png'), dpi=300, bbox_inches='tight', pad_inches=0.25); plt.close(fig)

# --- Figure 2: bounds for |d| when E(d) <= e
def bound(k, e):
    lo, hi = k * e, 10000.0
    f = lambda t: (t - math.log(4)) - k * e * math.log(2 * t)
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if f(mid) > 0 else (mid, hi)
    return hi / math.log(10)
es = np.arange(2, 11)
fig = plt.figure(figsize=(6.4, 4.6)); ax = fig.add_subplot(projection='3d'); pane(ax)
for j, (k, col, lab) in enumerate([(8, BLUE, r'$\mathrm{H}(7/8)$'), (12, ORANGE, r'$\mathrm{H}(11/12)$')]):
    hs = [bound(k, e) for e in es]
    ax.bar3d(es - 0.3, np.full(len(es), j) - 0.3, np.zeros(len(es)), 0.6, 0.6, hs, color=col,
             edgecolor='#fcfcfb', linewidth=0.6, shade=True, label=lab)
ax.set_xticks(es); ax.set_yticks([0, 1]); ax.set_yticklabels([r'$7/8$', r'$11/12$'])
ax.set_xlabel(r'exponent $e$'); ax.set_ylabel(r'half-plane $B$'); ax.set_zlabel(r'$\log_{10}$ of bound for $|d|$')
ax.view_init(elev=22, azim=-58)
fig.tight_layout(); fig.savefig(os.path.join(out, 'fig_exponent.png'), dpi=300); plt.close(fig)

# --- Figure 3: the 101 discriminants of exponent at most two
rows = [tuple(map(int, l.split())) for l in open(os.path.join(anc, 'exp2_1e6.txt'))]
def sqfree(n):
    p = 2
    while p * p <= n:
        if n % (p * p) == 0: return False
        p += 1
    return True
def fundamental(D):
    N = -D
    if D % 4 == 1: return sqfree(N)
    m = N // 4
    return D % 4 == 0 and m % 4 in (1, 2) and sqfree(m)
fig = plt.figure(figsize=(6.4, 4.8)); ax = fig.add_subplot(projection='3d'); pane(ax)
for fund, col, lab in [(True, BLUE, 'fundamental (65)'), (False, ORANGE, 'non-fundamental (36)')]:
    pts = [(D, h) for D, h in rows if fundamental(D) == fund]
    x = np.log10([-D for D, h in pts]); yv = np.log2([h for D, h in pts])
    z = np.array([math.pi * h / math.sqrt(-D) for D, h in pts])
    for xi, yi, zi in zip(x, yv, z):
        ax.plot([xi, xi], [yi, yi], [0, zi], color=col, lw=0.7, alpha=0.55)
    ax.scatter(x, yv, z, s=26, color=col, edgecolor='#fcfcfb', linewidth=0.8, depthshade=False, label=lab)
ax.set_xlabel(r'$\log_{10}|D|$'); ax.set_ylabel(r'$\log_2 h(D)$'); ax.set_zlabel(r'$\pi h(D)/\sqrt{|D|}$', labelpad=8)
ax.legend(loc='upper left', frameon=False, fontsize=9); ax.view_init(elev=20, azim=-62)
fig.tight_layout(); fig.savefig(os.path.join(out, 'fig_exptwo.png'), dpi=300, bbox_inches='tight', pad_inches=0.25); plt.close(fig)
print('figures written to', out)
