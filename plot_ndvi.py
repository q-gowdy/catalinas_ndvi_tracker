"""
plot_ndvi.py
Plot mean NDVI by elevation/aspect zone for the Catalina Mountains monsoon analysis.

Inputs (exported from the Earth Engine script):
  catalinas_ndvi_by_zone.csv          one row per zone, one column per 10-day window
                                      (column label = window midpoint date)
  catalinas_ndvi_pixel_retention.csv  per zone: total pixels vs pixels kept
                                      (valid in every window)

Outputs (written to --out, default ./figures):
  ndvi_all_zones.png        raw NDVI and change from baseline, all six zones
  ndvi_north_vs_south.png   north vs south within each elevation band
  ndvi_summary.csv          baseline, peak, peak date, and gain per zone

Usage:
  python plot_ndvi.py
  python plot_ndvi.py --ndvi data/catalinas_ndvi_by_zone.csv \
                      --retention data/catalinas_ndvi_pixel_retention.csv
"""

import argparse
import re
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd

# ---------- Configuration ----------

# Zone IDs match the Earth Engine script (section 14)
ZONES = {
    1: ('North', 'Low'),
    2: ('North', 'Mid'),
    3: ('North', 'High'),
    4: ('South', 'Low'),
    5: ('South', 'Mid'),
    6: ('South', 'High'),
}

COLORS = {'Low': '#c9822b', 'Mid': '#6a9e3f', 'High': '#1f5f8b'}  # color = elevation
STYLES = {'North': '-', 'South': '--'}                             # line = aspect

N_BASELINE = 3         # baseline = mean of each zone's first N windows
MIN_RETENTION = 0.70   # warn if fewer than this fraction of a zone's pixels survive

DATE_RE = re.compile(r'^\d{4}-\d{2}-\d{2}$')


# ---------- Loading ----------

def load_ndvi(path):
    """Read the wide zone-by-window CSV and return a tidy long-format DataFrame."""
    raw = pd.read_csv(path)
    date_cols = [c for c in raw.columns if DATE_RE.match(str(c))]
    if not date_cols:
        raise ValueError(f'No date columns found in {path}. Columns: {list(raw.columns)}')

    df = raw.melt(id_vars='zone', value_vars=date_cols,
                  var_name='date', value_name='ndvi')
    df['zone'] = df['zone'].astype(int)
    df['date'] = pd.to_datetime(df['date'])
    df = df.dropna(subset=['ndvi']).sort_values(['zone', 'date']).reset_index(drop=True)

    df['aspect'] = df['zone'].map(lambda z: ZONES[z][0])
    df['elev'] = df['zone'].map(lambda z: ZONES[z][1])
    df['label'] = df['aspect'] + ' ' + df['elev']

    # Change from each zone's early-season baseline
    df['baseline'] = df.groupby('zone')['ndvi'].transform(
        lambda s: s.iloc[:N_BASELINE].mean())
    df['ndvi_change'] = df['ndvi'] - df['baseline']
    return df


def check_retention(path):
    """Print how many pixels each zone kept after the same-pixels-in-every-window rule."""
    r = pd.read_csv(path)
    r['zone'] = r['zone'].astype(int)
    r['kept_frac'] = r['kept'] / r['total']
    r['label'] = r['zone'].map(lambda z: f'{ZONES[z][0]} {ZONES[z][1]}')

    print('\nPixel retention (pixels valid in every window):')
    print(r[['zone', 'label', 'total', 'kept', 'kept_frac']]
          .round({'kept_frac': 3}).to_string(index=False))

    low = r[r['kept_frac'] < MIN_RETENTION]
    if not low.empty:
        names = ', '.join(low['label'])
        print(f'\nWARNING: under {MIN_RETENTION:.0%} of pixels retained for: {names}.\n'
              'Those zones are under-sampled. Consider longer windows in the GEE script.')
    return r


# ---------- Plotting ----------

def plot_all_zones(df, out_path):
    fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=True)

    # Shade the baseline period
    base_end = df.groupby('zone')['date'].apply(lambda s: s.iloc[N_BASELINE - 1]).max()
    for ax in axes:
        ax.axvspan(df['date'].min(), base_end, color='gray', alpha=0.12,
                   label='baseline period' if ax is axes[0] else None)

    for (aspect, elev), g in df.groupby(['aspect', 'elev']):
        kw = dict(color=COLORS[elev], linestyle=STYLES[aspect],
                  marker='o', markersize=4, label=f'{aspect} {elev}')
        axes[0].plot(g['date'], g['ndvi'], **kw)
        axes[1].plot(g['date'], g['ndvi_change'], **kw)

    axes[0].set_ylabel('Mean NDVI (10-day median composite)')
    axes[0].set_title('Catalina Mountains NDVI by zone, monsoon season')
    axes[0].legend(ncol=4, fontsize=9)
    axes[1].set_ylabel('Change in NDVI from baseline')
    axes[1].axhline(0, color='gray', linewidth=0.8)
    axes[1].xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))

    for ax in axes:
        ax.grid(alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200)


def plot_north_vs_south(df, out_path):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    for ax, elev in zip(axes, ['Low', 'Mid', 'High']):
        for aspect, g in df[df['elev'] == elev].groupby('aspect'):
            ax.plot(g['date'], g['ndvi'], linestyle=STYLES[aspect], color=COLORS[elev],
                    marker='o', markersize=4, label=aspect)
        ax.set_title(f'{elev} elevation')
        ax.grid(alpha=0.3)
        ax.legend()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
        ax.tick_params(axis='x', rotation=45)
    axes[0].set_ylabel('Mean NDVI')
    plt.tight_layout()
    fig.savefig(out_path, dpi=200)


# ---------- Summary ----------

def build_summary(df):
    rows = []
    for zone, g in df.groupby('zone'):
        base = g['baseline'].iloc[0]
        peak_idx = g['ndvi'].idxmax()
        rows.append({
            'zone': zone,
            'label': g['label'].iloc[0],
            'baseline_ndvi': base,
            'peak_ndvi': g.loc[peak_idx, 'ndvi'],
            'peak_date': g.loc[peak_idx, 'date'].date(),
            'gain': g.loc[peak_idx, 'ndvi'] - base,
            'pct_gain': 100 * (g.loc[peak_idx, 'ndvi'] - base) / base,
            'n_windows': len(g),
        })
    return pd.DataFrame(rows).sort_values('gain', ascending=False)


# ---------- Main ----------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--ndvi', default='catalinas_ndvi_by_zone.csv')
    p.add_argument('--retention', default='catalinas_ndvi_pixel_retention.csv')
    p.add_argument('--out', default='figures')
    p.add_argument('--show', action='store_true', help='open the plots in a window')
    args = p.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if Path(args.retention).exists():
        check_retention(args.retention)
    else:
        print(f'(No retention file at {args.retention}; skipping retention check.)')

    df = load_ndvi(args.ndvi)

    plot_all_zones(df, out / 'ndvi_all_zones.png')
    plot_north_vs_south(df, out / 'ndvi_north_vs_south.png')

    summary = build_summary(df)
    summary.to_csv(out / 'ndvi_summary.csv', index=False)
    print('\nSummary (sorted by gain):')
    print(summary.round(3).to_string(index=False))
    print(f'\nSaved figures and summary to {out.resolve()}')

    if args.show:
        plt.show()


if __name__ == '__main__':
    main()
