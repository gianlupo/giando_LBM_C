# -------------------------------------------------------------------
# timeplots.py
# 
# Reads all '.time' files from ./data with format:
# 
#     iteration, time, value
# 
# and generates plots of value vs time, saving PNG files into ./post.
# -------------------------------------------------------------------

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

SRC_DIR  = Path('./src')
DATA_DIR = Path('./data')
POST_DIR = Path('./post')

POST_DIR.mkdir(parents=True, exist_ok=True)

def read_iout():

  # -------------------------------------------------------------------
  # Read IOUT value from src/param.h
  # -------------------------------------------------------------------

  param_file = SRC_DIR / 'param.h'

  with open(param_file, "r") as f:
    for line in f:
      line = line.strip()

      if "IOUT" in line:
        tokens = line.replace("#define", "").split()

        for i, tok in enumerate(tokens):
          if tok == "IOUT":
            return int(tokens[i+1])

  raise ValueError("IOUT not found in src/param.h")

# -------------------------------------------------------------------
# Fixed output image size (same for every PNG)
# -------------------------------------------------------------------
WIDTH_PX =  1000
HEIGHT_PX = 600
DPI = 300

FIGSIZE = (WIDTH_PX / DPI, HEIGHT_PX / DPI)

def load_time_file(filepath):

  # -----------------------------------------------------------------
  # Load a .time file with columns:
  #     iteration, time, value
  # -----------------------------------------------------------------

  try:
    data = np.genfromtxt(filepath, delimiter=",")
    if data.ndim == 1 or np.isnan(data).all():
      raise ValueError
  except Exception:
    # fallback: whitespace-separated
    data = np.genfromtxt(filepath)

  if data.ndim == 1:
    data = data.reshape(1, -1)

  if data.shape[1] < 3:
    raise ValueError(f"{filepath} does not contain 3 columns")

  iteration = data[:, 0]
  time = data[:, 1]
  value = data[:, 2]

  return iteration, time, value

def plot_time_series(filepath,IOUT,factor,col):

  # -----------------------------------------------------------------
  # Generate and save plot for a single .time file.
  # -----------------------------------------------------------------

  _, time, value = load_time_file(filepath)

  value = value/value[0]

  slope = (value[1]-value[0])/(time[1]-time[0])

  nu = 1./6.
  l = 100.
  k = 2.*np.pi/l
  om = nu*k*k

  indices = np.arange(0, len(time), IOUT)

  for idx in indices:

    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)

    ax.plot(time[idx], value[idx], 'o', color=col, linewidth=1.5,label=r'LBM')
    ax.plot(time, np.exp(-factor*om*time), 'k', linewidth=0.5,label=r'analytic')
    ax.plot([0.,-value[0]/slope],[value[0],0.], '--k', alpha=0.25, linewidth=1.5)
    print(slope,-value[0]/slope)


    # Custom x ticks
    xticks = [0, -1./slope, np.max(time)/2., np.max(time)]
    
    # Custom labels
    xticklabels = [
        '0',
        r'$1/\omega$' if factor == 1 else rf'$1/{int(factor)}\omega$',
        f'{np.max(time)/2:.0f}',
        f'{np.max(time):.0f}',
    ]
    
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)

    ax.set_xlabel("Time")
    ax.set_ylabel(filepath.stem+'/'+filepath.stem+'$_0$')

    # Show only gridlines through the origin
    ax.axhline(0, color="gray", linewidth=0.8, alpha=0.5)
    ax.axvline(0, color="gray", linewidth=0.8, alpha=0.5)

    # Legend
    ax.legend(frameon=False)

    output_path = POST_DIR / f'{filepath.stem}_{idx:06d}.png'

    fig.tight_layout()

    # Force exact pixel dimensions
    fig.savefig(
        output_path,
        dpi=DPI,
        bbox_inches=None,
        pad_inches=0,
    )

    plt.close(fig)

    print(f"Saved: {output_path} ({WIDTH_PX}x{HEIGHT_PX}px)")

def main():

  IOUT = read_iout()
  print('IOUT = ',IOUT)

  files = sorted(DATA_DIR.glob('*.time'))
  
  if not files:
    print(f"No .time files found in {DATA_DIR}")
    return
  
  for filepath in files:
    if filepath.stem == 'u_peak':
      factor = 1
      col = 'r'
    if filepath.stem == 'K_peak':
      factor = 2
      col = 'g'
    try:
      plot_time_series(filepath,IOUT,factor,col)
    except Exception as e:
      print(f"Error processing {filepath}: {e}")

if __name__ == "__main__":
  main()
