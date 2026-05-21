from pathlib import Path
import gc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

DATA_DIR = Path("../../data")
POST_DIR = Path("../../post")

POST_DIR.mkdir(exist_ok=True)

def load_field(filename):

  # ---------------------------
  # Load field stored as:
  #     ix iy iz value
  # 
  # and reconstruct a 2D array.
  # ---------------------------

  data = np.loadtxt(filename)

  ix = data[:, 0].astype(int)
  iy = data[:, 1].astype(int)
  iz = data[:, 2].astype(int)
  val = data[:, 3]

  nx = ix.max()
  ny = iy.max()
  nz = iz.max()

  field = np.zeros((nz, nx))

  # iz -> rows, ix -> columns
  field[iz-1, ix-1] = val

  return nx, nz, field

def save_profile_plot(
    sizex,
    sizez,
    field,
    title,
    outfile,
    vmin,
    vmax
):

  z = np.arange(0,sizez) + 0.5
  za = np.insert(np.append(z,sizez),0,0.)
  u0 = 0.01

  WIDTH_PX =  1200
  HEIGHT_PX = 900
  DPI = 300
  
  FIGSIZE = (WIDTH_PX / DPI, HEIGHT_PX / DPI)

  # Fixed figure size
  fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)

  ax.plot(field[:,int(sizex/2)], z, '.', color='r', linewidth=1.5,label=r'LBM')
  ax.plot(4.*u0*(za/sizez)*(1-za/sizez), za, color='k', linewidth=1.5,label=r'analytic')

  # Show only gridlines representing the walls
  ax.axhline(0, color="gray", linewidth=0.8, alpha=0.5)
  ax.axhline(100, color="gray", linewidth=0.8, alpha=0.5)

  # Limits
  ax.set_xlim([0-0.5*u0,u0+0.2*u0])

  ax.set_title(title)
  ax.set_xlabel("u")
  ax.set_ylabel("z")

  # Legend
  ax.legend(frameon=False, loc=0)

  # Fixed layout margins
  fig.tight_layout()

  plt.savefig(outfile, dpi=300)

  plt.close(fig)


def main():

  u_files = sorted(DATA_DIR.glob("u_slice_*.dat"))

  if not u_files:
    print("No u_slice_*.dat files found in data/")
    return

  steps = [f.stem.split("_slice_")[1] for f in u_files]

  # --------------------------------------------------
  # FIRST PASS:
  # compute global min/max
  # --------------------------------------------------

  u_min = np.inf
  u_max = -np.inf

  print("Computing global limits...")

  for step in steps:
      
    _, _, u = load_field(DATA_DIR / f"u_slice_{step}.dat")

    u_min = min(u_min, u.min())
    u_max = max(u_max, u.max())

    del u
    gc.collect()

  print(f"Velocity limits : {u_min:.6e} -> {u_max:.6e}")

  # --------------------------------------------------
  # SECOND PASS:
  # generate plots with fixed limits
  # --------------------------------------------------

  print("Generating plots...")

  for step in steps:

    time = int(step)

    print(f"Processing timestep {step}")

    nx, nz, u = load_field(DATA_DIR / f"u_slice_{step}.dat")

    # u profile plot
    save_profile_plot(
        nx,
        nz,
        u,
        title=f"u (t = {time})",
        outfile=POST_DIR / f"u_profile_{step}.png",
        vmin=u_min,
        vmax=u_max
    )

    del u
    gc.collect()

  print("All plots saved in post/")

if __name__ == "__main__":
  main()
