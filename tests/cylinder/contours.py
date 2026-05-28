from pathlib import Path
import gc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as colors
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

  return field

def save_plot(
    field,
    cmap,
    title,
    outfile,
    vmin,
    vmax,
    norm=None,
    tickformat=None
):

  # Fixed figure size
  fig, ax = plt.subplots(figsize=(8, 6))

  if norm is not None:
    im = ax.imshow(
        field,
        origin="lower",
        cmap=cmap,
        aspect="equal",
        norm=norm
    )
  else:
    im = ax.imshow(
        field,
        origin="lower",
        cmap=cmap,
        aspect="equal",
        vmin=vmin,
        vmax=vmax
    )

  circle = plt.Circle((400, 412), 20., color='k', alpha=0.6)
  ax.add_patch(circle)

  ax.set_title(title)
  ax.set_xlabel("x")
  ax.set_ylabel("z")

  # Horizontal colorbar on top
  cbar = fig.colorbar(
      im,
      ax=ax,
      orientation="horizontal",
      location="top",
      pad=0.08
  )

  #cbar.ax.set_yscale('linear')
  ax.set(yscale='linear')

  cbar.ax.xaxis.set_ticks_position("top")
  cbar.ax.xaxis.set_label_position("top")

  # Optional tick formatting
  if tickformat is not None:
    ticks = np.linspace(vmin, vmax, 11)
    cbar.set_ticks(ticks)
    cbar.ax.xaxis.set_major_formatter(
        FormatStrFormatter(tickformat)
    )

  # Fixed layout margins
  fig.subplots_adjust(
      left=0.10,
      right=0.95,
      bottom=0.10,
      top=0.88
  )

  # Fixed pixel dimensions:
  # width  = 8 * 300 = 2400 px
  # height = 6 * 300 = 1800 px
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
  # compute global min/max for velocity and density
  # --------------------------------------------------

  u_min = np.inf
  u_max = -np.inf

  v_min = np.inf
  v_max = -np.inf

  w_min = np.inf
  w_max = -np.inf

  vel_min = np.inf
  vel_max = -np.inf

  vor_min = np.inf
  vor_max = -np.inf

  K_min = np.inf
  K_max = -np.inf

  rho_min = np.inf#0.
  rho_max = -np.inf

  print("Computing global limits...")

  for step in steps:
      
    u = load_field(DATA_DIR / f"u_slice_{step}.dat")
    v = load_field(DATA_DIR / f"v_slice_{step}.dat")
    w = load_field(DATA_DIR / f"w_slice_{step}.dat")
    K = load_field(DATA_DIR / f"K_slice_{step}.dat")
    rho = load_field(DATA_DIR / f"rho_slice_{step}.dat")

    vel = np.sqrt(u**2 + v**2 + w**2)

    dudz = np.diff(u,axis=0,prepend=0,append=0)
    dwdx = np.diff(w,axis=1,prepend=0,append=0)

    dudz = 0.5*(dudz[0:-1,:]+dudz[1:,:])
    dwdx = 0.5*(dwdx[:,0:-1]+dwdx[:,1:])

    vor = dudz - dwdx

    u_min = min(u_min, u.min())
    u_max = max(u_max, u.max())

    v_min = min(v_min, v.min())
    v_max = max(v_max, v.max())

    w_min = min(w_min, w.min())
    w_max = max(w_max, w.max())

    vel_min = min(vel_min, vel.min())
    vel_max = max(vel_max, vel.max())

    vor_min = min(vor_min, vor.min())
    vor_max = max(vor_max, vor.max())

    K_min = min(K_min, K.min())
    K_max = max(K_max, K.max())

    rho_min = min(rho_min, rho.min())
    rho_max = max(rho_max, rho.max())

    del u, v, w, K, rho, vel, vor
    gc.collect()

  print(f"Velocity limits : {vel_min:.6e} -> {vel_max:.6e}")
  print(f"Density limits  : {rho_min:.6e} -> {rho_max:.6e}")

  # --------------------------------------------------
  # SECOND PASS:
  # generate plots with fixed limits
  # --------------------------------------------------

  print("Generating plots...")

  for step in steps:

    time = int(step)

    print(f"Processing timestep {step}")

    u = load_field(DATA_DIR / f"u_slice_{step}.dat")
    v = load_field(DATA_DIR / f"v_slice_{step}.dat")
    w = load_field(DATA_DIR / f"w_slice_{step}.dat")
    K = load_field(DATA_DIR / f"K_slice_{step}.dat")
    rho = load_field(DATA_DIR / f"rho_slice_{step}.dat")

    vel = np.sqrt(u**2 + v**2 + w**2)

    dudz = np.diff(u,axis=0,prepend=0,append=0)
    dwdx = np.diff(w,axis=1,prepend=0,append=0)

    dudz = 0.5*(dudz[0:-1,:]+dudz[1:,:])
    dwdx = 0.5*(dwdx[:,0:-1]+dwdx[:,1:])

    vor = dudz - dwdx

    if (time == 0):
      angle = 0.
    else:
      angle = -180.*theta[time-1]/np.pi 

    # vor plot

    ncolors = 11
    seismic_colors = matplotlib.colormaps['seismic'](np.linspace(0,1,ncolors))
    a = np.flip(np.linspace(0,1,int((ncolors-1)/2+1)) * vor_min)
    b = np.linspace(0,1,int((ncolors-1)/2+1)) * vor_max
    boundaries = np.concatenate((a,b[1:])) 
    data_range = np.array(boundaries)
    norm_data = (data_range - data_range.min()) / (data_range.max() - data_range.min())
    mycolors = seismic_colors
    cmap_list = list(zip(norm_data, mycolors))
    mycmap = colors.LinearSegmentedColormap.from_list('custom_div', cmap_list)
    normalizer = colors.Normalize(vmin=vor_min, vmax=vor_max)

    save_plot(
        vor,
        cmap=mycmap,
        title=f"vorticity (t = {time})",
        outfile=POST_DIR / f"vor_{step}.png",
        vmin=vor_min,
        vmax=vor_max,
        norm=normalizer
    )

    # u plot

    ncolors = 11
    seismic_colors = matplotlib.colormaps['seismic'](np.linspace(0,1,ncolors))
    a = np.flip(np.linspace(0,1,int((ncolors-1)/2+1)) * u_min)
    b = np.linspace(0,1,int((ncolors-1)/2+1)) * u_max
    boundaries = np.concatenate((a,b[1:])) 
    data_range = np.array(boundaries)
    norm_data = (data_range - data_range.min()) / (data_range.max() - data_range.min())
    mycolors = seismic_colors
    cmap_list = list(zip(norm_data, mycolors))
    mycmap = colors.LinearSegmentedColormap.from_list('custom_div', cmap_list)
    normalizer = colors.Normalize(vmin=u_min, vmax=u_max)

    save_plot(
        u,
        cmap=mycmap,
        title=f"u (t = {time})",
        outfile=POST_DIR / f"u_{step}.png",
        vmin=u_min,
        vmax=u_max,
        norm=normalizer
    )

    # v plot
    save_plot(
        v,
        cmap="seismic",
        title=f"v (t = {time})",
        outfile=POST_DIR / f"v_{step}.png",
        vmin=v_min,
        vmax=v_max
    )

    # w plot
    save_plot(
        w,
        cmap="seismic",
        title=f"w (t = {time})",
        outfile=POST_DIR / f"w_{step}.png",
        vmin=w_min,
        vmax=w_max
    )

    # Velocity magnitude plot
    save_plot(
        vel,
        cmap="Reds",
        title=f"Velocity Magnitude (t = {time})",
        outfile=POST_DIR / f"vel_{step}.png",
        vmin=vel_min,
        vmax=vel_max
    )

    # Kinetic energy plot
    save_plot(
        K,
        cmap="Greens",
        title=f"Kinetic Energy (t = {time})",
        outfile=POST_DIR / f"K_{step}.png",
        vmin=K_min,
        vmax=K_max
    )

    # Density plot
    save_plot(
        rho,
        cmap="Blues",
        title=f"Density (t = {time})",
        outfile=POST_DIR / f"rho_{step}.png",
        vmin=rho_min,
        vmax=rho_max,
        tickformat="%.1f"
    )

    del u, v, w, K, rho, vel
    gc.collect()

  print("All plots saved in post/")

if __name__ == "__main__":
  main()
