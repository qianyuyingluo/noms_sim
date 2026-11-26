from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
import numpy as np


def read_grid(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    rows = int(data["row"].max()) + 1
    columns = int(data["col"].max()) + 1
    cells = data.reshape(rows, columns)
    x = np.append(cells["x_left_um"][0], cells["x_right_um"][0, -1])
    y = np.append(cells["y_top_um"][:, 0], cells["y_bottom_um"][-1, 0])
    return x, y, cells


def label_grid(ax, title):
    ax.axhline(0, color="black", linestyle="--", linewidth=0.75)
    ax.invert_yaxis()
    ax.set(xlabel="Position x (µm)", ylabel="Depth y (µm)", title=title)


def save_electrostatic_example(output):
    x, y, cells = read_grid(output / "potential.csv")
    fig, ax = plt.subplots(figsize=(7, 4.3), layout="constrained")
    colors = ax.pcolormesh(x, y, cells["potential_V"], cmap="viridis")
    label_grid(ax, "Electrostatic example (zero volume charge)")
    fig.colorbar(colors, ax=ax, label="Potential (V)")
    fig.savefig(output / "potential.png", dpi=200)
    plt.close(fig)


def save_field_plots(output):
    off = read_grid(output / "fields_off.csv")
    on = read_grid(output / "fields_on.csv")
    fields = [
        ("potential_V", "dc_potential", "Potential", "V", "linear"),
        ("field_V_m", "electric_field", "Electric field magnitude", "V/m", "linear"),
        ("electrons_m3", "electron_density", "Electron density", "m⁻³", "log"),
        ("holes_m3", "hole_density", "Hole density", "m⁻³", "log"),
        ("current_A_m2", "current_density", "Current density magnitude", "A/m²", "log"),
        ("net_doping_m3", "net_doping", "Net doping", "m⁻³", "signed"),
    ]
    for column, name, label, unit, scale in fields:
        values = np.concatenate([off[2][column].ravel(), on[2][column].ravel()])
        values = values[np.isfinite(values)]
        low, high = values.min(), values.max()
        cmap = "viridis"
        if scale == "log":
            positive = values[values > 0]
            if positive.size:
                norm = LogNorm(positive.min(), positive.max())
            else:
                norm = Normalize(0, 1)
        elif scale == "signed":
            limit = max(abs(low), abs(high))
            norm = SymLogNorm(linthresh=limit / 1000, vmin=-limit, vmax=limit)
            cmap = "RdBu_r"
        else:
            norm = Normalize(low, high)
        fig, axes = plt.subplots(1, 2, figsize=(9.1, 3.6), layout="constrained")
        for ax, (x, y, cells), state in zip(axes, (off, on), ("Off", "On")):
            colors = ax.pcolormesh(x, y, cells[column], cmap=cmap, norm=norm)
            vgs = cells["vgs_V"][0, 0]
            vds = cells["vds_V"][0, 0]
            label_grid(ax, f"{state} · VGS = {vgs:g} V, VDS = {vds:g} V")
        fig.suptitle(label)
        fig.colorbar(colors, ax=axes, label=unit)
        fig.savefig(output / f"{name}.png", dpi=180)
        plt.close(fig)


def save_iv_plots(output):
    data = np.atleast_1d(np.genfromtxt(output / "output_iv.csv", delimiter=",", names=True))
    fig, ax = plt.subplots(figsize=(6.2, 4), layout="constrained")
    for gate in np.unique(data["Vgs_V"]):
        curve = data[data["Vgs_V"] == gate]
        ax.plot(curve["Vds_V"], curve["Id_A"] * 1e6, "o-", markersize=3,
                label=f"VGS = {gate:g} V")
    ax.set(xlabel="VDS (V)", ylabel="ID (µA)", title="NMOS output characteristics")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(output / "output_iv.png", dpi=180)
    plt.close(fig)

    transfer = np.atleast_1d(np.genfromtxt(output / "transfer_iv.csv", delimiter=",", names=True))
    for logarithmic in (False, True):
        fig, ax = plt.subplots(figsize=(6.2, 4), layout="constrained")
        if logarithmic:
            positive = transfer["Id_A"] > 0
            ax.semilogy(transfer["Vgs_V"][positive], transfer["Id_A"][positive], "o-", markersize=3)
            filename = "transfer_log.png"
            unit = "A"
        else:
            ax.plot(transfer["Vgs_V"], transfer["Id_A"] * 1e6, "o-", markersize=3)
            filename = "transfer_iv.png"
            unit = "µA"
        ax.set(xlabel="VGS (V)", ylabel=f"ID ({unit})",
               title=f"NMOS transfer · VDS = {transfer['Vds_V'][0]:g} V")
        ax.grid(alpha=0.25)
        fig.savefig(output / filename, dpi=180)
        plt.close(fig)


def save_gm_id_plot(output):
    data = np.atleast_1d(np.genfromtxt(output / "gm_id.csv", delimiter=",", names=True))
    valid = (data["Id_A"] > 0) & np.isfinite(data["gm_over_Id_per_V"])
    data = data[valid]
    if not data.size:
        return
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.7), layout="constrained", sharey=True)
    ratio = data["gm_over_Id_per_V"]
    axes[0].plot(data["Vgs_V"], ratio, "o-", markersize=3)
    axes[0].set(xlabel="VGS (V)", ylabel="gm / ID (1/V)")
    axes[1].semilogx(data["Id_A"] * 1e6, ratio, "o-", markersize=3)
    axes[1].set(xlabel="ID (µA)")
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.suptitle(f"NMOS gm / ID · VDS = {data['Vds_V'][0]:g} V")
    fig.savefig(output / "gm_id.png", dpi=180)
    plt.close(fig)


def main():
    output = Path(__file__).resolve().parent.parent / "output"
    save_electrostatic_example(output)
    save_field_plots(output)
    save_iv_plots(output)
    if (output / "gm_id.csv").exists():
        save_gm_id_plot(output)


if __name__ == "__main__":
    main()
