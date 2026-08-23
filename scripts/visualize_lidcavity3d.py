"""
Visualizzazione 3D animata dell'output LBM lid-driven cavity.

Il file binario contiene:
  - header: nx, ny, nz (int32)
  - N frame di |u| (float32)
  - ordine: idx = x + nx*(y + ny*z)

Uso:

    python visualize_lidcavity3d.py output/lidcavity3d.bin \
        --html lidcavity3d.html

Il risultato è un HTML 3D interattivo e animato:
  - Play / Pause
  - slider dei frame
  - rotazione con mouse
  - zoom
  - pan
  - colore = |u|

Nota:
Il file contiene solo |u|. Non contiene ux, uy, uz,
quindi non è possibile mostrare vettori o streamlines.
"""

import argparse
import numpy as np


def load(path):
    """Carica il file binario LBM."""

    with open(path, "rb") as f:

        header = np.frombuffer(
            f.read(3 * 4),
            dtype=np.int32
        )

        if header.size != 3:
            raise RuntimeError(
                "Header non valido."
            )

        nx, ny, nz = (
            int(v) for v in header
        )

        if nx <= 0 or ny <= 0 or nz <= 0:
            raise RuntimeError(
                f"Dimensioni non valide: "
                f"{nx}x{ny}x{nz}"
            )

        frame_size = nx * ny * nz

        raw = np.frombuffer(
            f.read(),
            dtype=np.float32
        )

    n_frames = raw.size // frame_size

    if n_frames == 0:
        raise RuntimeError(
            "Nessun frame trovato."
        )

    if raw.size % frame_size != 0:
        print(
            f"[warn] {raw.size % frame_size} "
            f"float in coda scartati."
        )

    frames = raw[
        :n_frames * frame_size
    ].reshape(
        n_frames,
        nz,
        ny,
        nx
    )

    print(
        f"Griglia {nx}x{ny}x{nz}, "
        f"{n_frames} frame letti."
    )

    print(
        f"|u| min = {frames.min():.6e}"
    )

    print(
        f"|u| max = {frames.max():.6e}"
    )

    return frames, (nx, ny, nz)


def save_html_animation(
    frames,
    size,
    html_path,
    max_points=50000
):
    """
    Crea HTML Plotly 3D animato.

    Tutti i frame vengono inseriti nell'HTML.
    """

    try:
        import plotly.graph_objects as go

    except ImportError:
        raise RuntimeError(
            "\nPlotly non è installato.\n\n"
            "Installa con:\n"
            "    pip install plotly\n"
        )

    nx, ny, nz = size

    n_frames = frames.shape[0]

    print(
        f"Preparazione animazione 3D "
        f"con {n_frames} frame..."
    )

    # ------------------------------------------------------------
    # Coordinate della griglia
    # ------------------------------------------------------------

    z, y, x = np.indices(
        (nz, ny, nx)
    )

    x = x.ravel()
    y = y.ravel()
    z = z.ravel()

    n_points = x.size

    # ------------------------------------------------------------
    # Eventuale downsampling
    #
    # Per 32^3 = 32768 punti, NON viene fatto.
    # ------------------------------------------------------------

    if n_points > max_points:

        step = int(
            np.ceil(
                n_points / max_points
            )
        )

        x = x[::step]
        y = y[::step]
        z = z[::step]

        print(
            f"Downsampling coordinate: "
            f"{n_points} -> {x.size} punti"
        )

    # ------------------------------------------------------------
    # Valore globale per mantenere la stessa scala colori
    # durante tutta l'animazione.
    # ------------------------------------------------------------

    global_min = 0.0
    global_max = 0.1
    print(
        f"Scala colori globale: "
        f"{global_min:.6e} -> "
        f"{global_max:.6e}"
    )

    # ------------------------------------------------------------
    # Crea il primo frame
    # ------------------------------------------------------------

    values0 = frames[0].ravel()

    if n_points < values0.size:
        # stesso downsampling usato per le coordinate
        step = int(
            np.ceil(
                values0.size / max_points
            )
        )
        values0 = values0[::step]

    trace = go.Scatter3d(

        x=x,
        y=y,
        z=z,

        mode="markers",

        marker=dict(
            size=3,
            color=values0,
            colorscale="Viridis",
            cmin=global_min,
            cmax=global_max,
            opacity=0.70,

            colorbar=dict(
                title="|u|"
            )
        ),

        customdata=values0,

        hovertemplate=(
            "x=%{x}<br>"
            "y=%{y}<br>"
            "z=%{z}<br>"
            "|u|=%{customdata:.5e}"
            "<extra></extra>"
        )
    )

    # ------------------------------------------------------------
    # Costruzione dei frame Plotly
    # ------------------------------------------------------------

    plotly_frames = []

    for i in range(n_frames):

        values = frames[i].ravel()

        if n_points < values.size:
            step = int(
                np.ceil(
                    values.size / max_points
                )
            )
            values = values[::step]

        plotly_frames.append(
            go.Frame(

                name=str(i),

                data=[
                    go.Scatter3d(

                        x=x,
                        y=y,
                        z=z,

                        mode="markers",

                        marker=dict(
                            size=3,
                            color=values,
                            colorscale="Viridis",
                            cmin=global_min,
                            cmax=global_max,
                            opacity=0.70
                        ),

                        customdata=values,

                        hovertemplate=(
                            "x=%{x}<br>"
                            "y=%{y}<br>"
                            "z=%{z}<br>"
                            "|u|=%{customdata:.5e}"
                            "<extra></extra>"
                        )
                    )
                ]
            )
        )

        if i % 10 == 0 or i == n_frames - 1:
            print(
                f"  frame {i + 1}/{n_frames}"
            )

    # ------------------------------------------------------------
    # Slider
    # ------------------------------------------------------------

    slider_steps = []

    for i in range(n_frames):

        slider_steps.append(
            dict(

                args=[
                    [str(i)],

                    dict(
                        frame=dict(
                            duration=0,
                            redraw=True
                        ),

                        mode="immediate",

                        transition=dict(
                            duration=0
                        )
                    )
                ],

                label=str(i),

                method="animate"
            )
        )

    # ------------------------------------------------------------
    # Play / Pause
    # ------------------------------------------------------------

    play_button = dict(

        label="▶ Play",

        method="animate",

        args=[
            None,

            dict(

                frame=dict(
                    duration=80,
                    redraw=True
                ),

                transition=dict(
                    duration=0
                ),

                fromcurrent=True,

                mode="immediate"
            )
        ]
    )

    pause_button = dict(

        label="⏸ Pause",

        method="animate",

        args=[
            [None],

            dict(

                frame=dict(
                    duration=0,
                    redraw=False
                ),

                mode="immediate"
            )
        ]
    )

    # ------------------------------------------------------------
    # Figura
    # ------------------------------------------------------------

    fig = go.Figure(

        data=[trace],

        frames=plotly_frames
    )

    # ------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------

    fig.update_layout(

        title=(
            "Lid-driven cavity 3D — "
            "frame 0"
        ),

        scene=dict(

            xaxis=dict(
                title="X",
                range=[0, nx - 1]
            ),

            yaxis=dict(
                title="Y",
                range=[0, ny - 1]
            ),

            zaxis=dict(
                title="Z",
                range=[0, nz - 1]
            ),

            aspectmode="cube",

            camera=dict(
                eye=dict(
                    x=1.5,
                    y=1.5,
                    z=1.5
                )
            )
        ),

        updatemenus=[

            dict(

                type="buttons",

                showactive=False,

                x=0.05,
                y=1.12,

                buttons=[
                    play_button,
                    pause_button
                ]
            )
        ],

        sliders=[

            dict(

                active=0,

                x=0.10,
                y=0.02,

                len=0.85,

                currentvalue=dict(
                    prefix="Frame: "
                ),

                steps=slider_steps
            )
        ],

        margin=dict(
            l=0,
            r=0,
            b=0,
            t=80
        )
    )

    # ------------------------------------------------------------
    # Salvataggio
    # ------------------------------------------------------------

    print(
        "\nSalvataggio HTML..."
    )

    fig.write_html(
        html_path,
        include_plotlyjs=True,
        auto_open=False
    )

    print(
        f"\nHTML salvato in:\n"
        f"  {html_path}"
    )

    print(
        "\nPuoi aprirlo con:"
    )

    print(
        f"  explorer.exe {html_path}"
    )


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Visualizzazione 3D animata "
            "LBM lid-driven cavity"
        )
    )

    parser.add_argument(
        "path",
        help="file binario LBM"
    )

    parser.add_argument(
        "--html",
        metavar="PATH",
        required=True,
        help=(
            "file HTML 3D animato "
            "da generare"
        )
    )

    parser.add_argument(
        "--max-points",
        type=int,
        default=50000,
        help=(
            "numero massimo di punti "
            "visualizzati"
        )
    )

    args = parser.parse_args()

    frames, size = load(
        args.path
    )

    save_html_animation(
        frames,
        size,
        args.html,
        max_points=args.max_points
    )


if __name__ == "__main__":
    main()