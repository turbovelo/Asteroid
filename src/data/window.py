import tkinter as tk
from tkinter import ttk
from datetime import timedelta

import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.time import Time
from poliastro.bodies import Sun, Earth, Mars
from poliastro.ephem import Ephem
from poliastro.twobody import Orbit
from poliastro.plotting.static import StaticOrbitPlotter

# --- Create epoch in correct time scale ---
epoch = Time.now().tdb

# --- Load planetary data ---
earth_ephem = Ephem.from_body(Earth, epoch)
mars_ephem = Ephem.from_body(Mars, epoch)

earth = Orbit.from_ephem(Sun, earth_ephem, epoch)
mars = Orbit.from_ephem(Sun, mars_ephem, epoch)

# --- Tkinter setup ---
root = tk.Tk()
root.title("Planet Motion Viewer (Poliastro + Tkinter)")
root.geometry("900x700")

# --- Matplotlib figure setup ---
fig, ax = plt.subplots(figsize=(6, 6))
plotter = StaticOrbitPlotter(ax)
ax.set_title("Planet Orbits")
ax.set_xlabel("x [AU]")
ax.set_ylabel("y [AU]")

canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)


# --- Update function linked to slider ---
def update_plot(val):
    ax.clear()
    plotter = StaticOrbitPlotter(ax)

    # Plot orbits
    plotter.plot(earth, label="Earth", color="#2E8B57")
    plotter.plot(mars, label="Mars", color="#C1440E")

    # Time advance (days)
    t = epoch + timedelta(days=float(val))

    # Propagate orbits
    earth_pos = earth.propagate(t - earth.epoch)
    mars_pos = mars.propagate(t - mars.epoch)

    # Scatter positions
    ax.scatter(*earth_pos.r.to(u.AU).value[:2], color="#2E8B57", s=50)
    ax.scatter(*mars_pos.r.to(u.AU).value[:2], color="#C1440E", s=50)

    # Appearance
    ax.legend()
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True)
    canvas.draw()


# --- Slider (Tkinter Scale) ---
slider_frame = ttk.Frame(root)
slider_frame.pack(fill=tk.X, padx=20, pady=10)

slider_label = ttk.Label(slider_frame, text="Time offset (days):")
slider_label.pack(side=tk.LEFT)

slider = ttk.Scale(
    slider_frame, from_=0, to=800, orient="horizontal", command=update_plot
)
slider.set(0)
slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=10)

# --- Initial plot ---
update_plot(0)

root.mainloop()
