import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy import units as u
from poliastro.bodies import Earth
from poliastro.twobody import Orbit

# --- Configuración de la Página ---
st.set_page_config(
    page_title="Corrección Orbital (Hohmann)",
    page_icon="🛰️",
    layout="wide"
)

# Configurar matplotlib para usar backend Agg (mejor para entornos virtuales)
plt.switch_backend('Agg')

st.title("🛰️ Calculadora de Corrección Orbital")
st.write("""
Esta aplicación calcula la maniobra (Transferencia de Hohmann) necesaria para 
pasar de una órbita inicial a una órbita objetivo, cumpliendo con tus requisitos.
""")

# --- Constante Radio Terrestre ---
R_EARTH_KM = Earth.R.to(u.km).value

# --- Columna de Entradas (Inputs) ---
st.sidebar.header("Parámetros de Entrada")

st.sidebar.subheader("1. Órbita Inicial (Actual)")
st.sidebar.write("Define la órbita actual del satélite.")
alt_p_inicial = st.sidebar.slider(
    "Altitud Perigeo Inicial (km)", 150.0, 1000.0, 400.0, 50.0
)
alt_a_inicial = st.sidebar.slider(
    "Altitud Apogeo Inicial (km)", 150.0, 2000.0, 600.0, 50.0
)

st.sidebar.subheader("2. Órbita Óptima (Objetivo)")
st.sidebar.write("Define la órbita circular a la que quieres llegar.")
alt_target = st.sidebar.slider(
    "Altitud Circular Objetivo (km)", 500.0, 3000.0, 800.0, 50.0
)

st.sidebar.subheader("3. Info del Satélite")
sat_mass = st.sidebar.number_input("Masa del Satélite (kg)", 10.0, 5000.0, 100.0)

# --- Validaciones y Creación de Órbitas ---
if alt_p_inicial > alt_a_inicial:
    st.error("Error: La altitud del perigeo inicial no puede ser mayor que la del apogeo.")
    st.stop()

try:
    # Creando la órbita inicial
    r_p_inicial = (R_EARTH_KM + alt_p_inicial) * u.km
    r_a_inicial = (R_EARTH_KM + alt_a_inicial) * u.km
    
    # Calcular semieje mayor y excentricidad
    a_inicial = (r_p_inicial + r_a_inicial) / 2
    ecc_inicial = (r_a_inicial - r_p_inicial) / (r_a_inicial + r_p_inicial)
    
    orb_inicial = Orbit.from_classical(
        Earth,
        a_inicial,
        ecc_inicial,
        0 * u.deg,
        0 * u.deg,
        0 * u.deg,
        0 * u.deg
    )

    # Creando la órbita objetivo (circular)
    r_target = (R_EARTH_KM + alt_target) * u.km
    orb_target = Orbit.from_classical(
        Earth,
        r_target,
        0 * u.one,
        0 * u.deg,
        0 * u.deg,
        0 * u.deg,
        0 * u.deg
    )

except Exception as e:
    st.error(f"Error al crear las órbitas: {e}")
    st.stop()

# --- Cálculos de Transferencia de Hohmann ---
st.header("Resultados de la Corrección")

try:
    mu = Earth.k.to(u.km**3 / u.s**2).value
    
    r1 = orb_inicial.r_p.to(u.km).value
    r2 = orb_target.a.to(u.km).value
    
    a_transfer = (r1 + r2) / 2
    
    v1_inicial = np.sqrt(mu * (2/r1 - 1/orb_inicial.a.to(u.km).value))
    v1_transfer = np.sqrt(mu * (2/r1 - 1/a_transfer))
    dv1_mag = abs(v1_transfer - v1_inicial)
    
    v2_transfer = np.sqrt(mu * (2/r2 - 1/a_transfer))
    v2_circular = np.sqrt(mu / r2)
    dv2_mag = abs(v2_circular - v2_transfer)
    
    dv_total = dv1_mag + dv2_mag
    
    ecc_transfer = (r2 - r1) / (r2 + r1)
    orb_transfer = Orbit.from_classical(
        Earth,
        a_transfer * u.km,
        ecc_transfer * u.one,
        0 * u.deg,
        0 * u.deg,
        0 * u.deg,
        0 * u.deg
    )
    
    T_transfer = np.pi * np.sqrt(a_transfer**3 / mu)

except Exception as e:
    st.error(f"Error al calcular la maniobra de Hohmann: {e}")
    st.write(f"Detalles del error: {type(e).__name__}: {str(e)}")
    st.stop()

# --- Sección de Salidas (Outputs) ---

st.subheader(r"🚀 $\Delta V$ (Delta-V) Requerido") 
st.write("Se asume una Transferencia de Hohmann. Ambos impulsos son **Progrados** (en dirección de la velocidad).")

col1, col2, col3 = st.columns(3)
col1.metric("Impulso 1 (en Perigeo inicial)", f"{dv1_mag*1000:.2f} m/s")
col2.metric("Impulso 2 (en Apogeo de transf.)", f"{dv2_mag*1000:.2f} m/s")
col3.metric("🔥 Delta-V Total", f"{dv_total*1000:.2f} m/s")

st.info(f"⏱️ Tiempo de transferencia: {T_transfer/60:.1f} minutos ({T_transfer/3600:.2f} horas)")

st.subheader("📉 Error Orbital (vs. Objetivo)")
col1, col2 = st.columns(2)
with col1:
    st.write("**Antes de la Corrección (Inicial):**")
    err_p_antes = (orb_inicial.r_p.to(u.km).value - orb_target.r_p.to(u.km).value)
    err_a_antes = (orb_inicial.r_a.to(u.km).value - orb_target.r_a.to(u.km).value)
    st.metric("Error en Perigeo", f"{err_p_antes:.1f} km")
    st.metric("Error en Apogeo", f"{err_a_antes:.1f} km")
with col2:
    st.write("**Después de la Corrección (Final):**")
    st.metric("Error en Perigeo", "0.0 km")
    st.metric("Error en Apogeo", "0.0 km")

st.subheader("🪐 Comparativa de Parámetros Orbitales")
data = {
    "Parámetro": ["Semieje Mayor (a)", "Excentricidad (e)", "Altitud Perigeo", "Altitud Apogeo", "Periodo"],
    "Inicial": [
        f"{orb_inicial.a.to(u.km).value:.1f} km",
        f"{orb_inicial.ecc:.4f}",
        f"{(orb_inicial.r_p.to(u.km).value - R_EARTH_KM):.1f} km",
        f"{(orb_inicial.r_a.to(u.km).value - R_EARTH_KM):.1f} km",
        f"{orb_inicial.period.to(u.min).value:.1f} min"
    ],
    "Transferencia": [
        f"{orb_transfer.a.to(u.km).value:.1f} km",
        f"{orb_transfer.ecc:.4f}",
        f"{(orb_transfer.r_p.to(u.km).value - R_EARTH_KM):.1f} km",
        f"{(orb_transfer.r_a.to(u.km).value - R_EARTH_KM):.1f} km",
        f"{orb_transfer.period.to(u.min).value:.1f} min"
    ],
    "Final (Objetivo)": [
        f"{orb_target.a.to(u.km).value:.1f} km",
        f"{orb_target.ecc:.4f}",
        f"{(orb_target.r_p.to(u.km).value - R_EARTH_KM):.1f} km",
        f"{(orb_target.r_a.to(u.km).value - R_EARTH_KM):.1f} km",
        f"{orb_target.period.to(u.min).value:.1f} min"
    ]
}
st.dataframe(pd.DataFrame(data).set_index("Parámetro"), use_container_width=True)

st.subheader("📊 Gráficas")

col_graf1, col_graf2 = st.columns(2)

with col_graf1:
    st.write("**Gráfica Orbital (2D)**")
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    def plot_orbit(orbit, label, color, linestyle='-'):
        nus = np.linspace(0, 2*np.pi, 200)
        a_val = orbit.a.to(u.km).value
        e_val = orbit.ecc
        rs = a_val * (1 - e_val**2) / (1 + e_val * np.cos(nus))
        xs = rs * np.cos(nus)
        ys = rs * np.sin(nus)
        ax.plot(xs, ys, label=label, color=color, linestyle=linestyle, linewidth=2)
    
    plot_orbit(orb_inicial, "Órbita Inicial", "blue")
    plot_orbit(orb_target, "Órbita Objetivo", "green", "--")
    plot_orbit(orb_transfer, "Órbita de Transferencia", "red", ":")
    
    circle = plt.Circle((0, 0), R_EARTH_KM, color='lightblue', label='Tierra', zorder=10)
    ax.add_patch(circle)
    ax.plot(r1, 0, 'ro', markersize=10, label='ΔV1 (Perigeo)', zorder=11)
    ax.plot(-r2, 0, 'go', markersize=10, label='ΔV2 (Apogeo)', zorder=11)
    
    ax.set_xlabel('X (km)', fontsize=12)
    ax.set_ylabel('Y (km)', fontsize=12)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')
    ax.set_title('Vista Orbital 2D', fontsize=14, fontweight='bold')
    
    # Cerrar la figura después de mostrarla para liberar memoria
    st.pyplot(fig)
    plt.close(fig)

with col_graf2:
    st.write("**Altitud vs. Anomalía Verdadera**")
    
    anomalies_deg = np.linspace(0, 360, 400)
    anomalies_rad = np.radians(anomalies_deg)
    
    r_inicial = orb_inicial.a.to(u.km).value * (1 - orb_inicial.ecc**2) / (1 + orb_inicial.ecc * np.cos(anomalies_rad))
    r_target = orb_target.a.to(u.km).value * (1 - orb_target.ecc**2) / (1 + orb_target.ecc * np.cos(anomalies_rad))
    
    alt_inicial_plot = r_inicial - R_EARTH_KM
    alt_target_plot = r_target - R_EARTH_KM
    
    fig_params, ax_params = plt.subplots(figsize=(8, 6))
    ax_params.plot(anomalies_deg, alt_inicial_plot, label="Altitud Inicial", linewidth=2, color='blue')
    ax_params.plot(anomalies_deg, alt_target_plot, label="Altitud Objetivo", linestyle="--", linewidth=2, color='green')
    ax_params.axhline(y=alt_p_inicial, color='red', linestyle=':', alpha=0.5, label=f'Perigeo inicial ({alt_p_inicial:.0f} km)')
    ax_params.axhline(y=alt_a_inicial, color='orange', linestyle=':', alpha=0.5, label=f'Apogeo inicial ({alt_a_inicial:.0f} km)')
    
    ax_params.set_xlabel("Anomalía Verdadera (grados)", fontsize=12)
    ax_params.set_ylabel("Altitud (km)", fontsize=12)
    ax_params.legend()
    ax_params.grid(True, alpha=0.3)
    ax_params.set_title('Variación de Altitud en la Órbita', fontsize=14, fontweight='bold')
    
    # Cerrar la figura después de mostrarla para liberar memoria
    st.pyplot(fig_params)
    plt.close(fig_params)

st.subheader("ℹ️ Información Adicional")
st.write(f"""
**Requisitos de combustible estimado:**
- Asumiendo ISP específico de 300s (típico para propulsores químicos)
- Masa del satélite: {sat_mass:.1f} kg
- Masa de combustible requerida: {sat_mass * (np.exp(dv_total*1000/(300*9.81)) - 1):.2f} kg
- Proporción combustible/satélite: {(np.exp(dv_total*1000/(300*9.81)) - 1)*100:.1f}%
""")