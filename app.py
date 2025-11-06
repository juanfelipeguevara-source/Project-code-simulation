import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u
from astropy import constants as const

# --- Configuración de la Página ---
st.set_page_config(
    page_title="Corrección Orbital (Hohmann)",
    page_icon="🛰️",
    layout="wide"
)

# Configurar matplotlib para usar backend Agg
plt.switch_backend('Agg')

st.title("🛰️ Calculadora de Corrección Orbital")
st.write("""
Esta aplicación calcula la maniobra (Transferencia de Hohmann) necesaria para 
pasar de una órbita inicial a una órbita objetivo, cumpliendo con tus requisitos.
""")

# --- Constantes ---
R_EARTH_KM = 6378.137  # Radio terrestre en km
MU_EARTH = 398600.4418  # Parámetro gravitacional terrestre en km³/s²

# --- Funciones auxiliares ---
def calc_period(a_km):
    """Calcula el periodo orbital en segundos"""
    return 2 * np.pi * np.sqrt(a_km**3 / MU_EARTH)

def calc_velocity(r_km, a_km):
    """Calcula la velocidad en una órbita en km/s"""
    return np.sqrt(MU_EARTH * (2/r_km - 1/a_km))

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

# --- Validaciones ---
if alt_p_inicial > alt_a_inicial:
    st.error("Error: La altitud del perigeo inicial no puede ser mayor que la del apogeo.")
    st.stop()

# --- Calcular parámetros orbitales ---
try:
    # Órbita inicial
    r_p_inicial = R_EARTH_KM + alt_p_inicial
    r_a_inicial = R_EARTH_KM + alt_a_inicial
    a_inicial = (r_p_inicial + r_a_inicial) / 2
    ecc_inicial = (r_a_inicial - r_p_inicial) / (r_a_inicial + r_p_inicial)
    period_inicial = calc_period(a_inicial)
    
    # Órbita objetivo (circular)
    r_target = R_EARTH_KM + alt_target
    a_target = r_target
    ecc_target = 0.0
    period_target = calc_period(a_target)
    
    # Órbita de transferencia
    r1 = r_p_inicial
    r2 = r_target
    a_transfer = (r1 + r2) / 2
    ecc_transfer = (r2 - r1) / (r2 + r1)
    period_transfer = calc_period(a_transfer)
    
    # Cálculo de Delta-V
    v1_inicial = calc_velocity(r1, a_inicial)
    v1_transfer = calc_velocity(r1, a_transfer)
    dv1_mag = abs(v1_transfer - v1_inicial)
    
    v2_transfer = calc_velocity(r2, a_transfer)
    v2_circular = calc_velocity(r2, a_target)
    dv2_mag = abs(v2_circular - v2_transfer)
    
    dv_total = dv1_mag + dv2_mag
    
    # Tiempo de transferencia (medio periodo de la órbita de transferencia)
    T_transfer = period_transfer / 2

except Exception as e:
    st.error(f"Error al calcular: {e}")
    st.stop()

# --- Sección de Salidas (Outputs) ---
st.header("Resultados de la Corrección")

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
    err_p_antes = r_p_inicial - r_target
    err_a_antes = r_a_inicial - r_target
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
        f"{a_inicial:.1f} km",
        f"{ecc_inicial:.4f}",
        f"{alt_p_inicial:.1f} km",
        f"{alt_a_inicial:.1f} km",
        f"{period_inicial/60:.1f} min"
    ],
    "Transferencia": [
        f"{a_transfer:.1f} km",
        f"{ecc_transfer:.4f}",
        f"{alt_p_inicial:.1f} km",
        f"{alt_target:.1f} km",
        f"{period_transfer/60:.1f} min"
    ],
    "Final (Objetivo)": [
        f"{a_target:.1f} km",
        f"{ecc_target:.4f}",
        f"{alt_target:.1f} km",
        f"{alt_target:.1f} km",
        f"{period_target/60:.1f} min"
    ]
}
st.dataframe(pd.DataFrame(data).set_index("Parámetro"), use_container_width=True)

st.subheader("📊 Gráficas")

col_graf1, col_graf2 = st.columns(2)

with col_graf1:
    st.write("**Gráfica Orbital (2D)**")
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    def plot_orbit(a, e, label, color, linestyle='-'):
        nus = np.linspace(0, 2*np.pi, 200)
        rs = a * (1 - e**2) / (1 + e * np.cos(nus))
        xs = rs * np.cos(nus)
        ys = rs * np.sin(nus)
        ax.plot(xs, ys, label=label, color=color, linestyle=linestyle, linewidth=2)
    
    plot_orbit(a_inicial, ecc_inicial, "Órbita Inicial", "blue")
    plot_orbit(a_target, ecc_target, "Órbita Objetivo", "green", "--")
    plot_orbit(a_transfer, ecc_transfer, "Órbita de Transferencia", "red", ":")
    
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
    
    st.pyplot(fig)
    plt.close(fig)

with col_graf2:
    st.write("**Altitud vs. Anomalía Verdadera**")
    
    anomalies_deg = np.linspace(0, 360, 400)
    anomalies_rad = np.radians(anomalies_deg)
    
    r_inicial_plot = a_inicial * (1 - ecc_inicial**2) / (1 + ecc_inicial * np.cos(anomalies_rad))
    r_target_plot = a_target * (1 - ecc_target**2) / (1 + ecc_target * np.cos(anomalies_rad))
    
    alt_inicial_plot = r_inicial_plot - R_EARTH_KM
    alt_target_plot = r_target_plot - R_EARTH_KM
    
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
