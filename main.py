import numpy as np
from IPython.display import display, Math, Latex, Markdown
import sympy as sp
# Condition and geometry information
rho = 0.002377  # slugs/cubic feet
# rho = rho * 32.174 # to lb / cubic feet
rho = 0.0765  # lb/ft^3
U0 = 165  # Knots
U0 = U0 * 1.68781  # to feet/second
S = 5500  # square feet
b = 195.68  # feet
c = 27.31  # feet
theta_0 = 0
m = 564000  # lb

Iyy = 32.3e6  # slug-ft^2
Ixx = 14.3e6
Izz = 45.3e6
Ixz = -2.23e6

SIunits = False

if SIunits:
    # Conversion factors
    ft_to_metre = 0.3048
    lb_to_kg = 0.453592
    slug_to_kg = 14.5939

    MOIconvert = slug_to_kg * ft_to_metre ** 2

    # Constants
    a = 340  # sonic velocity in m/s
    g = 9.80665  # acceleration due to gravity
    rho = 1.225  # density in kg/m^3

    # Convert
    U0 = U0 * ft_to_metre
    S = S * ft_to_metre ** 2
    b = b * ft_to_metre
    c = c * ft_to_metre
    m = m * lb_to_kg
    Ixx = Ixx * MOIconvert
    Iyy = Iyy * MOIconvert
    Izz = Izz * MOIconvert
    Ixz = Ixz * MOIconvert

else:
    a = 1125.33  # sonic velocity in ft/s
    g = 32.174

    # Convert mass moments of inertia to consistent units --> INTO lb-ft^2 FROM slug-ft^2
    Ixx = Ixx * g
    Iyy = Iyy * g
    Izz = Izz * g
    Ixz = Ixz * g

# Store the derivatives in two dictionaries
B747_lon_ders = {'C_L': 1.11, 'C_D': 0.102, 'C_L_a': 5.7, 'C_D_a': 0.66, 'C_m_a': -1.26,
                 'C_L_da': -6.7, 'C_m_da': -3.2, 'C_L_hq': 5.4, 'C_m_hq': -20.8, 'C_L_M': -0.81,
                 'C_m_M': 0.27, 'C_L_de': 0.338, 'C_m_de': -1.34}

B747_lat_ders = {'C_y_b': -0.96, 'C_l_b': -0.221, 'C_n_b': 0.150, 'C_l_hp': -0.45, 'C_n_hp': -0.121,
                 'C_l_hr': 0.101, 'C_n_hr': -0.30, 'C_l_da': 0.0461, 'C_n_da': 0.0064, 'C_y_dr': 0.175,
                 'C_l_dr': 0.007, 'C_n_dr': -0.109}

# Put the dictionaries into the local namespace
# This might seem a bit convoluted but it enables us to store values of derivatives in the dicionary, above,
# and then put them all into the global namespace.
#
# It's fairly easy to get (key: value) pairs into a dictionary from a text file or xls so this will be handy
# later.
locals().update(B747_lon_ders)

# Note that without doing the dictionary --> local namespace, the we'd have to write:
# Xu = q * S / m / U0 * (2 * B747_lon_ders["C_D"] + M * B747_lon_ders["C_D_M"])

# Convert to dimensional form
q = 0.5 * rho * U0 ** 2
M = U0 / a

Xu = -q * S / m / U0 * (2 * C_D)  # No C_D_M term so assumed zero
Xw = q * S / m / U0 * (C_L - C_D_a)
Xq = 0  # No CDq term given in the table so assumed zero
Zu = -q * S / m / U0 * (2 * C_L + M * C_L_M)
Zw = -q * S / m / U0 * (C_D + C_L_a)
Zdw = q * S * c / m / 2 / U0 ** 2 * C_L_da  # This is a NEW term for us, but since it was given as C_L_da, must be included
Zq = -q * S * c / 2 / m / U0 * C_L_hq
Mu = q * S * c / Iyy / U0 * M * C_m_M
Mw = q * S * c / Iyy / U0 * C_m_a
Mdw = q * S * c ** 2 / 2 / Iyy / U0 ** 2 * C_m_da
Mq = q * S * c ** 2 / 2 / Iyy / U0 * C_m_hq
Mq = q * S * c ** 2 / 2 / Iyy / U0 * C_m_hq
Zde = -q * S / m * C_L_de
Mde = q * S * c / Iyy * C_m_de


Mustar = Mu + Mdw * Zu
Mwstar = Mw + Mdw * Zw
Mqstar = Mq + Mdw * Zq
Mthetastar = -Mdw * g * np.sin(theta_0)
Mdestar = Mde + Mdw * Zde

print(f"Xu = {Xu:1.4f}")
print(f"Xw = {Xw:1.4f}")
print(f"Zu = {Zu:1.4f}")
print(f"Zw = {Zw:1.4f}")
print(f"Zdw = {Zdw:1.4f}")
print(f"Zq = {Zq:1.4f}")
print(f"Mu = {Mu:1.4f}")
print(f"Mw = {Mw:1.4f}")
print(f"Mw = {Mw:1.4f}")
print(f"Mdw = {Mdw:1.4f}")
print(f"Mq = {Mq:1.4f}")

print(f"Zde = {Zde:1.4f}")
print(f"Mde = {Mde:1.4f}")


Alon = np.matrix([[Xu, Xw, 0, -g*np.cos(theta_0)],
               [Zu, Zw, U0 + Zq, -g*np.sin(theta_0)],
               [Mustar, Mwstar, Mqstar, Mthetastar],
               [0, 0, 1, 0]])

print("The system matrix for longitudinal motion is ")
display( Math('[A_{lon}] = ' + sp.latex(sp.Matrix(Alon).evalf(5))))

# Get the eigenvalues - this is my own checking I got the matrix correct. It'll make sense later.
eigs, _ = np.linalg.eig(Alon)

print("The eigenvalues are:", eigs)

print("The coefficients of the CE are:", np.poly(Alon))
