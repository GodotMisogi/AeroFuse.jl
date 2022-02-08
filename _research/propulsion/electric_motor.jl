"""
Reference: Drela's QPROP DC Electric Motor Model - http://web.mit.edu/drela/Public/web/qprop/motor1_theory.pdf
"""
## Model relations
shaft_torque(i, i₀, K_Q) = (i - i₀) / K_Q
back_EMF(Ω, K_v) = Ω / K_v
motor_terminal_voltage(i, Ω, R, K_v) = back_EMF(Ω, K_v) + i * R

# Inverse relations
current(v, Ω, R, K_v) = (v - Ω / K_v) / R
motor_efficiency(Ω, v)

