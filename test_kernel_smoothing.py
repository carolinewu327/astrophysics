#!/usr/bin/env python
"""
Quick test of the new kernel smoothing functionality.
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, 'analysis/boss/scripts')

from stack_single import extract_radial_profile_kde, extract_radial_profile_binned
from constants import GRID_SIZE, BOX_SIZE_HMPC

# Create a synthetic stacked map with a Gaussian signal
x = np.linspace(-50, 50, GRID_SIZE)
y = np.linspace(-50, 50, GRID_SIZE)
X, Y = np.meshgrid(x, y)
R = np.sqrt(X**2 + Y**2)

# Synthetic signal: Gaussian peak + noise
signal = 0.02 * np.exp(-R**2 / (2 * 10**2))  # Peak at center, sigma=10 Mpc/h
noise = np.random.normal(0, 0.005, (GRID_SIZE, GRID_SIZE))
kappa_map = signal + noise

# Test kernel smoothing with different kernel widths
kernel_widths = [2.0, 3.0, 5.0, 8.0]
plt.figure(figsize=(12, 8))

# Plot original map
plt.subplot(2, 3, 1)
plt.imshow(kappa_map, extent=[-50, 50, -50, 50], origin='lower', cmap='RdBu_r')
plt.colorbar(label='κ')
plt.title('Synthetic κ Map')
plt.xlabel('x [h⁻¹ Mpc]')
plt.ylabel('y [h⁻¹ Mpc]')

# Traditional binned profile
r_binned, kappa_binned, _ = extract_radial_profile_binned(kappa_map, n_bins=50)
plt.subplot(2, 3, 2)
plt.plot(r_binned, kappa_binned, 'o-', label='Binned', alpha=0.7)
plt.axhline(0, color='k', linestyle='--', alpha=0.3)
plt.xlabel('r [h⁻¹ Mpc]')
plt.ylabel('κ(r)')
plt.title('Traditional Binned Profile')
plt.legend()
plt.grid(True, alpha=0.3)

# KDE profiles with different kernel widths
plt.subplot(2, 3, 3)
for kw in kernel_widths:
    r_kde, kappa_kde, _ = extract_radial_profile_kde(
        kappa_map, kernel_width_hmpc=kw, n_radial_bins=500
    )
    plt.plot(r_kde, kappa_kde, label=f'σ={kw} Mpc/h', alpha=0.7)

plt.axhline(0, color='k', linestyle='--', alpha=0.3)
plt.xlabel('r [h⁻¹ Mpc]')
plt.ylabel('κ(r)')
plt.title('KDE Profiles (varying kernel width)')
plt.legend()
plt.grid(True, alpha=0.3)

# Comparison: binned vs KDE
plt.subplot(2, 3, 4)
plt.plot(r_binned, kappa_binned, 'o', label='Binned (50 bins)', alpha=0.6, markersize=4)

r_kde, kappa_kde, _ = extract_radial_profile_kde(kappa_map, kernel_width_hmpc=3.0)
plt.plot(r_kde, kappa_kde, '-', label='KDE (σ=3 Mpc/h)', linewidth=2)

# True signal
r_true = np.linspace(0, 50, 500)
signal_true = 0.02 * np.exp(-r_true**2 / (2 * 10**2))
plt.plot(r_true, signal_true, 'k--', label='True signal', linewidth=1.5)

plt.axhline(0, color='k', linestyle='--', alpha=0.3)
plt.xlabel('r [h⁻¹ Mpc]')
plt.ylabel('κ(r)')
plt.title('Binned vs KDE vs Truth')
plt.legend()
plt.grid(True, alpha=0.3)

# Residuals
plt.subplot(2, 3, 5)
# Interpolate true signal to KDE radii
signal_interp = 0.02 * np.exp(-r_kde**2 / (2 * 10**2))
residuals = kappa_kde - signal_interp
plt.plot(r_kde, residuals, label='KDE - True')
plt.axhline(0, color='k', linestyle='--', alpha=0.3)
plt.xlabel('r [h⁻¹ Mpc]')
plt.ylabel('Residual κ')
plt.title('KDE Recovery (σ=3 Mpc/h)')
plt.legend()
plt.grid(True, alpha=0.3)

# Effect of kernel width on smoothing
plt.subplot(2, 3, 6)
for kw in kernel_widths:
    r_kde, kappa_kde, _ = extract_radial_profile_kde(
        kappa_map, kernel_width_hmpc=kw, n_radial_bins=500
    )
    signal_interp = 0.02 * np.exp(-r_kde**2 / (2 * 10**2))
    rms = np.sqrt(np.mean((kappa_kde - signal_interp)**2))
    plt.scatter(kw, rms, s=100)

plt.xlabel('Kernel width σ [h⁻¹ Mpc]')
plt.ylabel('RMS (KDE - True)')
plt.title('Recovery vs Kernel Width')
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('test_kernel_smoothing.png', dpi=150)
print("✓ Test complete! Saved test_kernel_smoothing.png")
print("\nKey findings:")
print(f"  - Binned profile has {len(r_binned)} points")
print(f"  - KDE profile has {len(r_kde)} points (smooth)")
print(f"  - Kernel widths tested: {kernel_widths} Mpc/h")
print("  - Larger kernels → smoother but potentially oversmoothed")
print("  - Recommended: kernel_width ≈ Planck resolution (3-4 Mpc/h)")
