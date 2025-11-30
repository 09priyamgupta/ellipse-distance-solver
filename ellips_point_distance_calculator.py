import numpy as np
import matplotlib.pyplot as plt
import random

# --------------------------- Ellipse Point Distance Class ---------------------------
class EllipsePointDistance:
    '''
        Computes closest/farthest points on a general conic ellipse to a given point (u0, v0).
        
        The ellipse is defined by the implicit equation:
                        A*u² + B*v² + C*u*v + D*u + E*v + F = 0

        The algorithm follows the analytic quartic approach:
        - Parameterizes constrained optimization using Lagrange multiplier λ.
        - Builds quartic polynomial in λ.
        - Computes all real roots.
        - Calculates corresponding extremal points.
    '''
    def __init__(self, u0, v0):
        # Fixed ellipse coefficients
        self.A = 2
        self.B = 1
        self.C = 0.5
        self.D = -1
        self.E = 2
        self.F = -3

        # Input point (u0, v0)
        self.u0 = u0
        self.v0 = v0

        # Determinant check for ellipse validity
        det = self.A * self.B - (self.C/2)**2
        if det <= 0:
            print("Invalid ellipse.")
            print(f"Determinant = {det:.3f} <= 0")

        # Precompute matrix form of quadratic terms
        self.M_tilde = np.array([[self.A, self.C/2],
                                 [self.C/2, self.B]])
        
        self.a = np.array([self.D/2, self.E/2])

        print(f"Ellipse: {self.A:.3f}u² + {self.B:.3f}v² + {self.C:.3f}uv + {self.D:.3f}u + {self.E:.3f}v + {self.F:.3f} = 0")
        print(f"Input point: u0 = {self.u0:.3f}, v0 = {self.v0:.3f}")

    def get_ellipse_parameters(self):
        '''
            Return the parameters (A, B, C, D, E, F) of the ellipse equation.
        '''
        return self.A, self.B, self.C, self.D, self.E, self.F

    def find_extremal_points(self):
        '''
            Find the closest and farthest points on the ellipse from the given input point (u0, v0).

            Inputs:
                u0, v0 : Coordinates of the 2D input point.

            Returns:
                closest_point, farthest_point, all_real_points
        '''
        # Extract input point
        u0, v0 = self.u0, self.v0

        # Define the polynomial coefficients P1, P2, P3 for the lagrange multiplier λ
        # P1^2(λ) = (EC/4 - DB/2)λ^2 + (B*u0 - D/2 - C*v0/2)λ + u0
        P1_coeffs = np.array([
            (self.E*self.C/4 - self.D*self.B/2),
            (self.B*u0 - self.D/2 - self.C*v0/2),
            u0
        ])
        
        # P2^2(λ) = (DC/4 - EA/2)λ^2 + (A*v0 - E/2 - C*u0/2)λ + v0
        P2_coeffs = np.array([
            (self.D*self.C/4 - self.E*self.A/2),
            (self.A*v0 - self.E/2 - self.C*u0/2),
            v0
        ])
        
        # P3^2(λ) = (AB - C^2/4)λ^2 + (A+B)λ + 1
        P3_coeffs = np.array([
            (self.A*self.B - self.C**2/4),
            (self.A + self.B),
            1
        ])

        # Construct polynomial terms: P^T M P, 2P3 a^T P, F P3^2

        # Term 1: P^T M P = A*P1^2 + B*P2^2 + C*P1*P2
        P1_sq = np.polymul(P1_coeffs, P1_coeffs)
        P2_sq = np.polymul(P2_coeffs, P2_coeffs)
        P1_P2 = np.polymul(P1_coeffs, P2_coeffs)
        term1 = self.A * P1_sq + self.B * P2_sq + self.C * P1_P2
        
        # Term 2: 2P3 a^T P = P3*(D*P1 + E*P2)
        P3_P1 = np.polymul(P3_coeffs, P1_coeffs)
        P3_P2 = np.polymul(P3_coeffs, P2_coeffs)
        term2 = self.D * P3_P1 + self.E * P3_P2
        
        # Term 3: F * P3^2
        P3_sq = np.polymul(P3_coeffs, P3_coeffs)
        term3 = self.F * P3_sq
        
        # Final quartic polynomial whose roots give extremal points
        quartic_poly = term1 + term2 + term3
        
        # Find roots of the quartic polynomial
        roots = np.roots(quartic_poly)
        
        real_points = []

        # Filter real roots and compute corresponding points
        for lam in roots:
            if np.isreal(lam):
                lam_real = np.real(lam)
                
                # Compute denominator P3^2(λ)
                denom = (lam_real*self.A + 1)*(lam_real*self.B + 1) - (lam_real*self.C/2)**2
                
                if abs(denom) > 1e-10:
                    # Compute the point using the matrix formula
                    inv_matrix = np.linalg.inv([
                        [lam_real*self.A + 1, lam_real*self.C/2],
                        [lam_real*self.C/2, lam_real*self.B + 1]
                    ])
                    
                    right_vector = np.array([
                        u0 - lam_real*self.D/2,
                        v0 - lam_real*self.E/2
                    ])
                    
                    point = inv_matrix @ right_vector
                    u_pt, v_pt = point
                    
                    # Verify the point is on the ellipse (within tolerance)
                    ellipse_eq = (self.A * u_pt**2 + self.B * v_pt**2 + self.C * u_pt * v_pt + 
                                 self.D * u_pt + self.E * v_pt + self.F)
                    
                    if abs(ellipse_eq) < 1e-6:
                        real_points.append((u_pt, v_pt))
        
        # Find closest and farthest points
        if real_points:
            distances = [np.sqrt((p[0] - u0)**2 + (p[1] - v0)**2) for p in real_points]
            closest_idx = np.argmin(distances)
            farthest_idx = np.argmax(distances)
            
            return {
                        "closest_point": real_points[closest_idx],
                        "closest_distance": distances[closest_idx],
                        "farthest_point": real_points[farthest_idx],
                        "farthest_distance": distances[farthest_idx],
                        "all_points": real_points,
                        "all_distances": distances,
                        "roots": roots
                    }
        
        return None
    
# --------------------------- Helper functions ---------------------------
def get_user_input():
    '''
        Get 2D point coordinates (u0, v0) from the user via terminal.
    '''
    print("\nEnter 2D point coordinates (u, v):")

    while True:
        try:
            u0 = float(input("Enter u coordinate: "))
            v0 = float(input("Enter v coordinate: "))
            return u0, v0
            
        except ValueError:
            print("Error: Please enter valid numeric values.")


def generate_random_2D_point(range=5):
    '''
        Generate a random 2D input point (u0, v0) within [-r, r] as fallback.
    '''
    u0 = random.uniform(-range, range)
    v0 = random.uniform(-range, range)

    return u0, v0


def draw_ellipse(ax, A, B, C, D, E, F, r=5, points=1000):
    '''
        Draw the ellipse defined by the quadratic form on the given axes.
    '''
    u_vals = np.linspace(-r, r, points)
    v_vals = np.linspace(-r, r, points)

    U, V = np.meshgrid(u_vals, v_vals)
    Z = (A*U**2 + B*V**2 + C*U*V + D*U + E*V + F)

    contour = ax.contour(U, V, Z, levels=[0], colors='blue', linewidths=2)
    return contour


# --------------------------- Main pipeline ---------------------------
def main():
    print("\nChoose input method:")
    print("1. Enter ellipse points (u0, v0) manually")
    print("2. Use randomly generated input point (u0, v0)")

    while True:
        choice = input("Enter your choice (1 or 2): ").strip()

        if choice == '1':
            u0, v0 = get_user_input()
            break

        elif choice == '2':
            u0, v0 = generate_random_2D_point()
            print(f"Random point selected: ({u0:.3f}, {v0:.3f})")
            break

        else:
            print("Invalid choice. Please enter 1 or 2.")

    print("\n" + "="*50)
    print("Ellipse Point Distance Calculator:")
    print("="*50)

    # Initialize the ellipse distance calculator
    ellipse_solver = EllipsePointDistance(u0, v0)
    A, B, C, D, E, F = ellipse_solver.get_ellipse_parameters()

    # Create visualization
    fig, ax = plt.subplots(figsize=(10, 8))
    draw_ellipse(ax, A, B, C, D, E, F, r=10)

    # Plot input point
    ax.plot(u0, v0, 'ro', markersize=8, label=f'Input Point ({u0:.2f}, {v0:.2f})')
    
    # Compute and display results immediately
    print("\nComputing extremal points...")
    result = ellipse_solver.find_extremal_points()

    if result is None:
        print("No valid extremal points found.")
        return

    # Extract results
    closest = result["closest_point"]
    closest_dist = result["closest_distance"]
    farthest = result["farthest_point"]
    farthest_dist = result["farthest_distance"]
    all_points = result["all_points"]
    all_dists = result["all_distances"]
    roots = result["roots"]
    
    print(f"Number of real extremal points found: {len(all_points)}\n")

    print("Quartic roots λ:")
    for r in roots:
        print(f"   λ = {r}")

    print("\nCandidate extremal points:")
    for p, d in zip(all_points, all_dists):
        print(f"   Point: (u={p[0]:.6f}, v={p[1]:.6f})   Distance = {d:.6f}")

    print("\nClosest Point:")
    print(f"   Coordinates: (u={closest[0]:.6f}, v={closest[1]:.6f})")
    print(f"   Distance: {closest_dist:.6f}")

    print("\nFarthest Point:")
    print(f"   Coordinates: (u={farthest[0]:.6f}, v={farthest[1]:.6f})")
    print(f"   Distance: {farthest_dist:.6f}")

    # Plot all real solution points
    for pt in all_points:
        ax.plot(pt[0], pt[1], 'go', alpha=0.6)
        ax.plot([u0, pt[0]], [v0, pt[1]], 'gray', linestyle='--', alpha=0.4)

    if closest:
        ax.plot(closest[0], closest[1], 'g*', markersize=15, label='Closest')

    if farthest:
        ax.plot(farthest[0], farthest[1], 'r*', markersize=15, label='Farthest')

    ax.set_aspect('equal')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.show()
    
# --------------------------- Run Main ---------------------------
if __name__ == "__main__":
    main()