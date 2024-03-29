def damperposition(theta):
    # End-point
    H = L_shortbeam + L_shortbeam * np.cos(theta)
    V = L_shortbeam * np.sin(theta) - L_base
    R = np.sqrt(H ** 2 + V ** 2)
    # Blue beam
    H_blue = L_shortbeam + 0.5 * L_shortbeam * np.cos(theta)
    V_blue = 0.5 * L_shortbeam * np.sin(theta) - L_base
    R_blue = np.sqrt(H_blue ** 2 + V_blue ** 2)
    # Green beam
    R_green = np.sqrt(
        (0.5 * L_longbeam) ** 2 + (0.25 * L_longbeam) ** 2 - 0.25 * (L_longbeam ** 2) * np.cos(np.pi - phi(theta)))
    H_green = 0.5 * L_longbeam * np.cos(theta) + 0.25 * L_longbeam * np.cos(phi(theta) - theta)
    V_green = np.sqrt(R_green ** 2 - H_green ** 2)

    # Calculate the difference vectors
    H_diff_green = H - H_green
    V_diff_green = V - V_green

    H_diff_blue = H - H_blue
    V_diff_blue = V - V_blue

    # Calculate the dot product of the difference vectors and the end-point vector
    dot_product_green = H * H_diff_green + V * V_diff_green
    dot_product_blue = H * H_diff_blue + V * V_diff_blue
    dot_product_both = H_diff_green * H_diff_blue + V_diff_green * V_diff_blue

    # Calculate the magnitudes of the end-point vector and the difference vectors
    magnitude_diff_green = np.sqrt(H_diff_green ** 2 + V_diff_green ** 2)
    magnitude_diff_blue = np.sqrt(H_diff_blue ** 2 + V_diff_blue ** 2)

    # Calculate the angles between the difference vectors and the end-point vector
    angle_green = np.arccos(dot_product_green / (R * magnitude_diff_green))
    angle_blue = np.arccos(dot_product_blue / (R * magnitude_diff_blue))

    return angle_green, angle_blue, H, V
def dampertorque(c, theta, omega):
    angle_a, angle_b, H, V = damperposition(theta)
    #Calculate speed and torque
    H_prime = -omega * L_shortbeam * np.sin(theta)
    V_prime = omega * L_shortbeam * np.cos(theta)
    speed = np.sqrt(H_prime ** 2 + V_prime ** 2)
    F_magnitude = c * speed
    F_green = F_magnitude * np.cos(angle_a)
    Torque = F_green * L_shortbeam * np.sin(phi(theta))
    return Torque

def springtorque(k, theta, theta_initial):
    #initial position of the spring
    A1, A2, H_initial, V_initial = damperposition(theta_initial)
    angle_a, angle_b, H, V = damperposition(theta)
    R_initial = np.sqrt(H_initial ** 2 + V_initial ** 2)
    R = np.sqrt(H ** 2 + V ** 2)
    # Calculate displacement and torque
    displacement = abs(R_initial - R)
    F_magnitude = k * displacement
    F_green = F_magnitude * np.cos(angle_a)
    Torque = F_green * L_shortbeam * np.sin(phi(theta))
    return Torque

def springtorqueconstant(k, theta, theta_initial):
    angle_a, angle_b, H, V = damperposition(theta)
    # Calculate displacement and torque
    F_magnitude = k
    F_green = F_magnitude * np.cos(angle_a)
    Torque = F_green * L_shortbeam * np.sin(phi(theta))
    return Torque
