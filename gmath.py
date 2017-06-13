import math

e = 120
V = [0, 0, 1]

def calculate_normal(polygons, i):

    A = [0, 0, 0]
    B = [0, 0, 0]
    N = [0, 0, 0]

    A[0] = polygons[i+1][0] - polygons[i][0]
    A[1] = polygons[i+1][1] - polygons[i][1]
    A[2] = polygons[i+1][2] - polygons[i][2]

    B[0] = polygons[i+2][0] - polygons[i][0]
    B[1] = polygons[i+2][1] - polygons[i][1]
    B[2] = polygons[i+2][2] - polygons[i][2]

    N[0] = A[1] * B[2] - A[2] * B[1]
    N[1] = A[2] * B[0] - A[0] * B[2]
    N[2] = A[0] * B[1] - A[1] * B[0]

    return N

def calculate_mag(vector):
    return math.sqrt( vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2] )

def scalar_mult(vector, scalar):
    return [ scalar * i for i in vector ]

def vector_subtraction(a, b):
    return [ (a[0] - b[0]), (a[1] - b[1]), (a[2] - b[2]) ]

def to_unit_vector(vector):
    v = [0, 0, 0]
    mag = calculate_mag(vector)

    v[0] = (vector[0])/mag
    v[1] = (vector[1])/mag
    v[2] = (vector[2])/mag
    return v

def calculate_dot_product(a, b):
    a = to_unit_vector(a)
    b = to_unit_vector(b)

    p = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    return p if p >= 0 else 0

def i_ambient(L_c, K_a):
    return L_c * K_a

def i_diffuse(L_c, K_d, N, L):
    p = calculate_dot_product(N, L)
    return L_c * K_d * p

def i_specular(L_c, K_s, N, L):
    N_L = calculate_dot_product(N, L)
    p = scalar_mult(scalar_mult(N, 2), N_L)
    d = vector_subtraction(p, L)
    p = calculate_dot_product(d, V)

    return L_c * K_s * (p ** e)

def intensity(i_ambient, i_diffuse, i_specular):
    s = i_ambient + i_diffuse + i_specular
    return int(round(s)) if s <= 255 else 255

def calc_light(N, constants, light_source):
    L = light_source['location']

    K_a_r = constants["red"][0]
    K_d_r = constants["red"][1]
    K_s_r = constants["red"][2]

    K_a_g = constants["green"][0]
    K_d_g = constants["green"][1]
    K_s_g = constants["green"][2]

    K_a_b = constants["blue"][0]
    K_d_b = constants["blue"][1]
    K_s_b = constants["blue"][2]

    L_r = light_source['color'][0]
    L_g = light_source['color'][1]
    L_b = light_source['color'][2]

    red = intensity(i_ambient(L_r, K_a_r), i_diffuse(L_r, K_d_r, N, L), i_specular(L_r, K_s_r, N, L))
    green = intensity(i_ambient(L_g, K_a_g), i_diffuse(L_g, K_d_g, N, L), i_specular(L_g, K_s_g, N, L))
    blue = intensity(i_ambient(L_b, K_a_b), i_diffuse(L_b, K_d_b, N, L), i_specular(L_b, K_s_b, N, L))

    return [red, green, blue]

def calc_total_light(N, constants, light_sources):
    colour = [0, 0, 0]

    for light_source in light_sources:
        light = calc_light(N, constants, light_source)
        colour[0] += light[0]
        colour[1] += light[1]
        colour[2] += light[2]

    colour[0] = colour[0] if colour[0] <= 255 else 255
    colour[1] = colour[1] if colour[1] <= 255 else 255
    colour[2] = colour[2] if colour[2] <= 255 else 255
    return colour

def calc_avg_normal(normals):
    x_sum = 0.0
    y_sum = 0.0
    z_sum = 0.0

    for i in normals:
        x_sum += i[0]
        y_sum += i[1]
        z_sum += i[2]

    length = len(normals)
    return [ x_sum/length, y_sum/length, z_sum/length ]

def list_vertex_normals(matrix):
    v_n = {}

    point = 0
    while point < len(matrix) - 2:
        normal = calculate_normal(matrix, point)[:]

        if (int(matrix[point][0]), int(matrix[point][1]), matrix[point][2]) in v_n:
            v_n[(int(matrix[point][0]), int(matrix[point][1]), matrix[point][2])].append(normal)
        else:
            v_n[(int(matrix[point][0]), int(matrix[point][1]), matrix[point][2])] = [ normal ]

        if (int(matrix[point+1][0]), int(matrix[point+1][1]), matrix[point+1][2]) in v_n:
            v_n[(int(matrix[point+1][0]), int(matrix[point+1][1]), matrix[point+1][2])].append(normal)
        else:
            v_n[(int(matrix[point+1][0]), int(matrix[point+1][1]), matrix[point+1][2])] = [ normal ]

        if (int(matrix[point+2][0]), int(matrix[point+2][1]), matrix[point+2][2]) in v_n:
            v_n[(int(matrix[point+2][0]), int(matrix[point+2][1]), matrix[point+2][2])].append(normal)
        else:
            v_n[(int(matrix[point+2][0]), int(matrix[point+2][1]), matrix[point+2][2])] = [ normal ]

        point += 3

    for key in v_n:
        v_n[key] = calc_avg_normal(v_n[key])

    return v_n
