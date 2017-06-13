from display import *
from matrix import *
from math import *
from gmath import *

import random

def scanline_convert(polygons, i, screen, zbuffer, color, shading_type, intensities, constants, light_sources):
    y_list = [ polygons[i][1], polygons[i+1][1], polygons[i+2][1] ]

    b = y_list.index(min(y_list))
    Bpt = polygons[i+b]
    t = y_list.index(max(y_list))
    Tpt = polygons[i+t]
    m = 3 - (t + b)
    Mpt = polygons[i+m]

    y_0 = int(Bpt[1])
    y_1 = int(Tpt[1])
    x_0 = x_1 = int(Bpt[0])
    z_0 = z_1 = Bpt[2]

    x0_coords = scanline_helper(x_0, y_0, z_0, int(Tpt[0]), y_1, Tpt[2])[::-1] if (x_0 > int(Tpt[0])) else scanline_helper(x_0, y_0, z_0, int(Tpt[0]), y_1, Tpt[2])
    x1_coords = scanline_helper(x_0, y_0, z_0, int(Mpt[0]), int(Mpt[1]), Mpt[2])[::-1] if (x_0 > int(Mpt[0])) else scanline_helper(x_0, y_0, z_0, int(Mpt[0]), int(Mpt[1]), Mpt[2])
    i_0 = i_1 = 0

    if shading_type == "goroud":
        I_1 = intensities[t]
        I_2 = intensities[b]
        I_3 = intensities[m]

        I_a = intensities[b][:]
        I_b = intensities[b][:]
        (y_3, y_2) = (int(Mpt[1]), int(Bpt[1]))

        swap = True if (int(Mpt[0]) <= int(Tpt[0])) and (int(Mpt[0]) <= int(Bpt[0])) else False
        if int(Mpt[0]) != x_0 and int(Tpt[0]) != x_0:
            swap = True if (int(Mpt[0]) <= int(Tpt[0])) and (int(Mpt[0]) >= int(Bpt[0])) and (float(y_1-y_0)/(int(Tpt[0])-x_0) <= float(int(Mpt[1])-y_0)/(int(Mpt[0])-x_0)) else swap
            swap = True if (int(Bpt[0]) >= int(Tpt[0])) and (int(Mpt[0]) <= int(Bpt[0])) and (float(y_1-y_0)/(x_0-int(Tpt[0])) >= float(int(Mpt[1])-y_0)/(x_0-int(Mpt[0]))) else swap

    elif shading_type == "phong":
        N_Bpt = intensities[b]
        N_Mpt = intensities[m]
        N_Tpt = intensities[t]

        if y_1 - y_0 != 0:
            d_N_0 = [ float(N_Tpt[0] - N_Bpt[0])/(y_1 - y_0), float(N_Tpt[1] - N_Bpt[1])/(y_1 - y_0), float(N_Tpt[2] - N_Bpt[2])/(y_1 - y_0) ]
            if int(Mpt[1])- y_0 != 0:
                d_N_1 = [ float(N_Mpt[0] - N_Bpt[0])/(int(Mpt[1])-y_0), float(N_Mpt[1] - N_Bpt[1])/(int(Mpt[1])-y_0), float(N_Mpt[2] - N_Bpt[2])/(int(Mpt[1])-y_0) ]

        N_a = N_Bpt[:]
        N_b = N_Bpt[:]

        swap = True if (int(Mpt[0]) <= int(Tpt[0])) and (int(Mpt[0]) <= int(Bpt[0])) else False
        if int(Mpt[0]) != x_0 and int(Tpt[0]) != x_0:
            swap = True if (int(Mpt[0]) <= int(Tpt[0])) and (int(Mpt[0]) >= int(Bpt[0])) and (float(y_1-y_0)/(int(Tpt[0])-x_0) <= float(int(Mpt[1])-y_0)/(int(Mpt[0])-x_0)) else swap
            swap = True if (int(Bpt[0]) >= int(Tpt[0])) and (int(Mpt[0]) <= int(Bpt[0])) and (float(y_1-y_0)/(x_0-int(Tpt[0])) >= float(int(Mpt[1])-y_0)/(x_0-int(Mpt[0]))) else swap

    while y_0 < y_1:
        if y_0 == int(Mpt[1]):
            x_1 = int(Mpt[0])
            z_1 = Mpt[2]
            x1_coords = scanline_helper(x_1, y_0, z_1, int(Tpt[0]), y_1, Tpt[2])[::-1] if (x_1 > int(Tpt[0])) else scanline_helper(x_1, y_0, z_1, int(Tpt[0]), y_1, Tpt[2])
            i_1 = 0

            if shading_type == "goroud" and int(Mpt[1]) != y_1:
                I_b = intensities[m][:]
                I_2 = intensities[m]
                I_3 = intensities[t]
                (y_3, y_2) = (y_1, int(Mpt[1])
                              )

            elif shading_type == "phong" and int(Mpt[1]) != y_1:
                d_N_1 = [ float(N_Tpt[0] - N_Mpt[0])/(y_1 - int(Mpt[1])), float(N_Tpt[1] - N_Mpt[1])/(y_1 - int(Mpt[1])), float(N_Tpt[2] - N_Mpt[2])/(y_1 - int(Mpt[1])) ]
                N_b = N_Mpt[:]

        if shading_type == "goroud":
            if swap:
                draw_line( x_0, y_0, z_0, x_1, y_0, z_1, screen, zbuffer, color, I_b, I_a, shading_type )
            else:
                draw_line( x_0, y_0, z_0, x_1, y_0, z_1, screen, zbuffer, color, I_a, I_b, shading_type )
        elif shading_type == "phong":
            # print N_a
            # print N_b
            if swap:
                draw_line( x_0, y_0, z_0, x_1, y_0, z_1, screen, zbuffer, color, N_b, N_a, shading_type, constants, light_sources )
            else:
                draw_line( x_0, y_0, z_0, x_1, y_0, z_1, screen, zbuffer, color, N_a, N_b, shading_type, constants, light_sources )
        else:
            draw_line( x_0, y_0, z_0, x_1, y_0, z_1, screen, zbuffer, color )

        y_0 += 1
        i_0 += 1
        i_1 += 1
        x_0 = x0_coords[i_0][0]
        z_0 = x0_coords[i_0][1]
        x_1 = x1_coords[i_1][0]
        z_1 = x1_coords[i_1][1]

        if shading_type == "goroud":
            for j in range(len(color)):
                I_a[j] = int(round((float(y_0 - int(Bpt[1]))/(y_1 - int(Bpt[1]))) * I_1[j] + (float(y_1 - y_0)/(y_1 - int(Bpt[1]))) * intensities[b][j]))
                I_b[j] = int(round((float(y_0 - y_2)/(y_3 - y_2)) * I_3[j] + (float(y_3 - y_0)/(y_3 - y_2)) * I_2[j]))

        elif shading_type == "phong":
            for j in range(len(color)):
                N_a[j] += d_N_0[j]
                N_b[j] += d_N_1[j]

    if shading_type == "goroud":
        if swap:
            draw_line( x_0, y_0, z_0, x_1, y_0, z_1, screen, zbuffer, color, I_b, I_1, shading_type )
        else:
            draw_line( x_0, y_0, z_0, x_1, y_0, z_1, screen, zbuffer, color, I_1, I_b, shading_type )
    elif shading_type == "phong":
        if swap:
            draw_line( x_0, y_0, z_0, x_1, y_0, z_1, screen, zbuffer, color, N_b, N_Tpt, shading_type, constants, light_sources )
        else:
            draw_line( x_0, y_0, z_0, x_1, y_0, z_1, screen, zbuffer, color, N_Tpt, N_b, shading_type, constants, light_sources )
    else:
        draw_line( x_0, y_0, z_0, x_1, y_0, z_1, screen, zbuffer, color )

def add_polygon( polygons, x0, y0, z0, x1, y1, z1, x2, y2, z2 ):
    add_point(polygons, x0, y0, z0);
    add_point(polygons, x1, y1, z1);
    add_point(polygons, x2, y2, z2);

def draw_polygons( matrix, screen, zbuffer, color, constants = [], light_sources = [], shading_type = "" ):
    if len(matrix) < 2:
        print 'Need at least 3 points to draw'
        return

    if shading_type == "goroud" or shading_type == "phong":
        v_n = list_vertex_normals(matrix)

    intensities = []
    point = 0
    while point < len(matrix) - 2:

        normal = calculate_normal(matrix, point)[:]
        #print normal
        if normal[2] > 0:
            if shading_type == "flat":
                color = calc_total_light(normal, constants, light_sources)

            elif shading_type == "goroud":
                intensities = [ calc_total_light(v_n[(int(matrix[point][0]), int(matrix[point][1]), matrix[point][2])], constants, light_sources),
                                calc_total_light(v_n[(int(matrix[point+1][0]), int(matrix[point+1][1]), matrix[point+1][2])], constants, light_sources),
                                calc_total_light(v_n[(int(matrix[point+2][0]), int(matrix[point+2][1]), matrix[point+2][2])], constants, light_sources) ]

            elif shading_type == "phong":
                intensities = [ v_n[(int(matrix[point][0]), int(matrix[point][1]), matrix[point][2])],
                                v_n[(int(matrix[point+1][0]), int(matrix[point+1][1]), matrix[point+1][2])],
                                v_n[(int(matrix[point+2][0]), int(matrix[point+2][1]), matrix[point+2][2])] ]
            # print intensities
            # color = [ random.randrange(256), random.randrange(256), random.randrange(256) ]
            scanline_convert(matrix, point, screen, zbuffer, color, shading_type, intensities, constants, light_sources )
            # draw_line( int(matrix[point][0]),
            #            int(matrix[point][1]),
            #            matrix[point][2],
            #            int(matrix[point+1][0]),
            #            int(matrix[point+1][1]),
            #            matrix[point+1][2],
            #            screen, zbuffer, color)
            # draw_line( int(matrix[point+2][0]),
            #            int(matrix[point+2][1]),
            #            matrix[point+2][2],
            #            int(matrix[point+1][0]),
            #            int(matrix[point+1][1]),
            #            matrix[point+1][2],
            #            screen, zbuffer, color)
            # draw_line( int(matrix[point][0]),
            #            int(matrix[point][1]),
            #            matrix[point][2],
            #            int(matrix[point+2][0]),
            #            int(matrix[point+2][1]),
            #            matrix[point+2][2],
            #            screen, zbuffer, color)
        point += 3


def add_box( polygons, x, y, z, width, height, depth ):
    x1 = x + width
    y1 = y - height
    z1 = z - depth

    #front
    add_polygon(polygons, x, y, z, x1, y1, z, x1, y, z);
    add_polygon(polygons, x, y, z, x, y1, z, x1, y1, z);

    #back
    add_polygon(polygons, x1, y, z1, x, y1, z1, x, y, z1);
    add_polygon(polygons, x1, y, z1, x1, y1, z1, x, y1, z1);

    #right side
    add_polygon(polygons, x1, y, z, x1, y1, z1, x1, y, z1);
    add_polygon(polygons, x1, y, z, x1, y1, z, x1, y1, z1);
    #left side
    add_polygon(polygons, x, y, z1, x, y1, z, x, y, z);
    add_polygon(polygons, x, y, z1, x, y1, z1, x, y1, z);

    #top
    add_polygon(polygons, x, y, z1, x1, y, z, x1, y, z1);
    add_polygon(polygons, x, y, z1, x, y, z, x1, y, z);
    #bottom
    add_polygon(polygons, x, y1, z, x1, y1, z1, x1, y1, z);
    add_polygon(polygons, x, y1, z, x, y1, z1, x1, y1, z1);

def add_sphere( edges, cx, cy, cz, r, step ):
    points = generate_sphere(cx, cy, cz, r, step)
    num_steps = int(1/step+0.1)

    lat_start = 0
    lat_stop = num_steps
    longt_start = 0
    longt_stop = num_steps

    num_steps+= 1
    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):

            p0 = lat * (num_steps) + longt
            p1 = p0+1
            p2 = (p1+num_steps) % (num_steps * (num_steps-1))
            p3 = (p0+num_steps) % (num_steps * (num_steps-1))

            if longt != num_steps - 2:
	        add_polygon( edges, points[p0][0],
		             points[p0][1],
		             points[p0][2],
		             points[p1][0],
		             points[p1][1],
		             points[p1][2],
		             points[p2][0],
		             points[p2][1],
		             points[p2][2])
            if longt != 0:
	        add_polygon( edges, points[p0][0],
		             points[p0][1],
		             points[p0][2],
		             points[p2][0],
		             points[p2][1],
		             points[p2][2],
		             points[p3][0],
		             points[p3][1],
		             points[p3][2])

def generate_sphere( cx, cy, cz, r, step ):
    points = []
    num_steps = int(1/step+0.1)

    rot_start = 0
    rot_stop = num_steps
    circ_start = 0
    circ_stop = num_steps

    for rotation in range(rot_start, rot_stop):
        rot = step * rotation
        for circle in range(circ_start, circ_stop+1):
            circ = step * circle

            x = r * math.cos(math.pi * circ) + cx
            y = r * math.sin(math.pi * circ) * math.cos(2*math.pi * rot) + cy
            z = r * math.sin(math.pi * circ) * math.sin(2*math.pi * rot) + cz

            points.append([x, y, z])
            #print 'rotation: %d\tcircle%d'%(rotation, circle)
    return points

def add_torus( edges, cx, cy, cz, r0, r1, step ):
    points = generate_torus(cx, cy, cz, r0, r1, step)
    num_steps = int(1/step+0.1)

    lat_start = 0
    lat_stop = num_steps
    longt_start = 0
    longt_stop = num_steps

    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):

            p0 = lat * (num_steps) + longt;
            if (longt == num_steps - 1):
	        p1 = p0 - longt;
            else:
	        p1 = p0 + 1;
            p2 = (p1 + num_steps) % (num_steps * num_steps);
            p3 = (p0 + num_steps) % (num_steps * num_steps);

            add_polygon(edges,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p3][0],
                        points[p3][1],
                        points[p3][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2] )
            add_polygon(edges,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2],
                        points[p1][0],
                        points[p1][1],
                        points[p1][2] )

def generate_torus( cx, cy, cz, r0, r1, step ):
    points = []
    num_steps = int(1/step+0.1)

    rot_start = 0
    rot_stop = num_steps
    circ_start = 0
    circ_stop = num_steps

    for rotation in range(rot_start, rot_stop):
        rot = step * rotation
        for circle in range(circ_start, circ_stop):
            circ = step * circle

            x = math.cos(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cx;
            y = r0 * math.sin(2*math.pi * circ) + cy;
            z = -1*math.sin(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cz;

            points.append([x, y, z])
    return points

def add_circle( points, cx, cy, cz, r, step ):
    x0 = r + cx
    y0 = cy
    t = step

    while t <= 1.00001:
        x1 = r * math.cos(2*math.pi * t) + cx;
        y1 = r * math.sin(2*math.pi * t) + cy;

        add_edge(points, x0, y0, cz, x1, y1, cz)
        x0 = x1
        y0 = y1
        t+= step

def add_curve( points, x0, y0, x1, y1, x2, y2, x3, y3, step, curve_type ):

    xcoefs = generate_curve_coefs(x0, x1, x2, x3, curve_type)[0]
    ycoefs = generate_curve_coefs(y0, y1, y2, y3, curve_type)[0]

    t = step
    while t <= 1.00001:
        x = xcoefs[0] * t*t*t + xcoefs[1] * t*t + xcoefs[2] * t + xcoefs[3]
        y = ycoefs[0] * t*t*t + ycoefs[1] * t*t + ycoefs[2] * t + ycoefs[3]

        add_edge(points, x0, y0, 0, x, y, 0)
        x0 = x
        y0 = y
        t+= step

def draw_lines( matrix, screen, zbuffer, color ):
    if len(matrix) < 2:
        print 'Need at least 2 points to draw'
        return

    point = 0
    while point < len(matrix) - 1:
        draw_line( int(matrix[point][0]),
                   int(matrix[point][1]),
                   matrix[point][2],
                   int(matrix[point+1][0]),
                   int(matrix[point+1][1]),
                   matrix[point+1][2],
                   screen, zbuffer, color)
        point+= 2

def add_edge( matrix, x0, y0, z0, x1, y1, z1 ):
    add_point(matrix, x0, y0, z0)
    add_point(matrix, x1, y1, z1)

def add_point( matrix, x, y, z=0 ):
    matrix.append( [x, y, z, 1] )

def scanline_helper(x0, y0, z0, x1, y1, z1):
    #swap points if going right -> left
    if x0 > x1:
        xt = x0
        yt = y0
        zt = z0
        x0 = x1
        y0 = y1
        z0 = z1
        x1 = xt
        y1 = yt
        z1 = zt

    x = x0
    y = y0
    z = z0
    A = 2 * (y1 - y0)
    B = -2 * (x1 - x0)
    wide = False
    tall = False

    if ( abs(x1-x0) >= abs(y1 - y0) ): #octants 1/8
        wide = True
        loop_start = x
        loop_end = x1
        dx_east = dx_northeast = 1
        dy_east = 0
        d_east = A
        distance = x1 - x
        if ( A > 0 ): #octant 1
            d = A + B/2
            dy_northeast = 1
            d_northeast = A + B
        else: #octant 8
            d = A - B/2
            dy_northeast = -1
            d_northeast = A - B

    else: #octants 2/7
        tall = True
        dx_east = 0
        dx_northeast = 1
        distance = abs(y1 - y)
        if ( A > 0 ): #octant 2
            d = A/2 + B
            dy_east = dy_northeast = 1
            d_northeast = A + B
            d_east = B
            loop_start = y
            loop_end = y1
        else: #octant 7
            d = A/2 - B
            dy_east = dy_northeast = -1
            d_northeast = A - B
            d_east = -1 * B
            loop_start = y1
            loop_end = y

    delta_z = float(z1 - z)/distance if distance != 0 else 0
    old_y = y
    l = [ [x, z] ]
    while ( loop_start < loop_end ):
        if y != old_y:
            l.append([x, z])
            old_y = y
        if ( (wide and ((A > 0 and d > 0) or (A < 0 and d < 0))) or
             (tall and ((A > 0 and d < 0) or (A < 0 and d > 0 )))):
            x+= dx_northeast
            y+= dy_northeast
            d+= d_northeast
        else:
            x+= dx_east
            y+= dy_east
            d+= d_east
        z += delta_z
        loop_start+= 1

    if y != old_y:
        l.append([x, z])
    return l



def draw_line( x0, y0, z0, x1, y1, z1, screen, zbuffer, color, left = [], right = [], shading_type = "", constants = [], light_sources = [] ):

    #swap points if going right -> left
    if x0 > x1:
        xt = x0
        yt = y0
        zt = z0
        x0 = x1
        y0 = y1
        z0 = z1
        x1 = xt
        y1 = yt
        z1 = zt

    x = x0
    y = y0
    z = z0
    A = 2 * (y1 - y0)
    B = -2 * (x1 - x0)
    wide = False
    tall = False

    if ( abs(x1-x0) >= abs(y1 - y0) ): #octants 1/8
        wide = True
        loop_start = x
        loop_end = x1
        dx_east = dx_northeast = 1
        dy_east = 0
        d_east = A
        distance = x1 - x
        if ( A > 0 ): #octant 1
            d = A + B/2
            dy_northeast = 1
            d_northeast = A + B
        else: #octant 8
            d = A - B/2
            dy_northeast = -1
            d_northeast = A - B

    else: #octants 2/7
        tall = True
        dx_east = 0
        dx_northeast = 1
        distance = abs(y1 - y)
        if ( A > 0 ): #octant 2
            d = A/2 + B
            dy_east = dy_northeast = 1
            d_northeast = A + B
            d_east = B
            loop_start = y
            loop_end = y1
        else: #octant 7
            d = A/2 - B
            dy_east = dy_northeast = -1
            d_northeast = A - B
            d_east = -1 * B
            loop_start = y1
            loop_end = y

    delta_z = float(z1 - z)/distance if distance != 0 else 0

    if shading_type == "goroud":
        color = left[:]

    elif shading_type == "phong":
        # print A
        # print B
        color = calc_total_light(left, constants, light_sources)
        if distance != 0:
            d_N = [ float(right[0] - left[0])/(distance), float(right[1] - left[1])/(distance), float(right[2] - left[2])/(distance) ]
            N = left[:]

    while ( loop_start < loop_end ):
        plot( screen, zbuffer, color, x, y, z )

        if ( (wide and ((A > 0 and d > 0) or (A < 0 and d < 0))) or
             (tall and ((A > 0 and d < 0) or (A < 0 and d > 0 )))):
            x+= dx_northeast
            y+= dy_northeast
            d+= d_northeast
        else:
            x+= dx_east
            y+= dy_east
            d+= d_east

        if shading_type == "goroud":
            for i in range(len(color)):
                I = int(round((float(x - x0)/(x1 - x0)) * right[i] + (float(x1 - x)/(x1 - x0)) * left[i]))
                color[i] = I if I <= 255 else 255
        elif shading_type == "phong":
            for i in range(len(color)):
                N[i] += d_N[i]
            color = calc_total_light(N, constants, light_sources)
                

        # print color
        z += delta_z
        loop_start+= 1

    plot( screen, zbuffer, color, x1, y1, z1 )
