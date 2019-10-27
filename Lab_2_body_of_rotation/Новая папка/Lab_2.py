import pylab
import frechet
import matplotlib.pyplot as plt
# noinspection PyInterpreter
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


class TorPoint:
    def __init__(self, z_, r1_, r2_):
        self.z = z_
        self.r1 = r1_
        self.r2 = r2_


class Point:
    def __init__(self, x_, y_, z_):
        self.x = x_
        self.y = y_
        self.z = z_

class Tor:
    def __init__(self, z_, r_rotate, r_circle):
        self.z = z_
        self.R = r_rotate
        self.r = r_circle
        self.points = []

    def generate_points(self, num_points):
        self.points = []
        n = num_points
        step = np.pi / n
        fi = -np.pi / 2.0

        for i in range(0, n):
            z = self.r * np.sin(fi)
            dr = self.r * np.cos(fi)
            point = TorPoint(z, self.R - dr, self.R + dr)
            self.points.append(point)
            fi += step


    def get_tor_points(self):
        return self.points


def points_on_plane(tor_point, x_plane):
    z = tor_point.z
    r1 = tor_point.r1
    r2 = tor_point.r2

    point1 = Point(x_plane, np.sqrt(r1**2 - x_plane**2), z)
    point2 = Point(x_plane, np.sqrt(r2**2 - x_plane**2), z)

    point3 = Point(x_plane, -np.sqrt(r1 ** 2 - x_plane ** 2), z)
    point4 = Point(x_plane, -np.sqrt(r2 ** 2 - x_plane ** 2), z)
    return point1, point2, point3, point4


def get_y_array(points):
    data = []
    for point in points:
        data.append(point.y)
    return data

def get_z_array(points):
    data = []
    for point in points:
        data.append(point.z)
    return data


def get_intersection_tor_plane(tor, x_plane):
    tor_points = tor.get_tor_points()

    result = []
    left1 = []
    left2 = []
    left3 = []
    left4 = []

    left_jump_2_3_flag = False
    left_jump_1_4_flag = False
    right_jump_2_3_flag = False
    right_jump_1_4_flag = False

    right1 = []
    right2 = []
    right3 = []
    right4 = []
    for point in tor_points:
        tmp = points_on_plane(point, x_plane)

        if not np.isnan(tmp[0].y):
            if not right_jump_1_4_flag:
                right1.append(tmp[0])
            else:
                right4.append(tmp[0])
        else:
            right_jump_1_4_flag = True

        if not np.isnan(tmp[1].y):
            if not right_jump_2_3_flag:
                right2.append(tmp[1])
            else:
                right3.append(tmp[1])
        else:
            right_jump_2_3_flag = True

        if not np.isnan(tmp[2].y):
            if not left_jump_1_4_flag:
                left1.append(tmp[2])
            else:
                left4.append(tmp[2])
        else:
            left_jump_1_4_flag = True

        if not np.isnan(tmp[3].y):
            if not left_jump_2_3_flag:
                left2.append(tmp[3])
            else:
                left3.append(tmp[3])
        else:
            left_jump_2_3_flag = True

    right1.reverse()
    right4.reverse()
    left1.reverse()
    left4.reverse()

    left = left2 + left3 + left4 + left1
    right = right2 + right3 + right4 + right1
    right.reverse()

    return left, right


def print_points(points):
    print()
    for point in points:
        print(point.x, point.y, point.z)


def plot_3D(plane_cut, name, filename):

    data = []

    fig = pylab.figure()
    axes = Axes3D(fig)

    i = 0
    number = len(plane_cut)
    cmap =plt.get_cmap('gnuplot')
    colors = [cmap(i) for i in np.linspace(0, 1, number)]
    for points_arr in plane_cut:
        for points in points_arr:
            x = []
            y = []
            z = []
            x1 = []
            y1 = []
            z1 = []
            for elem in points:
                x.append([elem.x])
                y.append([elem.y])
                z.append([elem.z])
                x1.append(elem.x)
                y1.append(elem.y)
                z1.append(elem.z)

            x = np.matrix(x)
            y = np.matrix(y)
            z = np.matrix(z)
            axes.plot(x1, y1, z1, ".", color=colors[i])
            #axes.plot_wireframe(x, y, z, color=colors[i], linewidth=5, linestyle='-')

        i += 1

    fig.savefig(filename, dpi=300, format='png', bbox_inches='tight')
    fig.show()
    plt.close(fig)


def draw_cut(points_arr, name, filename):
    plt.axis('scaled')
    fig, ax = plt.subplots(nrows=1, ncols=1, sharey=True)
    ax.set_xlabel("y")
    ax.set_ylabel("z")
    ax.set_title(name)

    for points in points_arr:
        x = []
        y = []
        z = []
        for elem in points:
            x.append(elem.x)
            y.append(elem.y)
            z.append(elem.z)
        ax.plot(y, z, ".", color="orange")
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])
    #fig.legend()
    fig.savefig(filename, dpi=300, format='png', bbox_inches='tight')

    fig.show()
    plt.close(fig)

def plane_research(tor, planes):
    data = []
    for plane in planes:
        left, right = get_intersection_tor_plane(tor, plane)
        draw_cut([left, right], "H = %i" % plane, "Tor_section(H = %i).png" % plane)
        data.append([left, right])

    plot_3D(data, "H = %i" % len(planes), "Tor_cuts(H = %i).png" % len(planes))


def lemniscate(c):
    result = []
    n = 360
    fi = 0
    step = 2 * np.pi / n
    left = []
    right = []
    for i in range(0, n//4):
        t = np.tan(fi)
        x = c * (t + t ** 3) / (1 + t ** 4)
        y = c * (t - t ** 3) / (1 + t ** 4)
        right.append([x, y])
        fi += step

    for i in range(n//4, n // 2):
        t = np.tan(fi)
        x = c * (t + t ** 3) / (1 + t ** 4)
        y = c * (t - t ** 3) / (1 + t ** 4)
        left.append([x, y])
        fi += step
    '''
    names = ["left", "right"]
    fig, ax = plt.subplots(nrows=1, ncols=1, sharey=True)
    data = [left, right]
    i = 0
    for points in data:
        x = []
        y = []
        for elem in points:
            x.append(elem[0])
            y.append(elem[1])
        ax.plot(x, y, label=names[i])
        i+=1
    fig.legend()
    fig.show()
    '''
    return left, right



def draw_line_on_plot(line, ax, color):
    data = [line]
    for i in range(0, len(data)):
        x = []
        y = []
        for elem in data[i]:
            x.append(elem[0])
            y.append(elem[1])
        ax.plot(x, y, ":", label="lemniscate", color=color)


def process_lemniscate(tor, plane, c):
    name = "lemniscate"
    filename = "lemniscate.png"


    left, right = get_intersection_tor_plane(tor, plane)
    points_arr = [left, right]
    l, r = lemniscate(c)


    plt.axis('scaled')
    fig, ax = plt.subplots(nrows=1, ncols=1, sharey=True)
    ax.set_xlabel("y")
    ax.set_ylabel("z")
    ax.set_title(name + "\n")

    data = []
    data_arr = []
    draw_flag = True
    for points in points_arr:
        x = []
        y = []
        z = []
        tmp = []
        for elem in points:
            x.append(elem.x)
            y.append(elem.y)
            z.append(elem.z)
            data.append([elem.y, elem.z])
            tmp.append([elem.y, elem.z])
        data_arr.append(tmp)
        if(draw_flag):
            #ax.plot(y, z, "-", color="red", label="dist")
            ax.plot(y, z, ".", color="orange", label="section")
            draw_flag = False
        else:
            ax.plot(y, z, ".", color="orange")

    draw_line_on_plot(l + r, ax, "b")

    temp_lemn = [r, l]
    temp_data = [data_arr[0], data_arr[1]]
    dist1, ind1 = frechet.d_Frechet(temp_lemn[0], temp_data[1])
    dist2, ind2 = frechet.d_Frechet(temp_lemn[0], temp_data[1])

    if dist1 > dist2 + 100:
        if ind1[0] > -1 and ind1[1] > -1:
            ax.plot([temp_lemn[0][ind1[0]][0], temp_data[0][ind1[1]][0]], [temp_lemn[0][ind1[0]][1], temp_data[0][ind1[1]][1]], label="dist")
        print("frechet distance = ", dist1)
    else:
        if ind2[0] > -1 and ind2[1] > -1:
            ax.plot([temp_lemn[1][ind2[0]][0], temp_data[1][ind2[1]][0]], [temp_lemn[1][ind2[0]][1], temp_data[1][ind2[1]][1]], label="dist")
        print("frechet distance = ", dist1)

    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    leg = fig.legend()
    fig.legend()
    fig.savefig(filename, dpi=300, format='png', bbox_inches='tight')

    fig.show()
    plt.close(fig)






def main():
    R = 20
    r1 = 10
    numpoints = 360

    tor = Tor(0, R, r1)
    tor.generate_points(100)
    H = np.array([i for i in range(0, R + r1)])
    #plane_research(tor, H)
    H = [0, 5, 10, 15, 20, 30]
    #plane_research(tor, H)
    process_lemniscate(tor, R - r1, R + r1 - r1 / 6)

    print("plane = ", R - r1, "c = ", 28)



if __name__ == "__main__":
    main()
