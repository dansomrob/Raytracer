import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import raytracer

def create_display(element, spot):
    if (type(element) is raytracer.OutputPlane):
        global fig, ax
        fig, ax = plt.subplots(1, 1)
        ax.set_ylabel("x /mm", rotation = 0)
        if (spot is False):       
            ax.set_yticks(np.linspace(element.xlower, element.xupper, 11))
            ax.set_ylim(top = element.xupper, bottom = element.xlower) 
            ax.set_xticks(np.linspace(element.zlower, element.zupper, 11))
            ax.set_xlim(left = element.zlower, right = element.zupper)
            ax.set_xlabel("z /mm")
        if (spot is True):
            ax.set_xlabel("y /mm")
    else:
        raise Exception("need an output plane")

def plot_element(element):
    if (type(element) is raytracer.SphericalRefraction):
        if (element.r == 0): #The optical element is a flat plane
            ax.vlines(x = element.z, ymin = -1 * element.max, ymax = element.max, colors = 'black', lw = 1.5)
        elif (element.r != 0): #the optical element is a spherical lens
            if (element.inverse is False):
                ax.add_patch(Arc((element.z + element.r, 0), width = 2 * element.r, height \
                    = 2 * element.r, angle = 180, theta1 = -np.rad2deg(np.arcsin(element.max / \
                    element.r)), theta2 = np.rad2deg(np.arcsin(element.max / element.r)), color = 'black', lw = 1.5))
            if (element.inverse is True):
                ax.add_patch(Arc((element.z - element.r, 0), width = 2 * element.r, height \
                    = 2 * element.r, angle = 0, theta1 = -np.rad2deg(np.arcsin(element.max / \
                    element.r)), theta2 = np.rad2deg(np.arcsin(element.max / element.r)), color = 'black', lw = 1.5))
    elif (type(element) is raytracer.OutputPlane):
        ax.vlines(x = element.zlower, ymin = element.xlower, ymax = element.xupper,\
            ls = 'dashed', color = 'gray')
        ax.vlines(x = element.zupper, ymin = element.xlower, ymax = element.xupper,\
            ls = 'dashed', color = 'gray')
        ax.hlines(y = element.xlower, xmin = element.zlower, xmax = element.zupper,\
            ls = 'dashed', color = 'gray')
        ax.hlines(y = element.xupper, xmin = element.zlower, xmax = element.zupper,\
            ls = 'dashed', color = 'gray')
    elif (type(element) is raytracer.Plano_Convex):
        plot_element(element.spherical)
        plot_element(element.planar)
        if (element.inverse is False):
            ax.hlines(y = element.max, xmin = element.planar.z - element.d +\
                element.spherical.r - np.sqrt(element.spherical.r**2 - element.max**2),\
                    xmax = element.planar.z, color = 'black', lw = 1.5)
            ax.hlines(y = -element.max, xmin = element.planar.z - element.d +\
                element.spherical.r - np.sqrt(element.spherical.r**2 - element.max**2),\
                    xmax = element.planar.z, color = 'black', lw = 1.5)
        if (element.inverse is True):
            ax.hlines(y = element.max, xmin = element.planar.z, xmax = element.z +\
                element.d - element.spherical.r + np.sqrt(element.spherical.r**2 -\
                    element.max**2), color = 'black', lw = 1.5)
            ax.hlines(y = -element.max, xmin = element.planar.z, xmax = element.z +\
                element.d - element.spherical.r + np.sqrt(element.spherical.r**2 -\
                    element.max**2), color = 'black', lw = 1.5)
    elif (type(element) is raytracer.Bi_Convex):
        plot_element(element.spherical1)
        plot_element(element.spherical2)
        ax.hlines(y = element.max, xmin = element.z + element.r1 - np.sqrt(element.r1**2\
            - element.max**2), xmax = element.z + element.d - element.r2 + np.sqrt(element.r2**2\
                - element.max**2), color = 'black', lw = 1.5)
        ax.hlines(y = -element.max, xmin = element.z + element.r1 - np.sqrt(element.r1**2\
            - element.max**2), xmax = element.z + element.d - element.r2 + np.sqrt(element.r2**2\
                - element.max**2), color = 'black', lw = 1.5)
    else:
        raise Exception("This optical element is not of form 'SphericalRefraction', \
            check input")

def plot_ray(ray):
    if (type(ray) is raytracer.Ray):
        plt.plot(ray.vertices()[2::3], ray.vertices()[::3], color = "#e2ab27")
    elif (type(ray) is raytracer.Ray_Bundle):
        if (ray.type == '2d'):
            number = raytracer.hex_numbers(ray.number)
        if (ray.type == '1d'):
            number = ray.number
        for i in range(number):
            plot_ray(ray.rayname[i])
    else:
        raise Exception("This object is not of form 'Ray' or 'Ray_Bundle, check input")

def show(name, title):
    fig.suptitle(f"{title}")
    fig.savefig(f"{name}.png", dpi = 300)
    #fig.show()

def plot_surface_normals(ray):
    if (type(ray) is raytracer.Ray):
        for i in range(int(len(ray.surfaces) / 3)):
            plt.quiver(ray.position[5 + 3 * i], ray.position[3 + 3 * i],\
                ray.surfaces[2 + 3 * i], ray.surfaces[3 * i], color = '#f6546a', width = 0.004)
    if (type(ray) is raytracer.Ray_Bundle):
        if (ray.type == '2d'):
            number = raytracer.hex_numbers(ray.number)
        else:
            number = ray.number
        for i in range(number):
            plot_surface_normals(ray.rayname[i])

def plot_focal_length(input):
    plt.axvline(x = input, ls = 'dashed', color = '#8a2be2')

def plot_spot(ray, z):
    if (type(ray) is raytracer.Ray):
        if (ray.p()[2] - z <= 0):
            pass
        else:
            if (ray.direction[-4] == 0):
                newx, newy = ray.position[-6:-4]
            else:
                newx = ray.p()[0] - ((ray.direction[-6] / ray.direction[-4]) * (ray.p()[2] - z))
                newy = ray.p()[1] - ((ray.direction[-5] / ray.direction[-4]) * (ray.p()[2] - z))
            plt.scatter(newy, newx, marker = 'x', color = '#e2ab27')
    if (type(ray) is raytracer.Ray_Bundle):
        if (ray.type == '2d'):
            number = raytracer.hex_numbers(ray.number)
        else:
            number = ray.number
        for i in range(number):
            plot_spot(ray.rayname[i], z)

def plot_spot_in(ray):
    if (type(ray) is raytracer.Ray):
        plt.scatter(ray.position[1], ray.position[0], marker = 'x', color = '#e2ab27')
    if (type(ray) is raytracer.Ray_Bundle):
        if (ray.type == '2d'):
            number = raytracer.hex_numbers(ray.number)
        else:
            number = ray.number
        for i in range(number):
            plot_spot_in(ray.rayname[i])