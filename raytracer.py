import numpy as np

def lensmaker(r1, r2, n, d):
    if (r2 != 0 and r1 != 0):
        focal_length = (n - 1) * ((1 / r1) - (1 / r2) + (((n - 1) * d) / (n * r1 * r2)))
        focal_length = 1 / focal_length
    elif (r2 == 0 and r1 != 0):
        focal_length = (r1 * n) / (n - 1)
    elif (r2 ==0 and r1 == 0):
        focal_length = np.nan
    return focal_length

hex_numbers = lambda n: (3 * n * n) - (3 * n) + 1
"""
Class for an optical ray. Position (p) and direction (k) are stored as 3-element 
numpy arrays in a list.
"""
class Ray:
    def __init__(self, starting_direction, \
        starting_position = np.array([0, 0, 0])):
        starting_direction = np.divide(starting_direction, np.linalg.norm(starting_direction))
        self.direction = list(starting_direction)
        self.position = list(starting_position)
        self.surfaces = list()

    def k(self):
        # this should (always) be a unit vector!
        return([self.direction[-3],self.direction[-2],self.direction[-1]])

    def p(self):
        return([self.position[-3],self.position[-2],self.position[-1]])

    def append(self, p, k):
        for i in p:
            self.position.append(i)
        for j in k:
            self.direction.append(j)

    def vertices(self):
        return self.position

    # This function only works on rays restricted to the xz plane. All y components should be zero.
    def paraxial_focal_length(self):
        deltax = self.p()[0] * self.direction[-4] / self.direction[-6]
        return self.p()[2] - deltax

"""
The base class 'OpticalElement', main purpose is 'propogate_ray' function which
propogates ray through system.
"""

class OpticalElement():
    def __init__(self, pos):
        self.position = pos
    def propagate_ray(self, ray):
        "propagate a ray through the optical element"
        raise NotImplementedError()
    def end_of_ray(self, ray):
        ray.append('End', 'End')


"""
An inherited class, 
"""

class SphericalRefraction(OpticalElement):
    """This function .intercept finds the closest intercept of the ray with optical \
    element. It assumes the optical element is spherical but can also deal with \
    strainght lines (constant x).
    z - the intercept of the sphere with the z axis
    r - 1/(radius of curvature)
    n1 - refractive index outside the lens
    n2 - refractive index inside the lens
    max - the apature radius 
    """
    def __init__(self, z, r, n1, n2, max, inverse):
        self.z = z
        self.r = r
        self.n1 = n1
        self.n2 = n2
        self.max = max
        self.inverse = inverse
        if (inverse is False):
            self.position = (0, 0, self.z + self.r)
        if (inverse is True):
            self.position = (0, 0, self.z - self.r)
        #self.lensmaker_focal_length = lensmaker(r, 0, n2, 0)

    def intercept(self, ray):
        vector_diff = np.subtract(ray.p(), self.position[-3:])
        if (self.r != 0):
            l1 = np.abs(-1 * np.dot(vector_diff, ray.k()) + np.sqrt((np.dot(vector_diff, \
                ray.k())**2 - (np.linalg.norm(vector_diff)**2 - self.r**2))))
            l2 = np.abs(-1 * np.dot(vector_diff, ray.k()) - np.sqrt((np.dot(vector_diff, \
                ray.k())**2 - (np.linalg.norm(vector_diff)**2 - self.r**2))))
            out = min(l1, l2)
        if (self.r == 0):
            delta_z = self.position[2] - ray.p()[2]
            delta_x = ((ray.k()[0] / ray.k()[2]) * (delta_z))
            delta_y = ((ray.k()[1] / ray.k()[2]) * (delta_z))
            l = np.sqrt(delta_z**2 + delta_x**2 + delta_y**2)
            out = l
        return out
    
    def Snells_Law(self, ray, distance): # ray and surface vectors should be unit vectors!
        inverse_ray_k = [-x for x in ray.k()]
        if (self.r != 0):
            vector_diff = np.subtract(self.position, ray.p())
            surface = -1 * np.divide((vector_diff - np.multiply(distance, ray.k())), self.r)
            inverse_surface = [-x for x in surface]

            surface = np.divide(surface, np.linalg.norm(surface))
            
            if (np.linalg.norm(surface) < 0.99 or np.linalg.norm(surface) > 1.01):
                raise Exception(f"surface vector not normalised {np.linalg.norm(surface)} {surface}")
            
            theta1 = np.arccos(np.dot(ray.k(), surface))
            theta2 = np.arcsin(self.n1 * np.sin(theta1) / self.n2)
            if (np.isnan(theta1)):
                theta1, theta2 = 0.0, 0.0
            #print(f"theta1 = {np.rad2deg(theta1)} and theta2 = {np.rad2deg(theta2)}")

            """
            To perform the mapping from k1 to k2 we will rotate vector n by angle theta2.

            To do this we use the Euler-Rodrigues formula
            (https://en.wikipedia.org/wiki/Euler-Rodrigues_formula).

            First we need the surface normal vector - found by the cross product of the ray
            vector and the surface normal vector. We need it to be a unit vector.

            """
            k = np.cross(surface, ray.k())
            if (np.linalg.norm(k) == 0.):
                k[1] = -1
            else:
                k = np.divide(k, np.linalg.norm(k)) 

            if (self.inverse is True):
                out = np.multiply(surface, np.cos(theta2)) + np.multiply(np.cross(k,\
                     surface), np.sin(theta2)) + (k * np.dot(k, surface) * (1 - \
                        np.cos(theta2)))
            else:
                out = np.multiply(inverse_surface, np.cos(-theta2)) + np.multiply\
                    (np.cross(k, inverse_surface), np.sin(-theta2)) + (k * np.dot\
                        (k, inverse_surface) * (1 - np.cos(-theta2)))

        if (self.r == 0):
            vector_diff = np.subtract(self.position, ray.p())
            surface = [0, 0, -1]
            inverse_surface = [0, 0, 1]

            theta1 = np.arccos(np.dot(inverse_ray_k, surface))
            theta2 = np.arcsin(self.n1 * np.sin(theta1) / self.n2)
            if (np.isnan(theta1)):
                theta1, theta2 = 0.0, 0.0
            #print(f"theta1 = {np.rad2deg(theta1)} and theta2 = {np.rad2deg(theta2)}")

            k = np.cross(inverse_surface, ray.k())
            if (np.linalg.norm(k) == 0.):
                k[1] = -1
            else:
                k = np.divide(k, np.linalg.norm(k))

            omega = np.multiply(np.sin(theta2 / 2), k)
            a = np.cos(theta2 / 2)
            out = inverse_surface + (2 * a * (np.cross(omega, inverse_surface)))\
                + (2 * np.cross(omega, (np.cross(omega, inverse_surface))))

            #k = np.cross(surface, ray.k())
            #k = k / np.linalg.norm(k)
            #out = np.multiply(inverse_surface, np.cos(-theta2)) + np.multiply(np.cross(k, inverse_surface), np.sin(-theta2)) + (k * np.dot(k, inverse_surface) * (1 - np.cos(-theta2)))

        deviation = self.n1 * np.sin(np.arccos(np.dot(inverse_ray_k, surface))) -\
            self.n2 * np.sin(np.arccos(np.dot(out, inverse_surface)))
        
        if (np.abs(deviation) >= 0.01):
            raise Exception(f"snells laws error: theta1 = {np.rad2deg(np.arccos(np.dot(inverse_ray_k, surface)))}\
                theta2 = {np.rad2deg(np.arccos(np.dot(out, inverse_surface)))} deviation = {deviation}")

        ray.surfaces.append(surface[0])
        ray.surfaces.append(surface[1])
        ray.surfaces.append(surface[2])

        if (np.sin(theta1) > (self.n2 / self.n1)):
            return [0, 0, 0]
        else:
            return out

    def propagate_ray(self, ray):
        if (type(ray) is Ray):
            if (ray.k() == [0, 0, 0]):
                pass
            else:
                # find intercept of ray with surface
                distance_to_intercept = self.intercept(ray)
                distance_from_centre = np.sqrt(np.add(ray.p(), np.multiply(distance_to_intercept,\
                    ray.k()))[0]**2 + np.add(ray.p(), np.multiply(distance_to_intercept, ray.k()))[1]**2)
                if (distance_from_centre >= self.max):
                    print(f"outside of lens apature - radius from centre = {distance_from_centre}")
                    ray.append(ray.p(), [0, 0, 0])
                else:
                    # get new position of the ray
                    new_position_of_ray = np.add(ray.p(), np.multiply(distance_to_intercept, ray.k()))
                    # find direction of refracted ray
                    new_direction_of_ray = self.Snells_Law(ray, distance_to_intercept)
                    # update with new position and direction
                    ray.append(new_position_of_ray, new_direction_of_ray)
                    #print(f"position = {new_position_of_ray} direction = {new_direction_of_ray}")
        if (type(ray) is Ray_Bundle):
            if (ray.type == '2d'):
                number = hex_numbers(ray.number)
            if (ray.type == '1d'):
                number = ray.number
            for i in range(number):
                self.propagate_ray(ray.rayname[i]) # a bit of cheeky recursion!


class OutputPlane(OpticalElement):
    def __init__(self, zlower, zupper, xlower, xupper, ylower, yupper):
        self.zlower = zlower 
        self.zupper = zupper 
        self.xlower = xlower 
        self.xupper = xupper
        self.ylower = ylower
        self.yupper = yupper

    def intercept(self, ray):
        zdifZ = self.zupper - ray.p()[2]
        xdifZ = ((ray.k()[0] / ray.k()[2]) * (zdifZ))
        ydifZ = ((ray.k()[1] / ray.k()[2]) * (zdifZ))
        lenZ = np.sqrt(xdifZ**2 + zdifZ**2 + ydifZ**2)

        if (ray.k()[0] == 0):
            return lenZ

        if (ray.k()[0] < 0):
            xdifX = ray.p()[0] - self.xlower
        else:
            xdifX = self.xupper - ray.p()[0]
        zdifX = ((ray.k()[2] / ray.k()[0]) * (xdifX))
        ydifX = ((ray.k()[1] / ray.k()[0]) * (xdifX))
        lenX = np.sqrt(xdifX**2 + zdifX**2 + ydifX**2)

        if (ray.k()[1] < 0):
            ydifY = ray.p()[1] - self.ylower
        else:
            ydifY = self.yupper - ray.p()[1]
        zdifY = ((ray.k()[2] / ray.k()[1]) * (ydifY))
        xdifY = ((ray.k()[0] / ray.k()[1]) * (ydifY))
        lenY = np.sqrt(xdifY**2 + zdifY**2 + ydifY**2)
        return min(lenX, lenY, lenZ)

    def propagate_ray(self, ray):
        if (type(ray) is Ray):
            if (ray.k() == [0, 0, 0]):
                pass
            else:
                # find intercept of ray with surface
                distance_to_intercept = self.intercept(ray)
                #get new position of ray
                new_position_of_ray = np.add(ray.p(), np.multiply(distance_to_intercept, ray.k()))
                # update with new position
                ray.append(new_position_of_ray, [0, 0, 0])
            #print(ray.direction, ray.vertices())
        elif (type(ray) is Ray_Bundle):
            if (ray.type == '2d'):
                number = hex_numbers(ray.number)
            if (ray.type == '1d'):
                number = ray.number
            for i in range(number):
                self.propagate_ray(ray.rayname[i]) # a bit of cheeky recursion!

class Plano_Convex(SphericalRefraction):
    """
    d is the distance between the two z axis intercepts
    """
    def __init__(self, z, r, n1, n2, max, d, inverse):
        self.max = max
        self.d = d
        self.z = z
        self.inverse = inverse
        self.lensmaker_focal_length = lensmaker(r, 0, n2, d) / n2

        if (inverse is False):
            self.spherical = SphericalRefraction(z, r, n1, n2, max, False)
            self.planar = SphericalRefraction(z + d, 0, n2, n1, max, False)
        elif (inverse is True):
            self.planar = SphericalRefraction(z, 0, n1, n2, max, False)       
            self.spherical = SphericalRefraction(z + d, r, n2, n1, max, True)

        if (self.spherical.r - self.d - np.sqrt(self.spherical.r**2 - self.max**2) > 0):
            raise Exception(f"increase d - this value is not possible")

    def propagate_ray(self, ray):
        if (type(ray) is Ray):
            if (ray.k() == [0, 0, 0]):
                pass
            else:
                if (self.inverse is False):
                    self.spherical.propagate_ray(ray)
                    self.planar.propagate_ray(ray)
                elif (self.inverse is True):
                    self.planar.propagate_ray(ray)
                    self.spherical.propagate_ray(ray)
        elif (type(ray) is Ray_Bundle):
            if (ray.type == '2d'):
                number = hex_numbers(ray.number)
            if (ray.type == '1d'):
                number = ray.number
            for i in range(number):
                self.propagate_ray(ray.rayname[i])
        

class Bi_Convex(SphericalRefraction):
    """
    r1 and r2 are the respective curavtures of the two spherical lenses
    d is the distance between the two z axis intercepts
    """
    def __init__(self, z, r1, r2, n1, n2, max, d):
        self.d = d
        self.z = z
        self.r1 = r1
        self.r2 = r2
        self.max = max
        self.spherical1 = SphericalRefraction(z, r1, n1, n2, max, False)
        self.spherical2 = SphericalRefraction(z + d, r2, n2, n1, max, True)
        self.lensmaker_focal_length = lensmaker(r1, r2, n2, d)

    def propagate_ray(self, ray):
        if (type(ray) is Ray):
            if (ray.k() == [0, 0, 0]):
                pass
            else:
                self.spherical1.propagate_ray(ray)
                self.spherical2.propagate_ray(ray)
        elif (type(ray) is Ray_Bundle):
            if (ray.type == '2d'):
                number = hex_numbers(ray.number)
            if (ray.type == '1d'):
                number = ray.number
            for i in range(number):
                self.propagate_ray(ray.rayname[i])

class Ray_Bundle(Ray):
    def __init__(self, type, number, central_position, direction, width):
        direction = np.divide(direction, np.linalg.norm(direction))
        self.number =  number
        self.rayname = dict()
        self.type = type
        self.pos_dif = width / (number - 1)
        self.width = width

        if (type == '1d'):
            rotation_matrix = [[0, 0, -1], [0, 1, 0], [1, 0, 0]]
            new_vec = np.dot(rotation_matrix, direction)

            for i in range(self.number):
                self.rayname[i] = Ray(direction, [central_position[0] + (width / 2),\
                    central_position[1], central_position[2]] + np.multiply(self.pos_dif\
                        * i, new_vec))

        elif (type == '2d'):
            pos = central_position
            move1 = [[1, 0, central_position[0]], [0, 1, central_position[1]], [0, 0, 1]]
            move2 = [[1, 0, -central_position[0]], [0, 1, -central_position[1]], [0, 0, 1]] 
            radius = [0, self.pos_dif, 0]
            i = 1
            for n in range(hex_numbers(number)):
                self.rayname[n] = Ray(direction, pos)
                if ((n + 1) % hex_numbers(i) == 0):
                    i += 1
                    if (n == 0):
                        pos = np.add(pos, radius)
                        continue
                    else:
                        theta = (2 * np.pi) / (6 * (i - 2))
                        rotation_matrix = [[np.cos(theta), -np.sin(theta), 0],\
                            [np.sin(theta), np.cos(theta), 0], [0, 0, 1]]
                        pos = np.matmul(move1, np.matmul(rotation_matrix, np.matmul(move2,\
                            [pos[0], pos[1], 1])))
                        pos = [pos[0], pos[1], 0]
                        pos = np.add(pos, radius)
                else:
                    theta = (2 * np.pi) / (6 * (i - 1))
                    rotation_matrix = [[np.cos(theta), -np.sin(theta), 0], [np.sin(theta),\
                        np.cos(theta), 0], [0, 0, 1]]
                    pos = np.matmul(move1, np.matmul(rotation_matrix, np.matmul(move2,\
                        [pos[0], pos[1], 1])))
                    pos = [pos[0], pos[1], 0]
        else:
            raise Exception("please specify if ray bundle is either 1d or 2d")
    
    def RMS(self, z):
        distance = list()
        for i in range(hex_numbers(self.number)):
            newx = self.rayname[i].p()[0] - ((self.rayname[i].direction[-6] / \
                self.rayname[i].direction[-4]) * (self.rayname[i].p()[2] - z))
            newy = self.rayname[i].p()[1] - ((self.rayname[i].direction[-5] / \
                self.rayname[i].direction[-4]) * (self.rayname[i].p()[2] - z))
            distance.append(np.linalg.norm(np.subtract([newx, newy], [0, 0]))**2)
        MS = np.average(distance)
        RMS = np.sqrt(MS)
        return RMS

    def inital_RMS(self):
        distance = list()
        for i in range(hex_numbers(self.number)):
            x = self.rayname[i].position[0]
            y = self.rayname[i].position[1]
            distance.append(np.linalg.norm([x, y])**2)
        MS = np.average(distance)
        RMS = np.sqrt(MS)
        return RMS
