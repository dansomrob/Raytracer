README

The main simulation code is stored in the module 'raytracer.py' with the
visulisation code stored in 'visualisation.py'. The visualisation code
has a dependancy on the 'raytracer.py' module.

To run a simulation, first we need to create some optical elements.
These include: ray(s), lens(es) and an output plane. Lenses are centred
on the z axis (the direction of propogation), however rays can be created
with a displacement from the z axis.

Rays

raytracer.Ray(direction, position)
raytracer.Ray_Bundle(type, number, position, direction, width)

type = '1d' or '2d'
number = integer
position, direction = [a, b, c]
width = float (really the starting radius)

Rays can either be created as a singular ray, or as a bundle. A bundle
of rays can either be 1 dimensional (with a constant y value) or 2
dimensional. In the 2d case, the ray objects are stored in a
dictionary '.rayname[i]'.

A 2d ray bundle is made using hexagonal numbers, you can use the
initial_RMS method to find the starting RMS for a 2d bundle.

The Ray object has the method paraxial_focal_length which will
determine the focal length of a lens for a paraxial ray. For this to
work, y = 0 and x should be small and non-zero. It determines the
focal length by finding the intercept with the y axis. If calling
this method only create a single Ray object, not as part of a
Ray_Bundle.

Lenses

raytracer.SphericalRefraction(z, r, n1, n2, max, inverse)
raytracer.Plano_Convex(z, r1, r2, n1, n2, max, d, inverse)
raytracer.Bi_Convex(z, r1, r2, n1, n2, max, d)

z = float (z axis intercept of the lens)
r(1), r2 = float (radius of curvature of the lens)
n1, n2 = float, (n1 is the refractive index outside the lens,
	n2 is the refractive index in the lens)
max = float (the maxiumum radius of the lens
d = float (distance between the front and back lenses)
inverse = True/False (controls direction of lens)

The SphericalRefraction object creates a spherical lens, with
radius of curvature r. When r = 0, it creates a straight vertical
boundary. The Plano- and Bi-Convex lenses inherit their behaviour 
from SphericalRefraction.

When inverse is False, the Plano-Convex lens has its spherical side
facing the incoming light rays. When inverse is True, the planar
side faces the incoming light rays.

Output Plane

raytracer.OutputPlane(zlower, zupper, xlower, xupper, ylower, yupper) 

The six floats determine the size of the cuboid in which the simulation
occurs. Generally zlower is 0 and zupper around 250, all other values
should take -20 or 20.

After creating our optical elements we use the method 'propogate_ray'
to see how the ray evolves. This method can have both a single ray
and ray bundles as inputs. Generally we propogate through the lenses
first and finish with the optical element.

Propgation works by finding the distance to intercept for the ray,
propogating the ray to that position, and using snell's law to
determine the new direction of the ray.

Once a ray hits the OutputPlane, or total internal refraction occurs,
or is outside of the lens apature it is terminated. This is done by
setting the direction vector to [0, 0, 0].

Visualisation

All plots must begin with the command 

visualisation.create_display(outputplane, spot)

outputplane = the optical element 'OutputPlane' used in the simulation
spot = True/False (if False plots the xz graph of the rays propgating
	through the system, if True plots a spot diagram in xy)

If spot is False, then run the commands

visualisation.plot_element(opticalelement)
visualisation.plot_ray(ray)

which will plot the optical elements and the ray(s) in the xz plane.
These commands will works for all optical elements and all rays/ray
bundles.

There are additional commands

visualisation.plot_surface_normals(ray)
visualisation.plot_focal_length(in)

The first of which will plot the surface normals for all refracting
surfaces that the ray encounters.

plot_focal_lengths(in) takes as its input a float, this was done as
there are many ways of estimating the focal length.

If spot is True, run the command

visualisation.plot_spot(ray, z)

Where z is the location of the xy plane. Due to the way this function
works, z must be greater then an optical element.

The function visualisation.plot_spot_in(ray) can also be called which
plots the Ray_Bundle at z = 0 before any propogation. This only works
for 2d ray bundles.

Finally the command visualisation.show("name", "graph_title") will
save the plot as name.png and add 'graph_title' as a title to the
plot. 