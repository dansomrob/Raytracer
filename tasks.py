import raytracer
import visualisation as vis

from copy import deepcopy

# Task 9

T9_ray = raytracer.Ray_Bundle('2d', 5, [0, 0, 0], [0, 0, 1], 10)
T9_lens = raytracer.SphericalRefraction(100, 1/0.03, 1, 1.5, 11, False)
T9_plane = raytracer.OutputPlane(0, 250, -20, 20, -20, 20)

T9_lens.propagate_ray(T9_ray)
T9_plane.propagate_ray(T9_ray)

vis.create_display(T9_plane, False)

vis.plot_element(T9_lens)
vis.plot_element(T9_plane)
vis.plot_ray(T9_ray)

vis.plot_surface_normals(T9_ray)
vis.plot_focal_length(199.99980003015014)

vis.show("T9", "paraxial focal length: f = 99.99mm")

# Task 10

T10_ray = raytracer.Ray([0, 0, 1], [-0.1, 0, 0])
T10_lens = T9_lens
T10_plane = T9_plane

T10_lens.propagate_ray(T10_ray)
T10_plane.propagate_ray(T10_ray)

T10_focal_length = T10_ray.paraxial_focal_length()

print(f" the focal length on the lens is {T10_focal_length - T10_lens.z}")

# Task 11

# Task 12

T11_ray = raytracer.Ray_Bundle('2d', 5, [0, 0, 0], [0, 0, 1], 2.5)
T11_lens = T9_lens
T11_plane = T9_plane

T11_lens.propagate_ray(T11_ray)
T11_plane.propagate_ray(T11_ray)

vis.create_display(T11_plane, False)

vis.plot_ray(T11_ray)
vis.plot_element(T11_lens)
vis.plot_element(T11_plane)

vis.plot_surface_normals(T11_ray)
vis.plot_focal_length(T10_focal_length)

vis.show("T12", "f = 99.99mm")

# Task 13

vis.create_display(T11_plane, True)

vis.plot_spot(T11_ray, T10_focal_length)
vis.show("T13", "f = 99.99mm")

T13_RMS = T11_ray.RMS(T10_focal_length)
print(f"Task 13, the RMS value is {T13_RMS}")
print(f"inital RMS = {T11_ray.inital_RMS()}")

# hence the focal size is ...

# Task 15

T15_lens_convexfacing = raytracer.Plano_Convex(50, 1/0.02, 1, 1.5168, 15, 5, False)
T15_lens_planarfacing = raytracer.Plano_Convex(50, 1/0.02, 1, 1.5168, 15, 5, True)
T15_plane = T9_plane

T15_ray_10mm = raytracer.Ray_Bundle('2d', 5, [0, 0, 0], [0, 0, 1], 10)
T15_ray_10mm_1 = deepcopy(T15_ray_10mm)

T15_rays_initialRMS = T15_ray_10mm.inital_RMS()

# Finding the paraxial focus of the lens

T15_ray_focal = raytracer.Ray([0, 0, 1], [0.1, 0, 0])
T15_ray_focal_1 = deepcopy(T15_ray_focal)

T15_lens_convexfacing.propagate_ray(T15_ray_focal)
T15_plane.propagate_ray(T15_ray_focal)

T15_lens_planarfacing.propagate_ray(T15_ray_focal_1)
T15_plane.propagate_ray(T15_ray_focal_1)

T15_focal_length_convexfacing = T15_ray_focal.paraxial_focal_length()
T15_focal_length_planarfacing = T15_ray_focal_1.paraxial_focal_length()

print(f"The focal length of the convex facing lens is {T15_focal_length_convexfacing} and the focal length of the planar facing lens is {T15_focal_length_planarfacing}")

# 10mm width beam, convex facing input
T15_lens_convexfacing.propagate_ray(T15_ray_10mm)
T15_plane.propagate_ray(T15_ray_10mm)

vis.create_display(T15_plane, False)

vis.plot_focal_length(T15_focal_length_convexfacing)

vis.plot_element(T15_lens_convexfacing)
vis.plot_element(T15_plane)
vis.plot_ray(T15_ray_10mm)

vis.show("T15_convexfacing_10mm", f"f = {(T15_focal_length_convexfacing - T15_lens_convexfacing.z):.2f} mm")

# 10mm width beam, planar facing input
T15_lens_planarfacing.propagate_ray(T15_ray_10mm_1)
T15_plane.propagate_ray(T15_ray_10mm_1)

vis.create_display(T15_plane, False)

vis.plot_focal_length(T15_focal_length_planarfacing)

vis.plot_element(T15_lens_planarfacing)
vis.plot_element(T15_plane)
vis.plot_ray(T15_ray_10mm_1)

vis.show("T15_planarfacing_10mm", f"f = {(T15_focal_length_planarfacing - T15_lens_planarfacing.z):.1f} mm")

# Investigating RMS with beam width, convex facing input
T15_rays_convexfacing = dict()
RMS_values_convexfacing = list()
beam_widths_convexfacing = list()
for i in range(10):
    T15_rays_convexfacing[i] = raytracer.Ray_Bundle('2d', 5, [0, 0, 0], [0, 0, 1], i + 1)
    T15_lens_convexfacing.propagate_ray(T15_rays_convexfacing[i])
    T15_plane.propagate_ray(T15_rays_convexfacing[i])
    RMS = T15_rays_convexfacing[i].RMS(T15_focal_length_convexfacing)
    RMS_values_convexfacing.append(RMS)
    beam_widths_convexfacing.append(T15_rays_convexfacing[i].width)

# Investigating RMS with beam width, planar facing input
T15_rays_planarfacing = dict()
RMS_values_planarfacing = list()
beam_widths_planarfacing = list()
for i in range(10):
    T15_rays_planarfacing[i] = raytracer.Ray_Bundle('2d', 5, [0, 0, 0], [0, 0, 1], i + 1)
    T15_lens_planarfacing.propagate_ray(T15_rays_planarfacing[i])
    T15_plane.propagate_ray(T15_rays_planarfacing[i])
    RMS = T15_rays_planarfacing[i].RMS(T15_focal_length_planarfacing)
    RMS_values_planarfacing.append(RMS)
    beam_widths_planarfacing.append(T15_rays_planarfacing[i].width)

# Plotting the variation of RMS with beam width
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1)
plt.plot(beam_widths_convexfacing, RMS_values_convexfacing, color = '#006666', label = 'Convex Facing')
plt.plot(beam_widths_planarfacing, RMS_values_planarfacing, color = '#ef2f2f', label = 'Planar Facing')
plt.legend(loc = 'best')
plt.xlabel("Beam Radius /mm")
plt.ylabel("RMS", rotation = 0)

plt.savefig("T15_RMS", dpi = 300)

print(f"Task 15, The initial RMS was {T15_rays_initialRMS}")
print(f"convex rms = {RMS_values_convexfacing[-1]} planar rms = {RMS_values_planarfacing}")