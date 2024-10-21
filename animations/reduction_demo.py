from manim import *
import numpy as np

class ThreeDData(ThreeDScene):
    def construct(self):
        # Set up 3D axes
        axes = ThreeDAxes()
        self.set_camera_orientation(phi=70 * DEGREES, theta=45 * DEGREES)
        self.begin_ambient_camera_rotation(rate=0.1)  # Start rotation
        self.add(axes)

        # Define five cluster centers
        cluster_centers = [
            np.array([2, 2, 2]),
            np.array([-2, -2, 2]),
            np.array([-2, 2, -2]),
            np.array([2, -2, -2]),
            np.array([0, 0, 0]),
        ]

        # Define colors for clusters
        cluster_colors = [RED, GREEN, BLUE, YELLOW, PURPLE]

        # Create a group to hold all points
        points = VGroup()

        # Set a random seed for reproducibility
        np.random.seed(42)

        # Generate points for each cluster
        for center, color in zip(cluster_centers, cluster_colors):
            for _ in range(20):  # 20 points per cluster
                # Generate a point around the cluster center
                point = center + np.random.normal(scale=0.5, size=3)
                # Create a 3D dot at the generated point
                dot = Dot3D(point=point, radius=0.05, color=color)
                points.add(dot)

        # Animate the points appearing
        self.play(LaggedStart(*[FadeIn(point, scale=0.5) for point in points], lag_ratio=0.01), run_time=2)
        self.wait(1)
        self.next_section()

        # Define the plane
        # Plane normal vector
        n = np.array([1, 1, 1])
        n = n / np.linalg.norm(n)  # Normalize n

        # Point on the plane
        P0 = np.array([0, 0, 0])

        # Create two vectors perpendicular to n
        if not np.allclose(n, [1, 0, 0]):
            v1 = np.cross(n, [1, 0, 0])
        else:
            v1 = np.cross(n, [0, 1, 0])
        v1 = v1 / np.linalg.norm(v1)
        v2 = np.cross(n, v1)
        v2 = v2 / np.linalg.norm(v2)

        # Define the size of the plane
        plane_size = 5

        # Calculate the four corners of the plane
        corners = [
            P0 - plane_size * v1 - plane_size * v2,
            P0 + plane_size * v1 - plane_size * v2,
            P0 + plane_size * v1 + plane_size * v2,
            P0 - plane_size * v1 + plane_size * v2,
        ]

        # Create the plane as a Polygon
        plane = Polygon(
            *corners,
            color=WHITE,
            fill_color=WHITE,
            fill_opacity=0.4,
            stroke_opacity=0,  # Hide edges
        )
        plane.set_shade_in_3d(True)

        # Animate the plane appearing
        self.play(Create(plane, scale=0.5), run_time=1)
        self.wait(1)
        self.next_section()

        # Project each point onto the plane and create lines
        projection_lines = VGroup()
        projected_points = VGroup()
        for dot in points:
            P = dot.get_center()
            # Compute the projection of P onto the plane
            projection = P - n * np.dot(P - P0, n)
            # Create a line from P to its projection
            line = Line3D(P, projection, color=YELLOW, stroke_width=0.5)
            projection_lines.add(line)
            # Create a dot at the projection
            proj_dot = Dot3D(point=projection, radius=0.05, color=dot.get_color())
            projected_points.add(proj_dot)

        # Animate the projection lines appearing
        self.play(
            LaggedStart(*[Create(line) for line in projection_lines], lag_ratio=0.01),
            run_time=3
        )
        self.wait(1)
        self.next_section()

        # Remove the original points and projection lines, keep the projected points
        self.play(
            FadeOut(projection_lines),
            ReplacementTransform(points, projected_points),
            run_time=2
        )
        self.wait(1)

        # Rotate the plane and projected points onto the xy-plane
        # Compute the rotation matrix to align n with the z-axis
        # We need a rotation that maps n to [0, 0, 1]
        # The axis of rotation is the cross product of n and [0, 0, 1]
        k = np.array([0, 0, 1])
        v = np.cross(n, k)
        s = np.linalg.norm(v)
        c = np.dot(n, k)
        if s == 0:
            # n is already aligned with k
            R = np.identity(3)
        else:
            vx = np.array([
                [0, -v[2], v[1]],
                [v[2], 0, -v[0]],
                [-v[1], v[0], 0]
            ])
            R = np.identity(3) + vx + np.matmul(vx, vx) * ((1 - c) / (s ** 2))

        # Create a rotation function for Manim
        def rotate_point(p):
            return np.dot(R, p - P0) + P0
        
        self.next_section()

        # Apply rotation to the plane and projected points, and remove plane
        self.play(
            plane.animate.apply_function(rotate_point),
            projected_points.animate.apply_function(rotate_point),
            FadeOut(plane),
            run_time=3
        )
        self.wait(1)
        self.stop_ambient_camera_rotation()

        # Add 2D axes to the scene
        axes2d = Axes(
            x_range=[-5, 5, 1],
            y_range=[-5, 5, 1],
            tips=False
        )

        # Transform the projected points into 2D dots
        projected_points_2d = VGroup()
        for proj_dot in projected_points:
            proj_point = proj_dot.get_center()
            # Since the plane is now in the xy-plane, we can ignore the z-coordinate
            x, y, z = proj_point
            dot2d = Dot(point=axes2d.coords_to_point(x, y), color=proj_dot.get_color())
            projected_points_2d.add(dot2d)

        # Replace the 3D projected points with the 2D projected points
        # and the 3D axes with the 2D axes
        self.move_camera(phi=0 * DEGREES, theta=90 * DEGREES)
        self.play(ReplacementTransform(axes, axes2d), ReplacementTransform(projected_points, projected_points_2d), run_time=2)
        self.wait()