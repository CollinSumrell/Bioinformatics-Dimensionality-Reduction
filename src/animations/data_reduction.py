"""
This is all mostly the same as the demo from the presentation. But instead of randomly generating
data, we'll be importing it from a CSV that has had PCA applied to it.
"""

from manim import *
import numpy as np
import pandas as pd
import random

class PCA3DVisualization(ThreeDScene):
    def construct(self):
        # Import data
        datapoints_df = pd.read_csv('results/wine_3_PCs.csv')

        # Sample 200 points for visualization.
        # we don't want to plot every single point because it would be too cluttered
        sampled_df = datapoints_df.sample(n=200, random_state=42)

        # Extract the proper columns
        data = sampled_df[["PC1", "PC2", "PC3"]].values
        clusters = sampled_df["Cluster"].values

        # Define colors for clusters
        cluster_colors = [
            RED, GREEN, BLUE, YELLOW, PURPLE, ORANGE, PINK, TEAL, GOLD,
            MAROON, GRAY, DARK_BLUE, WHITE
        ]
        # map each cluster to a color
        cluster_colors = cluster_colors[:len(np.unique(clusters))]
        color_map = {cluster: cluster_colors[i] for i, cluster in enumerate(np.unique(clusters))}

        # Set up 3D axes and camera
        axes = ThreeDAxes()
        self.set_camera_orientation(phi=70 * DEGREES, theta=45 * DEGREES)
        self.begin_ambient_camera_rotation(rate=0.1)
        self.add(axes)

        # Create a group to hold all points
        points = VGroup()
        for point, cluster in zip(data, clusters):
            dot = Dot3D(point=point, radius=0.05, color=color_map[cluster])
            points.add(dot)

        # Animate the points appearing
        self.play(LaggedStart(*[FadeIn(point, scale=0.5) for point in points], lag_ratio=0.01), run_time=2)
        self.wait(1)
        self.next_section()

        # Define the plane. The demo presentation just used a predefined plane to project onto.
        # That plane was chosen because it was easy to visualize and maintained most clustering
        # that the data had there. Here, instead, since we've already used PCA to reduce the data,
        # we'll just project onto the PC1-PC2 plane.
        P0 = np.array([0, 0, 0])
        plane_size = 5
        corners = [
            P0 + plane_size * np.array([1, 0, 0]) + plane_size * np.array([0, 1, 0]),
            P0 - plane_size * np.array([1, 0, 0]) + plane_size * np.array([0, 1, 0]),
            P0 - plane_size * np.array([1, 0, 0]) - plane_size * np.array([0, 1, 0]),
            P0 + plane_size * np.array([1, 0, 0]) - plane_size * np.array([0, 1, 0]),
        ]

        plane = Polygon(
            *corners,
            color=WHITE,
            fill_color=WHITE,
            fill_opacity=0.4,
            stroke_opacity=0
        )
        plane.set_shade_in_3d(True)

        self.play(Create(plane, scale=0.5), run_time = 1)
        self.wait(1)
        self.next_section()

        # Project the points onto the plane. Draw a line from the point to the plane.
        projection_lines = VGroup()
        projected_points = VGroup()
        for dot in points:
            P = dot.get_center()
            projection = np.array([P[0], P[1], 0])
            line = Line(P, projection, color=YELLOW, stroke_width=0.5)
            projection_lines.add(line)
            proj_dot = Dot3D(point=projection, radius=0.05, color=dot.get_color())
            projected_points.add(proj_dot)

        self.play(
            LaggedStart(*[Create(line) for line in projection_lines], lag_ratio=0.01),
            run_time=3
        )
        self.wait(1)
        self.next_section()

        # Get rid of the plane, have the points move to their projections
        self.play(
            FadeOut(plane),
            FadeOut(projection_lines),
            ReplacementTransform(points, projected_points),
            run_time=2
        )
        self.wait(1)
        self.stop_ambient_camera_rotation()

        # Make 2D axes and projected points
        axes2d = Axes(
            x_range=[-5, 5, 1],
            y_range=[-5, 5, 1],
            tips=False
        )

        projected_points_2d = VGroup()
        for proj_dot in projected_points:
            proj_point = proj_dot.get_center()
            x, y, _ = proj_point
            dot2d = Dot(axes2d.coords_to_point(x, y), color=proj_dot.get_color())
            projected_points_2d.add(dot2d)

        # Move the camera to a 2D view
        self.move_camera(phi=0 * DEGREES, theta = 90 * DEGREES)
        self.play(
            ReplacementTransform(axes,axes2d),
            ReplacementTransform(projected_points, projected_points_2d),
            run_time=2
        )
        self.wait(1)
