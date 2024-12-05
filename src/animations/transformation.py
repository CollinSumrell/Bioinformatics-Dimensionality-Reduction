from manim import *
import numpy as np

class LinearTransformationSceneExample(LinearTransformationScene):
    def __init__(self, **kwargs):
        LinearTransformationScene.__init__(
            self,
            show_coordinates=True,
            leave_ghost_vectors=False,
            **kwargs
        )

    # def setup(self):
    #     super().setup()
    #     self.remove(self.i_hat, self.j_hat)

    def construct(self):
        # Generate clusters
        cluster1 = np.random.multivariate_normal(
            mean=[1, 0], cov=[[0.5, 0], [0, 0.5]], size=100
        )
        cluster2 = np.random.multivariate_normal(
            mean=[3, 3], cov=[[0.5, 0], [0, 0.5]], size=100
        )

        # Combine clusters
        points = np.vstack([cluster1, cluster2])

        # Create dots for the initial positions
        initial_dots = VGroup(
            *[
                Dot(self.plane.coords_to_point(x, y))
                for x, y in points
            ]
        )
        original_dots = initial_dots.copy()
        self.add(initial_dots)
        self.wait()

        # Project onto x-axis
        x_projected_points = np.array([[x, 0] for x, y in points])
        x_projected_dots = VGroup(
            *[
                Dot(self.plane.coords_to_point(x, y), color=ORANGE)
                for x, y in x_projected_points
            ]
        )
        self.play(Transform(initial_dots, x_projected_dots))
        self.wait()

        # Undo x-axis projection
        self.play(Transform(initial_dots, original_dots))
        self.wait()

        # Project onto y-axis
        y_projected_points = np.array([[0, y] for x, y in points])
        y_projected_dots = VGroup(
            *[
                Dot(self.plane.coords_to_point(x, y), color=ORANGE)
                for x, y in y_projected_points
            ]
        )
        self.play(Transform(initial_dots, y_projected_dots))
        self.wait()

        # Undo y-axis projection
        self.play(Transform(initial_dots, original_dots))
        self.wait()

        # Calculate the centroid
        centroid = np.mean(points, axis=0)

        # Center the data
        centered_points = points - centroid
        centered_dots = VGroup(
            *[
                Dot(self.plane.coords_to_point(x, y), color=YELLOW)
                for x, y in centered_points
            ]
        )

        # Animate from initial to centered positions
        self.play(Transform(initial_dots, centered_dots))
        self.wait()

        # Calculate covariance matrix and eigenvalues/eigenvectors
        covariance_matrix = np.cov(centered_points, rowvar=False)
        eigenvalues, eigenvectors = np.linalg.eig(covariance_matrix)

        # Display eigenvectors
        vectors = VGroup()
        for eigenvalue, eigenvector in zip(eigenvalues, eigenvectors.T):
            direction = eigenvector * np.sqrt(eigenvalue)
            vector = Arrow(
                self.plane.coords_to_point(0, 0),
                self.plane.coords_to_point(*direction),
                buff=0,
                color=RED,
            )
            vectors.add(vector)

        self.play(AnimationGroup(*[GrowArrow(vector) for vector in vectors]))
        self.wait()

        for dot in initial_dots:
            self.add_moving_mobject(dot)

        # Apply linear transformation (eigenbasis)
        transformation_matrix = eigenvectors * np.sqrt(eigenvalues)
        self.apply_matrix(transformation_matrix)

        # Keep everything visible after transformation
        self.wait()

        self.apply_matrix((transformation_matrix / eigenvalues).T)

        self.wait()

        # Collapse points vertically to the x-axis
        collapsed_dots = VGroup(
            *[
                Dot(self.plane.coords_to_point(self.plane.point_to_coords(dot.get_center())[0], 0), color=YELLOW)
                for dot in initial_dots
            ]
        )

        # Animate points collapsing to the x-axis
        self.play(Transform(initial_dots, collapsed_dots))
        self.wait()