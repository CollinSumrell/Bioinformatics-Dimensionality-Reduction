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

    def setup(self):
        super().setup()
        self.remove(self.i_hat, self.j_hat)

    def add_text_with_background(self, text, position=UP, scale=0.8, bg_color=BLACK, bg_opacity=0.7):
        """
        Add text with a background rectangle.
        """
        text.scale(scale)
        text.to_edge(position)
        bg=Rectangle(width=text.get_width()+0.2, height=text.get_height()+0.2, color=bg_color, fill_opacity=bg_opacity, stroke_opacity=0).move_to(text)
        group = VGroup(bg, text)
        return(group)

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

        text1 = self.add_text_with_background(Tex("Project data onto x-axis"))
        self.play(FadeIn(text1))

        # Project onto x-axis
        x_projected_points = np.array([[x, 0] for x, _ in points])
        x_projected_dots = VGroup(
            *[
                Dot(self.plane.coords_to_point(x, y), color=ORANGE)
                for x, y in x_projected_points
            ]
        )
        self.play(Transform(initial_dots, x_projected_dots))
        self.wait()

        text2 = self.add_text_with_background(Tex("Can't distinguish clusters"))
        self.play(Transform(text1, text2))
        self.wait(3)

        # Undo x-axis projection
        self.play(Transform(initial_dots, original_dots))
        self.wait()

        text3 = self.add_text_with_background(Tex("Project data onto y-axis"))
        self.play(Transform(text1, text3))

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

        text4 = self.add_text_with_background(Tex("Still can't distinguish clusters"))
        self.play(Transform(text1, text4))
        self.wait(3)

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