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
        bg = Rectangle(
            width=text.get_width() + 0.2,
            height=text.get_height() + 0.2,
            color=bg_color,
            fill_opacity=bg_opacity,
            stroke_opacity=0
        ).move_to(text)
        group = VGroup(bg, text)
        return group

    def add_mathtex_with_background(self, mathtex, position=DOWN, scale=0.8, bg_color=BLACK, bg_opacity=0.7):
        """
        Add MathTex with a background rectangle.
        """
        mathtex.scale(scale)
        mathtex.to_edge(position)
        bg = Rectangle(
            width=mathtex.get_width() + 0.2,
            height=mathtex.get_height() + 0.2,
            color=bg_color,
            fill_opacity=bg_opacity,
            stroke_opacity=0
        ).move_to(mathtex)
        group = VGroup(bg, mathtex)
        return group

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

        # Step 1: Project data onto x-axis
        text1 = self.add_text_with_background(Tex("Step 1: Project data onto x-axis"))
        self.play(FadeIn(text1))
        equation1 = self.add_mathtex_with_background(
            MathTex("X_{\\text{projected}} = [x, 0]")
        )
        self.play(FadeIn(equation1))

        x_projected_points = np.array([[x, 0] for x, _ in points])
        x_projected_dots = VGroup(
            *[
                Dot(self.plane.coords_to_point(x, y), color=ORANGE)
                for x, y in x_projected_points
            ]
        )
        self.play(Transform(initial_dots, x_projected_dots))
        self.wait()

        # Step 2: Explain lack of distinguishability
        text2 = self.add_text_with_background(Tex("Clusters are indistinguishable"))
        self.play(Transform(text1, text2), FadeOut(equation1))
        self.wait(3)

        # Undo x-axis projection
        self.play(Transform(initial_dots, original_dots))
        self.wait()

        # Step 3: Project data onto y-axis
        text3 = self.add_text_with_background(Tex("Step 3: Project data onto y-axis"))
        self.play(Transform(text1, text3))
        equation2 = self.add_mathtex_with_background(
            MathTex("Y_{\\text{projected}} = [0, y]")
        )
        self.play(FadeIn(equation2))

        y_projected_points = np.array([[0, y] for x, y in points])
        y_projected_dots = VGroup(
            *[
                Dot(self.plane.coords_to_point(x, y), color=ORANGE)
                for x, y in y_projected_points
            ]
        )
        self.play(Transform(initial_dots, y_projected_dots))
        self.wait()

        # Explain lack of distinguishability
        text4 = self.add_text_with_background(Tex("Clusters are still indistinguishable"))
        self.play(Transform(text1, text4), FadeOut(equation2))
        self.wait(3)

        # Undo y-axis projection
        self.play(Transform(initial_dots, original_dots))
        self.wait()

        # Step 4: Center the data
        text5 = self.add_text_with_background(Tex("Step 4: Center the data"))
        self.play(Transform(text1, text5))

        centroid = np.mean(points, axis=0)
        centered_points = points - centroid
        equation3 = self.add_mathtex_with_background(
            MathTex("X_{\\text{centered}} = X - \\bar{X}")
        )
        self.play(FadeIn(equation3))

        centered_dots = VGroup(
            *[
                Dot(self.plane.coords_to_point(x, y), color=YELLOW)
                for x, y in centered_points
            ]
        )
        self.play(Transform(initial_dots, centered_dots))
        self.wait()

        # Step 5: Compute covariance matrix
        text6 = self.add_text_with_background(Tex("Step 5: Compute covariance matrix"))
        self.play(Transform(text1, text6), FadeOut(equation3))
        equation4 = self.add_mathtex_with_background(
            MathTex("C = \\frac{1}{n-1} X^T X")
        )
        self.play(FadeIn(equation4))

        covariance_matrix = np.cov(centered_points, rowvar=False)
        eigenvalues, eigenvectors = np.linalg.eig(covariance_matrix)
        self.wait()

        # Step 6: Rotate data to PC space
        text7 = self.add_text_with_background(Tex("Step 6: Rotate data to PC space"))
        self.play(Transform(text1, text7), FadeOut(equation4))
        self.wait()

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

        transformation_matrix = eigenvectors * np.sqrt(eigenvalues)
        self.apply_matrix(transformation_matrix)
        self.wait()

        # Step 7: Collapse to original x-axis
        text8 = self.add_text_with_background(Tex("Step 7: Collapse data to original x-axis"))
        self.play(Transform(text1, text8))
        self.wait()

        self.apply_matrix((transformation_matrix / eigenvalues).T)
        self.wait()

        collapsed_dots = VGroup(
            *[
                Dot(self.plane.coords_to_point(self.plane.point_to_coords(dot.get_center())[0], 0), color=YELLOW)
                for dot in initial_dots
            ]
        )
        self.play(Transform(initial_dots, collapsed_dots))
        self.wait()

        text9 = self.add_text_with_background(Tex("Clusters are easily distinguishable"))
        self.play(Transform(text1, text9))