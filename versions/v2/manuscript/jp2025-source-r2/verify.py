#!/usr/bin/env python3
"""Exact rational certificates for the Hadamard-closed polygons and cones."""

import json
from fractions import Fraction
from itertools import product


def scalar(value, denominator=1):
    return Fraction(value, denominator)


def dot(left, right):
    return sum((first * second for first, second in zip(left, right)), scalar(0))


def orientation(origin, left, right):
    return ((left[0] - origin[0]) * (right[1] - origin[1])
            - (left[1] - origin[1]) * (right[0] - origin[0]))


def convex_hull(points):
    ordered = sorted(set(points))
    chains = []
    for sequence in (ordered, list(reversed(ordered))):
        chain = []
        for point in sequence:
            while len(chain) >= 2 and orientation(chain[-2], chain[-1], point) <= 0:
                chain.pop()
            chain.append(point)
        chains.extend(chain[:-1])
    return chains


def construction(parameter):
    arc = [(scalar(index * index - parameter * parameter,
                   4 * (index * index + parameter * parameter)),
            scalar(-parameter * index, 2 * (index * index + parameter * parameter)))
           for index in range(parameter + 1)]
    vertices = arc + [(scalar(1, 4), scalar(0)), (scalar(1), scalar(1)),
                      (scalar(0), scalar(1, 4))]
    normals = [(scalar(1),
                scalar(4 * (parameter * parameter - index * (index + 1)),
                       parameter * parameter + index * (index + 1)),
                scalar(4 * parameter * (2 * index + 1),
                       parameter * parameter + index * (index + 1)))
               for index in range(parameter)]
    normals.extend(tuple(map(scalar, normal))
                   for normal in [(1, -4, 4), (1, -4, 3), (1, 3, -4), (1, 4, -4)])
    return vertices, normals


def diamond_certificate(point):
    horizontal, vertical = point
    slack = 1 - 4 * (abs(horizontal) + abs(vertical))
    assert slack >= 0
    weights = (4 * max(horizontal, scalar(0)) + slack / 2,
               4 * max(vertical, scalar(0)),
               4 * max(-horizontal, scalar(0)) + slack / 2,
               4 * max(-vertical, scalar(0)))
    axes = [(scalar(1, 4), scalar(0)), (scalar(0), scalar(1, 4)),
            (scalar(-1, 4), scalar(0)), (scalar(0), scalar(-1, 4))]
    assert sum(weights) == 1 and min(weights) >= 0
    assert tuple(sum(weight * axis[coordinate] for weight, axis in zip(weights, axes))
                 for coordinate in range(2)) == point
    return weights


def stochastic_matrix(point):
    horizontal, vertical = point
    unit = (1, 1, 1)
    first = (1, -1, 0)
    second = (1, 1, -2)
    return [[scalar(unit[row] * unit[column], 3)
             + horizontal * scalar(first[row] * first[column], 2)
             + vertical * scalar(second[row] * second[column], 6)
             for column in range(3)] for row in range(3)]


def check_family(parameter):
    vertices, normals = construction(parameter)
    lifted = [(scalar(1), *vertex) for vertex in vertices]
    assert len(convex_hull(vertices)) == len(vertices) == parameter + 4
    for index, normal in enumerate(normals):
        slacks = [dot(normal, vertex) for vertex in lifted]
        assert min(slacks) == 0
        assert {position for position, value in enumerate(slacks) if value == 0} == {
            index, (index + 1) % len(vertices)}
        if index < parameter:
            for position in range(parameter + 1):
                denominator = ((parameter * parameter + index * (index + 1))
                               * (position * position + parameter * parameter))
                assert slacks[position] == scalar(
                    2 * parameter * parameter * (position - index) * (position - index - 1),
                    denominator)
    unit = (scalar(1), scalar(1))
    product_count = 0
    for left, right in product(vertices, repeat=2):
        result = tuple(first * second for first, second in zip(left, right))
        assert all(dot(normal, (scalar(1), *result)) >= 0 for normal in normals)
        if left != unit and right != unit:
            diamond_certificate(result)
        else:
            assert result in vertices
        product_count += 1
    basis = [(1, 1, 1), (1, -1, 0), (1, 1, -2)]
    for point in vertices:
        matrix = stochastic_matrix(point)
        assert min(value for row in matrix for value in row) >= 0
        if point != unit:
            assert min(value for row in matrix for value in row) >= scalar(1, 6)
        assert [sum(row) for row in matrix] == [1, 1, 1]
        assert matrix == [list(column) for column in zip(*matrix)]
        for eigenvalue, vector in zip((1, *point), basis):
            assert [dot(row, vector) for row in matrix] == [eigenvalue * value for value in vector]
    for corner in product((scalar(-1, 4), scalar(1, 4)), repeat=2):
        assert min(value for row in stochastic_matrix(corner) for value in row) >= scalar(1, 6)
    assert stochastic_matrix(unit) == [[int(row == column) for column in range(3)] for row in range(3)]
    identity = (scalar(1), scalar(1), scalar(1))
    away = [normal for normal in normals if dot(normal, identity) != 0]
    assert len(away) == parameter + 2
    if parameter == 1:
        normalized = [tuple(value / dot(normal, identity) for value in normal)
                      for normal in away]
        assert sum(normal[0] for normal in normalized) == scalar(19, 9)
        assert sum(normal[0] for normal in normalized) != 1
    return {"m": parameter, "vertices": len(vertices), "facets": len(normals),
            "ordered_generator_products": product_count,
            "facets_omitting_identity": len(away),
            "nonrepresentation_certificate": "normalized diagonal trace" if parameter == 1
                                               else "more than three facets omit the identity"}


def check_product_dimension(parameter, dimension):
    vertices, normals = construction(parameter)
    cube_vertices = list(product((scalar(0), scalar(1)), repeat=dimension - 3))
    basis = [tuple(scalar(1) for coordinate in range(dimension)),
             (scalar(1), scalar(-1), scalar(0)) + (scalar(0),) * (dimension - 3),
             (scalar(1), scalar(1), scalar(-2)) + (scalar(0),) * (dimension - 3)]
    basis.extend(tuple(scalar(int(coordinate == index)) for coordinate in range(dimension))
                 for index in range(3, dimension))
    checked = 0
    for point, extra in product(vertices, cube_vertices):
        top = [row + [scalar(0)] * (dimension - 3) for row in stochastic_matrix(point)]
        bottom = [[(1 - value) / 3] * 3
                  + [value if column == index else scalar(0)
                     for column in range(dimension - 3)]
                  for index, value in enumerate(extra)]
        matrix = top + bottom
        assert all(value >= 0 for row in matrix for value in row)
        assert all(sum(row) == 1 for row in matrix)
        eigenvalues = (scalar(1), *point, *extra)
        for eigenvalue, vector in zip(eigenvalues, basis):
            assert tuple(dot(row, vector) for row in matrix) == tuple(
                eigenvalue * value for value in vector)
        checked += 1
    away = [normal for normal in normals if dot(normal, (1, 1, 1)) != 0]
    away_count = len(away) + dimension - 3
    if parameter == 1:
        assert away_count == dimension
        assert sum(normal[0] / dot(normal, (1, 1, 1)) for normal in away) == scalar(19, 9)
    else:
        assert away_count > dimension
    return {"m": parameter, "dimension": dimension,
            "facets": parameter + 2 * dimension - 2,
            "facets_omitting_identity": away_count,
            "stochastic_eigensystems": checked}


def canonical_ten():
    vertices, normals = construction(10)
    expected_arc = [(-scalar(1, 4), scalar(0)), (-scalar(99, 404), -scalar(5, 101)),
                    (-scalar(3, 13), -scalar(5, 52)), (-scalar(91, 436), -scalar(15, 109)),
                    (-scalar(21, 116), -scalar(5, 29)), (-scalar(3, 20), -scalar(1, 5)),
                    (-scalar(2, 17), -scalar(15, 68)), (-scalar(51, 596), -scalar(35, 149)),
                    (-scalar(9, 164), -scalar(10, 41)), (-scalar(19, 724), -scalar(45, 181)),
                    (scalar(0), -scalar(1, 4))]
    expected_normals = [(4, scalar(2, 5)), (scalar(196, 51), scalar(20, 17)),
                        (scalar(188, 53), scalar(100, 53)), (scalar(22, 7), scalar(5, 2)),
                        (scalar(8, 3), 3), (scalar(28, 13), scalar(44, 13)),
                        (scalar(116, 71), scalar(260, 71)), (scalar(44, 39), scalar(50, 13)),
                        (scalar(28, 43), scalar(170, 43)), (scalar(4, 19), 4),
                        (-4, 4), (-4, 3), (3, -4), (4, -4)]
    assert vertices[:11] == expected_arc
    assert [normal[1:] for normal in normals] == expected_normals
    axis_lifts = [(scalar(1), scalar(1, 4), scalar(0)),
                  (scalar(1), scalar(0), scalar(1, 4)),
                  (scalar(1), scalar(-1, 4), scalar(0)),
                  (scalar(1), scalar(0), scalar(-1, 4))]
    assert tuple(sum(vertex[coordinate] for vertex in axis_lifts) / 4
                 for coordinate in range(3)) == (1, 0, 0)
    return {"m": 10, "lifted_vertices": [[str(value) for value in (scalar(1), *point)]
                                         for point in vertices],
            "facet_normals": [[str(value) for value in normal] for normal in normals]}


def check_boundary_spectrum(dimension):
    matrix = [[scalar(1, dimension) for column in range(dimension)]
              for row in range(dimension)]
    assert all(sum(row) == 1 for row in matrix)
    unit = tuple(scalar(1) for coordinate in range(dimension))
    assert tuple(dot(row, unit) for row in matrix) == unit
    for index in range(1, dimension):
        vector = tuple(scalar(int(coordinate == index) - int(coordinate == 0))
                       for coordinate in range(dimension))
        assert all(dot(row, vector) == 0 for row in matrix)
    tail = (scalar(0),) * (dimension - 1)
    assert tuple(2 * value for value in tail) == tail
    imaginary_parts = (scalar(0), scalar(1, 10)) + (scalar(0),) * (dimension - 2)
    assert sum(imaginary_parts) == scalar(1, 10) != 0
    return {"dimension": dimension, "spectrum": ["1"] + ["0"] * (dimension - 1),
            "radial_factor": "2", "tail_unchanged": True,
            "perturbation_trace_imaginary_part": str(sum(imaginary_parts))}


def main():
    parameters = list(range(1, 17)) + [32, 64]
    families = [check_family(parameter) for parameter in parameters]
    higher_dimensions = [check_product_dimension(parameter, dimension)
                         for parameter in (1, 2) for dimension in range(3, 9)]
    result = {"arithmetic": "exact rational", "status": "PASS", "families": families,
              "cube_extensions": higher_dimensions,
              "pentagon_trace": {"normalized_diagonal_sum_at_100": "19/9",
                                 "required_trace_at_100": "1"},
              "ambient_boundary_examples": [check_boundary_spectrum(dimension)
                                            for dimension in range(2, 11)],
              "canonical_example": canonical_ten(),
              "scope": "Finite generator, facet, convex-product, stochastic-eigenbasis and normalized-trace identities; the unrestricted parameter families and identity-facet obstruction have written proofs in note.pdf."}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
