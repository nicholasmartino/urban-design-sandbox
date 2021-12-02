import itertools
import gc
import itertools
import math
import sys
import time

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.affinity import scale, translate
from shapely.geometry import Point, Polygon, LineString, MultiPoint, MultiLineString, MultiPolygon
from shapely.ops import split, unary_union, nearest_points, linemerge, snap
from tqdm import tqdm
import matplotlib
matplotlib.use('TkAgg')

plt.show()


def azimuth(ln):
    """
	Returns the azimuth of a LineString.
	:param ln: LineString
	:return: number
	"""
    return math.degrees(math.atan2((ln.xy[0][1] - ln.xy[0][0]), (ln.xy[1][1] - ln.xy[1][0])))


def relative_max_min(df):
    """
	Returns difference between maximum and minimum X and Y coordinates.
	:param df: GeoDataFrame with LineStrings as geometries
	:return: GeoDataFrame
	"""

    df['maxX'] = [max(ln.xy[0]) for ln in df.geometry]
    df['minX'] = [min(ln.xy[0]) for ln in df.geometry]
    df['maxY'] = [max(ln.xy[1]) for ln in df.geometry]
    df['minY'] = [min(ln.xy[1]) for ln in df.geometry]

    # Get X coordinates of line vertices that have the maximum and minimum Y coordinates
    df['X_of_maxY'] = [ln.xy[0][ln.xy[1].index(max(ln.xy[1]))] for ln in df.geometry]
    df['X_of_minY'] = [ln.xy[0][ln.xy[1].index(min(ln.xy[1]))] for ln in df.geometry]

    # Get Y coordinates of line vertices that have the maximum and minimum X coordinates
    df['Y_of_maxX'] = [ln.xy[1][ln.xy[0].index(max(ln.xy[0]))] for ln in df.geometry]
    df['Y_of_minX'] = [ln.xy[1][ln.xy[0].index(min(ln.xy[0]))] for ln in df.geometry]

    # Get difference between X and Y coordinates normalized by length
    df['X_x'] = (df['maxX'] - df['minX']) / df.length
    df['Y_y'] = (df['maxY'] - df['minY']) / df.length

    return df


def get_gtypes(gdf):
    return [geom.__class__.__name__ for geom in gdf['geometry']]


def drop_index(gdf2drop):
    for col in ['index_left', 'index_right', 'join_uid']:
        if col in gdf2drop.columns:
            gdf2drop = gdf2drop.drop(col, axis=1)
    return gdf2drop


def bivariate_map(gdf, columns, classes=3):
    """
    Adds a new column in the gdf with classes to produce a bivariate map
    :param gdf: GeoDataFrame
    :param columns: tuple or list with the name of two columns from the gdf
    :param classes: number of percentile divisions
    :return:
    """

    for col in columns:
        assert col in gdf.columns, KeyError(f'{col} column not found in GeoDataFrame')


    return gdf


# Class Polyline (?)

def segments(curve):
    return list(map(LineString, zip(curve.coords[:-1], curve.coords[1:])))

def divide_line_by_length(gdf, length):
    """
	Returns a list of MultiPoints separated according to a defined length.
	:param gdf: GeoDataFrame
	:param length: number
	:return: list
	"""

    if type(length) in [int, float]:
        length = [length]

    length = [-lgth if lgth < 0 else lgth for lgth in length]

    mpts = []
    for geom, lgth in zip(gdf['geometry'], length):
        if geom.__class__.__name__ in ['MultiPolygon', 'Polygon']:
            geom = geom.boundary
        pts = []
        for distance in np.arange(0, geom.length, lgth):
            pts.append(unary_union(geom.interpolate(distance)))
        mpts.append(MultiPoint(pts))
    return mpts

def divide_line_by_count(line, n_segments):
    """
	Returns a list whose length is equal to n_segments of LineStrings divided from the line input
	:param line: LineString
	:param n_segments: int
	:return:
	"""
    assert len(line) == len(n_segments), AssertionError("Arrays of different sizes")

    splitters = [[ln.interpolate((j/i), normalized=True) for j in range(1, i)] for ln, i in zip(line, n_segments)]
    splitters = [MultiPoint([snap(pt, ln, tolerance=0) for pt in points]).buffer(0.001) for points, ln in zip(splitters, line)]
    return [MultiLineString([segment for segment in split(ln, splitter) if segment.length > 0.1]) for ln, splitter in zip(line, splitters)]

def get_point(geom, parameter, buffer=1):
    if parameter == 'start':
        return Point(geom.xy[0][0], geom.xy[1][0]).buffer(buffer)
    elif parameter == 'end':
        return Point(geom.xy[1][0], geom.xy[0][0]).buffer(buffer)
    else:
        return geom.centroid


class Shape:
    def __init__(self, gdf, verbose=False, reset_index=False):
        """
		Shape Set / Vector Layer
		:param gdf: GeoDataFrame
		"""

        if 'shapely' in str(type(gdf)): gdf = gpd.GeoDataFrame({'geometry': [gdf.copy()]})
        if reset_index:
            gdf = gdf.copy().reset_index(drop=True)

        else:
            assert type(gdf).__name__ == 'GeoDataFrame', AssertionError("Input is not a valid GeoDataFrame")
            assert len(list(gdf.index)) == len(list(gdf.index.unique())), IndexError("Duplicated indices found in reference geometry, set reset_index to True")

        gdf = gdf.copy()
        gdf['id'] = gdf.index
        gdf['area'] = gdf.area
        self.gdf = gdf
        self.splitters = []  # List of shapely LineStrings
        self.meridian = []
        self.verbose = verbose
        return

    def is_convex(self):
        """
        Adds a column in the gdf testing whether the shapes are convex or not
        :return:
        """

        self.gdf['convex'] = [False if geom.minimum_rotated_rectangle.area > 1.02 * (geom.area) else True for geom in self.gdf.geometry]
        return self.gdf

    def assemble_to_area(self, area, tolerance=0.2):
        """
        Recursivelly assemble adjacent geometries until reaches an approximate area
        :return:
        """
        gdf = self.gdf.copy()
        assert len(gdf) > 0, AssertionError("Input GeoDataFrame is empty")
        
        # Find clusters of adjacent geometries by dissolving them and assigning a cluster index
        dis = self.dissolve().reset_index(drop=True)
        dis['cluster'] = dis.index
        dis['n_parts'] = ((dis.area/area).apply(np.floor)).astype(int)
        gdf = gpd.sjoin(gdf, dis.loc[:, ['cluster', 'n_parts', 'geometry']])

        # Divide pieces to assemble and pieces to dissolve
        condition = (gdf['n_parts'] <= 1)
        if len(gdf[condition]) > 0:
            no_assemble = Shape(gdf[condition], reset_index=True).dissolve()
        else:
            no_assemble = gpd.GeoDataFrame()
        assemble = gdf[~condition].copy()

        # Sort pieces to assemble clockwise
        assembled = gpd.GeoDataFrame({'geometry': []}, crs=gdf.crs)
        for i in assemble['cluster'].unique():
            to_assemble = assemble[assemble['cluster'] == i]
            target_area = (dis[dis['cluster'] == i].area.sum() / to_assemble['n_parts']).mean()

            while len(to_assemble) > 0:

                # Start with furtherst from centroid
                to_assemble = Shape(to_assemble).get_distances_from_centroid()

                # Pick closest geometry in each iteration
                if len(to_assemble) > 1:
                    initial = to_assemble['d2centroid'] == max(to_assemble['d2centroid'])

                    shape_area = to_assemble[initial].area.sum()

                    while shape_area < (target_area * (1 - tolerance)):
                        closest_i = list(Shape(to_assemble[initial]).get_closest(to_assemble[~initial])['ref_id'])[0]
                        initial.loc[closest_i] = True
                        shape_area = to_assemble[initial].area.sum()

                        if (shape_area >= (target_area * (1 - tolerance))) or (initial.sum() == len(to_assemble)):
                            assembled = pd.concat([assembled, Shape(to_assemble[initial]).dissolve()])
                            to_assemble = to_assemble[~initial]
                            break

                else:
                    assembled = pd.concat([assembled, to_assemble])
                    break

        return pd.concat([no_assemble, assembled])

    def get_rectangularity(self):
        """
        http://www.cyto.purdue.edu/cdroms/micro2/content/education/wirth10.pdf page 45
        :return:
        """
        gdf = self.gdf.copy()
        gdf['rectangularity'] = [geom.area/geom.minimum_rotated_rectangle.area for geom in self.gdf['geometry']]
        return gdf

    def get_distances_from_centroid(self):
        gdf = self.gdf.copy()
        centroid = gdf.unary_union.centroid

        for i in gdf.index:
            geom = gdf.loc[i, 'geometry']
            gdf.loc[i, 'd2centroid'] = geom.distance(centroid)
        return gdf

    def export(self, name='shape'):
        return self.gdf.to_file(f'{name}.geojson', driver='GeoJSON')

    def loc_id(self, geom, round=5):
        return f'{round(geom.x, round)}_{round(geom.y, round)}'

    def simplify(self, tolerance=0.1):
        """

        :param tolerance:
        :return: MultiPolygon
        """
        gdf = self.gdf.copy()
        gdf['geometry'] = gdf.simplify(tolerance, preserve_topology=True).buffer(1).buffer(-1).simplify(tolerance,
                                                                                                        preserve_topology=True)
        return gdf

    def divorce(self):
        """
		Multi-part to single-part geometry
		:return:
		"""
        gdf = self.gdf.copy()
        gdf['geom_type'] = self.get_geom_types()
        multi_gdf = gdf[gdf['geom_type'].isin(['MultiPolygon', 'MultiLineString'])]
        multi_gdf = multi_gdf.explode().reset_index(drop=True).rename(columns={0: 'geometry'}).set_geometry('geometry')
        single_gdf = gdf[~gdf['geom_type'].isin(['MultiPolygon', 'MultiLineString'])]
        if len(single_gdf) > 0:
            gdf = pd.concat([single_gdf, multi_gdf])
        else:
            gdf = multi_gdf

        if 'id' in gdf.columns:
            gdf['old_id'] = gdf['id']
        gdf['id'] = gdf.index
        return gdf

    def explode(self):
        """
		MultiPolygons to Polygons. LineString polylines to LineString segments.
		:return: GeoDataFrame with exploded elements as MultiGeometries in the Geometry column
		"""

        gdf = self.gdf.copy()

        # Explode multi geometries
        gdf = Shape(gdf).divorce()

        # Get geometry types
        gdf['geom_type'] = Shape(gdf, reset_index=True).get_geom_types()
        geom_types = get_gtypes(gdf)
        unique_types = list(set(geom_types))
        cols = [col for col in gdf.columns if col != 'geometry']

        # Transform polygons into linestrings and clean them
        gdf.loc[gdf['geom_type'] == 'Polygon', 'geometry'] = gdf[gdf['geom_type'] == 'Polygon'].boundary
        gdf = gdf.dropna(subset=['geometry']).reset_index(drop=True)

        # # Explode LineStrings into segments
        # gdf['geom_type'] = Shape(gdf).get_geom_types()
        # gdf['is_multi'] = [True  if 'Multi' in geom.geom_type else False for geom in gdf['geometry']]
        # m_gdf = Shape(gdf.loc[gdf['is_multi'], :].copy()).divorce()
        # m_gdf.index = m_gdf['old_id']
        # gdf = pd.concat([gdf[~gdf['is_multi']], m_gdf])
        try:
            gdf['segments'] = [segments(geom) for geom in gdf['geometry']]
        except:
            pass
        gdf = gdf.reset_index(drop=True)

        exploded = gpd.GeoDataFrame()
        exploded['parent_id'] = [gdf.loc[i, 'old_id'] for i in gdf.index for line in gdf.loc[i, 'segments']]
        exploded['geometry'] = [line for i in gdf.index for line in gdf.loc[i, 'segments']]
        exploded = exploded.set_geometry('geometry')

        # for i in list(gdf.index):
        #     geom = gdf.loc[i, 'geometry']
        #     if geom.geom_type == 'LineString':
        #         segs = segments(geom)
        #     else:
        #         if len(geom) > 0:
        #             segs = []
        #             for sub_geom in geom:
        #                 segs = segs + segments(sub_geom)
        #
        #     try:
        #         if len(segs) > 0:
        #             for seg in segs:
        #                 k = len(exploded)
        #                 exploded.loc[k, cols] = gdf.loc[i, cols]
        #                 exploded.loc[k, 'geometry'] = seg
        #         else:
        #             k = len(exploded)
        #             exploded.loc[k, cols] = gdf.loc[i, cols]
        #             exploded.loc[k, 'geometry'] = seg
        #     except:
        #         pass

        # exploded = exploded.set_geometry('geometry')
        # exploded['sid'] = exploded.reset_index(drop=True).index
        # exploded = exploded.apply(pd.to_numeric, errors='ignore').set_geometry('geometry')

        exploded.crs = self.gdf.crs
        exploded['length'] = exploded.length
        exploded['sid'] = exploded.index
        return exploded
        """
		if (('MultiPolygon' in geom_types) or ('MultiLineString' in geom_types)) and len(unique_types) < 2:
			gdf = gdf.explode()
			gdf = gdf.reset_index(drop=True).rename(columns={0: 'geometry'})
			gdf = gdf.set_geometry('geometry')
			gdf.crs = self.gdf.crs
			# return gdf

		else:
			exploded = gpd.GeoDataFrame({'geometry': []}, geometry='geometry')
			for i in list(gdf.index):
				geom = gdf.loc[i, 'geometry'] 
				geom_type = geom.geom_type

				if geom_type == 'LineString':
					segs = MultiLineString(list(map(LineString, zip(geom.coords[:-1], geom.coords[1:]))))

				elif geom_type == 'Polygon':
					segs = segments(geom)

				elif geom_type == 'MultiLineString':
					for line in geom:
						segs = MultiLineString(list(map(LineString, zip(line.coords[:-1], line.coords[1:]))))
						# for seg in segs:
						#     k = len(exploded)
						#     exploded.loc[k, cols] = gdf.loc[i, cols]
						#     exploded.loc[k, 'geometry'] = seg

				if len(segs) > 0:
					for seg in segs:
						k = len(exploded)
						exploded.loc[k, cols] = gdf.loc[i, cols]
						exploded.loc[k, 'geometry'] = seg
				else:
					k = len(exploded)
					exploded.loc[k, cols] = gdf.loc[i, cols]
					exploded.loc[k, 'geometry'] = seg
					
			exploded['sid'] = exploded.index
			return exploded
		"""

    def reduce(self, buffer=5):
        """
		Scale down using a negative buffer
		:return: GeoDataFrame
		"""

        uu = self.gdf.copy().buffer(1).unary_union.buffer(-buffer)
        if 'Multi' in uu.geom_type:
            return gpd.GeoDataFrame({'geometry': [pol for pol in uu]}, crs=self.gdf.crs)
        else:
            return gpd.GeoDataFrame({'geometry': [self.gdf.copy().buffer(1).unary_union.buffer(-buffer)]}, crs=self.gdf.crs)

    def dissolve(self, tolerance=0.1, join_operations=None, min_area=1):
        """
		Merge adjacent polygons with buffer and unary_union
		:return: GeoDataFrame
		"""
    
        uu = self.gdf.copy().buffer(tolerance, join_style=2, cap_style=2).unary_union.buffer(-tolerance, join_style=2, cap_style=2)

        if 'Multi' in uu.geom_type:
            gdf = gpd.GeoDataFrame({'geometry': [Polygon(geom.exterior) for geom in uu]}, crs=self.gdf.crs)
        else:
            try:
                gdf = gpd.GeoDataFrame({'geometry': [Polygon(uu.exterior)]}, crs=self.gdf.crs)
            except:
                gdf is None

        if join_operations is not None:
            gdf = Analyst(gdf, self.gdf.copy()).spatial_join(operations=join_operations)
        
        raw_gdf = self.gdf.copy()
        raw_gdf['id'] = raw_gdf.index
        gdf['centroid'] = gdf.centroid.buffer(1)
        gdf = gdf[gdf.area > min_area]
        gdf = gdf.set_geometry('centroid')
        gdf['parent_id'] = [int(i) for i in Analyst(gdf, self.gdf.copy().loc[:, ['id', 'geometry']]).spatial_join(operations=['max'])['id_max']]
        gdf = gdf.set_geometry('geometry')
        return gdf.drop('centroid', axis=1)

    def extract_open_boundaries(self, inner_rings=True):
        """
		Dissolve adjacent polygons and extract segments from original shape located on the outer boundary of dissolved geometry

		:return: GeoDataFrame with MultiLineStrings
		"""

        # Explode shape to get segments
        segments = self.explode()
        segments['pid'] = list(self.gdf.loc[segments['parent_id'], 'pid'])
        segments['length'] = segments.length

        # Intersect with reduced polygons
        reduced = self.reduce()
        overlay = gpd.overlay(segments, reduced)
        segments.crs = self.gdf.crs

        segments = segments[~segments['sid'].isin(overlay['sid'])]
        if inner_rings:
            return segments

        else:
            blocks = self.dissolve()
            blocks['geometry'] = blocks.buffer(-2)
            overlay = gpd.overlay(segments, blocks)

            return segments[~segments['sid'].isin(overlay['sid'])]

    def extract_inner_rings(self, multi=False):
        """
		Extract inner rings from polygons if exists
		:return:
		"""

        # Explode shape to get segments
        segments = self.explode()
        segments['length'] = segments.length

        # Intersect with reduced polygons
        dissolved = self.dissolve()
        reduced = Shape(dissolved).reduce()

        overlay = gpd.overlay(segments, reduced)

        return segments[segments['sid'].isin(overlay['sid'])]

    def extract_vertices(self):

        gdf = self.explode()
        vertices = gpd.GeoDataFrame()

        for i in gdf.index:
            coord = gdf.loc[i, 'geometry'].coords
            for pt in [Point(coord[0]), Point(coord[len(coord) - 1])]:
                j = len(vertices)
                vertices.loc[j, gdf.columns] = gdf.loc[i, :]
                vertices.loc[j, 'geometry'] = pt

        vertices = vertices.rename({'id': 'line_id'}, axis=1)
        vertices['vid'] = vertices.index
        vertices = vertices.drop_duplicates('geometry')
        return vertices

    def get_geom_types(self):
        return [geom.geom_type for geom in self.gdf['geometry']]

    def get_indicators(self, explode=False):
        if explode: self.gdf = self.explode()
        gdf = self.gdf.copy()

        gdf['n_vertices'] = [len(geom.exterior.coords) for geom in gdf['geometry']]
        gdf['area'] = [geom.area for geom in gdf['geometry']]
        gdf['perimeter'] = [geom.length for geom in gdf['geometry']]
        gdf['ch_perimeter'] = [geom.convex_hull.length for geom in gdf['geometry']]

        # Calculate elongation based on minimum rotated rectangles
        min_rot_rec = self.min_rot_rec()
        gdf['vector'] = [geom[0] if geom is not None else None for geom in min_rot_rec['segments']]
        gdf['azimuth'] = [azimuth(geom) if geom is not None else None for geom in gdf['vector']]

        elongation = []
        width = []
        length = []
        for lengths in min_rot_rec['lengths']:
            if lengths is not None:
                if max(lengths) != 0:
                    elongation.append(1 - (min(lengths) / max(lengths)))
                else:
                    elongation.append(None)
                width.append(min(lengths))
                length.append(max(lengths))
            else:
                elongation.append(None)
                width.append(None)
                length.append(None)

        gdf['elongation'] = elongation
        gdf['width'] = width
        gdf['length'] = length
        return gdf

    def get_closest(self, reference_geom, reset_index=True):
        """
		For each base geometry (left) in the object's gdf, return the closest element
		from a reference set (right) that is directly connected to the base geometry.

		:param reference_geom: Shapely geometry to search to
		:param reset_index: If true returns indexes different from output geometry
		:return:
		"""

        gdf = self.gdf.copy()
        reference_geom = reference_geom.copy()
        reference_geom['ref_id'] = reference_geom.index

        # Draw a line from centroid of each base element to closest point on reference layer
        st = time.time()
        gdf['nearest'] = [
            scale(LineString([ctr, nearest_points(ctr, reference_geom.unary_union)[1]]), 1.1, 1.1, 1.1, 'centroid')
            for ctr in gdf.centroid]
        nearest = gdf.copy()
        nearest = nearest.drop('geometry', axis=1).set_geometry('nearest')
        nearest['base_id'] = nearest['id']
        nearest['nearest'] = nearest.buffer(1)
        nearest.crs = gdf.crs
        if self.verbose: print(f"Nearest Points: {round((time.time() - st) / 60, 3)} minutes")

        # Identify base line directly connected to reference layer
        inters = gpd.overlay(reference_geom.loc[:, ['ref_id', 'geometry']], nearest.loc[:, ['base_id', 'nearest']])
        ref_outs = [list(inters[inters['base_id'] == b_id]['ref_id'])[0] for b_id in gdf['id'].unique()]

        out_gdf = reference_geom.loc[ref_outs, :].copy()
        out_gdf['base_id'] = list(gdf['id'])

        if reset_index:
            out_gdf.index = list(out_gdf['base_id'])
            out_gdf['id'] = out_gdf.index

        return out_gdf

    def get_largest_by_column(self, column):
        """
        For elements within the same group in a defined column, get largest Shape.
        :return: GeoDataFrame with one geometry per unique element of the input column
        """

        gdf = self.gdf.copy()
        assert column in gdf, KeyError(f"{column} not found in GeoDataFrame")

        gdf['length'] = gdf.length
        gdf = gdf.sort_values(by=[column, 'length'], ascending=False)
        return gpd.GeoDataFrame({
            column: [i  for i in gdf[column].unique()],
            'geometry': [list(gdf[gdf[column] == i]['geometry'])[0] for i in gdf[column].unique()]
        }, crs=gdf.crs)

    def move(self, reference_geom):
        """
        Move elements of the gdf altogether to the centroid of the reference_geom
        :return:
        """
        gdf = self.gdf.copy()
        gdf_ctr = self.gdf.unary_union.centroid
        ref_ctr = reference_geom.centroid
        gdf['geometry'] = [translate(geom, xoff=ref_ctr.x - gdf_ctr.x, yoff=ref_ctr.y - gdf_ctr.y) for geom in gdf['geometry']]
        return gdf

    def flip(self, x, y, reference_geom=None):
        """
        Flip Shapes according to a reference geometry
        :param x: 
        :param y: 
        :param reference_geom: 
        :return: 
        """

        gdf = self.gdf.copy()
        gdf['geometry'] = [scale(geom, xfact=x, yfact=y, origin=reference_geom.unary_union.centroid) for geom in gdf['geometry']]
        assert len(gdf['geometry']) == len(self.gdf['geometry'])
        return gdf

    def min_rot_rec(self, quick_run=False):
        gdf = gpd.GeoDataFrame()
        gdf['min_rot_rec'] = [geom.minimum_rotated_rectangle.boundary if geom.length > 1 else None for geom in
                              self.gdf['geometry']]
        gdf['id'] = self.gdf.index
        gdf['segments'] = [list(map(LineString, zip(ln.coords[:-1], ln.coords[1:]))) if ln is not None else None for ln
                           in gdf['min_rot_rec']]
        gdf['lengths'] = [[int(ln.length) for ln in geom] if geom is not None else None for geom in gdf['segments']]

        # Check if shape is a square (all sides of bounding box are equal)
        result = all(element == list(gdf['lengths'])[0][0] for element in list(gdf['lengths'])[0])
        if result:
            largest = [[0, 2] for lengths in gdf['lengths']]
            smallest = [[1, 3] for lengths in gdf['lengths']]

        else:
            largest = [[i for i, length in enumerate(lengths) if length >= max(lengths)]
                       if lengths is not None else None for lengths in gdf['lengths']]
            smallest = [[i for i, length in enumerate(lengths) if length <= min(lengths)]
                        if lengths is not None else None for lengths in gdf['lengths']]

        gdf['largest_segment'] = [geom[i[0]] if geom is not None else None for i, geom in zip(largest, gdf['segments'])]
        gdf['largest_parallel'] = [geom[i[1]] if geom is not None else None for i, geom in zip(largest, gdf['segments'])]
        gdf['smallest_segment'] = [geom[i[0]] if geom is not None else None for i, geom in zip(smallest, gdf['segments'])]
        gdf['smallest_parallel'] = [geom[i[1]] if geom is not None else None for i, geom in zip(smallest, gdf['segments'])]

        if not quick_run:
            # Get line connecting the centroids of the smaller sides of bounding rectangle
            gdf['meridian'] = [LineString([smallest.centroid, parallel.centroid]) if smallest is not None else None
                               for smallest, parallel in zip(gdf['smallest_segment'], gdf['smallest_parallel'])]
            gdf['equator'] = [LineString([largest.centroid, parallel.centroid]) if largest is not None else None
                              for largest, parallel in zip(gdf['largest_segment'], gdf['largest_parallel'])]

            # Get diagonal line
            gdf['diagonal1'] = [LineString([Point(ln.coords[0]), Point(ln.coords[2])]) if ln is not None else None for ln in gdf['min_rot_rec']]
            gdf['diagonal2'] = [LineString([Point(ln.coords[1]), Point(ln.coords[3])]) if ln is not None else None for ln in gdf['min_rot_rec']]

        gdf['geometry'] = gdf['min_rot_rec']
        return gdf

    def min_bbox(self):
        gdf = self.gdf.copy()
        gdf['geometry'] = [Polygon(
            [Point(geom.bounds[0], geom.bounds[1]), Point(geom.bounds[2], geom.bounds[1]),
             Point(geom.bounds[2], geom.bounds[3]), Point(geom.bounds[0], geom.bounds[3])]
        ) for geom in gdf['geometry']]
        return gdf

    def offset_in(self, distance, polygon):
        """
		Offsets a line string and keep only shapes located inside a polygon.

		:param distance: float
		:param polygon: Polygon
		:return: GeoDataFrame
		"""

        shapes = gpd.GeoDataFrame(columns=self.gdf.columns)
        for i in list(self.gdf.index):
            shape = self.gdf.loc[i, 'geometry'].simplify(1, preserve_topology=True)
            left_offset = shape.parallel_offset(distance=distance, side='left')
            right_offset = shape.parallel_offset(distance=distance, side='right')
            if left_offset.centroid.intersects(polygon):
                l = len(shapes)
                shapes.loc[l, :] = self.gdf.loc[i, :]
                shapes.loc[l, 'geometry'] = left_offset
            if right_offset.centroid.intersects(polygon):
                l = len(shapes)
                shapes.loc[l, :] = self.gdf.loc[i, :]
                shapes.loc[l, 'geometry'] = right_offset

        shapes = shapes.apply(pd.to_numeric, errors='ignore')
        shapes = shapes.set_geometry('geometry')
        shapes.crs = self.gdf.crs
        return shapes

    def connect_parallels(self, spines, width=20, tolerance=0.2):
        """

		:param spines:
		:param width:
		:param tolerance:
		:return:
		"""
        spines = relative_max_min(spines)
        spines.crs = self.gdf.crs

        # Load shape boundaries
        boundaries = self.gdf.copy()
        boundaries['geometry'] = boundaries.boundary
        bnd_segments = Shape(boundaries).divorce()
        bnd_segments['bound_id'] = bnd_segments.reset_index(drop=True).index

        # Get relative XY coordinates for boundary line segments
        bnd_segments = relative_max_min(bnd_segments)
        bnd_segments['length'] = bnd_segments.length

        # Extract boundary segments with XY coordinates similar to boundary XY coordinates (parallel)
        for i in bnd_segments['bound_id'].unique():

            sdf = spines[spines['id'] == i]

            max_Yy = max(sdf['Y_y'], default=0) * (1 + tolerance)
            min_Yy = min(sdf['Y_y'], default=0) * (1 - tolerance)
            max_Xx = max(sdf['X_x'], default=0) * (1 + tolerance)
            min_Xx = min(sdf['X_x'], default=0) * (1 - tolerance)

            # Get parallel based on directions
            parallel = bnd_segments[
                (bnd_segments['bound_id'] == i) &
                (((bnd_segments['X_x'] < max_Xx) & (bnd_segments['X_x'] > min_Xx)) |
                 ((bnd_segments['Y_y'] < max_Yy) & (bnd_segments['Y_y'] > min_Yy)))
                ]
            parallel = linemerge(list(parallel['geometry']))
            if parallel.__class__.__name__ == 'LineString': parallel = [parallel]

            # Divide spine and connect divisors with closest point on parallel segments
            divisions = divide_line_by_length(sdf, length=width)
            for points in divisions:
                for pt in points:
                    for segment in parallel:
                        if segment.length > (width * 2):
                            j = len(spines)
                            spines.at[j, 'pol_id'] = i
                            spines.at[j, 'geometry'] = LineString([nearest_points(pt, segment)[1], pt])

        return spines

    def create_court(self, depth):
        """
		Creates a hole inside polygon according to defined depth of border
		:return: GeoDataFrame
		"""
        if type(depth) == int:
            depth = [depth for i in self.gdf.index]

        # Offset shape boundary according to defined depth to find inner boundary
        inner = self.gdf.copy()
        inner['geometry'] = [geom.buffer(-float(d)) for d, geom in zip(depth, inner['geometry'])]

        # Subtract court from outer boundary
        outer = self.gdf.copy()
        return gpd.overlay(outer, inner, how="difference")

    def subdivide_obb(self, max_area=10000):
        """
        Recursively subdivide geometries until areas are smaller than maximum area
        :param max_area:
        :return:
        """

        gdf = self.gdf.copy()
        assert len(gdf) > 0, AssertionError("Empty input, subdivision not performed")
        all_geometries = gpd.GeoDataFrame()

        # Recursively break parcel until all geometries have areas smaller than maximum area
        for i in gdf.index:
            # Get original geometries
            geometries = gdf.loc[[i], :].reset_index(drop=True)
            geometries['parent_id'] = i
            max_geom_area = max([gdf.loc[i, :]['geometry'].area])

            while max_geom_area > max_area:
                geometry = geometries.at[0, 'geometry']
                if geometry.area > max_area:

                    # Get two largest line segments of minimum bounding rectangle
                    min_rot_rec = Shape(gpd.GeoDataFrame({'geometry': [geometry]}, geometry='geometry')).min_rot_rec(quick_run=True) # !!!
                    lar = min_rot_rec['largest_segment'].values[0]
                    par = min_rot_rec['largest_parallel'].values[0]

                    # Split geometries by line between the centroids of the two largest line segments
                    splitter = scale(LineString([lar.centroid, par.centroid]), 1.2, 1.2, 1.2, origin='centroid')
                    splits = gpd.GeoDataFrame({'geometry': split(geometry, splitter)}, geometry='geometry')

                    # Save splitter on global parameters
                    self.splitters.append(splitter)

                    # Delete original geometry and add splits to geometry list
                    geometries = geometries.reset_index(drop=True).drop(0)
                    geometries = pd.concat([geometries, splits], ignore_index=True)
                    geometries['area'] = geometries.area

                    # Filter geometries smaller and larger than maximum area
                    smaller = geometries[geometries['area'] < max_area]
                    geometries = geometries[geometries['area'] > max_area]

                    # Collect parcels smaller than maximum area
                    smaller['parent_id'] = i
                    all_geometries = pd.concat([all_geometries, smaller], ignore_index=True)

                    # Sort values by area
                    if len(geometries) > 0:
                        geometries = geometries.sort_values(by='area', ascending=False).reset_index(drop=True)
                    else:
                        break

                else:
                    all_geometries = pd.concat([all_geometries, geometries], ignore_index=True)
                    break
            else:
                all_geometries = pd.concat([all_geometries, geometries], ignore_index=True)
        assert len(all_geometries) > 0, AssertionError("Subdivision returned empty geometries")
        return all_geometries

    def subdivide_offset(self, width=30, depth=12, outwards=False):
        """
        Offset geometry according to a defined depth and subdivide it into pieces of an approximate width
        :param width:
        :param depth:
        :param outwards:
        :return:
        """

        if type(depth) in [int, float]:
            depth = [depth for geom in self.gdf['geometry']]

        if type(width) in [int, float]:
            width = [width for geom in self.gdf['geometry']]

        assert len(depth) == len(self.gdf['geometry'])

        border = self.create_court(depth=depth)
        inner = self.gdf.copy()
        outer = self.gdf.copy()

        inner['geometry'] = [geom.buffer(float(min([d, -d]))) for d, geom in zip(depth, inner['geometry'])]

        outer = outer[inner.area > 0]
        inner = inner[inner.area > 0]

        if len(inner) > 0:
            if not outwards:
                # Simplify geometries to extract important vertices
                inner_bounds = inner.copy()
                inner_bounds['geometry'] = inner_bounds['geometry'].simplify(0.1)
                inner_bounds = Shape(inner_bounds).explode()
                inner_vertices = Shape(inner_bounds).extract_vertices()

                # Subdivide inner boundary according to defined width
                inner['divider'] = divide_line_by_length(gdf=inner, length=width)

                # # Append vertices to inner boundary dividers
                # for i in inner.index:
                #     for geom in list(inner_vertices[inner_vertices['old_id'] == i]['geometry']):
                #         inner.loc[i, 'divider'].append(geom)

                # Get point on outer boundary closest to inner boundary dividers
                outer['geometry'] = outer.boundary
                outer['divider'] = [MultiPoint([nearest_points(pt, outer_ln)[1] for pt in mpt]) if outer_ln.length > 0 else None
                                    for mpt, outer_ln in zip(inner['divider'], outer['geometry'])]

                gdf = inner
                other = outer

            else:
                # Subdivide outer boundary according to defined width
                outer['divider'] = divide_line_by_length(gdf=outer, length=width)

                # Get point on inner boundary closest to outer boundary divider
                inner['geometry'] = inner.boundary
                inner['divider'] = [MultiPoint([nearest_points(pt, inner_ln)[1] for pt in mpt]) if inner_ln.length > 0 else None
                                    for mpt, inner_ln in zip(outer['divider'], inner['geometry'])]
                gdf = outer
                other = inner

            # Create the splitter lines
            gdf['splitters'] = [
                MultiLineString([LineString([pt0, pt1]) for pt0, pt1 in zip(mpt0, mpt1)]) if mpt1 is not None else None for
                mpt0, mpt1 in zip(gdf['divider'], other['divider'])]

            outer_bounds = outer.copy()
            outer_bounds['geometry'] = outer.simplify(1)
            outer_bounds = Shape(outer_bounds).explode()

            # Connect inner polygon vertices to closest point on the two closest edges of the outer polygon
            splitters = []
            for geom in list(inner_vertices['geometry']): # Iterate over vertices of inner boundary

                # Get closes points on two closest edges
                pt1 = nearest_points(geom, outer_bounds.unary_union)[1]
                other_bounds = outer_bounds.copy()[~outer_bounds.intersects(pt1.buffer(2))]
                pt2 = nearest_points(geom, other_bounds.unary_union)[1]

                # Create two splitter lines connecting the vertices to closest points and append to "splitters"
                gdf.loc[max(gdf.index)+1, 'splitters'] = LineString([geom, pt1])
                gdf.loc[max(gdf.index)+1, 'splitters'] = LineString([geom, pt2])

            gdf = gdf.set_geometry('splitters')
            gdf['splitters'] = gdf['splitters'].buffer(0.01)

            # Split by subtracting polygon by splitters
            diff = gpd.overlay(border, gdf, how="difference")
            diff.index = self.gdf.index
            diff['Subdivision'] = 'Offset'
            return diff
        else:
            print("ShapeTools warning: Inner court not found, offset subdivision not performed")
            return self.gdf.copy()

    def subdivide_skeleton(self, width=20, resolution=5, spine_tolerance=0.2):
        gdf = self.gdf.copy()

        # Generate skeletons and buffer to create polygons
        skeletons = Skeleton(gdf, crs=26910, resolution=resolution).generate()
        skeletons_pol = skeletons.copy()
        skeletons_pol['geometry'] = skeletons_pol.buffer(0.1)

        # Find connectivity of each line in the skeleton and filter centerlines (spine)
        skeletons_pol = get_connectivity(skeletons_pol)
        skeletons['connectivity'] = skeletons_pol['connectivity']
        skeletons = skeletons[skeletons['connectivity'] > 3]
        skeletons_pol = skeletons_pol[skeletons_pol['connectivity'] > 3]

        # Scale spine lines up intil it reaches the parcel boundary
        geometries = {}
        spines = gpd.GeoDataFrame({'id': [], 'geometry': []}, geometry='geometry')
        for pol, boundary in zip(skeletons['id'].unique(), gdf.boundary):
            spine = skeletons[skeletons['id'] == pol].unary_union

            # Get distance from endpoints to boundary vertices and divide to spine length to find scale factor
            e0 = Point(spine.bounds[0], spine.bounds[1])
            e1 = Point(spine.bounds[2], spine.bounds[3])
            max_d = max([Point(vertice).distance(ep) for vertice in boundary.coords for ep in [e0, e1]])
            fact = (max_d / spine.length) * 1.25

            # Scale up spines
            scaled = scale(spine, xfact=fact, yfact=fact, zfact=fact, origin='centroid')
            geometries[pol] = scaled
            j = len(spines)
            spines.at[j, 'id'] = pol
            spines.at[j, 'geometry'] = scaled

        # Get relative XY coordinates for spines
        if 'MultiLineString' in get_gtypes(spines):
            spines = Shape(spines).explode()
        spines = relative_max_min(spines)
        spines = self.connect_parallels(spines=spines)

        '''# Extract spine segments overlapping with parcels
		skeletons['geometry'] = [geometries[pol] for pol in skeletons['pol_id']]
		skeletons = gpd.overlay(skeletons, gdf, how="intersection")

		# Divide spines according to parcel width and connect
		skeletons = skeletons.reset_index(drop=True)
		divisions = divide_line(skeletons, length=width)
		for i, points in zip(skeletons['pol_id'], divisions):
			for pt in points:
				j = len(skeletons)
				skeletons.at[j, 'pol_id'] = i
				skeletons.at[j, 'geometry'] = LineString([nearest_points(pt, boundary)[1], pt])'''

        # Buffer edges and subtract to create polygons
        spines.crs = skeletons.crs
        spines['geometry'] = spines.buffer(0.01)
        diff = gpd.overlay(gdf, spines, how="difference")
        diff.index = gdf.index
        diff['Subdivision'] - 'Skeleton'
        return diff

    def subdivide_spine(self, width=20, spine='meridian', spine_splitter=True):
        """
		Divide geometry in two parts according to a spine dividor and connect edges back to spine according to width spacing

		:param width:
		:param spine: One of 'meridian' or 'equator'
		:return:
		"""

        gdf = self.gdf.copy()

        if len(gdf) > 0:
            # Get line in the middle of minimum rotated ractangle for each polygon
            bbox = self.min_rot_rec()
            bbox['geometry'] = bbox[spine]
            bbox.index = gdf.index
            merid = bbox.loc[:, ['id', 'geometry']]
            self.meridian.append(bbox[spine])

            # Connect middle line with edges from polygons that are relatively parallel to it
            gdf = gdf[merid.is_valid]
            merid = merid[merid.is_valid]
            merid = relative_max_min(merid)
            merid = Shape(gdf).connect_parallels(merid, width=width)

            # Include or not spine in the division
            if not spine_splitter:
                merid = merid.drop(0)

            # Buffer edges and subtract to create polygons
            merid['geometry'] = merid.buffer(0.01)
            diff = gpd.overlay(gdf, merid, how="difference")
            diff.index = gdf.index
            diff['Subdivision'] = 'BBox'
            return diff

        else:
            return gdf

    def subdivide_parallels(self, width=20, direction='uv'):
        """
		Subdivide polygons with lines parallel in a defined direction (meridian or equator)
		until the average area of parts is smaller than a defined max_area.

		:param width: Targeted area of polygons
		:param direction: One of 'u', 'v' or 'uv'
		:return:
		"""

        gdf = self.gdf.copy()

        # Get line in the middle of minimum rotated ractangle for each polygon
        bbox = self.min_rot_rec()

        overlays = gpd.GeoDataFrame(crs=gdf.crs)
        for i in gdf.index:
            geom_u = bbox.loc[i, 'meridian']
            geom_v = bbox.loc[i, 'equator']

            parts_u = gpd.GeoDataFrame({'geometry': divide_line_by_count(geom_u, int(geom_u.length / width))})
            parts_v = gpd.GeoDataFrame({'geometry': divide_line_by_count(geom_v, int(geom_v.length / width))})

            # Buffer divisions on the width of minimum rotated rectangle
            buffer_u = bbox.loc[i, 'equator'].length
            buffer_v = bbox.loc[i, 'meridian'].length

            parts_u['geometry'] = parts_u.buffer(buffer_u, join_style=2, cap_style=2)
            parts_v['geometry'] = parts_v.buffer(buffer_v, join_style=2, cap_style=2)
            parts = pd.concat([parts_u, parts_v])

            # Overlay with main shape to extract intersection
            overlays = pd.concat([overlays, gpd.overlay(parts, gdf)])

        overlays = Shape(overlays).divorce()
        overlays['geometry'] = overlays.boundary.buffer(0.01)
        return Shape(gpd.overlay(gdf, overlays, how="difference")).divorce()

    def sort_clockwise(self, pivot='centroid'):
        """
        Sort geometries clockwise around their center point
        :return:
        """
        gdf = self.gdf.copy()
        if pivot == 'centroid': pivot = gdf.unary_union.centroid
        elif pivot == 'left': pivot = gdf.unary_union

        for i in gdf.index:
            geom_ctr = gdf.loc[i, 'geometry'].centroid
            gdf.loc[i, 'angle'] = math.atan2(geom_ctr.y - pivot.y, geom_ctr.x - pivot.x)

        return gdf.sort_values('angle')

    def split_diagonal(self):
        """
        Split geometries based on diagonals of minimum rotated rectangle
        :return:
        """
        mrr = self.min_rot_rec()
        gdf = self.gdf.copy()

        for i, geom, diagonal1, diagonal2 in zip(gdf.index, gdf.geometry, mrr['diagonal1'], mrr['diagonal2']):
            d1_spt = [i.area for i in split(geom, diagonal1)]
            d2_spt = [i.area for i in split(geom, diagonal2)]
            if (max(d1_spt) - min(d1_spt)) < (max(d2_spt) - min(d2_spt)):
                gdf.at[i, 'geometry'] = MultiPolygon(split(geom, diagonal1))
            else:
                gdf.at[i, 'geometry'] = MultiPolygon(split(geom, diagonal2))
        return gdf

    def stack(self, spacing=3, n_stacks=2):
        """
		Replicates geometries from the Shape gdf according to the n_stacks and add the spacing parameter in a 'height' column

		:return:
		"""

        if type(n_stacks) in [int, float, np.float64]:
            n_stacks = [int(n_stacks) for i in self.gdf.index]

        assert len(n_stacks) == len(self.gdf), AssertionError("Size of n_stacks different from size of Shape gdf")

        gdf = self.gdf.copy()
        if True in [True if 'Multi' in geom.type else False for geom in gdf['geometry']]:
            gdf = self.divorce()
            gdf = gdf.reset_index(drop=True)

        gdf['height'] = spacing
        gdf_bounds = gdf.copy()
        gdf_bounds['geometry'] = gdf.boundary

        stacked = gpd.GeoDataFrame({
            'geometry': [geom for geom, ns in zip(gdf['geometry'], n_stacks) for j in range(ns)],
            'top_hgt': [h * (j + 1) for h, ns in zip(gdf['height'], n_stacks) for j in range(ns)],
            'parent_id': [i for i, ns in zip(gdf.index, n_stacks) for j in range(ns)]
        }, crs=self.gdf.crs)
        stacked['bottom_hgt'] = stacked['top_hgt'] - spacing
        stacked['height'] = spacing
        stacked['geometry'] = [Polygon([tuple(list(coords) + [bottom]) for coords in geom.exterior.coords])
                               for bottom, geom in zip(stacked['bottom_hgt'], stacked['geometry'])]

        """
        for i in gdf.index:
            for j in range(int(n_stacks)):
                top = gdf.loc[i, 'height'] * (j + 1)
                bottom = top - spacing

                if bottom != 0:
                    l = max(gdf.index) + 1
                    gdf.loc[l, 'bottom_height'] = bottom
                    gdf.loc[l, :] = gdf.loc[i, :]
                    # gdf.loc[l, 'height'] = top
                    # gdf.loc[l, 'geometry'] = gdf.loc[l, 'geometry'].buffer(-0.01 * (j + 1))
                    gdf.loc[l, 'geometry'] = Polygon([tuple(list(coords) + [bottom]) for coords in gdf.loc[l, 'geometry'].exterior.coords])

                else:
                    gdf.loc[i, 'bottom_height'] = bottom
        """
        return stacked


class Analyst:
    def __init__(self, left_gdf, right_gdf, network=None):
        self.left_gdf = left_gdf
        self.right_gdf = right_gdf
        self.network = network
        return

    def check_crs(self):
        assert self.left_gdf.crs == self.right_gdf.crs, \
            AttributeError(f'Left ({left_gdf.crs}) and right ({right_gdf.crs}) coordinate systems are different')
        return

    def build_network(self, rebuild=False):
        assert self.network is not None, AttributeError("Analyst object has no attribute 'network'")
        if (self.network.pdn_net is None) or rebuild:
            self.network.pdn_net = self.network.build()
        return self.network.pdn_net

    def spatial_join(self, right_suffix='', operations=None):
        """
		Join and aggregate attributes from right_gdf to gdf
		:param right_suffix:
		:param operations: List of aggregation operations, one or more of: ['sum', 'median', 'mean', 'count', 'min', 'max', 'range']
		:return:
		"""

        right_gdf = self.right_gdf.copy()
        left_gdf = drop_index(self.left_gdf).reset_index(drop=True)
        right_gdf = drop_index(right_gdf)
        if left_gdf.crs != right_gdf.crs: right_gdf = right_gdf.to_crs(left_gdf.crs)
        left_gdf['join_uid'] = list(left_gdf.reset_index(drop=True).index)

        assert left_gdf.sindex is not None and right_gdf.sindex is not None
        join = gpd.sjoin(left_gdf.loc[:, ['join_uid', 'geometry']], right_gdf, how='left').groupby('join_uid',
                                                                                                   as_index=False)
        old_columns = list(left_gdf.columns)

        # Join numerical variables
        if operations is None: operations = ['mean']
        for operation in operations:
            if operation == 'sum':
                left_joined = drop_index(join.sum())
            elif operation == 'median':
                left_joined = drop_index(join.median())
            elif operation == 'count':
                left_joined = drop_index(join.count())
            elif operation == 'min':
                left_joined = drop_index(join.min())
            elif operation == 'max':
                left_joined = drop_index(join.max())
            elif operation == 'range':
                left_joined = pd.DataFrame(drop_index(join.max())) - pd.DataFrame(drop_index(join.min())).drop(
                    'join_uid', axis=1)
            else:
                left_joined = drop_index(join.mean())
            left_joined.columns = [f'{col}_{operation}{right_suffix}'
                                   if (col in right_gdf.columns) and (col not in ['join_uid', 'geometry'])
                                   else col for col in left_joined.columns]
            left_gdf = pd.concat(
                [left_gdf, left_joined.loc[:, [col for col in left_joined.columns if col not in old_columns]]], axis=1)

        # # Join categorical variables
        # overlay = gpd.overlay(left_gdf, right_gdf)
        # for col in right_gdf.columns:
        #     if (is_string_dtype(right_gdf[col])) and (col != 'geometry'):
        #         left_gdf.loc[set(overlay['join_uid']), col] = overlay.loc[set(overlay['join_uid']), col]

        return left_gdf.sort_values('join_uid')

    def get_closest_node_id(self, left=True):
        """
		Gets the closest node of left or right GeoDataFrame to network
		:param left:
		:return:
		"""

        net = self.build_network()

        if left:
            gdf = self.left_gdf
        else:
            gdf = self.right_gdf

        gdf['node'] = net.get_node_ids(gdf.geometry.centroid.x, gdf.geometry.centroid.y)
        return gdf

    def get_distance_to_nearest(self, max_dist=1000, max_items=10, suffix='nearest'):
        """
		Gets the distance from every element in the left gdf to closest element on right gdf
		:return:
		"""

        net = self.build_network()
        net.precompute(max_dist)

        net.set_pois(suffix, maxdist=max_dist, maxitems=max_items, x_col=self.right_gdf.geometry.x, y_col=self.right_gdf.geometry.y)
        dist = net.nearest_pois(distance=max_dist, category=suffix, num_pois=max_items).drop([1], axis=1)
        dist.columns = [f'd2_{suffix}_{col}' for col in dist.columns]

        gdf = self.get_closest_node_id()
        gdf.index = gdf['node']
        gdf.loc[:, list(dist.columns)] = dist.loc[list(gdf.index), :]

        gc.collect()
        return gdf

    def convert_categorical_to_numerical(self, left=True):

        if left:
            gdf = self.left_gdf
        else:
            gdf = self.right_gdf

        for col in [col for col in gdf.columns if col != 'geometry']:
            try:
                gdf[col] = pd.to_numeric(gdf[col])
            except:
                for item in gdf[col].unique():
                    gdf.loc[gdf[col] == item, item] = 1
                    gdf.loc[gdf[col] != item, item] = 0
                gdf = gdf.drop(col, axis=1)
        return gdf

    def buffer_join_aerial(self, radius):
        """
		Spatial join features from the right_gdf located on the surroundings of the Shape
		:param radius:
		:return:
		"""

        right_gdf = self.right_gdf.copy()

        old_gdf = self.left_gdf.copy()
        self.left_gdf['geometry'] = self.left_gdf.centroid.buffer(radius)
        right_gdf.columns = [f'{col}_r{radius}' if col != 'geometry' else col for col in right_gdf.columns]
        self.left_gdf = self.spatial_join(right_gdf, operations=['count', 'mean', 'sum', 'max', 'min', 'range'])
        gdf = self.left_gdf.copy()
        self.left_gdf = old_gdf

        return gdf

    def buffer_join_network(self, radius, decay='flat', operations=None):
        """
		Aggregates data from the right_gdf located on the surroundings of the gdf elements according to a street network
		:return:
		"""
        self.check_crs()
        self.right_gdf = self.get_closest_node_id(left=False)
        right_gdf = self.convert_categorical_to_numerical(left=False)
        left_gdf = self.get_closest_node_id()

        self.build_network()
        net = self.network.pdn_net
        net.precompute(radius)

        if operations is None: operations = ['ave']

        for col in [col for col in right_gdf.columns if col != 'geometry']:
            net.set(node_ids=right_gdf["node"], variable=right_gdf[col])
            for operation in operations:
                agg = net.aggregate(distance=radius, type=operation, decay=decay)
                left_gdf[f"{col}_r{radius}_{operation}_{decay}"] = list(agg.loc[left_gdf["node"]])
            gc.collect()

        """
		# Calculate diversity index for categorical variables
		if diversity:
			for k, v in uniques.items():
				simpson = True
				shannon = True
				print(f"> Calculating diversity index for {k}")
				for decay in decays:
					categories = [f"{category}_r{radius}_sum_{decay[0]}" for category in v if
					              category != "other"]
					for category in categories:
						try:
							del sample_gdf[f"{category}_r{radius}_ave_{decay[0]}"]
						except:
							pass
						try:
							del sample_gdf[f"{category}_r{radius}_rng_{decay[0]}"]
						except:
							pass
					gc.collect()
					cat_gdf = sample_gdf.loc[:, categories]
					if simpson: sample_gdf[
						f"{k}_r{radius}_si_div_{decay[0]}"] = div.alpha_diversity('simpson',
					                                                                           cat_gdf)
					if shannon: sample_gdf[
						f"{k}_r{radius}_sh_div_{decay[0]}"] = div.alpha_diversity('shannon',
					                                                                           cat_gdf)
				gc.collect()
		gc.collect()

		# Clean count columns
		for col in sample_gdf.columns:
			if ('_ct_' in col) & ('_cnt' in col): sample_gdf = sample_gdf.drop([col], axis=1)
			if ('_ct_' in col) & ('_ave' in col): sample_gdf = sample_gdf.drop([col], axis=1)
			gc.collect()
		"""

        return left_gdf


class Generator:
    def __init__(self, gdf, spacing=10, direction=0):
        self.gdf = gdf
        self.spacing = spacing
        self.direction = direction
        return

    def open_yards(self):
        gdf = self.gdf.copy()
        min_rot_rec = Shape(gdf).min_rot_rec()
        min_rot_rec['vertices'] = [[Point(coord) for coord in geom.coords] for geom in min_rot_rec['min_rot_rec']]
        min_rot_rec['mid_points'] = [[ln.centroid for ln in geom] for geom in min_rot_rec['segments']]

        yards = gpd.GeoDataFrame()
        for i in tqdm(gdf.index):
            for j in range(1, 3):

                # Iterate over possible combinations given the list of vertices and the number of vertices in the combination
                for combination in itertools.combinations(
                        min_rot_rec.iloc[i]['vertices'] + min_rot_rec.iloc[i]['mid_points'], j):
                    for factor in [0.2, 0.4, 0.6, 0.8]:
                        parcel_geom = gdf.loc[[i], ['geometry']]
                        for vertex in combination:
                            scaled = gpd.GeoDataFrame({'geometry': [
                                scale(gdf.at[i, 'geometry'], factor, factor, factor, origin=vertex)]},
                                geometry='geometry')
                            scaled.crs = gdf.crs
                            parcel_geom = gpd.overlay(parcel_geom, scaled, how="difference")

                        if len(parcel_geom['geometry']) > 0:
                            k = len(yards)
                            yards.at[k, 'parcel_id'] = i
                            yards.at[k, 'origins'] = j
                            yards.at[k, 'scale'] = factor
                            yards.at[k, 'geometry'] = parcel_geom['geometry'].values[0]
        yards.crs = gdf.crs
        return yards


if __name__ == '__main__':
    bivariate_map(gpd.read_file('/Volumes/ELabs/50_projects/19_SSHRC_Insight/02_shading_model/Results_shp/roofs/roofs.shp'),
        columns=('height_agl', 'delta_pct'))
