import geopandas as gpd
import pandas as pd

from rtree import index
from shapely.affinity import scale
from shapely.geometry import LineString, Point
from shapely.geometry import Point, LineString
from shapely.ops import snap

pd.options.display.max_columns = 15
pd.options.display.max_rows = 6
pd.options.display.max_colwidth = 40
pd.options.display.width = 300

def get_connectivity(gdf):
    gdf['id'] = gdf.index

    connectivity = []
    idx = index.Index()

    # Populate R-tree index with bounds of grid cells
    for pos, cell in enumerate(gdf.geometry):
        # assuming cell is a shapely object
        idx.insert(pos, cell.bounds)

    for i, pol in enumerate(gdf.geometry):

        # Merge cells that have overlapping bounding boxes
        potential_conn = [gdf.id[pos] for pos in idx.intersection(pol.bounds)]

        # Now do actual intersection
        conn = []
        for j in potential_conn:
            if gdf.loc[j, 'geometry'].intersects(pol):
                conn.append(j)
        connectivity.append(len(conn))

    gdf['connectivity'] = connectivity
    return gdf


class Streets:
    def __init__(self, gdf, buildings=None, crs=26910, widths=None, trees=None, intersections=None, verbose=False):
        self.gdf = gdf.to_crs(crs)
        self.trees = trees
        self.widths = widths
        self.crs = crs
        if buildings is not None: self.barriers = buildings.gdf
        if verbose: print("Streets object created")
        self.verbose = verbose
        self.intersections = intersections
        return

    def extract_intersections(self):
        gdf = self.gdf.copy().reset_index(drop=True)
        gdf['geometry'] = [geom.simplify(0.2) for geom in gdf['geometry']]
        vertices = Shape(gdf).extract_vertices()

        # For each line check intersection with all other lines
        for i in gdf.index:
            all_net = gdf.loc[[j for j in gdf.index if j != i], :].unary_union
            line = gdf.loc[i, 'geometry']
            inters = line.intersection(all_net)
            if inters.type == 'MultiPoint':
                for pt in inters:
                    k = len(vertices)
                    vertices.at[k, 'geometry'] = pt
            elif inters.type == 'Point':
                vertices.at[len(vertices), 'geometry'] = inters

        vertices.crs = gdf.crs
        return vertices.drop_duplicates('geometry')

    def buffer(self, width_column='width'):
        """
        Buffer street network according to a defined width_column
        :param width_column:
        :return:
        """
        gdf = self.gdf.copy().reset_index(drop=True)
        for w in gdf[width_column].unique():
            gdf.loc[gdf[width_column] == w, 'geometry'] = gdf.loc[gdf[width_column] == w, 'geometry'].buffer(w/2, cap_style=2)
        gdf['Area'] = gdf.area
        return gdf

    def dimension(self):
        gdf = self.gdf.copy().reset_index(drop=True)

        print("> Cleaning street widths")
        if self.widths is not None:
            widths = gpd.read_file(self.widths)
            widths = widths.to_crs(self.crs)

            # Clean data from City of Vancouver open data catalogue
            for i in widths.index:
                if '(m)' in str(widths.at[i, 'width']):
                    try: widths.at[i, 'width'] = float(widths.loc[i, 'width'].split('(m)')[0])
                    except: widths = widths.drop(i)
                elif 'D.L.' in str(widths.at[i, 'width']):
                    widths.at[i, 'width'] = float(widths.loc[i, 'width'].split('D.L.')[1])
                elif '±' in str(widths.at[i, 'width']):
                    widths.at[i, 'width'] = float(widths.loc[i, 'width'].split('±')[0])
                elif 'm' in str(widths.at[i, 'width']):
                    try: widths.at[i, 'width'] = float(widths.loc[i, 'width'].split('m')[0])
                    except: widths = widths.drop(i)
                elif 'R' in str(widths.at[i, 'width']):
                    try: widths.at[i, 'width'] = float(widths.loc[i, 'width'].split('R')[0])
                    except: widths = widths.drop(i)
                elif 'M' in str(widths.at[i, 'width']):
                    widths.at[i, 'width'] = float(widths.loc[i, 'width'].split('M')[0])
                elif widths.at[i, 'width'] in ['-', '+', 'CHD.', 'ST', 'ST.', None]: widths = widths.drop(i)
                else: widths.at[i, 'width'] = float(widths.loc[i, 'width'])

            # Buffer street segments based on average street width
            widths = widths[widths['width'] < 100]
            widths.crs = self.crs
            widths.to_file(f'{directory}/row.shp', driver='ESRI Shapefile')
            widths['geometry'] = widths.buffer(10)

            gdf['id'] = gdf.index
            joined = gpd.sjoin(gdf, widths, how='left')
            joined['width'] = pd.to_numeric(joined['width'])
            joined = joined.groupby('id', as_index=False).mean()
            joined = pd.merge(gdf, joined, on='id', copy=False)
            joined['geometry'] = list(gdf.loc[gdf['id'].isin(joined['id'])]['geometry'])

            # Replace NaN values
            print(f'Width information from {joined["width"].isna().sum()} features could not be joined')
            for use in joined['streetuse'].unique():
                joined.loc[(joined['streetuse'] == use) & (joined['width'].isna()), 'width'] = joined[joined['streetuse'] == use]['width'].mean()

            joined['length'] = [geom.length for geom in joined['geometry']]
            return joined
        else:
            print('Widths layer not found!')
            return gdf

    def direction(self):
        gdf = self.gdf.copy()

        # Calculate street segment in relation to shortest line from start to end of segment
        print("> Calculating shortest path")
        shortest = []
        for i, segment in zip(gdf.index, gdf['geometry']):
            if segment.__class__.__name__ == 'LineString':
                shortest.append(LineString([Point(segment.coords[0]), Point(segment.coords[1])]).length)
            elif segment.__class__.__name__ == 'MultiLineString':
                shortest.append(LineString([Point(segment[0].coords[0]), Point(segment[0].coords[1])]).length)
                gdf.at[i, 'geometry'] = segment[0]

        gdf['shortest'] = shortest

        # Calculate straightness
        print("> Calculating straightness")
        gdf['straight'] = gdf['shortest']/gdf['length']

        # Calculate azimuth
        print("> Calculating azimuth")
        gdf['azimuth'] = [azimuth(ln) for ln in gdf.geometry]

        # Get one way streets from OpenStreetMap

        return gdf

    def connectivity(self):
        gdf = self.gdf.copy().reset_index()
        gdf['id'] = gdf.index
        buffer = gpd.GeoDataFrame(gdf.buffer(1), columns=['geometry'], crs=26910)
        buffer.to_file(f'{directory}/street_segs_buffer.geojson', drive='GeoJSON')

        print(f"> Calculating connections between {buffer.sindex} and {gdf.sindex}")
        overlay = gpd.overlay(gdf, buffer, how='intersection')
        gdf['conn'] = [len(overlay[overlay['id'] == i]) for i in gdf['id'].unique()]
        gdf['deadend'] = [1 if gdf.loc[i, 'conn'] < 2 else 0 for i in gdf.index]
        return gdf

    def visibility(self):
        if self.barriers is not None:
            gdf = self.gdf.copy()
            c_gdf = self.gdf.copy()
            c_gdf.crs = self.crs

            print("> Generating isovists")
            c_gdf['geometry'] = c_gdf.centroid
            gdf['isovists'] = Isovist(origins=c_gdf, barriers=self.barriers).create()
            gdf['iso_area'] = [geom.area for geom in gdf['isovists']]
            gdf['iso_perim'] = [geom.length for geom in gdf['isovists']]
            gdf.drop(['geometry'], axis=1).set_geometry('isovists').to_file(f'{directory}/isovists.geojson', driver='GeoJSON')
            return gdf.drop(['isovists'])

    def greenery(self):
        gdf = self.gdf.copy()
        if self.trees is not None:
            trees = gpd.read_file(self.trees)
            trees = trees.dropna(subset=['geometry'])
            trees.crs = self.crs
            print(trees.sindex)

            # Buffer streets geometries
            if 'width' not in gdf.columns:
                buffer = gpd.GeoDataFrame(gdf.buffer(10), columns=['geometry'], crs=26910)
            else:
                buffer = gpd.GeoDataFrame([g.buffer(w/2) for g, w in zip(gdf['geometry'], gdf['width'])], columns=['geometry'], crs=26910)
            buffer['id'] = buffer.index
            buffer.crs = self.crs
            print(buffer.sindex)

            # Overlay trees and street buffers
            overlay = gpd.overlay(trees, buffer, how='intersection')
            gdf['n_trees'] = [len(overlay[overlay['id'] == i]) for i in buffer['id'].unique()]
            gdf.crs = self.crs
        else: print("Trees layer not found!")
        return gdf

    def all(self):
        print("> Calculating all indicators for Streets")
        self.gdf = self.dimension()
        self.gdf = self.direction()
        self.gdf = self.connectivity()

        try:
            iso_gdf = gpd.read_file(f'{directory}/isovists.geojson').loc[:, ['iso_area', 'iso_perim']]
            iso_gdf['id'] = iso_gdf.index
            self.gdf['id'] = self.gdf.index
            self.gdf = pd.merge(self.gdf, iso_gdf, how="left", copy=False, on='id')
        except:
            print("> Isovists not found")
            self.gdf = self.visibility()

        self.gdf = self.greenery()
        self.gdf.crs = self.crs
        return self.gdf

    def get_length(self):
        if 'length' not in self.gdf.columns:
            self.gdf['length'] = self.gdf.length
            return self.gdf

    def calculate_row(self, parcels_gdf, max_buffer=15):
        gdf = self.get_length().copy()
        gdf_buffer = gdf.copy()
        gdf_buffer['geometry'] = gdf_buffer.buffer(max_buffer, cap_style=2)
        gdf_buffer = gpd.overlay(gdf_buffer, parcels_gdf.loc[:, ['geometry']], how='difference')
        gdf['row'] = gdf_buffer.area/gdf_buffer['length']
        return gdf

    def segmentize(self, min_length=5):
        """
        Divide street network into segments based on cross-intersections of lines
        :return:
        """
        inters = self.extract_intersections()
        inters = Shape(gpd.GeoDataFrame({'geometry': [inters.buffer(0.1).unary_union]}, crs=self.crs)).divorce()
        inters['geometry'] = inters.centroid
        inters['geometry'] = inters.buffer(0.1)
        gdf = self.gdf.copy()
        gdf = Shape(gpd.overlay(gdf, inters.loc[:, ['geometry']], how='difference')).divorce()
        gdf = gdf[gdf.length > min_length]

        # Snap lines to intersections
        inters_centroids = inters.buffer(1).centroid.unary_union
        gdf['geometry'] = [snap(geom, inters_centroids, 1) for geom in gdf['geometry']]
        gdf['length'] = gdf.length
        return gdf


class Axial:
    def __init__(self, buffer=7, tolerance=15):
        self.buffer = buffer
        self.s_tol = tolerance
        self.connectivity = []
        self.gdf = gpd.GeoDataFrame()
        return

    def gdf_from_segments(self, links):
        """
        :param links: GeoDataFrame with LineString geometry
        :return: GeoDataFrame with axial polygons
        """

        # Simplify links geometry
        s_tol = self.s_tol
        links.loc[:, 'geometry'] = links.simplify(s_tol)

        print("> Processing axial map, exploding polylines")
        n_links = []
        for ln in links.geometry:
            if len(ln.coords) > 2:
                for i, coord in enumerate(ln.coords):
                    if i == len(ln.coords) - 1:
                        pass
                    else:
                        n_links.append(LineString([Point(coord), Point(ln.coords[i + 1])]))
            else:
                n_links.append(ln)

        print("> Lines exploded, calculating azimuth")
        links = gpd.GeoDataFrame(n_links, columns=['geometry'])
        links = links.set_geometry('geometry')
        links = azimuth(links)
        links = links.reset_index()

        n_clusters = 6
        print(f"> Clustering segments into {n_clusters} clusters based on azimuth")
        bgm = GaussianMixture(
            n_components=n_clusters,
            # weight_concentration_prior_type='dirichlet_process',
            # weight_concentration_prior=0.001
        )
        to_cluster = links.loc[:, ['azimuth', 'azimuth_n']]
        bgm.fit(to_cluster)
        links['axial_labels'] = bgm.predict(to_cluster)

        print("> Isolating geometries to create axial-like lines")
        clusters = [links.loc[links.axial_labels == i] for i in links['axial_labels'].unique()]

        print(f"> Buffering and iterating over multi geometries (!!!)")
        buffer_r = self.buffer
        m_pols = [df.buffer(buffer_r).unary_union for df in clusters]
        geoms = []
        for m_pol in m_pols:
            if m_pol.__class__.__name__ == 'Polygon':
                geoms.append(m_pol)
            else:
                for pol in m_pol:
                    geoms.append(pol)

        print("> Creating axial GeoDataFrame")
        axial_gdf = gpd.GeoDataFrame(geometry=geoms)
        axial_gdf = axial_gdf.reset_index()
        axial_gdf.to_feather('data/axial.feather')
        return axial_gdf

    def create_graph(self):
        axial_gdf = self.gdf

        g = Graph(directed=True)
        for i, pol in enumerate(axial_gdf.geometry):
            v = g.add_vertex()
            v.index = i

        axial_gdf['id'] = axial_gdf.index
        axial_gdf['axial_length'] = axial_gdf.area / self.buffer

        print("> Creating spatial index")
        idx = index.Index()
        # Populate R-tree index with bounds of grid cells
        for pos, cell in enumerate(axial_gdf.geometry):
            # assuming cell is a shapely object
            idx.insert(pos, cell.bounds)

        print("> Finding connected lines")
        g_edges = []
        # Loop through each Shapely polygon (axial line)
        for i, pol in enumerate(axial_gdf.geometry):

            # Merge cells that have overlapping bounding boxes
            potential_conn = [axial_gdf.id[pos] for pos in idx.intersection(pol.bounds)]

            # Now do actual intersection
            conn = []
            for j in potential_conn:
                if axial_gdf.loc[j, 'geometry'].intersects(pol):
                    conn.append(j)
                    g_edges.append([i, j, 1])  # axial_gdf.loc[i, 'length']/1000])
            self.connectivity.append(len(conn))
            print(f"> Finished adding edges of axial line {i + 1}/{len(axial_gdf)} to dual graph")
        return
