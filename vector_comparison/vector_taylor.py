## Author: Jens von der Linden 2017-01-26
## Started with Yannick Copin's public domain code: https://gist.github.com/ycopin/3342888


import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../vector_calculus')
from vector_calculus import vector_calculus as vc
from matplotlib.ticker import MultipleLocator, MaxNLocator

class TaylorDiagram(object):
    """Taylor diagram: plot model standard deviation and correlation
    to reference (data) sample in a single-quadrant polar plot, with
    r=stddev and theta=arccos(correlation).

    Example
    -------
    dia = TaylorDiagram(refstd, fig=fig, rect=111, label="Reference")

    colors = plt.matplotlib.cm.jet(np.linspace(0,1,len(samples)))

    # Add samples to Taylor diagram
    for i,(stddev,corrcoef) in enumerate(samples):
        dia.add_sample(stddev, corrcoef, marker='s', ls='', c=colors[i],
                       label="Model %d" % (i+1))

    # Add RMS contours, and label them
    contours = dia.add_contours(colors='0.5')
    plt.clabel(contours, inline=1, fontsize=10)
    """

    def __init__(self, refstd, fig=None, rect=111, label='_',
                 std_multiplier=1.5, axes=None,
                 multiple_locator=None, r_locations=None, r_label_pos=(.3, 0.02),
                 noref=False):
        """Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using mpl_toolkits.axisartist.floating_axes. refstd is
        the reference standard deviation to be compared to.
        """

        from matplotlib.projections import PolarAxes
        import mpl_toolkits.axisartist.floating_axes as fa
        import mpl_toolkits.axisartist.grid_finder as gf

        self.refstd = refstd            # Reference standard deviation

        tr = PolarAxes.PolarTransform()

        # Correlation labels
        rlocs = np.concatenate(([-0.99, -0.95, -0.9],
                                np.arange(-8, 9, 2)/10.,[0.9, 0.95, 0.99]))
        tlocs = np.arccos(rlocs)        # Conversion to polar angles
        gl1 = gf.FixedLocator(tlocs)    # Positions
        tf1 = gf.DictFormatter(dict(zip(tlocs, map(str,rlocs))))

        # Standard deviation axis extent
        self.smin = 0.
        self.smax = std_multiplier*self.refstd

        if r_locations is None:
            ghelper = fa.GridHelperCurveLinear(tr,
                                               extremes=(0, np.pi,
                                                         self.smin,self.smax),
                                               grid_locator1=gl1,
                                               tick_formatter1=tf1)
        else:
            r_grid_locator = gf.FixedLocator(r_locations)
            r_grid_labels = gf.DictFormatter(dict(zip(r_locations,
                                                      map(str, r_locations))))
            ghelper = fa.GridHelperCurveLinear(tr,
                                               extremes=(0, np.pi,
                                                         self.smin,self.smax),
                                               grid_locator1=gl1,
                                               tick_formatter1=tf1,
                                               grid_locator2=r_grid_locator,
                                               tick_formatter2=r_grid_labels)

        if (axes is None):
            if fig is None:
                fig = plt.figure()

            ax = fa.FloatingSubplot(fig, rect, grid_helper=ghelper)
            fig.add_subplot(ax)
        else:
            ax = fa.FloatingSubplot(fig, rect, grid_helper=ghelper)
            axes = ax


        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")  # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text("Normalized helicity density sum")

        ax.axis["left"].set_axis_direction("bottom") # "right X axis
        #ax.axis["left"].label.set_text("Standard deviation")
        #ax.axis["left"].major_ticklabels.set_axis_direction("left")

        ax.axis["right"].set_axis_direction("top")   # "Left X axis"
        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction("bottom")

        ax.axis["bottom"].set_visible(False) # Useless

        if multiple_locator:
            majorLocator  = MultipleLocator(multiple_locator)
            ax.axis["left"].set_major_locator(majorLocator)

        # Contours along standard deviations
        ax.grid(True, which='major', linestyle='-', alpha=0.3, color='grey')
        # Use this to edit gridlines:
        gridlines = ax.get_xgridlines() + ax.get_ygridlines()
        for gridline in gridlines:
            gridline.set_alpha=0.3

        self._ax = ax                   # Graphical axes
        self.ax = ax.get_aux_axes(tr)   # Polar coordinates

        # Add reference point and stddev contour
        if not noref:
            l, = self.ax.plot([0], self.refstd, 'k*',
                              ls='', ms=10, label=label)
            t = np.linspace(0, np.pi)
            r = np.zeros_like(t) + self.refstd
            self.ax.plot(t,r, 'k--', label='_')
            self.samplePoints = [l]
        else:
            self.samplePoints = []

        plt.figtext(r_label_pos[0], r_label_pos[1],
                    'Root mean square vector length')

    def add_sample(self, stddev, corrcoef, *args, **kwargs):
        """Add sample (stddev,corrcoeff) to the Taylor diagram. args
        and kwargs are directly propagated to the Figure.plot
        command."""
        angle = np.arccos(corrcoef)
        if np.allclose(corrcoef, 1.0):
            angle = 0.
        l, = self.ax.plot(angle, stddev,
                          *args, **kwargs) # (theta,radius)
        self.samplePoints.append(l)

        return l

    def add_contours(self, levels=5, **kwargs):
        """Add constant centered RMS difference contours."""

        rs,ts = np.meshgrid(np.linspace(self.smin,self.smax),
                            np.linspace(0,np.pi))
        # Compute centered RMS difference
        rms = np.sqrt(self.refstd**2 + rs**2 - 2*self.refstd*rs*np.cos(ts))
        contours = self.ax.contour(ts, rs, rms, levels, **kwargs)

        return contours

def vector_simililarity_coefficient(ref_field, field):
    r"""
    """
    assert ref_field[0].shape == field[0].shape, 'Vector fields do not have the same dimensions.'
    dot_product = vc.dot_product(ref_field, field)
    numerator = np.sum(dot_product)
    denominator = (np.sqrt(np.sum(vc.magnitude(ref_field)**2))*
                   np.sqrt(np.sum(vc.magnitude(field)**2)))
    return numerator / denominator


def root_mean_square_lenth(field):
    r"""
    """
    size = field[0].size
    return np.sqrt(1./size*np.sum(vc.magnitude(field)**2))


def root_mean_square_vector_difference(ref_field, field):
    r"""
    """
    assert ref_field[0].shape == field[0].shape, 'Vector fields do not have the same dimensions.'
    assert 4 > len(ref_field) > 1, 'Vectors should have at least 2 no more then 3 components.'
    difference_x = field[0] - ref_field[0]
    difference_y = field[1] - ref_field[1]
    if len(field) == 3:
        difference_z = field[2] - ref_field[2]
        difference_mag = vc.magnitude([difference_x, difference_y, difference_z])
    else:
        difference_mag = vc.magnitude([difference_x, difference_y])
    return np.sqrt(1/ref_field[0].size*np.sum(difference_mag**2))


def calc_and_plot(ref_field, fields, fig=None,
                  std_multiplier=1.5, labels=None,
                  colors=None, markers=None, axes=None,
                  bbox_to_anchor=[0.1, 1.2], multiple_locator=None,
                  r_locations=None, r_label_pos=(0.3, 0.02),
                  norm_rmsl=False, noref=False):
    r"""
    """
    ref_rmsl = root_mean_square_lenth(ref_field)
    if norm_rmsl:
        diagram = TaylorDiagram(1., fig=fig,
                                std_multiplier=std_multiplier, axes=axes,
                                multiple_locator=multiple_locator,
                                r_locations=r_locations, r_label_pos=r_label_pos,
                                noref=noref)
    else:
        diagram = TaylorDiagram(ref_rmsl, fig=fig,
                                std_multiplier=std_multiplier, axes=axes,
                                multiple_locator=multiple_locator,
                                r_locations=r_locations, r_label_pos=r_label_pos,
                                noref=noref)

    for i, field in enumerate(fields):
        similarity = vector_simililarity_coefficient(ref_field, field)
        rmsl = root_mean_square_lenth(field)
        if norm_rmsl:
            rmsl = rmsl/ref_rmsl
        rmsvds = root_mean_square_vector_difference(ref_field, field)
        if not colors:
            colors = ['red',] * len(fields)
        if not markers:
            markers = ['s',] * len(fields)
        if labels:
            diagram.add_sample(rmsl, similarity, marker=markers[i],
                               ls='', c=colors[i], label=labels[i])
        else:
            diagram.add_sample(rmsl, similarity, marker=markers[i],
                               ls='', c=colors[i])
    diagram.ax.legend(bbox_to_anchor=bbox_to_anchor,
                      loc='upper center',
                      borderaxespad=0.25)
    #diagram.add_contours(levels=10)
    return diagram
