# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 14:38:10 2016

@author: Jens von der Linden
"""

from datetime import date
from datetime import datetime
import numpy as np
import os
from scipy.constants import elementary_charge, proton_mass

tan = (209, 178, 111, 255)
olive = (110, 117, 14, 255)

def define_expressions(visit, alpha=8.1e5):
    r"""
    Define Visit expressions with repective alpha.
    """
    visit.DefineVectorExpression("B", "{B_x, B_y, B_z}")
    visit.DefineVectorExpression("B_norm", "{B_norm_x, B_norm_y, "
                                 "B_norm_z}")
    visit.DefineVectorExpression("B_perp", "{B_x, B_y, 0}")
    visit.DefineVectorExpression("B_para", "{0, 0, B_z}")
    visit.DefineScalarExpression("B_para_scalar", "B_z")

    visit.DefineVectorExpression("J", "{j_x, j_y, j_z}")
    visit.DefineScalarExpression("J_mag", "sqrt(j_x^2 + j_y^2 + j_z^2)")
    visit.DefineVectorExpression("J_perp", "{j_x, j_y, 0}")
    visit.DefineVectorExpression("J_para", "{0, 0, j_z}")
    visit.DefineScalarExpression("J_para_scalar", "j_z")

    visit.DefineVectorExpression("u_e_norm", "{u_e_norm_x, u_e_norm_y, "
                                 "u_e_norm_z}")
    visit.DefineVectorExpression("u_e", "%1.2e * u_e_norm" % alpha)
    visit.DefineVectorExpression("u_e_perp", "{dot(u_e, {1, 0, 0}), "
                                 "dot(u_e, {0, 1, 0}), 0}")
    visit.DefineScalarExpression("u_e_para_scalar", "dot(u_e, {0, 0, 1})")


    visit.DefineVectorExpression("u_i_term1", "{u_i_term1_x, u_i_term1_y, "
                                 "u_i_term1_z}")
    visit.DefineVectorExpression("u_i_term1_perp", "{u_i_term1_x, "
                                 "u_i_term1_y, 0}")
    visit.DefineVectorExpression("u_i_term1_para", "{0, 0, u_i_term1_z}")
    visit.DefineScalarExpression("u_i_term1_para_scalar", "u_i_term1_z")

    visit.DefineVectorExpression("omega_i_term1", "{w_i_term1_x, w_i_term1_y, "
                                 "w_i_term1_z}")
    visit.DefineVectorExpression("omega_i_term1_perp", "{w_i_term1_x, "
                                 "w_i_term1_y, 0}")
    visit.DefineVectorExpression("omega_i_term1_para", "{0, 0, w_i_term1_z}")
    visit.DefineScalarExpression("omega_i_term1_para_scalar", "w_i_term1_z")

    visit.DefineVectorExpression("omega_i_term2", "{w_i_term2_x, w_i_term2_y, "
                                 "w_i_term2_z} * %1.2e" % alpha)
    visit.DefineVectorExpression("omega_i_term2_perp", "{w_i_term2_x, "
                                 "w_i_term2_y, 0} * %1.2e" % alpha)
    visit.DefineVectorExpression("omega_i_term2_para", "{0, 0, "
                                 "w_i_term2_z} * %1.2e" % alpha)
    visit.DefineScalarExpression("omega_i_term2_para_scalar",
                                 "w_i_term2_z * %1.2e" % alpha)

    visit.DefineVectorExpression("u_i", "u_i_term1 + u_e")
    visit.DefineVectorExpression("u_i_perp", "{dot(u_i, {1, 0, 0}), dot(u_i, "
                                 "{0, 1, 0}), 0}")
    visit.DefineVectorExpression("u_i_para", "{0, 0, dot(u_i, {0, 0, 1})}")
    visit.DefineScalarExpression("u_i_para_scalar", "dot(u_i, {0, 0, 1})")

    visit.DefineVectorExpression("omega_i", "omega_i_term1 + "
                                 "omega_i_term2 * %1.2e" % alpha)
    visit.DefineVectorExpression("omega_i_perp", "{dot(omega_i, {1, 0, 0}), "
                                 "dot(omega_i, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("omega_i_para", "{0, 0, dot(omega_i, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("omega_i_para_scalar", "dot(omega_i, "
                                 "{0, 0, 1})")

    visit.DefineVectorExpression("Omega_i", "B +" + str(proton_mass) + "/" +
                                 str(elementary_charge) + "*omega_i")
    visit.DefineVectorExpression("Omega_i_perp", "{dot(Omega_i, {1, 0, 0}), "
                                 "dot(Omega_i, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("Omega_i_para", "{0, 0, dot(Omega_i, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("Omega_i_para_scalar",
                                 "dot(Omega_i, {0, 0, 1})")


def launch_points(center, plane=0.249, num_inner=80, num_outer=60):
    r"""
    Calculate points on a circle outline for a given center point.
    """
    thetas = np.linspace(0, 3./4.*np.pi, num_outer)
    thetas = np.concatenate((thetas, np.linspace(5./4.*np.pi, 2.*np.pi,
                                                 num_outer)))
    x_outer = 0.01 * np.cos(thetas) + center[0]
    y_outer = 0.01 * np.sin(thetas) + center[1]
    z_outer = plane * np.ones((x_outer.size))
    points_outer = np.empty((x_outer.size + y_outer.size + z_outer.size))
    points_outer[0::3] = x_outer
    points_outer[1::3] = y_outer
    points_outer[2::3] = z_outer
    points_outer = tuple(points_outer)

    thetas = np.linspace(0, 2.*np.pi, num_inner)
    x_inner = 0.005 * np.cos(thetas) + center[0]
    y_inner = 0.005 * np.sin(thetas) + center[1]
    z_inner = plane * np.ones(x_inner.size)
    points_inner = np.empty((x_inner.size + y_inner.size + z_inner.size))
    points_inner[0::3] = x_inner
    points_inner[1::3] = y_inner
    points_inner[2::3] = z_inner
    points_inner = tuple(points_inner)

    return points_outer, points_inner


def setup_current_pseudocolor(visit, colortable="Greens"):
    r"""
    Setup pseudocolor current plot.
    """
    visit.AddPlot("Pseudocolor", "J_mag", 1, 0)
    PseudocolorAtts = visit.PseudocolorAttributes()
    PseudocolorAtts.scaling = PseudocolorAtts.Linear
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData
    PseudocolorAtts.colorTableName = colortable
    visit.SetPlotOptions(PseudocolorAtts)
    visit.AddOperator("Slice", 0)
    SliceAtts = visit.SliceAttributes()
    SliceAtts.originType = SliceAtts.Intercept
    SliceAtts.originIntercept = 0.249
    SliceAtts.axisType = SliceAtts.ZAxis
    SliceAtts.project2d = 0
    visit.SetOperatorOptions(SliceAtts, 0)
    return PseudocolorAtts, SliceAtts


def setup_massless_electron_canonical_flux_tubes(visit, points_outer,
                                                 points_inner):
    r"""
    Setup two massless electron canonical flux tubes i.e. magnetic flux tubes.
    Inteteded to be inner and outer flux tubes.
    """
    visit.AddPlot("Streamline", "B", 1, 0)
    StreamlineAtts_outer = visit.StreamlineAttributes()
    StreamlineAtts_outer.sourceType = StreamlineAtts_outer.SpecifiedPointList
    StreamlineAtts_outer.SetPointList(points_outer)
    StreamlineAtts_outer.coloringMethod = StreamlineAtts_outer.Solid
    StreamlineAtts_outer.colorTableName = "Default"
    StreamlineAtts_outer.singleColor = (255, 0, 0, 255)
    StreamlineAtts_outer.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_outer)

    visit.AddPlot("Streamline", "B", 1, 0)
    StreamlineAtts_inner = visit.StreamlineAttributes()
    StreamlineAtts_inner.sourceType = StreamlineAtts_inner.SpecifiedPointList
    StreamlineAtts_inner.SetPointList(points_inner)
    StreamlineAtts_inner.coloringMethod = StreamlineAtts_inner.Solid
    StreamlineAtts_inner.colorTableName = "Default"
    StreamlineAtts_inner.singleColor = (190, 64, 0, 255)
    StreamlineAtts_inner.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_inner)

    return StreamlineAtts_outer, StreamlineAtts_inner


def setup_inner_outer_ion_canonical_flux_tubes(visit, points_outer,
                                               points_inner,
                                               outer_color=tan,
                                               inner_color=olive):
    r"""
    Setup two ion canonical flux tubes.
    Inteteded to be inner and outer flux tubes.
    """
    visit.AddPlot("Streamline", "Omega_i", 1, 0)
    StreamlineAtts_outer = visit.StreamlineAttributes()
    StreamlineAtts_outer.sourceType = StreamlineAtts_outer.SpecifiedPointList
    StreamlineAtts_outer.SetPointList(points_outer)
    StreamlineAtts_outer.coloringMethod = StreamlineAtts_outer.Solid
    StreamlineAtts_outer.colorTableName = "Default"
    StreamlineAtts_outer.singleColor = outer_color
    StreamlineAtts_outer.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_outer)

    visit.AddPlot("Streamline", "Omega_i", 1, 0)
    StreamlineAtts_inner = visit.StreamlineAttributes()
    StreamlineAtts_inner.sourceType = StreamlineAtts_inner.SpecifiedPointList
    StreamlineAtts_inner.SetPointList(points_inner)
    StreamlineAtts_inner.coloringMethod = StreamlineAtts_inner.Solid
    StreamlineAtts_inner.colorTableName = "Default"
    StreamlineAtts_inner.singleColor = inner_color
    StreamlineAtts_inner.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_inner)

    return StreamlineAtts_outer, StreamlineAtts_inner


def setup_forward_backward_ion_canonical_flux_tubes(visit, points_foward,
                                                    points_backward,
                                                    forward_color=tan,
                                                    backward_color=olive):
    r"""
    Setup two ion canonical flux tubes, one integrating in the forward
    direction, one integrating in the backward direction.
    """
    visit.AddPlot("Streamline", "Omega_i", 1, 0)

    StreamlineAtts_forward = visit.StreamlineAttributes()
    StreamlineAtts_forward.sourceType = StreamlineAtts_forward.SpecifiedPointList
    StreamlineAtts_forward.SetPointList(points_foward)
    StreamlineAtts_forward.coloringMethod = StreamlineAtts_forward.Solid
    StreamlineAtts_forward.colorTableName = "Default"
    StreamlineAtts_forward.singleColor = forward_color
    StreamlineAtts_forward.integrationDirection = StreamlineAtts_forward.Forward
    StreamlineAtts_forward.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_forward)

    visit.AddPlot("Streamline", "Omega_i", 1, 0)
    StreamlineAtts_backward = visit.StreamlineAttributes()
    StreamlineAtts_backward.sourceType = StreamlineAtts_backward.SpecifiedPointList
    StreamlineAtts_backward.SetPointList(points_backward)
    StreamlineAtts_backward.coloringMethod = StreamlineAtts_backward.Solid
    StreamlineAtts_backward.colorTableName = "Default"
    StreamlineAtts_backward.singleColor = backward_color
    StreamlineAtts_backward.integrationDirection = StreamlineAtts_backward.Backward
    StreamlineAtts_backward.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_backward)

    return StreamlineAtts_forward, StreamlineAtts_backward


def setup_field_line(visit, center=(0.01, 0.01, 0.249),
                     outer_radius=0.01):
    visit.AddPlot("Streamline", "B", 1, 0)
    StreamlineAtts_line = visit.StreamlineAttributes()
    StreamlineAtts_line.sourceType = StreamlineAtts_line.SpecifiedPoint
    StreamlineAtts_line.pointSource = (center[0], center[1] + outer_radius, center[2])
    StreamlineAtts_line.coloringMethod = StreamlineAtts_line.Solid
    StreamlineAtts_line.singleColor = (255, 255, 153, 255)
    StreamlineAtts_line.legendFlag = 0
    StreamlineAtts_line.showSeeds = 0
    visit.SetPlotOptions(StreamlineAtts_line)
    return StreamlineAtts_line


def setup_annotations(visit, time_scale=1):
    r"""
    """
    AnnotationAtts = visit.AnnotationAttributes()
    AnnotationAtts.axes3D.autoSetScaling = 0
    AnnotationAtts.axes3D.xAxis.title.visible = 0
    AnnotationAtts.axes3D.yAxis.title.visible = 0
    AnnotationAtts.axes3D.zAxis.title.visible = 0
    AnnotationAtts.axes3D.xAxis.label.font.scale = 3
    AnnotationAtts.axes3D.xAxis.label.scaling = -2
    AnnotationAtts.axes3D.yAxis.label.font.scale = 3
    AnnotationAtts.axes3D.yAxis.label.scaling = -2
    AnnotationAtts.axes3D.zAxis.label.font.scale = 3
    AnnotationAtts.axes3D.zAxis.label.scaling = -2
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.databaseInfoFlag = 0
    AnnotationAtts.databaseInfoTimeScale = time_scale
    AnnotationAtts.databaseInfoTimeOffset = 0
    visit.SetAnnotationAttributes(AnnotationAtts)
    return AnnotationAtts


def set_default_view(visit):
    r"""
    Set default view for viewing fluxtubes.
    If view needs to be modified it is best to align visit with gui and save
    parameters from visit.GetView3D().
    """
    view = visit.GetView3D()
    view.SetViewNormal((-0.731293, 0.40847, 0.546227))
    view.SetFocus((0.00202222, 0.000976744, 0.331997))
    view.SetViewUp((0.322268, 0.91274, -0.251095))
    view.SetViewAngle(30)
    view.SetParallelScale(0.088383)
    view.SetNearPlane(-0.176766)
    view.SetImagePan((0, 0))
    view.SetImageZoom(1)
    view.SetPerspective(1)
    view.SetEyeAngle(2)
    view.SetCenterOfRotationSet(0)
    view.SetCenterOfRotation((0.00202222, 0.000976744, 0.331997))
    view.SetAxis3DScaleFlag(0)
    view.SetAxis3DScales((1,1,1))
    view.SetShear((0,0,1))
    view.SetWindowValid(0)
    visit.SetView3D(view)