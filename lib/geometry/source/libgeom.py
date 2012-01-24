# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _libgeom

def _swig_setattr(self,class_type,name,value):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    self.__dict__[name] = value

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


POINT_3D_LABEL_SIZE = _libgeom.POINT_3D_LABEL_SIZE
class Vector2D(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vector2D, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Vector2D, name)
    def __repr__(self):
        return "<C Vector2D instance at %s>" % (self.this,)
    __swig_setmethods__["x"] = _libgeom.Vector2D_x_set
    __swig_getmethods__["x"] = _libgeom.Vector2D_x_get
    if _newclass:x = property(_libgeom.Vector2D_x_get, _libgeom.Vector2D_x_set)
    __swig_setmethods__["y"] = _libgeom.Vector2D_y_set
    __swig_getmethods__["y"] = _libgeom.Vector2D_y_get
    if _newclass:y = property(_libgeom.Vector2D_y_get, _libgeom.Vector2D_y_set)
    __swig_setmethods__["code"] = _libgeom.Vector2D_code_set
    __swig_getmethods__["code"] = _libgeom.Vector2D_code_get
    if _newclass:code = property(_libgeom.Vector2D_code_get, _libgeom.Vector2D_code_set)
    def __init__(self, *args):
        _swig_setattr(self, Vector2D, 'this', _libgeom.new_Vector2D(*args))
        _swig_setattr(self, Vector2D, 'thisown', 1)
    def __del__(self, destroy=_libgeom.delete_Vector2D):
        try:
            if self.thisown: destroy(self)
        except: pass

class Vector2DPtr(Vector2D):
    def __init__(self, this):
        _swig_setattr(self, Vector2D, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Vector2D, 'thisown', 0)
        _swig_setattr(self, Vector2D,self.__class__,Vector2D)
_libgeom.Vector2D_swigregister(Vector2DPtr)

class Vector3D(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vector3D, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Vector3D, name)
    def __repr__(self):
        return "<C Vector3D instance at %s>" % (self.this,)
    __swig_setmethods__["x"] = _libgeom.Vector3D_x_set
    __swig_getmethods__["x"] = _libgeom.Vector3D_x_get
    if _newclass:x = property(_libgeom.Vector3D_x_get, _libgeom.Vector3D_x_set)
    __swig_setmethods__["y"] = _libgeom.Vector3D_y_set
    __swig_getmethods__["y"] = _libgeom.Vector3D_y_get
    if _newclass:y = property(_libgeom.Vector3D_y_get, _libgeom.Vector3D_y_set)
    __swig_setmethods__["z"] = _libgeom.Vector3D_z_set
    __swig_getmethods__["z"] = _libgeom.Vector3D_z_get
    if _newclass:z = property(_libgeom.Vector3D_z_get, _libgeom.Vector3D_z_set)
    __swig_setmethods__["label"] = _libgeom.Vector3D_label_set
    __swig_getmethods__["label"] = _libgeom.Vector3D_label_get
    if _newclass:label = property(_libgeom.Vector3D_label_get, _libgeom.Vector3D_label_set)
    __swig_setmethods__["code"] = _libgeom.Vector3D_code_set
    __swig_getmethods__["code"] = _libgeom.Vector3D_code_get
    if _newclass:code = property(_libgeom.Vector3D_code_get, _libgeom.Vector3D_code_set)
    def __init__(self, *args):
        _swig_setattr(self, Vector3D, 'this', _libgeom.new_Vector3D(*args))
        _swig_setattr(self, Vector3D, 'thisown', 1)
    def __del__(self, destroy=_libgeom.delete_Vector3D):
        try:
            if self.thisown: destroy(self)
        except: pass

class Vector3DPtr(Vector3D):
    def __init__(self, this):
        _swig_setattr(self, Vector3D, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Vector3D, 'thisown', 0)
        _swig_setattr(self, Vector3D,self.__class__,Vector3D)
_libgeom.Vector3D_swigregister(Vector3DPtr)


create_point_3D = _libgeom.create_point_3D

free_point_3D = _libgeom.free_point_3D

point_3D_set_label = _libgeom.point_3D_set_label

point_3D_get_label = _libgeom.point_3D_get_label

create_point_3D_array = _libgeom.create_point_3D_array

get_point_3D_ptr = _libgeom.get_point_3D_ptr

vector_copy_3D = _libgeom.vector_copy_3D

magnitude_3D = _libgeom.magnitude_3D

normalize_3D = _libgeom.normalize_3D

vector_scale_3D = _libgeom.vector_scale_3D

scalar_prod_3D = _libgeom.scalar_prod_3D

vector_prod_3D = _libgeom.vector_prod_3D

triple_prod_3D = _libgeom.triple_prod_3D

vector_sum_3D = _libgeom.vector_sum_3D

point_3D_translate = _libgeom.point_3D_translate

vector_diff_3D = _libgeom.vector_diff_3D

quad_properties_3D = _libgeom.quad_properties_3D

tetrahedron_properties_3D = _libgeom.tetrahedron_properties_3D

wedge_properties_3D = _libgeom.wedge_properties_3D

hexahedron_properties_3D = _libgeom.hexahedron_properties_3D

quad_area = _libgeom.quad_area

quad_cent = _libgeom.quad_cent

local_frame_3D = _libgeom.local_frame_3D

xyz_frame_3D = _libgeom.xyz_frame_3D
GEOM_H = _libgeom.GEOM_H
GPATH_NONE = _libgeom.GPATH_NONE
GPATH_LINE = _libgeom.GPATH_LINE
GPATH_ARC = _libgeom.GPATH_ARC
GPATH_BEZIER = _libgeom.GPATH_BEZIER
class GPathElement(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, GPathElement, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, GPathElement, name)
    def __repr__(self):
        return "<C GPathElement instance at %s>" % (self.this,)
    __swig_setmethods__["type"] = _libgeom.GPathElement_type_set
    __swig_getmethods__["type"] = _libgeom.GPathElement_type_get
    if _newclass:type = property(_libgeom.GPathElement_type_get, _libgeom.GPathElement_type_set)
    __swig_setmethods__["np"] = _libgeom.GPathElement_np_set
    __swig_getmethods__["np"] = _libgeom.GPathElement_np_get
    if _newclass:np = property(_libgeom.GPathElement_np_get, _libgeom.GPathElement_np_set)
    __swig_setmethods__["p"] = _libgeom.GPathElement_p_set
    __swig_getmethods__["p"] = _libgeom.GPathElement_p_get
    if _newclass:p = property(_libgeom.GPathElement_p_get, _libgeom.GPathElement_p_set)
    def __init__(self, *args):
        _swig_setattr(self, GPathElement, 'this', _libgeom.new_GPathElement(*args))
        _swig_setattr(self, GPathElement, 'thisown', 1)
    def __del__(self, destroy=_libgeom.delete_GPathElement):
        try:
            if self.thisown: destroy(self)
        except: pass

class GPathElementPtr(GPathElement):
    def __init__(self, this):
        _swig_setattr(self, GPathElement, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, GPathElement, 'thisown', 0)
        _swig_setattr(self, GPathElement,self.__class__,GPathElement)
_libgeom.GPathElement_swigregister(GPathElementPtr)

GPATH_MAX_ELEMENTS = _libgeom.GPATH_MAX_ELEMENTS
class GPathPolyLine(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, GPathPolyLine, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, GPathPolyLine, name)
    def __repr__(self):
        return "<C GPathPolyLine instance at %s>" % (self.this,)
    __swig_setmethods__["ne"] = _libgeom.GPathPolyLine_ne_set
    __swig_getmethods__["ne"] = _libgeom.GPathPolyLine_ne_get
    if _newclass:ne = property(_libgeom.GPathPolyLine_ne_get, _libgeom.GPathPolyLine_ne_set)
    __swig_setmethods__["pe"] = _libgeom.GPathPolyLine_pe_set
    __swig_getmethods__["pe"] = _libgeom.GPathPolyLine_pe_get
    if _newclass:pe = property(_libgeom.GPathPolyLine_pe_get, _libgeom.GPathPolyLine_pe_set)
    __swig_setmethods__["el"] = _libgeom.GPathPolyLine_el_set
    __swig_getmethods__["el"] = _libgeom.GPathPolyLine_el_get
    if _newclass:el = property(_libgeom.GPathPolyLine_el_get, _libgeom.GPathPolyLine_el_set)
    __swig_setmethods__["length"] = _libgeom.GPathPolyLine_length_set
    __swig_getmethods__["length"] = _libgeom.GPathPolyLine_length_get
    if _newclass:length = property(_libgeom.GPathPolyLine_length_get, _libgeom.GPathPolyLine_length_set)
    __swig_setmethods__["t_star"] = _libgeom.GPathPolyLine_t_star_set
    __swig_getmethods__["t_star"] = _libgeom.GPathPolyLine_t_star_get
    if _newclass:t_star = property(_libgeom.GPathPolyLine_t_star_get, _libgeom.GPathPolyLine_t_star_set)
    __swig_setmethods__["t0"] = _libgeom.GPathPolyLine_t0_set
    __swig_getmethods__["t0"] = _libgeom.GPathPolyLine_t0_get
    if _newclass:t0 = property(_libgeom.GPathPolyLine_t0_get, _libgeom.GPathPolyLine_t0_set)
    __swig_setmethods__["t1"] = _libgeom.GPathPolyLine_t1_set
    __swig_getmethods__["t1"] = _libgeom.GPathPolyLine_t1_get
    if _newclass:t1 = property(_libgeom.GPathPolyLine_t1_get, _libgeom.GPathPolyLine_t1_set)
    __swig_setmethods__["closed"] = _libgeom.GPathPolyLine_closed_set
    __swig_getmethods__["closed"] = _libgeom.GPathPolyLine_closed_get
    if _newclass:closed = property(_libgeom.GPathPolyLine_closed_get, _libgeom.GPathPolyLine_closed_set)
    def __init__(self, *args):
        _swig_setattr(self, GPathPolyLine, 'this', _libgeom.new_GPathPolyLine(*args))
        _swig_setattr(self, GPathPolyLine, 'thisown', 1)
    def __del__(self, destroy=_libgeom.delete_GPathPolyLine):
        try:
            if self.thisown: destroy(self)
        except: pass

class GPathPolyLinePtr(GPathPolyLine):
    def __init__(self, this):
        _swig_setattr(self, GPathPolyLine, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, GPathPolyLine, 'thisown', 0)
        _swig_setattr(self, GPathPolyLine,self.__class__,GPathPolyLine)
_libgeom.GPathPolyLine_swigregister(GPathPolyLinePtr)


gpath_init = _libgeom.gpath_init

gpath_add_element = _libgeom.gpath_add_element

gpath_element_length = _libgeom.gpath_element_length

gpath_write_element_to_file = _libgeom.gpath_write_element_to_file

gpath_write_all_elements_to_file = _libgeom.gpath_write_all_elements_to_file

gpath_scan_from_string_and_add_element = _libgeom.gpath_scan_from_string_and_add_element

gpath_append_polyline = _libgeom.gpath_append_polyline

gpath_get_first_point_on_path = _libgeom.gpath_get_first_point_on_path

gpath_get_last_point_on_path = _libgeom.gpath_get_last_point_on_path

gpath_polyline_set_subrange = _libgeom.gpath_polyline_set_subrange

gpath_polyline_eval = _libgeom.gpath_polyline_eval

gpath_polyline_translate = _libgeom.gpath_polyline_translate

coons_patch = _libgeom.coons_patch

TFI_3D_no_array = _libgeom.TFI_3D_no_array

TFI_3D = _libgeom.TFI_3D

line_translate = _libgeom.line_translate

line_eval = _libgeom.line_eval

line_length = _libgeom.line_length

arc_translate = _libgeom.arc_translate

arc_eval_position_and_length = _libgeom.arc_eval_position_and_length

arc_eval = _libgeom.arc_eval

arc_length = _libgeom.arc_length
GPATH_H_ALREADY_INCLUDED = _libgeom.GPATH_H_ALREADY_INCLUDED
BEZIER_H = _libgeom.BEZIER_H

bezier_3_eval = _libgeom.bezier_3_eval

bezier_3_deriv = _libgeom.bezier_3_deriv
MAX_BEZ_DEGREE = _libgeom.MAX_BEZ_DEGREE

bezier_eval = _libgeom.bezier_eval

bezier_translate = _libgeom.bezier_translate

bezier_length = _libgeom.bezier_length

bezier_3_spline = _libgeom.bezier_3_spline

