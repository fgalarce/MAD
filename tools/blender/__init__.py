bl_info = {
    "name": "INRIA mesh tools",
    "author": "Mauricio Ivan Portilla",
    "version": (0, 1),
    "blender": (2, 80, 0),
    "location": "View3D > Add > Mesh > New Object",
    "description": "Importer and exporter for INRIA's' <<.mesh>> format",
    "warning": "The current version cannot export 3D meshes properly",
    "doc_url": "",
    "category": "Import-Export",
}

import bpy
#from .import_inria import ImportSomeData
#from .DefBorders_inria import MESH_OT_define_edge_group
#from .export_inria import ExportSomeData
from .import_inria import *
from .DefBorders_inria import *
from .export_inria import *


def register():
    bpy.utils.register_class(VIEW3D_PT_cust1)
    bpy.utils.register_class(MESH_OT_define_edge_group)
    bpy.utils.register_class(MESH_OT_exp_edges)
    bpy.utils.register_class(ExportSomeData)
    bpy.types.TOPBAR_MT_file_export.append(menu_func_export)
    bpy.utils.register_class(ImportSomeData)
    bpy.types.TOPBAR_MT_file_import.append(menu_func_import)


def unregister():
    bpy.utils.unregister_class(VIEW3D_PT_cust1)
    bpy.utils.unregister_class(MESH_OT_define_edge_group)
    bpy.utils.unregister_class(MESH_OT_exp_edges)
    bpy.utils.unregister_class(ExportSomeData)
    bpy.types.TOPBAR_MT_file_export.remove(menu_func_export)
    bpy.utils.unregister_class(ImportSomeData)
    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import)
