import bpy
from bpy import context
import os

#REVISAR CÓMO SE EXPORTAN LOS EDGES
tittle='INRIA'
dsp = 3*'\n'  # Default new lines
def autowrite_elem(file, list_tuples):
    try:
        n = len(list_tuples[0])
        for elem in list_tuples:
            file.write('\n')
            for i in range(0,n):
                if elem[i] != 0.0:
                    file.write(str(elem[i]))
                else:
                    file.write(str(0))
                if i < n-1:
                    file.write(' ')
            if elem == list_tuples[-1]:
                file.write(dsp)
    except:
        for elem in list_tuples:
            file.write('\n'+str(elem))
        if elem == list_tuples[-1]:
            file.write(dsp)
    


def write_some_data(context, filepath, use_some_setting):
    print("running write_some_data...")
    
    ### PARÁMETROS (AÚN NECESITO INVESTIGAR SU IMPORTANCIA)
    dimension = 2
    #if use_some_setting:
    #    dimension = 3 
    
    corners = [(1), (2), (3), (4)]  # ???
    req_vertices = [(1), (2), (3), (4)]  # ???
    ###

    #Se selecciona el objeto y se extrae su información geométrica
    obj = context.active_object
    dat = obj.data

    #vert-> Vértices
    #nvert -> Num. de vértices
    #edge-> Aristas
    #nedge -> Num. de aristas
    #face-> Caras
    #nface -> Num. de caras

    #Vertices
    vert = [tuple(i.co) for i in dat.vertices]
    nvert = len(vert)

    #Edges
    #nedge = len(dat.edges)
    #edge = []
    #for edges in dat.edges:
    #    edge.append( tuple( list(edges.vertices[:])+[1] ) )
    #    edge[-1] = tuple([e1 + e2 for e1, e2 in zip(edge[-1],[1,1,0])])
    
    #Faces
    nface = len(dat.polygons)
    face = []
    for faces in dat.polygons:
        face.append( tuple( list(faces.vertices[:])+[5] ) )
        face[-1] = tuple([e1 + e2 for e1, e2 in zip(face[-1],[1,1,1,0])])
    
    
    f = open(filepath, 'w', encoding='utf-8')
    f.write('MeshVersionFormatted 2'+dsp+'Dimension '+str(dimension)+dsp)
    
    f.write('Vertices\n'+str(nvert))
    autowrite_elem(f, vert)
    
    f.write('Corners\n'+str(len(corners)))
    autowrite_elem(f, corners)
        
    f.write('RequiredVertices\n'+str(len(req_vertices)))
    autowrite_elem(f, req_vertices)
    
    f.write('Edges\n')
    with open(os.path.join(os.path.dirname(bpy.data.filepath), "edges.txt"), 'r') as file:
        #SE CUENTAN LAS LÍNEAS
        # n -> número de líneas
        for n, line in enumerate(file):
            pass
        n = n + 1
        file.seek(0)
        lines = file.readlines()
        for i in range(0,n):
            f.write(lines[i])
    f.write(2*'\n')  # El mejor parche que se me ocurre xD
    
    f.write('Triangles\n'+str(nface))
    autowrite_elem(f, face)
    
    f.write('End\n')
    f.close()

    return {'FINISHED'}


# ExportHelper is a helper class, defines filename and
# invoke() function which calls the file selector.
from bpy_extras.io_utils import ExportHelper
from bpy.props import StringProperty, BoolProperty, EnumProperty
from bpy.types import Operator


class ExportSomeData(Operator, ExportHelper):
    """Export Ansys .mesh (Not sure if works in other software)"""
    bl_idname = "export_test.some_data"  # important since its how bpy.ops.import_test.some_data is constructed
    bl_label = "Export Some Data"

    # ExportHelper mix-in class uses this.
    filename_ext = ".mesh"

    filter_glob: StringProperty(
        default="*.mesh",
        options={'HIDDEN'},
        maxlen=255,  # Max internal buffer length, longer would be clamped.
    )

    # List of operator properties, the attributes will be assigned
    # to the class instance from the operator settings before calling.

    use_setting: BoolProperty(
        name="3D?",
        description="Example Tooltip",
        default=False,
    )

    def execute(self, context):
        return write_some_data(context, self.filepath, self.use_setting)


# Only needed if you want to add into a dynamic menu
def menu_func_export(self, context):
    self.layout.operator(ExportSomeData.bl_idname, text=tittle+" .mesh")

