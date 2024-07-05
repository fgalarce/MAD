import bpy
tittle = 'INRIA'





def read_some_data(context, filepath, param):
    print("running read_some_data...")
    f = open(filepath, 'r', encoding='utf-8')
    data = f.read()
    f.close()

    with open(filepath, 'r') as file:
        #SE CUENTAN LAS LÍNEAS
        # n -> número de líneas
        for n, line in enumerate(file):
            pass
        n = n + 1
        file.seek(0)
        lines = file.readlines()

    #Dimension ???
    dimension = int()

    #Parámetros geométricos
    Vert_range = []
    RVert_range = []
    Edge_range = []
    Corner_range = []
    Vert_range = []
    Trian_range = []
    head = 2
    
    def get_range(i,head,num):
        # Dado un título, devuelve el rango en el que se encuentra la lista
        #'i' representa el número de una línea
        x_range = [0,0]
        x_range[0] = i+head
        x_range[1] = x_range[0] + num
        return x_range
    
    def line2list(str):
        # Dada una línea la convierte en una lista
        pass
    
    for i in range(0, n):
        title = lines[i].strip()
        
        try:
            num = int(lines[i+1].strip())
        except:
            num = 0
        
        match title:
            case "Vertices":
                Vert_range = get_range(i,head,num)
            case "Corners":
                Corner_range = get_range(i,head,num)
            case "RequiredVertices":
                RVert_range = get_range(i,head,num)
            case "Edges":
                Edge_range = get_range(i,head,num)
            case "Triangles":
                Trian_range = get_range(i,head,num)
            case "End":
                pass

    # Se guardan los vértices y las caras

    verts = []
    for i in range(Vert_range[0],Vert_range[1]):
        verts.append(tuple([float(x) for x in lines[i].strip().split(' ')][0:3]))

    edges = []
    for i in range(Edge_range[0],Edge_range[1]):
        # 1 -> IN
        # 2 -> OUT
        # 3 -> PAREDES
        # 4 -> BORDE DEL CILINDRO
        if [int(x) for x in lines[i].strip().split(' ')][-1] == 1:
            edges.append([int(x) for x in lines[i].strip().split(' ')][0:2])
            # A continuación se corrigen los índices para que empiecen en 0
            edges[-1] = tuple([e1 - e2 for e1, e2 in zip(edges[-1],[1,1])])

    faces = []
    print('VertRange:')
    print(Vert_range)
    for i in range(Trian_range[0],Trian_range[1]):
        # 5 -> ??? Todas las caras tienen este valor
        if [int(x) for x in lines[i].strip().split(' ')][-1] > 0:
            faces.append([int(x) for x in lines[i].strip().split(' ')][:])
            # A continuación se corrigen los índices para que empiecen en 0
            faces[-1] = tuple([e1 - e2 for e1, e2 in zip(faces[-1],[1,1,1])])

    name = "f1"
    mesh = bpy.data.meshes.new(name)
    object = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(object)
    mesh.from_pydata(verts,edges,faces)

    return {'FINISHED'}


# ImportHelper is a helper class, defines filename and
# invoke() function which calls the file selector.
from bpy_extras.io_utils import ImportHelper
from bpy.props import StringProperty, BoolProperty, EnumProperty
from bpy.types import Operator


class ImportSomeData(Operator, ImportHelper):
    """Import INRIA .mesh file"""
    bl_idname = "import_test.some_data"  # important since its how bpy.ops.import_test.some_data is constructed
    bl_label = "Import Some Data"

    # ImportHelper mix-in class uses this.
    filename_ext = ".mesh"

    filter_glob: StringProperty(
        default="*.mesh",
        options={'HIDDEN'},
        maxlen=255,  # Max internal buffer length, longer would be clamped.
    )

    # List of operator properties, the attributes will be assigned
    # to the class instance from the operator settings before calling.
    use_setting: BoolProperty(
        name="Example Boolean",
        description="Example Tooltip",
        default=True,
    )

    type: EnumProperty(
        name="Example Enum",
        description="Choose between two items",
        items=(
            ('OPT_A', "First Option", "Description one"),
            ('OPT_B', "Second Option", "Description two"),
        ),
        default='OPT_A',
    )

    def execute(self, context):
        return read_some_data(context, self.filepath, self.use_setting)


# Only needed if you want to add into a dynamic menu.
def menu_func_import(self, context):
    self.layout.operator(ImportSomeData.bl_idname, text=tittle+" .mesh")

