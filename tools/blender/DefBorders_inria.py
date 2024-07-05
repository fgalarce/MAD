import bpy
from bpy import context
import os
#RECORDATORIO: ELIMINAR EL ARCHIVO QUE GENERA PARA EVITAR PROBLEMAS
print('---STARTED---')
edges_to_export = []
selected_edge = []
obj = None
dat = None
tittle = 'INRIA'
#selected_edge = [ed.vertices[:] for ed in dat.edges if ed.select]
#print(selected_edge[0])

#selected_edge = [ed.vertices[:] for ed in dat.edges if ed.select]
#print(selected_edge)
#for i in range(0,len(selected_edge)):
#    selected_edge[i] = [e1 + e2 for e1, e2 in zip(list(selected_edge[i]),[1,1])]
#print(selected_edge)
#edges_to_export = []

#    #Edges
#    nedge = len(dat.edges)
#    edge = []
#    for edges in dat.edges:
#        edge.append( tuple( list(edges.vertices[:])+[1] ) )
#        edge[-1] = tuple([e1 + e2 for e1, e2 in zip(edge[-1],[1,1,0])])

class MESH_OT_define_edge_group(bpy.types.Operator):
    bl_idname = "mesh.define_edge_group"
    bl_label = "Define edge group"
    
    preset_enum : bpy.props.EnumProperty(
        name = "",
        description = "Select a group",
        items = [('OP1',"1","Flux in"),
                 ('OP2',"2","Flux out"),
                 ('OP3',"3","Wall"),
                 ('OP4',"4","Obstacle"),
                 ('OP5',"DELETE","Exclude edges")]
    )
    
    def invoke(self, context, event):
        wm = context.window_manager
        return wm.invoke_props_dialog(self)
    
    def draw(self, context):
        layout = self.layout
        layout.prop(self, "preset_enum")
    
    def execute(self, context):
        #y??? Si no reseteo el EDIT MODE, no funca.
        bpy.ops.object.editmode_toggle()
        bpy.ops.object.editmode_toggle()

        obj = context.active_object
        dat = obj.data
        selected_edge = [ed.vertices[:] for ed in dat.edges if ed.select]
        for i in range(0,len(selected_edge)):
            selected_edge[i] = [e1 + e2 for e1, e2 in zip(list(selected_edge[i]),[1,1])]
        print('---Size of selection:')
        print(len(selected_edge))
        
        match self.preset_enum:
            #Salen con los índices corregidos de aquí
            case 'OP1':
                for ed in selected_edge:
                    rec = list(ed)+[1]
                    edges_to_export.append(rec)
                    #edges_to_export[-1] = tuple([e1 + e2 for e1, e2 in zip(edges_to_export[-1],[1,1,0])])
                    
            case 'OP2':
                for ed in selected_edge:
                    rec = list(ed)+[2]
                    edges_to_export.append(rec)
                    #edges_to_export[-1] = tuple([e1 + e2 for e1, e2 in zip(edges_to_export[-1],[1,1,0])])
            case 'OP3':
                for ed in selected_edge:
                    rec = list(ed)+[3]
                    edges_to_export.append(rec)
                    #edges_to_export[-1] = tuple([e1 + e2 for e1, e2 in zip(edges_to_export[-1],[1,1,0])])
            case 'OP4':
                print('OP4!!!')
                for ed in selected_edge:
                    rec = list(ed)+[4]
                    edges_to_export.append(rec)
                    #edges_to_export[-1] = tuple([e1 + e2 for e1, e2 in zip(edges_to_export[-1],[1,1,0])])
            case 'OP5':
                #NO FUNCIONA BIEN
                print('OP5!!!')
                for ed in edges_to_export:
                    #Remueve las aristas seleccionadas de la lista para exportar
                    for sed in selected_edge:
                        if ed[0:2] == sed:
                            print('equal')
                            edges_to_export.remove(ed)
                            print(len(edges_to_export))
        print('---Current size of list of edges to export:')
        print(len(edges_to_export))
        return {"FINISHED"}

class MESH_OT_exp_edges(bpy.types.Operator):
    bl_idname = "mesh.exp_edges"
    bl_label = "Export edges"
    
    def execute(self, context):
        pass
        print('Edges to export')
        print(edges_to_export)
        print('Selected')
        print(selected_edge)
        with open(os.path.join(os.path.dirname(bpy.data.filepath), "edges.txt"),"w") as f:
            #FUNCIONA :)
            f.write(str(len(edges_to_export))+'\n')
            for i in range(0,len(edges_to_export)):
                f.write(str(edges_to_export[i][0])+' '+str(edges_to_export[i][1])+' '+str(edges_to_export[i][2])+'\n')
        edges_to_export = []
        selected_edge = []
        obj = None
        dat = None
        return {"FINISHED"}

class VIEW3D_PT_cust1(bpy.types.Panel):
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = tittle
    bl_label = tittle+" mesh tools"
    
    def draw(self, context):
        row = self.layout.row()
        row.operator("mesh.define_edge_group", text = "Save in group")
        row = self.layout.row()
        row.operator("mesh.exp_edges", text = "Export groups")

