import arcpy

arcpy.env.workspace = "C:/Users/g_filo01/sciebo/Scripts/Image of the City/Outputs/tmp"
env = arcpy.env.workspace
arcpy.CheckOutExtension("3D")
arcpy.env.overwriteOutput = True

# output folder from which loading
city_name = 'London'
load_folder = "C:/Users/g_filo01/sciebo/Scripts/Image of the City/Outputs/"+city_name
aprx = arcpy.mp.ArcGISProject("CURRENT")

buildings = env+"/"+city_name+"_buildings_sight.shp"
obstructions = env+"/"+city_name+"_obstructions.shp"

arcpy.MakeFeatureLayer_management(obstructions, "o_layer").getOutput(0)
arcpy.FeatureTo3DByAttribute_3d("o_layer", "o_layer3d", "height").getOutput(0)
arcpy.AddZInformation_3d("o_layer3d", 'Z_MAX', 'NO_FILTER')

# click on o_layer3d, properties, elevation, chooce relative to the ground (or absolute height), field: "base height"
# appeareance, extrusion, type = base_height, field = Z_Max or height

arcpy.Layer3DToFeatureClass_3d("o_layer3d", city_name+"_multiPatch", None, "ENABLE_COLORS_AND_TEXTURES").getOutput(0)

height_observer = "height"
heigth_buildings = "r_height"
direction = "NOT_OUTPUT_THE_DIRECTION"
observer_points =  env+"/"+city_name+"_nodes_cleaned_CCZ.shp"



geoDB = "C:/Users/g_filo01/sciebo/Scripts/ArcGis/GeoDB.gdb"
arcpy.FeatureClassToGeodatabase_conversion(city_name+"_multiPatch", geoDB)
obstructions = geoDB+"/"+city_name+"_multiPatch"
sightlines_file = geoDB+"/"+city_name+"_sightlines"
arcpy.ddd.ConstructSightLines(observer_points, buildings, sightlines_file, height_observer, heigth_buildings, None, 30, direction)

with arcpy.da.UpdateCursor(sightlines_file, 'SHAPE@LENGTH') as cursor:
    for row in cursor:
        if row[0] < 200:
            cursor.deleteRow()

maximum = int(arcpy.GetCount_management(sightlines_file).getOutput(0))
quantity = 1000000
cycles = int(maximum/quantity)+1
to_merge = []

for i in range(1,cycles):
	print("cycle nr "+str(i)+" in "+str(cycles))
	
	if i == 1: l = 1
	else: l = i*quantity+1
	r = i+quantity-1
	if r > maximum: r = maximum
		
	selection = geoDB+"/"+city_name+"_sight_selection"+str(i)
	to_merge.append(selection)
	arcpy.analysis.Select(sightlines_file, selection, "OID >="+str(l)+"And OID <="+str(r))
	arcpy.ddd.Intervisibility(selection, obstructions, visible_field = "Visible")
	with arcpy.da.UpdateCursor(selection, "Visible") as cursor:
		for row in cursor:
			if (row[0] == None): 
				cursor.deleteRow()
	

arcpy.Merge_management(to_merge, geoDB+"/"+city_name+"_visibile_sight_lines")
			
			
			
			
			
			
arcpy.ddd.Intervisibility(sightlines_file, obstructions, visible_field = "Visibility")

with arcpy.da.UpdateCursor(sightlines, "Visibility") as cursor:
    for row in cursor:
        if (row[0] == 0) | (row[0] == None): 
            cursor.deleteRow()
			
			
output = "C:/Users/g_filo01/sciebo/Scripts/Image of the City/Outputs/"+city_name

arcpy.JoinField_management(sightlines, "OID_OBSERV", observer_points, "FID", "nodeID")
arcpy.JoinField_management(sightlines, "OID_TARGET", buildings, "FID", "buildingID")
arcpy.DeleteField_management(sightlines, ["OID_OBSERV", "OID_TARGET"])
arcpy.FeatureClassToShapefile_conversion(sightlines, output)


