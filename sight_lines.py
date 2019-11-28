# to run in ArcGis Pro, Python Command Line

import arcpy
city_name = 'London'

# environment
arcpy.env.workspace = "C:/Users/g_filo01/sciebo/Scripts/Image of the City/Outputs/tmp/"+city_name
env = arcpy.env.workspace
arcpy.CheckOutExtension("3D")
arcpy.env.overwriteOutput = True

# folder from which loading
load_folder = "C:/Users/g_filo01/sciebo/Scripts/Image of the City/Outputs/"+city_name
aprx = arcpy.mp.ArcGISProject("CURRENT")

# buildings, obstructions files, observer points (nodes from where computing the sight_lines)
buildings = env+"/"+city_name+"_buildings_sight.shp"
obstructions = env+"/"+city_name+"_obstructions.shp"
observer_points =  load_folder+"/"+city_name+"_nodes.shp"

# transforming the obstructions layer in a 3d Layer and then multipatch
arcpy.MakeFeatureLayer_management(obstructions, "o_layer").getOutput(0)
arcpy.FeatureTo3DByAttribute_3d("o_layer", "o_layer3d", "height").getOutput(0)
arcpy.AddZInformation_3d("o_layer3d", 'Z_MAX', 'NO_FILTER')

''''
MANUAL PROCEDURE:
click on o_layer3d, properties, elevation, chooce relative to the ground (or absolute height), field: "base height"
appeareance, extrusion, type = base_height, field = Z_Max or height
''''
arcpy.Layer3DToFeatureClass_3d("o_layer3d", city_name+"_multiPatch", None, "ENABLE_COLORS_AND_TEXTURES").getOutput(0)

# sight_lines parameters
height_observer = "height"
heigth_buildings = "height"
direction = "NOT_OUTPUT_THE_DIRECTION"

# files
geoDB = "C:/Users/g_filo01/sciebo/Scripts/ArcGis/GeoDB.gdb"
arcpy.FeatureClassToGeodatabase_conversion(city_name+"_multiPatch", geoDB)
obstructions = geoDB+"/"+city_name+"_multiPatch"
sightlines_file = geoDB+"/"+city_name+"_sight_lines"

# construct sight-lines:
arcpy.ddd.ConstructSightLines(observer_points, buildings, sightlines_file, height_observer, heigth_buildings, None, 100, direction)

# remove line shortest than 200 m
with arcpy.da.UpdateCursor(sightlines_file, 'SHAPE@LENGTH') as cursor:
    for row in cursor:
        if row[0] < 300: cursor.deleteRow()
            
maximum = int(arcpy.GetCount_management(sightlines_file).getOutput(0))
quantity = 100000
cycles = int(maximum/quantity)+1
to_merge = []

# computing intervisibility sight-lines in blocks (to handle large files)

for i in range(1,cycles):

	print("cycle nr "+str(i)+" in "+str(cycles))
	if i == 1: l = 1
	else: l = i*quantity+1
	r = l+quantity-1
	if r > maximum: r = maximum
		
	selection = geoDB+"/"+city_name+"_sight_selection"+str(i)
	to_merge.append(selection)
	# selecting a limited amount of records, computing intervisibility and deleting non-visible lines
	arcpy.analysis.Select(sightlines_file, selection, "OID >="+str(l)+"And OID <="+str(r))
	arcpy.ddd.Intervisibility(selection, obstructions, visible_field = "visible")
	with arcpy.da.UpdateCursor(selection, "visible") as cursor:
		for row in cursor:
			if (row[0] == 0) | (row[0] == None): cursor.deleteRow()
	
# merging the file created above
arcpy.Merge_management(to_merge, geoDB+"/"+city_name+"_sightlines")
visible_sight_lines =  geoDB+"/"+city_name+"_sightlines"

# assigning nodeID-buildingID to the relative rows
arcpy.JoinField_management(visible_sight_lines , "OID_OBSERV", observer_points, "FID", "nodeID")
arcpy.JoinField_management(visible_sight_lines , "OID_TARGET", buildings, "FID", "buildingID")
arcpy.DeleteField_management(visible_sight_lines , ["OID_OBSERV", "OID_TARGET"])

# saving 
output = "C:/Users/g_filo01/sciebo/Scripts/Image of the City/Outputs/"+city_name
arcpy.FeatureClassToShapefile_conversion(visible_sight_lines , output)

# deleting temporary files
for i in to_merge:
    if arcpy.Exists(i): arcpy.Delete_management(i)
    
arcpy.Delete_management(sightlines_file)

