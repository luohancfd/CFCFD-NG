#!/usr/bin/env python
"""
Function to automatically convert e3prep output to OpenFoam mesh.
script performes the following tasks:
- combines the individual blocks into a single unstructured OpenFoam mesh
- stitch internal faces
- combine front and back faces (Top and Bottom in the e3prep blocks) to form an empty front and rear patch 
- combine individual boundaries into patches. The following boundaries are recognised and will be grouped automatically:
    - OF_inlet_1
    - OF_inlet_2
    - OF_outlet_1
    - OF_outlet_2
    - OF_wall_1
    - OF_wall_2
    - Anything else will retain name and be set as "patch"


Author: ingo jahn 04/09/2014
"""

import os as os
import numpy as np
import shutil as sh
from getopt import getopt
import sys as sys 

shortOptions = ""
longOptions = ["help", "job=", "create_0"]

def printUsage():
    print ""
    print "Usage: e3prepToFoam.py [--help] [--job=<jobFileName>] [--create_0]"
    return


def get_folders():
    print 'Obtaining directory from which code is executed'
    start_dir = os.getcwd()
    str2 = start_dir.split('/')
    n = len(str2)
    str3 = str2[n-1] 
    if str3 == 'e3prep':
        case_dir = os.path.dirname(start_dir)
        root_dir = os.path.dirname(case_dir)
        case_name = str2[n-2]
    else:
        case_dir = start_dir
        root_dir = os.path.dirname(case_dir) 
        case_name = str2[n-1]
    return root_dir, case_dir, start_dir, case_name 

def check_case_structure(case_dir):
    print 'Checking if correct OpenFOAM case structure exists'
    flag = 0
    if os.path.exists(case_dir+'/0') is not True:
        flag = 1
        print 'Error: Missing /0 directory'
    if os.path.exists(case_dir+'/system') is not True:
        flag = 1
        print 'Error: Missing /system directory'
    if os.path.exists(case_dir+'/constant') is not True:
        flag = 1
        print 'Error: Missing /constant directory'
    if os.path.exists(case_dir+'/constant/polyMesh') is not True:
        flag = 1
        print 'Error: Missing /constant/polyMesh directory'
    return flag


def face_index_to_string(ind):
    if ind == 0:
        return 'n'
    elif ind == 1: 
        return 'e'
    elif ind == 2:
        return 's'
    elif ind == 3: 
        return 'w'
    elif ind == 4:
        return 't'
    elif ind == 5:
        return 'b'
    else:
        print 'Error'
    return 

def get_job_config_data(job):
    print 'Extracting required variables from job.config'
    f = open((job + '.config'),'r')
    # find dimensions
    for line in f:
        if "dimensions" in line:
            temp = line.split()
            dimensions = int(temp[2])
            break   
    # find number of blocks
    for line in f:
        if "nblock" in line:
            temp = line.split()
            nblock = int(temp[2])
            break

    other_block = np.zeros((nblock,6))
    other_face = np.zeros((nblock,6))   
    # find connecting blocks
    f.seek(0,0)
    block = 0
    face = 0           
    for line in f:
        if "other_block" in line:
            temp = line.split()
            other_block[block,face] = int(temp[2])
            face = face + 1
            if dimensions == 2:
                if face == 4:
                    other_block[block,4] = -1
                    other_block[block,5] = -1
                    face = 0
                    block = block + 1
            elif dimensions == 3:
                if face == 6:
                    face = 0
                    block = block + 1
    # find connecting faces
    f.seek(0,0)
    block = 0
    face = 0           
    for line in f:
        if "other_face" in line:
            temp = line.split()
            if temp[2] == "none":
                other_face[block,face] = -1
            elif temp[2] == "north":
                other_face[block,face] = 0  
            elif temp[2] == "east":
                other_face[block,face] = 1 
            elif temp[2] == "south":
                other_face[block,face] = 2 
            elif temp[2] == "west":
                other_face[block,face] = 3 
            elif temp[2] == "top":
                other_face[block,face] = 4 
            elif temp[2] == "bottom":
                other_face[block,face] = 5 
            else:
                print "Error"
            face = face + 1
            if dimensions == 2:    
                if face == 4:
                    other_face[block,4] = -1
                    other_face[block,5] = -1
                    face = 0
                    block = block + 1
            elif dimensions == 3:
                if face == 6:
                    face = 0
                    block = block + 1   
    # find boundary labels
    f.seek(0,0)
    Label = [[None for i in range(6)] for i in range(nblock)]       
    for block in range(nblock):
        if dimensions == 2:
            for face in range(4):
                temp = find_boundary_info(f,block,face,"label")
                temp = temp.split()
                if len(temp) == 3:
                    Label[block][face] = temp[2]
                    # print "I was here", block, face, temp[2]
                else:
                    Label[block][face] = "EMPTY"
            Label[block][4] = "EMPTY"
            Label[block][5] = "EMPTY"
        if dimensions == 3:
            for face in range(6):
                temp = find_boundary_info(f,block,face,"label")
                temp = temp.split()
                if len(temp) == 3:
                    Label[block][face] = temp[2]
                    # print "I was here", block, face, temp[2]
                else:
                    Label[block][face] = "EMPTY"
    f.close()
    return (nblock, dimensions, other_block, other_face, Label)


def find_boundary_info(fp,block,face,lookup):
    if face == 0:
        phrase = ("[block/"+str(block)+"/face/north]")
    elif face == 1:
        phrase = ("[block/"+str(block)+"/face/east]")
    elif face == 2:
        phrase = ("[block/"+str(block)+"/face/south]")
    elif face == 3:
        phrase = ("[block/"+str(block)+"/face/west]")
    elif face == 4:
        phrase = ("[block/"+str(block)+"/face/top]")
    elif face == 5:
        phrase = ("[block/"+str(block)+"/face/bottom]")
    else:
        print "wrong face_index"
    fp.seek(0,0)
    for num, line in enumerate(fp, 1):
        if phrase in line:
            break
    for line in fp:
        if lookup in line:
            break  
    return line


def write_general_OpenFoam_header(fp):
    fp.write("/*--------------------------------*- C++ -*----------------------------------*\\\n")
    fp.write(" | ========                 |                                                 |\n")
    fp.write(" | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n")
    fp.write(" |  \\    /   O peration     | Version:  2.2.2                                 |\n")
    fp.write(" |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n")
    fp.write(" |    \\/     M anipulation  | This file generated by e3post.py                |\n")
    fp.write("\*---------------------------------------------------------------------------*/\n")
    fp.write("FoamFile\n")
    fp.write("{\n")
    fp.write("    version     2.0;\n")
    fp.write("    format      ascii;\n")
    return

def write_general_OpenFoam_bottom_round(fp):
    fp.write(");\n")
    fp.write("\n")
    fp.write("// ************************************************************************* //\n")
    return

def write_general_OpenFoam_bottom_curly(fp):
    fp.write("};\n")
    fp.write("\n")
    fp.write("// ************************************************************************* //\n")
    return

def write_createPatch_header(fp):
    #
    # ------------------- writing files now -----------------------------
    # points
    fp.write("    class       dictionary;\n")
    fp.write("    object      createPatchDict;\n")
    fp.write("}\n")
    fp.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
    fp.write("\n")
    fp.write("pointSync false;\n")
    fp.write("// Patches to create. \n")
    fp.write("patches \n")
    fp.write("(\n")
    return

def write_p_header(fp):
    #
    # ------------------- writing files now -----------------------------
    # points
    fp.write("    class       volScalarField;\n")
    fp.write("    location    \"0\";\n")
    fp.write("    object      p;\n")
    fp.write("}\n")
    fp.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
    fp.write("\n")
    fp.write("dimensions      [0 2 -2 0 0 0 0];\n")
    fp.write("\n")
    fp.write("internalField   uniform 0;\n")
    fp.write("\n")
    fp.write("boundaryField \n")
    fp.write("{ \n")
    return

def write_U_header(fp):
    #
    # ------------------- writing files now -----------------------------
    # points
    fp.write("    class       volVectorField;\n")
    fp.write("    location    \"0\";\n")
    fp.write("    object      U;\n")
    fp.write("}\n")
    fp.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
    fp.write("\n")
    fp.write("dimensions      [0 1 -1 0 0 0 0];\n")
    fp.write("\n")
    fp.write("internalField   uniform (0 0 0);\n")
    fp.write("\n")
    fp.write("boundaryField \n")
    fp.write("{ \n")
    return

def write_patches(fp,input_patch_str,output_name,output_type):
    fp.write("    {\n")
    fp.write(("        name " + output_name +";\n"));
    fp.write("        patchInfo\n")
    fp.write("        {\n")
    fp.write(("            type "+output_type+";\n"))
    fp.write("        }\n")
    fp.write("        constructFrom patches;\n")
    fp.write("        patches ("+input_patch_str+");\n")
    fp.write("    }\n")
    return 

def write_p_Boundary(fp,bname,btype):
    fp.write("    "+bname+"\n")
    fp.write("    {\n")
    if btype == "zeroGradient":
        fp.write("        type   zeroGradient; \n")
    elif btype == "empty":
        fp.write("        type   empty; \n")
    elif btype == "symmetry":
        fp.write("        type   symmetry; \n")
    elif btype == "fixedValue":
        fp.write("        type   fixedValue; \n")
        fp.write("        value   0; \n")
    else:
        print ("Boundary type, " + btype + " not recognised. Setting empty")
    fp.write("    }\n")
    return 

def write_U_Boundary(fp,bname,btype):
    fp.write("    "+bname+"\n")
    fp.write("    {\n")
    if btype == "zeroGradient":
        fp.write("        type   zeroGradient; \n")
    elif btype == "empty":
        fp.write("        type   empty; \n")
    elif btype == "symmetry":
        fp.write("        type   symmetry; \n")
    elif btype == "fixedValue":
        fp.write("        type   fixedValue; \n")
        fp.write("        value   uniform (0 0 0); \n")
    else:
        print ("Boundary type, " + btype + " not recognised. Setting empty")
    fp.write("    }\n")
    return 
    
def combine_faces(case_dir,start_dir,patch_str,patch_name,patch_type):
    file_createPatchDict = "createPatchDict"
    file_createPatchDict = os.path.join((case_dir+ '/system'), file_createPatchDict)
    OFFile0 = open(file_createPatchDict, "wb")
    
    write_general_OpenFoam_header(OFFile0)
    write_createPatch_header(OFFile0)
    write_patches(OFFile0,patch_str,patch_name,patch_type)
    write_general_OpenFoam_bottom_round(OFFile0)
    OFFile0.close()
    print "createPatchDict has been written. \n"
    # execute createPatch
    os.chdir(case_dir)
    os.system('createPatch -overwrite')
    # move back to starting_directory
    os.chdir(start_dir) 
    print ("The following boundaries" +patch_str+ " have been combined to form Patch: " +patch_name+ " with the type: " + patch_type)
    return 

def check_for_undefined_faces(case_dir,nblock):
    file_name = "boundary"
    file_name = os.path.join((case_dir+ '/constant/polyMesh/'), file_name)
    File = open(file_name, "r") 
    String = []   
    for n in range(nblock):
        for line in File:
            if ('n'+'%04d' % n) in line:
                String.append(line + '; ') 
            if ('e'+'%04d' % n) in line:
                String.append(line + '; ') 
            if ('s'+'%04d' % n) in line:
                String.append(line + '; ') 
            if ('w'+'%04d' % n) in line:
                String.append(line + '; ') 
        File.seek(0,0)       
    File.close()
    return String


def check_for_undefined_labels(patch_Label):
    A = [item for sublist in patch_Label for item in sublist]
    A = set(A)
    String = ['EMPTY']
    for i in range(10):
        String.append("OF_inlet_"+'%02d' % i)
        String.append("OF_outlet_"+'%02d' % i)
        String.append("OF_wall_"+'%02d' % i)  
        String.append("OF_symmetry_"+'%02d' % i) 
    String = set(String)

    return list(A.difference(String))


class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def main(uoDict):
    # main file to be executed 
    jobName = uoDict.get("--job", "test")

    # strip .py extension form jobName
    jobName = jobName.split('.')
    jobName = jobName[0]

    # establish case, root, and start directory
    root_dir, case_dir, start_dir, case_name  = get_folders()

    # check that correct directory structure exists
    dir_flag = check_case_structure(case_dir)
    if dir_flag == 1:
        raise MyError('ERROR: Incorrect Directory Structure. e3preToFoam must be run inside an OpenFoam case with appropriate sub-directories. \n See error message above and create missing folders or copy from existing case.')

   
    # change into e3prep directory    
    os.chdir((case_dir+'/e3prep'))

    # get data from job.config
    nblock, dimensions, other_block, other_face, patch_Label = get_job_config_data(jobName)

    # run e3post to generate /foam folder containing meshes for respective block
    os.system(("e3post.py --job=" + jobName + " --OpenFoam"))

    print 'e3post has been executed and individual foam meshes have been generated for each block \n \n '

    # merging individual blocks
    print ('Working on Case = '+case_name)

    ## move data currently in /polMesh
    #sh.move('polyMesh','polyMesh_old') 
    sh.rmtree(case_dir + '/constant/polyMesh')

    # create case file for slave_mesh
    sh.copytree(case_dir,(root_dir+'/slave_mesh'))       

    # copy Master mesh data into required folder
    sh.copytree((case_dir+'/e3prep/foam/b0000/'),case_dir+'/constant/polyMesh')
        
    for block in range(nblock-1):
        # copy correct slave_mesh into slabe_mesh case
        sh.copytree((case_dir+'/e3prep/foam/b'+ '%04d' % (block+1) + '/'),root_dir+'/slave_mesh/constant/polyMesh')

        # execute mergeMeshes command
        os.chdir(root_dir)
        flag = os.system('mergeMeshes -overwrite ' + case_name + ' slave_mesh') 
        if flag != 0:
            print  'Error with mergeMeshes. \n try running of230 to load OpenFOAM module'
            sh.rmtree(root_dir+'/slave_mesh') # removing slave_mesh directory before exiting
            os.chdir(start_dir)
            sys.exit(1)

        print ('Block ' + '%04d' % block + ' and ' + '%04d' % (block+1) + ' have been merged.')   
        # remove polyMesh from slave_mesh
        sh.rmtree(root_dir+'/slave_mesh/constant/polyMesh')

    # remove slave_mesh
    sh.rmtree(root_dir+'/slave_mesh')
    # move back to starting_directory
    os.chdir(start_dir)

    print "Merging of meshes complete. \n \n "


    #identify number of block connections
    interfaces = len(other_block[np.where(other_block != -1)]) # counts 2 x internal connections, as seen by other blocks
    if interfaces > 0:

        # move /0 directory
        sh.move((case_dir + '/0'), case_dir + '/temp')

        while True:
            (block,face) = np.where(other_block != -1)
            if len(block) == 0:
                break

            print (block,face)
            o_block = other_block[block[0],face[0]]
            o_face = other_face[block[0],face[0]]

            current_facename = (face_index_to_string(face[0]) + '%04d' % block[0])
            other_facename = (face_index_to_string(o_face) + '%04d' % o_block)
        
            # print (current_facename,other_facename)

            # overwrite matcing face in other block
            other_block[o_block,o_face] = -1
            other_face[o_block,o_face] = -1
            other_block[block[0],face[0]] = -1
    
            # execute stitchMesh command
            os.chdir(case_dir)
            os.system('stitchMesh -overwrite -perfect ' + current_facename + ' ' + other_facename) 

            print ('Face ' + current_facename + ' and ' + other_facename + ' have been stitched.')   

            # move back to starting_directory
            os.chdir(start_dir)        

            # print (other_block, other_face)

        # move /0 directory back
        sh.move((case_dir + '/temp'), case_dir + '/0')    

    print "Stitching of internal Faces complete. \n \n"

    # do automatic patch combination
    # top and bottom faces
    if dimensions == 2:
        patch_str = ' '
        for i in range(nblock):
            patch_str = (patch_str+' b'+ '%04d' % i + ' t' + '%04d' % i)
        patch_name = 'FrontBack'
        patch_type = 'empty'

        combine_faces(case_dir,start_dir,patch_str,patch_name,patch_type)

    # combine patches, based on block label.
    # Following labels are supported: 
    # OF_inlet_00, OF_inlet_01, OF_inlet_02 (up to 09)
    # OF_outlet_00, OF_outlet_01, OF_outlet_02 (up to 09)
    # OF_wall_00, OF_wall_01, OF_wall_02 (up to 09)
    # OF_symmetry_00, OF_symmetry_01, OF_symmetry_02 (up to 09)

    N_list_in = []
    N_list_out = []
    N_list_wall = []
    N_list_sym = []

    for i in range(10):
        in_n = ("OF_inlet_"+'%02d' % i)
        out_n = ("OF_outlet_"+'%02d' % i)
        wall_n = ("OF_wall_"+'%02d' % i)  
        sym_n = ("OF_symmetry_"+'%02d' % i)      
        
        inlet_str = ""
        outlet_str = ""
        wall_str = ""
        sym_str = ""
        for block in range(nblock):
            L_block = patch_Label[block]
            #print L_block
            i_ind = [n for n, s in enumerate(L_block) if in_n in s]
            o_ind = [n for n, s in enumerate(L_block) if out_n in s]
            w_ind = [n for n, s in enumerate(L_block) if wall_n in s]
            s_ind = [n for n, s in enumerate(L_block) if sym_n in s]
            #print i_ind != []
            #print o_ind != []

            if i_ind != []:
                for n in i_ind:
                    if n == 0:
                        inlet_str = (inlet_str + ' n' +'%04d' % block) 
                    if n == 1:
                        inlet_str = (inlet_str + ' e' +'%04d' % block) 
                    if n == 2:
                        inlet_str = (inlet_str + ' s' +'%04d' % block) 
                    if n == 3:
                        inlet_str = (inlet_str + ' w' +'%04d' % block) 
                    if n == 4:
                        inlet_str = (inlet_str + ' t' +'%04d' % block) 
                    if n == 5:
                        inlet_str = (inlet_str + ' b' +'%04d' % block) 
            if o_ind != []:
                for n in o_ind:
                    if n == 0:
                        outlet_str = (outlet_str + ' n' +'%04d' % block) 
                    if n == 1:
                        outlet_str = (outlet_str + ' e' +'%04d' % block) 
                    if n == 2:
                        outlet_str = (outlet_str + ' s' +'%04d' % block) 
                    if n == 3:
                        outlet_str = (outlet_str + ' w' +'%04d' % block) 
                    if n == 4:
                        outlet_str = (outlet_str + ' t' +'%04d' % block) 
                    if n == 5:
                        outlet_str = (outlet_str + ' b' +'%04d' % block)   
            if w_ind != []:
                for n in w_ind:
                    if n == 0:
                        wall_str = (wall_str + ' n' +'%04d' % block) 
                    if n == 1:
                        wall_str = (wall_str + ' e' +'%04d' % block) 
                    if n == 2:
                        wall_str = (wall_str + ' s' +'%04d' % block) 
                    if n == 3:
                        wall_str = (wall_str + ' w' +'%04d' % block) 
                    if n == 4:
                        wall_str = (wall_str + ' t' +'%04d' % block) 
                    if n == 5:
                        wall_str = (wall_str + ' b' +'%04d' % block) 
            if s_ind != []:
                for n in s_ind:
                    if n == 0:
                        sym_str = (sym_str + ' n' +'%04d' % block) 
                    if n == 1:
                        sym_str = (sym_str + ' e' +'%04d' % block) 
                    if n == 2:
                        sym_str = (sym_str + ' s' +'%04d' % block) 
                    if n == 3:
                        sym_str = (sym_str + ' w' +'%04d' % block) 
                    if n == 4:
                        sym_str = (sym_str + ' t' +'%04d' % block) 
                    if n == 5:
                        sym_str = (sym_str + ' b' +'%04d' % block)     


        print inlet_str, outlet_str, wall_str, sym_str
        if inlet_str != "":
            combine_faces(case_dir,start_dir,inlet_str,in_n,'patch')
            N_list_in.append(in_n)
        if outlet_str != "":
            combine_faces(case_dir,start_dir,outlet_str,out_n,'patch')
            N_list_out.append(out_n)
        if wall_str != "":
            combine_faces(case_dir,start_dir,wall_str,wall_n,'wall')
            N_list_wall.append(wall_n)
        if sym_str != "":
            combine_faces(case_dir,start_dir,sym_str,sym_n,'symmetry')
            N_list_sym.append(sym_n)


    # Re-order numbering of faces/cells for numerical efficiency
    # execute renumberMesh
    os.chdir(case_dir)
    os.system('renumberMesh -overwrite')
    # move back to starting_directory
    os.chdir(start_dir) 
    

    # check if there are patches remaining that havent been defined. 
    String1 = check_for_undefined_faces(case_dir,nblock)
    String2 = check_for_undefined_labels(patch_Label) 

    if not(String1 == []):
        print 'WARNING: Not all external boundaries were defined in e3prep'
        print 'Check these faces: ', String1    
        print "\n \n"

    if not(String2 == []):
        print 'WARNING: labels used to define boundary faces do not follow standard OF_names'
        print 'Check these labels: ', String2    
        print "\n \n"

    print "\n \n"
    print "The multi-block mesh created by e3prep.py has been converted into a single Polymesh for use with OpenFoam."
    print "\n \n"


    # Option to create template entries for /0.
    if uoDict.has_key("--create_0"): 

        # check if /0/p file exists
        if os.path.isfile(case_dir+'/0/'+'p') == 1:
            sh.copyfile(case_dir+'/0/'+'p', case_dir+'/0/'+'p.bak')
            print "WARNING: Existing copy of /0/p has been copied to /0/p.bak"
        # check if /0/U file exists
        if os.path.isfile(case_dir+'/0/'+'U') == 1:
            sh.copyfile(case_dir+'/0/'+'U', case_dir+'/0/'+'U.bak')
            print "WARNING: Existing copy of /0/U has been copied to /0/U.bak"
 
        # U and p template are created. The others can be duplicated form these
        file_name = "p"
        file_name = os.path.join((case_dir+ '/0/'), file_name)
        OFFile0 = open(file_name, "wb")
    
        write_general_OpenFoam_header(OFFile0)
        write_p_header(OFFile0)

        if dimensions == 2:
            write_p_Boundary(OFFile0,'FrontBack','empty')
        for n in range(len(N_list_in)):
            write_p_Boundary(OFFile0,N_list_in[n],'zeroGradient')
        for n in range(len(N_list_out)):
            write_p_Boundary(OFFile0,N_list_out[n],'zeroGradient')
        for n in range(len(N_list_wall)):
            write_p_Boundary(OFFile0,N_list_wall[n],'zeroGradient')
        for n in range(len(N_list_sym)):
            write_p_Boundary(OFFile0,N_list_sym[n],'symmetry')

        write_general_OpenFoam_bottom_curly(OFFile0)
        OFFile0.close()
        print "/0/p has been written. \n"

        file_name = "U"
        file_name = os.path.join((case_dir+ '/0/'), file_name)
        OFFile0 = open(file_name, "wb")
    
        write_general_OpenFoam_header(OFFile0)
        write_U_header(OFFile0)
        if dimensions == 2:
            write_U_Boundary(OFFile0,'FrontBack','empty')
        for n in range(len(N_list_in)):
            write_U_Boundary(OFFile0,N_list_in[n],'fixedValue')
        for n in range(len(N_list_out)):
            write_U_Boundary(OFFile0,N_list_out[n],'zeroGradient')
        for n in range(len(N_list_wall)):
            write_U_Boundary(OFFile0,N_list_wall[n],'zeroGradient')
        for n in range(len(N_list_sym)):
            write_U_Boundary(OFFile0,N_list_sym[n],'symmetry')

        write_general_OpenFoam_bottom_curly(OFFile0)
        OFFile0.close()
        print "/0/U has been written. \n"


if __name__ == "__main__":
    userOptions = getopt(sys.argv[1:], shortOptions, longOptions)
    uoDict = dict(userOptions[0])
    
    if len(userOptions[0]) == 0 or uoDict.has_key("--help"):
        printUsage()
        sys.exit(1)
    
    try:
        main(uoDict)
    except MyError as e:
        print "This run of e3prepToFoam.py has gone bad. \n"
        print e.value
