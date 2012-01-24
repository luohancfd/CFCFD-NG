# pyglet-bits.py
# Code extracted from e3prep.py 30-Sep-2008

try:
    import pyglet
    from pyglet import gl as opengl
except:
    print "Could not load pyglet (OpenGL) module."
    opengl = None

# --------------------------------------------------------------------

def render_to_opengl():
    window = pyglet.window.Window()
    batch = pyglet.graphics.Batch()
    @window.event
    def on_draw():
        opengl.glClear(opengl.GL_COLOR_BUFFER_BIT)
        opengl.glLoadIdentity()
        # camera_pos = Vector3(0.0, 0.0, 5.0)
        # look_at = Vector3(0.0, 0.0, -100.0)
        # up_vector = Vector3(0.0, 1.0, 0.0)
        #opengl.gluLookAt(camera_pos.x, camera_pos.y, camera_pos.z,
        #                 look_at.x, look_at.y, look_at.z,
        #                 up_vector.x, up_vector.y, up_vector.z)
        left = -2.0; right = 2.0; bot = -2.0; top = 2.0
        near = -50.0; far = 50.0
        opengl.glOrtho(left, right, bot, top, near, far)
        opengl.glScalef(200.0, 200.0, 200.0)
        opengl.glTranslatef(1.0, 1.0, 0.0)
        batch.draw()
    #
    for blk in Block.blockList:
        # Render bottom face as a mesh of line segments
        ni = blk.nni; nj = blk.nnj
        nvertices = (ni+1) * (nj+1)
        dr = 1.0/ni; ds = 1.0/nj
        # Assemble a list of vertex locations.
        vertex_pos_list = []
        for j in range(0,nj+1):
            for i in range(0,ni+1):
                r = dr * i; s = ds * j
                p = blk.parametric_volume.eval(r, s, 0.0)
                vertex_pos_list.append(p.x)
                vertex_pos_list.append(p.y)
                vertex_pos_list.append(p.z)
        vertex_colour_list = [255, 255, 255,] * nvertices
        # Now, assemble the list of vertex indices that define the line segments
        vertex_index_list = []
        # line segments in i-direction
        for i in range(0,ni):
            for j in range(0,nj+1):
                vertex_index_list.append((ni+1)*j+i)
                vertex_index_list.append((ni+1)*j+i+1)
        # line segments in j-direction
        for j in range(0,nj):
            for i in range(0,ni+1):
                vertex_index_list.append((ni+1)*j+i)
                vertex_index_list.append((ni+1)*(j+1)+i)
        # print "nvertices=", nvertices
        # print "vertex_pos_list=", vertex_pos_list, len(vertex_pos_list)
        # print "vertex_index_list=", vertex_index_list, len(vertex_index_list)
        # print "vertex_colour_list=", vertex_colour_list, len(vertex_colour_list)
        vl = batch.add_indexed(nvertices, opengl.GL_LINES, None,
                               vertex_index_list, ('v3f', vertex_pos_list), 
                               ('c3B', vertex_colour_list))
    pyglet.app.run()
    return

#------------------------------------------------------------------------------
        #
        if opengl and uoDict.has_key("--do-opengl") and gdata.dimensions == 3:
            render_to_opengl()
